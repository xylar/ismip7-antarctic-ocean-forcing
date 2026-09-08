"""
Generate and run the single-timestep MALI melt diagnostics.

The calibration needs the melt field MALI produces for a given ocean state and
parameter value, and nothing else: no ice dynamics, no thermal evolution, no
transient.  So each member is a MALI run with the velocity solver off that
takes one short timestep and writes the melt.

That is enough because melt is computed from the geometry and the forcing
alone.  MALI evaluates it inside the timestep rather than in the initial
diagnostic
solve, so a zero-length run would produce nothing -- hence one step rather than
none.

Melt is exactly linear in ``K`` (findings F1), so the parameter ensemble does
not need a run per parameter value: one run per ocean state at a reference
``K`` suffices and the rest follow by scaling.
:func:`linearity_check_runs` builds the handful of extra runs that verify that
rather than assuming it.
"""

from __future__ import annotations

import os
from collections import OrderedDict

import xarray as xr

#: namelist options set for every melt-diagnostic run
DIAGNOSTIC_NAMELIST = {
    # melt does not depend on velocity, and skipping it avoids needing Albany
    'config_velocity_solver': "'none'",
    # no thermal evolution, and no grounded basal melt
    'config_thermal_solver': "'none'",
    'config_thermal_calculate_bmb': '.false.',
    # nothing should evolve the geometry
    'config_calving': "'none'",
    'config_front_mass_bal_grounded': "'none'",
    # the ISMIP7 thermal forcing is already extrapolated into cavities on the
    # ISMIP grid, and is gap-free after remapping (findings F5)
    'config_ocean_data_extrapolation': '.false.',
    'config_basal_mass_bal_float': "'ismip7'",
}


def read_namelist(path):
    """
    Read an MPAS namelist into ``{record: {option: value}}``.

    Values are kept as their literal strings, so writing back an unmodified
    namelist reproduces it.

    Parameters
    ----------
    path : str
        Namelist file, normally MALI's generated
        ``src/default_inputs/namelist.landice``.

    Returns
    -------
    collections.OrderedDict
    """
    records = OrderedDict()
    current = None
    with open(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('&'):
                current = stripped[1:]
                records[current] = OrderedDict()
            elif stripped == '/':
                current = None
            elif '=' in stripped and current is not None:
                key, value = stripped.split('=', 1)
                records[current][key.strip()] = value.strip()
    return records


def write_namelist(records, path):
    """Write ``{record: {option: value}}`` back out in MPAS namelist form."""
    with open(path, 'w') as handle:
        for record, options in records.items():
            handle.write(f'&{record}\n')
            for key, value in options.items():
                handle.write(f'    {key} = {value}\n')
            handle.write('/\n')


def apply_overrides(records, overrides):
    """
    Set namelist options, finding the record each one belongs to.

    Raises
    ------
    KeyError
        If an option is not present in the namelist.  That is deliberate: a
        silently ignored option would mean the run quietly did the wrong thing.
    """
    for key, value in overrides.items():
        for options in records.values():
            if key in options:
                options[key] = value
                break
        else:
            raise KeyError(
                f"'{key}' is not in the namelist; check the spelling against "
                f'the Registry'
            )
    return records


STREAMS_TEMPLATE = """<streams>

<immutable_stream name="basicmesh"
                  filename_template="not-to-be-used.nc"
                  type="none"/>

<immutable_stream name="input"
                  filename_template="{mesh_file}"
                  input_interval="initial_only"
                  type="input"/>

<immutable_stream name="restart"
                  filename_template="rst.$Y-$M-$D.nc"
                  input_interval="initial_only"
                  output_interval="none"
                  type="input;output"/>

<stream name="ismip7_masks"
        type="input"
        filename_template="{masks_file}"
        input_interval="initial_only"
        runtime_format="single_file">
    <var name="ismip6shelfMelt_basin"/>
    <var name="ismip6shelfMelt_deltaT"/>
</stream>

<stream name="ismip7_ocean_forcing"
        type="input"
        filename_template="{forcing_file}"
        input_interval="initial_only"
        runtime_format="single_file">
    <var name="ismip6shelfMelt_3dThermalForcing"/>
    <var name="ismip6shelfMelt_zOcean"/>
    <var name="ismip7shelfMelt_3dSalinity"/>
</stream>

<stream name="melt_diagnostic"
        type="output"
        filename_template="output_melt.nc"
        output_interval="{output_interval}"
        clobber_mode="overwrite"
        precision="double"
        runtime_format="single_file">
    <stream name="basicmesh"/>
    <var name="xtime"/>
    <var name="floatingBasalMassBal"/>
    <var name="ismip6shelfMelt_TFdraft"/>
    <var name="ismip7shelfMelt_Sdraft"/>
    <var name="ismip6shelfMelt_basin"/>
    <var name="ismip6shelfMelt_deltaT"/>
    <var name="cellMask"/>
    <var name="thickness"/>
    <var name="lowerSurface"/>
    <var name="bedTopography"/>
    <var name="areaCell"/>
</stream>

</streams>
"""


def write_streams(path, mesh_file, masks_file, forcing_file, output_interval):
    """Write the streams file for a melt-diagnostic run."""
    with open(path, 'w') as handle:
        handle.write(
            STREAMS_TEMPLATE.format(
                mesh_file=mesh_file,
                masks_file=masks_file,
                forcing_file=forcing_file,
                output_interval=output_interval,
            )
        )


def setup_run(
    run_dir,
    default_namelist,
    mesh_file,
    masks_file,
    forcing_file,
    graph_file=None,
    ntasks=None,
    melt_k=8.5e-5,
    sin_slope=0.0051117,
    coriolis=1.4e-4,
    timestep='0000-00-01_00:00:00',
    extra_namelist=None,
):
    """
    Create one melt-diagnostic run directory.

    Parameters
    ----------
    run_dir : str
        Directory to create.
    default_namelist : str
        MALI's generated ``namelist.landice`` to use as the template.
    mesh_file, masks_file, forcing_file : str
        Inputs; symlinked into the run directory.
    graph_file : str, optional
        Graph partition file, for multi-task runs.  It is linked as
        ``graph.info.part.<ntasks>``, which is the name MALI looks for given
        the default ``config_block_decomp_file_prefix``.
    ntasks : int, optional
        Number of MPI tasks the partition file is for.  Required with
        ``graph_file``.
    melt_k, sin_slope, coriolis : float, optional
        ISMIP7 melt parameters.
    timestep : str, optional
        MPAS time string for the single step.
    extra_namelist : dict, optional
        Further namelist overrides.

    Returns
    -------
    str
        The run directory.
    """
    os.makedirs(run_dir, exist_ok=True)

    links = {
        'mesh.nc': mesh_file,
        'masks.nc': masks_file,
        'forcing.nc': forcing_file,
    }
    if graph_file is not None:
        if ntasks is None:
            raise ValueError('ntasks is required when graph_file is given')
        # MALI builds the name from config_block_decomp_file_prefix, which
        # defaults to 'graph.info.part.', plus the task count
        links[f'graph.info.part.{ntasks}'] = graph_file
    for name, target in links.items():
        link = os.path.join(run_dir, name)
        if os.path.lexists(link):
            os.remove(link)
        os.symlink(os.path.abspath(target), link)

    # MPAS takes dimension sizes from the *input* stream, and the mesh file
    # has no nISMIP6OceanLayers.  Without this the 3-D forcing fields are
    # allocated against a zero-length dimension and the run dies trying to
    # allocate hundreds of GB.  Read it from the forcing rather than hard-
    # coding 30, so it stays right if ISMIP7 changes the vertical grid.
    with xr.open_dataset(forcing_file) as ds_forcing:
        n_layers = ds_forcing.sizes['nISMIP6OceanLayers']

    overrides = dict(DIAGNOSTIC_NAMELIST)
    overrides.update(
        {
            'config_nISMIP6OceanLayers': str(n_layers),
            'config_ismip7_melt_K': repr(float(melt_k)),
            'config_ismip7_melt_sin_slope': repr(float(sin_slope)),
            'config_ismip7_melt_coriolis': repr(float(coriolis)),
            'config_dt': f"'{timestep}'",
            'config_run_duration': f"'{timestep}'",
            'config_stop_time': "'none'",
        }
    )
    if extra_namelist is not None:
        overrides.update(extra_namelist)

    records = apply_overrides(read_namelist(default_namelist), overrides)
    write_namelist(records, os.path.join(run_dir, 'namelist.landice'))

    write_streams(
        os.path.join(run_dir, 'streams.landice'),
        'mesh.nc',
        'masks.nc',
        'forcing.nc',
        timestep,
    )
    return run_dir


def setup_ensemble(
    base_dir,
    states,
    default_namelist,
    mesh_file,
    masks_file,
    forcing_dir,
    **kwargs,
):
    """
    Create one run directory per ocean state, all at the same reference ``K``.

    Because melt is linear in ``K``, a single reference value per state is
    enough; the parameter ensemble is formed by scaling afterwards.

    Parameters
    ----------
    base_dir : str
        Directory to hold the runs.
    states : list of mali_melt_calib.datasets.OceanState
        The ocean states to run.
    default_namelist, mesh_file, masks_file : str
        As for :func:`setup_run`.
    forcing_dir : str
        Directory holding ``ocean_forcing_<state>.nc``.
    **kwargs
        Passed to :func:`setup_run`.

    Returns
    -------
    dict
        Mapping from state name to run directory.
    """
    runs = {}
    for state in states:
        forcing = os.path.join(forcing_dir, f'ocean_forcing_{state.name}.nc')
        run_dir = os.path.join(base_dir, state.name)
        runs[state.name] = setup_run(
            run_dir,
            default_namelist,
            mesh_file,
            masks_file,
            forcing,
            **kwargs,
        )
    return runs


def linearity_check_runs(
    base_dir,
    state,
    default_namelist,
    mesh_file,
    masks_file,
    forcing_dir,
    k_values=(4.0e-5, 8.5e-5, 1.3e-4),
    **kwargs,
):
    """
    Build runs at several ``K`` for one ocean state, to test the linearity.

    Findings F1 argues melt is exactly proportional to ``K``; this is how that
    is checked rather than assumed, and it is what licenses one run per ocean
    state instead of one per parameter value.
    """
    forcing = os.path.join(forcing_dir, f'ocean_forcing_{state.name}.nc')
    runs = {}
    for melt_k in k_values:
        run_dir = os.path.join(base_dir, f'{state.name}_K{melt_k:.3e}')
        runs[melt_k] = setup_run(
            run_dir,
            default_namelist,
            mesh_file,
            masks_file,
            forcing,
            melt_k=melt_k,
            **kwargs,
        )
    return runs


def write_job_script(
    path,
    run_dirs,
    model_exe,
    load_script,
    ntasks=128,
    account=None,
    partition='compute',
    walltime='0:30:00',
):
    """
    Write a Slurm script that runs each directory in turn.

    The runs are small, so serialising them in one job is simpler than a job
    array and avoids queue churn.
    """
    lines = [
        '#!/bin/bash',
        '#SBATCH --job-name=mali-melt-calib',
        '#SBATCH --nodes=1',
        f'#SBATCH --time={walltime}',
        f'#SBATCH --partition={partition}',
        '#SBATCH --output=mali-melt-calib.o%j',
    ]
    if account is not None:
        lines.append(f'#SBATCH --account={account}')
    lines += [
        '',
        'set -e',
        '',
        '# The environment must be *activated*, not just referenced: MALI and',
        '# its tooling shell out to netcdf and PIO utilities that live in the',
        "# environment's bin.",
        f'source {load_script}',
        '',
    ]
    for run_dir in run_dirs:
        lines += [
            f'cd {os.path.abspath(run_dir)}',
            f'srun -n {ntasks} {os.path.abspath(model_exe)}',
            '',
        ]
    with open(path, 'w') as handle:
        handle.write('\n'.join(lines))
    os.chmod(path, 0o755)
    return path
