"""
Remap ISMIP7 thermal forcing and salinity onto the MALI mesh.

One file per ocean state, in the form MALI's input streams expect.  Thermal
forcing drives the existing ISMIP6 quadratic; salinity is additionally needed
by the Burgard et al. (2022) form on branch ``mali/add-burgard-melt-param``,
since protocol Eq. (1) contains ``beta_S S_loc``.  Both are remapped together
so that a late decision on the melt module costs nothing here.

The ISMIP7 vertical coordinate is already the one MALI expects -- 30 layers at
60 m spacing, centres from -30 m to -1770 m -- so there is no vertical
interpolation, only a horizontal remap.

Field naming and dimension order follow compass's
``landice/tests/ismip7_forcing/ocean_thermal/process_thermal_forcing.py`` so
that the files are interchangeable with the ones the production runs use.
"""

from __future__ import annotations

import os

import numpy as np
import xarray as xr
from pyremap import Remapper

from mali_melt_calib.datasets import (
    DEFAULT_DATA_ROOT,
    ISMIP_RESOLUTION_KM,
    climatology_files,
    ocean_states,
)

# EPSG:3031, the ISMIP Antarctic polar-stereographic projection
ISMIP_PROJ_STR = (
    '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 '
    '+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
)

#: half the ISMIP7 ocean layer thickness, m
HALF_LAYER = 30.0

TF_VAR = 'ismip6shelfMelt_3dThermalForcing'
SO_VAR = 'ismip7shelfMelt_3dSalinity'
Z_VAR = 'ismip6shelfMelt_zOcean'
ZBNDS_VAR = 'ismip6shelfMelt_zBndsOcean'


class ForcingRemapper:
    """
    Remap ISMIP-grid ocean fields onto a MALI mesh.

    The mapping weights depend only on the two grids, so the remapper is built
    once and reused for every ocean state.

    Parameters
    ----------
    mali_mesh : str
        MALI mesh file.
    map_filename : str, optional
        Where to write the mapping file.
    resolution : int, optional
        ISMIP grid resolution in km.
    method : str, optional
        Remapping method.  Bilinear matches compass's ocean-thermal step;
        the fields are continuous, so conservative would also be defensible.
    ntasks : int, optional
        MPI tasks for ``ESMF_RegridWeightGen``.
    mali_mesh_name : str, optional
        Name for the mesh in the mapping file; derived from the filename by
        default.
    """

    def __init__(
        self,
        mali_mesh,
        map_filename=None,
        resolution=ISMIP_RESOLUTION_KM,
        method='bilinear',
        ntasks=1,
        mali_mesh_name=None,
        logger=None,
    ):
        self.mali_mesh = mali_mesh
        self.resolution = resolution
        self.method = method
        self.ntasks = ntasks
        self.logger = logger

        if mali_mesh_name is None:
            mali_mesh_name = os.path.splitext(os.path.basename(mali_mesh))[0]
        self.mali_mesh_name = mali_mesh_name

        src_name = f'ismip{resolution}km'
        if map_filename is None:
            map_filename = f'map_{src_name}_to_{mali_mesh_name}_{method}.nc'
        self.map_filename = map_filename
        self._remapper = None
        self._src_name = src_name

    def _build(self, src_grid_file):
        """
        Construct the pyremap remapper and its weights.

        ``build_map`` is what creates the source and destination descriptors,
        so it has to be called even if the mapping file already exists;
        skipping it leaves ``remap_numpy`` with descriptors of None.
        """
        remapper = Remapper(
            ntasks=self.ntasks,
            map_filename=self.map_filename,
            method=self.method,
        )
        remapper.src_from_proj(
            src_grid_file, self._src_name, proj_str=ISMIP_PROJ_STR
        )
        remapper.dst_from_mpas(self.mali_mesh, self.mali_mesh_name)
        remapper.build_map(logger=self.logger)
        self._remapper = remapper
        return remapper

    def remap_state(self, state, out_file, data_root=DEFAULT_DATA_ROOT):
        """
        Remap one ocean state and write it in MALI's input form.

        Parameters
        ----------
        state : mali_melt_calib.datasets.OceanState
            The state to remap.
        out_file : str
            Output NetCDF file.
        data_root : str, optional
            Root of the ISMIP7 datasets, used to find the climatology when a
            partially covered state has to be filled.

        Returns
        -------
        xarray.Dataset
            What was written.
        """
        if self._remapper is None:
            self._build(state.tf_file)

        ds_in = xr.Dataset()
        ds_in['tf'] = xr.open_dataset(state.tf_file)['tf']
        ds_in['so'] = xr.open_dataset(state.so_file)['so']

        filled = _fill_from_climatology(ds_in, state, data_root, self.logger)

        z = ds_in['tf']['z']
        ds_out = self._remapper.remap_numpy(ds_in)
        return _to_mali_form(ds_out, z, state, out_file, filled)


def _fill_from_climatology(ds_in, state, data_root, logger=None):
    """
    Fill gaps in a partially covering ocean state with the climatology.

    The ISMIP7 near-ice-shelf observational datasets apply a single profile to
    the Amundsen basin only (protocol §A10) and are undefined elsewhere --
    about 6% of the ISMIP grid is valid.  The regional *model* datasets, by
    contrast, are distributed already filled with the ISMIP7 climatology
    outside their domains (protocol §A9), so they arrive complete.

    An ice-sheet model needs valid forcing everywhere it has ice, so the
    observational states are filled the same way the model states were.  This
    does not affect the calibration: J4 aggregates only over Pine Island and
    Dotson, both inside the covered basin.

    Filling happens *before* remapping, so that interpolation never blends a
    real value with a missing one at the edge of the covered region.

    Returns
    -------
    float
        Fraction of the source grid that was filled, for the record.
    """
    valid = ds_in['tf'].notnull()
    fraction = float((~valid).mean())
    if fraction == 0.0:
        return 0.0

    tf_file, so_file = climatology_files(data_root)
    clim_tf = xr.open_dataset(tf_file)['tf']
    clim_so = xr.open_dataset(so_file)['so']

    if logger is not None:
        logger.info(
            f'    filling {100.0 * fraction:.1f}% of the source grid from the '
            f'climatology'
        )

    ds_in['tf'] = ds_in['tf'].where(valid, clim_tf)
    ds_in['so'] = ds_in['so'].where(ds_in['so'].notnull(), clim_so)
    return fraction


def _to_mali_form(ds_remapped, z, state, out_file, filled=0.0):
    """
    Put remapped fields into the names, dimensions and order MALI expects.

    Registry declares the 3-D fields as ``nISMIP6OceanLayers nCells Time`` in
    Fortran order, which is ``Time, nCells, nISMIP6OceanLayers`` in the file.
    """
    ds = ds_remapped.rename({'z': 'nISMIP6OceanLayers'})
    ds = ds.rename({'tf': TF_VAR, 'so': SO_VAR})

    for var in (TF_VAR, SO_VAR):
        ds[var] = ds[var].expand_dims('Time', axis=0)
        ds[var] = ds[var].transpose('Time', 'nCells', 'nISMIP6OceanLayers')
        ds[var] = ds[var].astype(float)
        ds[var].encoding.clear()

    z_values = np.asarray(z.values, dtype=float)
    ds[Z_VAR] = ('nISMIP6OceanLayers', z_values)
    # the calibration files carry no z_bnds; layers are 60 m thick and
    # centred on z, matching the climatology's own bounds
    ds[ZBNDS_VAR] = (
        ('TWO', 'nISMIP6OceanLayers'),
        np.stack([z_values + HALF_LAYER, z_values - HALF_LAYER]),
    )

    ds[TF_VAR].attrs = {
        'long_name': 'thermal forcing for the ISMIP6 ice-shelf melting method',
        'units': 'degC',
    }
    ds[SO_VAR].attrs = {
        'long_name': 'practical salinity for the ISMIP7 quadratic melt '
        'parameterisation',
        'units': 'PSU',
    }
    ds[Z_VAR].attrs = {
        'long_name': 'depth coordinate for ocean thermal forcing',
        'units': 'm',
    }
    ds[ZBNDS_VAR].attrs = {
        'long_name': 'bounds of the ocean layers',
        'units': 'm',
    }

    if filled > 0.0:
        ds.attrs['climatology_filled_fraction'] = filled
        ds.attrs['climatology_fill_note'] = (
            'Cells outside this dataset coverage were filled with the ISMIP7 '
            'present-day climatology before remapping, matching how the '
            'regional ocean-model datasets are distributed. Only the covered '
            'basins are used in the objective function.'
        )
    ds.attrs['ismip7_ocean_state'] = state.name
    ds.attrs['ismip7_term'] = state.term
    ds.attrs['source_tf'] = state.tf_file
    ds.attrs['source_so'] = state.so_file
    if state.basins is not None:
        ds.attrs['constrains_ismip7_basins'] = ', '.join(
            str(b) for b in state.basins
        )

    # MPAS cannot read a file that carries a variable named after one of its
    # dimensions, and pyremap leaves several auxiliary coordinates behind.
    # compass's ismip7_forcing drops the same set; without this MALI builds a
    # bad decomposition and dies trying to allocate hundreds of GB.
    aux = [
        'lon',
        'lat',
        'lon_vertices',
        'lat_vertices',
        'lon_bnds',
        'lat_bnds',
        'lat_cell',
        'lon_cell',
        'area',
        'z_bnds',
        'time_bnds',
        'x_bnds',
        'y_bnds',
        'crs',
    ]
    drop = [name for name in aux if name in ds]
    if drop:
        ds = ds.drop_vars(drop)
    # the renamed vertical coordinate shares the dimension's name
    if 'nISMIP6OceanLayers' in ds.coords:
        ds = ds.drop_vars('nISMIP6OceanLayers')
    if 'Time' in ds.coords:
        ds = ds.drop_vars('Time')

    ds.to_netcdf(out_file, unlimited_dims=['Time'])
    return ds


def remap_all(
    mali_mesh,
    out_dir,
    data_root=DEFAULT_DATA_ROOT,
    subset='recommended',
    method='bilinear',
    ntasks=1,
    logger=None,
    overwrite=False,
):
    """
    Remap every ocean state in ``subset`` onto the MALI mesh.

    Returns
    -------
    dict
        Mapping from state name to the file written.
    """
    os.makedirs(out_dir, exist_ok=True)
    states = ocean_states(data_root=data_root, subset=subset)

    remapper = ForcingRemapper(
        mali_mesh,
        method=method,
        ntasks=ntasks,
        logger=logger,
        map_filename=os.path.join(
            out_dir, f'map_ismip8km_to_mali_{method}.nc'
        ),
    )

    written = {}
    for state in states:
        out_file = os.path.join(out_dir, f'ocean_forcing_{state.name}.nc')
        written[state.name] = out_file
        if os.path.exists(out_file) and not overwrite:
            if logger is not None:
                logger.info(f'  {state.name}: exists, skipping')
            continue
        if logger is not None:
            logger.info(f'  {state.name} ({state.term})')
        remapper.remap_state(state, out_file, data_root=data_root)

    return written
