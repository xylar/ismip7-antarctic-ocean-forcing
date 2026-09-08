"""
Evaluate the ISMIP7 objective-function terms on the MALI mesh.

:mod:`mali_melt_calib.replicate` reproduces the published calibration on the
8 km structured grid.  This does the same arithmetic against melt fields that
**MALI itself** produced, one run per ocean state, so the calibration is of
MALI's parameterisation on MALI's mesh rather than of a Python reimplementation
on the protocol's grid.

Two things make that possible without remapping any melt field:

* every target the protocol distributes is **already aggregated** -- melt per
  IMBIE2 basin, per BFRN buttressing bin, per ice-shelf region -- so the
  targets are mesh-independent scalars and only the *grouping* of cells has to
  be transferred to the MALI mesh, which the mask tool does with
  nearest-neighbour remapping of categorical fields;
* melt is exactly linear in the melt parameter (findings F1, measured in F9),
  so one MALI run per ocean state gives the whole parameter sweep by scaling.

Aggregation uses **MALI's own floating cells**, following the note the mask
tool writes into the mask file: the calibration holds the ice-sheet model
accountable for its own shelf extent.  That choice is not free -- MALI's
floating area differs from the observed extent by enough to matter, see
:func:`mask_disagreement` -- so the intersection with the observed mask is
available as a diagnostic to separate geometry error from parameter error.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
import xarray as xr

from mali_melt_calib.datasets import (
    DEFAULT_DATA_ROOT,
    OCEAN_MODELS,
    RECOMMENDED_MODELS,
    RECOMMENDED_OBS_YEARS,
)
from mali_melt_calib.quadratic import MALI
from mali_melt_calib.terms import average_by_group, integrate_by_group

#: cell dimensions of an MPAS mesh
CELL_DIMS = ('nCells',)

#: MALI cellMask bits, from mpas_li_mask.F
MASK_FLOATING = 4
MASK_ICE = 32

#: ismip7ShelfRegion codes written by interpolate_ismip7_masks_to_mali.py
REGION_LABELS = {0: '', 1: 'pig', 2: 'dotson'}


def melt_from_run(run_dir, var='floatingBasalMassBal'):
    """
    Melt rate in kg m-2 yr-1, positive for melting, from a MALI run.

    MALI stores ``floatingBasalMassBal`` in kg m-2 s-1 with melting negative,
    the opposite convention to the calibration toolbox.
    """
    path = os.path.join(run_dir, 'output_melt.nc')
    with xr.open_dataset(
        path, decode_times=False, decode_timedelta=False
    ) as ds:
        bmb = ds[var]
        if 'Time' in bmb.dims:
            bmb = bmb.isel(Time=-1)
        return (-bmb * MALI.seconds_per_year).compute()


def load_mali_static(masks_file, geometry_run):
    """
    Cell areas, masks and groupings on the MALI mesh.

    Parameters
    ----------
    masks_file : str
        Output of ``interpolate_ismip7_masks_to_mali.py``.
    geometry_run : str
        Any run directory; ``areaCell`` and ``cellMask`` are read from its
        output.  The geometry is identical across the ensemble.

    Returns
    -------
    dict
        ``area``, ``floating``, ``obs_floating``, ``basins``, ``bfrn``,
        ``region``.
    """
    ms = xr.open_dataset(masks_file)
    path = os.path.join(geometry_run, 'output_melt.nc')
    with xr.open_dataset(
        path, decode_times=False, decode_timedelta=False
    ) as ds:
        area = ds['areaCell']
        if 'Time' in area.dims:
            area = area.isel(Time=-1)
        cell_mask = ds['cellMask']
        if 'Time' in cell_mask.dims:
            cell_mask = cell_mask.isel(Time=-1)
        area = area.compute()
        cell_mask = cell_mask.compute()

    floating = ((cell_mask & MASK_FLOATING) > 0) & ((cell_mask & MASK_ICE) > 0)

    # invalid group ids become NaN so that _prepare drops those cells rather
    # than forming a spurious group
    basins = ms['ismip7BasinNumber'].astype(float)
    basins = basins.where(basins >= 0)
    bfrn = ms['ismip7BFRNBin'].astype(float)
    bfrn = bfrn.where(bfrn >= 0)

    codes = ms['ismip7ShelfRegion'].values
    region = xr.DataArray(
        np.array([REGION_LABELS.get(int(c), '') for c in codes]),
        dims=ms['ismip7ShelfRegion'].dims,
    )

    return dict(
        area=area,
        floating=floating,
        obs_floating=ms['ismip7FloatingMask'] > 0,
        basins=basins,
        bfrn=bfrn,
        region=region,
    )


def load_targets(data_root=DEFAULT_DATA_ROOT):
    """
    The observational targets and the BFRN weights.

    These are the same files :mod:`mali_melt_calib.replicate` uses.  They are
    aggregated quantities, so they apply unchanged on the MALI mesh.
    """
    param = os.path.join(data_root, 'parameterisations', 'ocean')

    melt_imbie = pd.read_csv(
        os.path.join(
            param, 'meltobs', 'Melt_Paolo_Davison_Adusumilli_imbie2.csv'
        ),
        index_col=0,
    )
    buttressing = xr.load_dataset(
        os.path.join(param, 'meltobs', 'melt_target_term2_v3.nc')
    )
    bfrn = xr.load_dataset(os.path.join(param, 'bfrns', 'BFRN_ismip8km_v2.nc'))
    cold = xr.load_dataset(
        os.path.join(
            param, 'ocean_modelling_data', 'melt_cold_target_term3_v2.nc'
        )
    )
    warm = xr.load_dataset(
        os.path.join(
            param, 'ocean_modelling_data', 'melt_warm_target_term3_v2.nc'
        )
    )
    t4 = xr.load_dataset(
        os.path.join(
            param,
            'ocean_observations_data',
            'melt_observations_target_term4.nc',
        )
    )

    return dict(
        t1_mean=xr.DataArray(
            melt_imbie['BMR (Gt/yr)'].values.astype(float),
            dims=['basins'],
            coords={'basins': np.arange(len(melt_imbie), dtype=float)},
        ),
        t1_sigma=xr.DataArray(
            melt_imbie['BMR uncert (Gt/yr)'].values.astype(float),
            dims=['basins'],
            coords={'basins': np.arange(len(melt_imbie), dtype=float)},
        ),
        t2_mean=buttressing['melt_mean'],
        t2_sigma=buttressing['melt_mean_err'],
        t2_weights=xr.DataArray(
            (bfrn['BFRN_medians'] / bfrn['BFRN_median']).values,
            dims=['BFRN_bins'],
            coords={'BFRN_bins': buttressing.BFRN_bins.values},
        ),
        t3_mean=warm['melt_rate'] - cold['melt_rate'],
        t3_sigma=np.sqrt(
            warm['melt_rate_uncert'] ** 2 + cold['melt_rate_uncert'] ** 2
        ),
        t4_mean=t4['melt_rate'],
        t4_sigma=t4['melt_rate_uncert'],
    )


def unit_aggregates(ensemble_dir, param_value, static, mask=None):
    """
    Aggregate every ocean state's melt to the four terms, at unit parameter.

    Dividing by ``param_value`` is exact, because melt is proportional to the
    parameter (F9), so these aggregates can be scaled to any parameter value
    afterwards rather than needing a run per value.

    Parameters
    ----------
    ensemble_dir : str
        Directory holding one run directory per ocean state.
    param_value : float
        The K or gamma0 the runs used.
    static : dict
        From :func:`load_mali_static`.
    mask : xarray.DataArray, optional
        Cells to aggregate over; defaults to MALI's own floating cells.

    Returns
    -------
    dict
        ``t1``, ``t2`` (per basin, per BFRN bin, Gt/yr) and ``t3``, ``t4``.
    """
    if mask is None:
        mask = static['floating']
    area, basins = static['area'], static['basins']

    def unit(name):
        return melt_from_run(os.path.join(ensemble_dir, name)) / param_value

    clim = unit('climatology')
    t1 = integrate_by_group(
        clim, area, mask, basins, CELL_DIMS, group_dim='basins'
    )
    t2 = integrate_by_group(
        clim, area, mask, static['bfrn'], CELL_DIMS, group_dim='BFRN_bins'
    )

    # J3: warm minus cold basin-mean melt, per ocean model.  Regional models
    # only constrain the basins their domain covers.
    covered_by = {label: covered for _, label, covered in OCEAN_MODELS}
    per_model = []
    for label in RECOMMENDED_MODELS:
        means = {}
        for state in ('cold', 'warm'):
            melt = unit(f'{label}_{state}')
            covered = covered_by[label]
            if covered is not None:
                melt = melt.where(basins.isin([float(b) for b in covered]))
            means[state] = average_by_group(
                melt, area, mask, basins, CELL_DIMS, group_dim='basins'
            )
        per_model.append(
            (means['warm'] - means['cold']).expand_dims({'model': [label]})
        )
    t3 = xr.concat(per_model, dim='model')

    # J4: PIG and Dotson integrated melt, per observation year
    per_year = []
    for year in RECOMMENDED_OBS_YEARS:
        agg = integrate_by_group(
            unit(f'obs_{year}'),
            area,
            mask,
            static['region'],
            CELL_DIMS,
            group_dim='region',
        )
        per_year.append(agg.expand_dims({'year': [year]}))
    t4 = xr.concat(per_year, dim='year')
    t4 = t4.sel(region=[r for r in t4.region.values if r])

    return dict(t1=t1, t2=t2, t3=t3, t4=t4)


def _mae(predicted, observed, weights, dims, skipna=True):
    """The toolbox's weighted mean absolute error, for one parameter value."""
    return abs(weights * (predicted - observed)).mean(dims, skipna=skipna)


def term_misfits(units, targets, parameters):
    """
    Weighted mean absolute error of each term, as a function of the parameter.

    ``parameters`` holds **actual parameter values** -- K for the local form,
    gamma0 for the non-local one -- not multipliers of the value the runs used.
    The aggregates in ``units`` are per unit parameter, so they are multiplied
    by these directly.

    Uses the targets at their means rather than sampling them, which makes the
    comparison between two parameterisations deterministic.  The toolbox
    samples the targets and the term weights instead; that spreads the
    resulting parameter distribution but does not change which parameterisation
    fits a given term better.

    Returns
    -------
    dict
        ``t1``..``t4``, each an array over ``scales``.
    """
    scale = xr.DataArray(
        parameters, dims=['parameter'], coords={'parameter': parameters}
    )

    t3_target = targets['t3_mean'].sel(model=units['t3'].model.values)
    # The protocol weights J4 to PIG alone (t4_weights zeroes Dotson), so the
    # misfit must be restricted the same way; averaging Dotson in as well
    # reverses which form appears to fit better.
    t4_model = units['t4'].sel(region='pig')
    t4_target = targets['t4_mean'].sel(
        region='pig', year=units['t4'].year.values
    )

    return {
        't1': _mae(
            units['t1'] * scale, targets['t1_mean'], 1.0, 'basins'
        ).values,
        't2': _mae(
            units['t2'] * scale,
            targets['t2_mean'],
            targets['t2_weights'],
            'BFRN_bins',
        ).values,
        't3': _mae(
            units['t3'] * scale, t3_target, 1.0, ['model', 'basins']
        ).values,
        't4': _mae(t4_model * scale, t4_target, 1.0, ['year']).values,
    }


def mask_disagreement(static):
    """
    How far MALI's floating extent differs from the observed ISMIP7 mask.

    The aggregation uses MALI's own floating cells, so any difference in shelf
    extent enters the misfit as if it were a melt-parameter error.  This
    quantifies that, in area and as a fraction of the cells involved.
    """
    area = static['area'].values
    mali = static['floating'].values
    obs = static['obs_floating'].values

    def km2(sel):
        return float(area[sel].sum()) / 1.0e9

    return {
        'mali_area': km2(mali),
        'obs_area': km2(obs),
        'intersection': km2(mali & obs),
        'mali_only': km2(mali & ~obs),
        'obs_only': km2(obs & ~mali),
    }


def toolbox_terms(units, targets, param_values):
    """
    Assemble the arguments ``calculate_objective_function`` expects.

    The toolbox wants each modelled term indexed by ``(p1, p2, ...)``, where
    ``p1`` is the melt parameter being selected and ``p2`` a second parameter
    the quadratic parameterisations do not use, so it is a singleton -- the
    same shape :mod:`mali_melt_calib.replicate` builds for the 8 km
    replication.  Scaling by ``p1`` is exact because melt is proportional to
    the parameter (F9).
    """
    p1 = xr.DataArray(param_values, dims=['p1'], coords={'p1': param_values})

    def scaled(unit):
        out = (unit * p1).expand_dims({'p2': np.ones(1)})
        return out.transpose('p1', 'p2', ...)

    t1_model = scaled(units['t1'])
    t2_model = scaled(units['t2'])
    t3_model = scaled(units['t3'])

    t4_model = scaled(units['t4'])
    t4_model = t4_model.where(targets['t4_mean'].notnull())
    t4_model = t4_model.reindex_like(targets['t4_mean'])

    t1_weights = xr.DataArray(
        np.ones(t1_model.sizes['basins']),
        dims=['basins'],
        coords={'basins': t1_model.basins.values},
    )
    t3_weights = xr.DataArray(
        np.ones((t3_model.sizes['model'], t3_model.sizes['basins'])),
        dims=['model', 'basins'],
        coords={
            'model': t3_model.model.values,
            'basins': t3_model.basins.values,
        },
    )
    t4_weights = xr.DataArray(
        np.ones(
            (
                targets['t4_mean'].sizes['region'],
                targets['t4_mean'].sizes['year'],
            )
        ),
        dims=['region', 'year'],
        coords={
            'region': targets['t4_mean'].region.values,
            'year': targets['t4_mean'].year.values,
        },
    )
    # only PIG, and only the years the ensemble actually ran
    t4_weights = t4_weights.where(t4_weights.region == 'pig', other=0)
    t4_weights = t4_weights.where(
        t4_weights.year.isin(list(units['t4'].year.values)), other=0
    )

    return dict(
        t1_model=t1_model,
        t1_obs_mean=targets['t1_mean'],
        t1_obs_sigma=targets['t1_sigma'],
        t1_weights=t1_weights,
        t2_model=t2_model,
        t2_obs_mean=targets['t2_mean'],
        t2_obs_sigma=targets['t2_sigma'],
        t2_weights=targets['t2_weights'],
        t3_model=t3_model,
        t3_obs_mean=targets['t3_mean'].sel(model=t3_model.model.values),
        t3_obs_sigma=targets['t3_sigma'].sel(model=t3_model.model.values),
        t3_weights=t3_weights,
        t4_model=t4_model,
        t4_obs_mean=targets['t4_mean'],
        t4_obs_sigma=targets['t4_sigma'],
        t4_weights=t4_weights,
    )


def run_optimisation(
    terms, param_values, resolution=8000.0, sample_size=100000, seed=None
):
    """
    Sample the objective function and return the parameter distribution.

    Parameters
    ----------
    terms : dict
        From :func:`toolbox_terms`.
    param_values : numpy.ndarray
        The ``p1`` grid the terms were built on.
    resolution : float, optional
        Passed through to the toolbox; it does not use it for these terms,
        which arrive already aggregated.
    sample_size : int, optional
        Number of random draws of the term weights and the targets.
    seed : int, optional
        Seed for reproducibility.

    Returns
    -------
    dict
        ``p5``, ``median``, ``p95``, ``mode`` and the raw ``min_p1``.
    """
    from mali_melt_calib.replicate import load_upstream_toolbox

    toolbox = load_upstream_toolbox()
    if seed is not None:
        np.random.seed(seed)

    min_p1, _ = toolbox.calculate_objective_function(
        sample_size,
        resolution,
        terms['t1_model'],
        terms['t1_obs_mean'],
        terms['t1_obs_sigma'],
        terms['t1_weights'],
        terms['t2_model'],
        terms['t2_obs_mean'],
        terms['t2_obs_sigma'],
        terms['t2_weights'],
        terms['t3_model'],
        terms['t3_obs_mean'],
        terms['t3_obs_sigma'],
        terms['t3_weights'],
        terms['t4_model'],
        terms['t4_obs_mean'],
        terms['t4_obs_sigma'],
        terms['t4_weights'],
    )
    min_p1 = np.asarray(min_p1, dtype=float)

    values = np.asarray(param_values, dtype=float)
    step = np.diff(values).min()
    edges = np.append(values[0] - 0.5 * step, values + 1.0e-7 * step)
    counts, _ = np.histogram(min_p1, bins=edges)

    return {
        'min_p1': min_p1,
        'p5': float(np.percentile(min_p1, 5)),
        'median': float(np.median(min_p1)),
        'p95': float(np.percentile(min_p1, 95)),
        'mode': float(values[int(np.argmax(counts))]),
    }
