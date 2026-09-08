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
    t4_target = targets['t4_mean'].sel(
        region=units['t4'].region.values, year=units['t4'].year.values
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
        't4': _mae(
            units['t4'] * scale, t4_target, 1.0, ['region', 'year']
        ).values,
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
