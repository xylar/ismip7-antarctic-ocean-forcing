"""
Replicate the published 8 km quadratic calibration, through our own code path.

The protocol reports, for the quadratic local parameterisation on the ISMIP
8 km grid (§4.3.1 and Fig. 5):

===========  ==========
percentile   K
===========  ==========
5th          4.75e-5
50th         8.5e-5
95th         13.75e-5
===========  ==========

Reproducing those numbers with :mod:`mali_melt_calib.terms` (area-weighted,
mesh-agnostic) driving the upstream objective function validates the whole of
calibration stage 2 before any MALI runs exist, and gives a regression test
that the unstructured generalisation must not break.

It also settles feedback item A3 empirically: the published percentiles came
out of the toolbox as it stands, so whichever sampling reproduces them is the
one that was used.

Run with ``python -m mali_melt_calib replicate``.
"""

from __future__ import annotations

import importlib.util
import os
import sys

import numpy as np
import pandas as pd
import xarray as xr

from mali_melt_calib.quadratic import local_quadratic_melt, mean_slope
from mali_melt_calib.terms import (
    average_by_group,
    integrate_by_group,
    uniform_area,
)

#: default root of the ISMIP7 AIS datasets on LCRC
DEFAULT_DATA_ROOT = '/lcrc/group/e3sm/ac.xylar/ismip7/forcing-data/data/AIS'

#: ISMIP grid resolution used for the published calibration, m
RESO = 8000.0

#: the 120 K values sampled by the protocol's worked example
K_VALUES = np.arange(0.25e-5, 3.025e-4, 0.25e-5)

CELL_DIMS = ('y', 'x')

#: ocean-model states, as file prefixes and the labels the toolbox expects
OCEAN_MODELS = [
    ('Mathiot_NEMO_{state}_v3_', 'mathiot'),
    ('Timmermann_FESOM_{state}_v3_', 'timmermann'),
    ('Naughten_FESOM_ACCESS_{state}_v2_', 'naughten_ais_1'),
    ('Naughten_FESOM_MMM_{state}_v2_', 'naughten_ais_2'),
    ('Jourdain-Naughten_NEMO-MITgcm_{state}_', 'jourdain_naughten'),
    ('Naughten_MITamu-MITwed_{state}_', 'naughten_naughten'),
]

#: models whose domain covers only the Weddell Sea (ISMIP7 basin 14)
WEDDELL_ONLY = ('timmermann', 'haid')

#: models covering the Amundsen (basin 9) and Weddell (basin 14) only
AMUNDSEN_WEDDELL_ONLY = ('jourdain_naughten', 'naughten_naughten')

#: years of Amundsen ocean observations
OBS_YEARS = [
    1994,
    2000,
    2006,
    2007,
    2009,
    2010,
    2011,
    2012,
    2014,
    2016,
    2018,
    2019,
    2020,
]

PIG_ID = 110
DOTSON_ID = 97


def load_upstream_toolbox():
    """Import ``parameter_selection_toolbox`` from the parent directory."""
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(
        here, os.pardir, os.pardir, 'parameter_selection_toolbox.py'
    )
    path = os.path.normpath(path)
    if not os.path.exists(path):
        raise FileNotFoundError(f'upstream toolbox not found at {path}')
    spec = importlib.util.spec_from_file_location(
        'parameter_selection_toolbox', path
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _paths(data_root):
    param = os.path.join(data_root, 'parameterisations', 'ocean')
    obs = os.path.join(data_root, 'obs', 'ocean')
    return param, obs


def load_static(data_root):
    """
    Load the topography, masks and targets that do not vary with ocean state.

    Returns
    -------
    dict
        Keys: ``draft``, ``floating``, ``basins``, ``bfrn``, ``mask``,
        ``area``, ``slope``, ``region_label`` and the four sets of targets.
    """
    param, obs = _paths(data_root)

    topo = xr.load_dataset(
        os.path.join(
            obs,
            'topography',
            'bedmap3',
            'v3',
            'bedmap3_AIS_obs_ocean_topography_v3.nc',
        )
    )
    basins = xr.load_dataset(
        os.path.join(param, 'imbie2', 'basin_numbers_ismip8km_v2.nc')
    ).basinNumber.rename('basins')
    bfrn = xr.load_dataset(os.path.join(param, 'bfrns', 'BFRN_ismip8km_v2.nc'))
    mask = xr.load_dataset(
        os.path.join(param, 'floatingmasks', 'floatingmask_ismip8km.nc')
    ).mask

    floating = topo['floating_frac'] > 0.5
    slope = mean_slope(topo['draft'], floating, RESO, RESO)

    shelves = xr.load_dataset(
        os.path.join(param, 'shelfmask', 'shelf_mask_ismip8km.nc')
    ).shelf_mask.isel(time=0)
    # restrict PIG to its main trunk, as the worked example does
    x = shelves['x'] if 'x' in shelves.coords else basins['x']
    pig = (shelves == PIG_ID) & (x > -1.625e6)
    dotson = shelves == DOTSON_ID
    region_label = xr.where(pig, 'pig', xr.where(dotson, 'dotson', ''))

    melt_imbie = pd.read_csv(
        os.path.join(
            param, 'meltobs', 'Melt_Paolo_Davison_Adusumilli_imbie2.csv'
        ),
        index_col=0,
    )
    buttressing_target = xr.load_dataset(
        os.path.join(param, 'meltobs', 'melt_target_term2_v3.nc')
    )
    cold_target = xr.load_dataset(
        os.path.join(
            param, 'ocean_modelling_data', 'melt_cold_target_term3_v2.nc'
        )
    )
    warm_target = xr.load_dataset(
        os.path.join(
            param, 'ocean_modelling_data', 'melt_warm_target_term3_v2.nc'
        )
    )
    t4_obs = xr.load_dataset(
        os.path.join(
            param,
            'ocean_observations_data',
            'melt_observations_target_term4.nc',
        )
    )

    area = uniform_area(topo['draft'], RESO, CELL_DIMS)

    return dict(
        topo=topo,
        draft=topo['draft'],
        floating=floating,
        basins=basins,
        bfrn=bfrn,
        mask=mask,
        area=area,
        slope=slope,
        region_label=region_label,
        melt_imbie=melt_imbie,
        buttressing_target=buttressing_target,
        cold_target=cold_target,
        warm_target=warm_target,
        t4_obs=t4_obs,
    )


def unit_melt(tf_file, so_file, static, tf_var='tf', so_var='so'):
    """
    Melt at ``K = 1`` for one ocean state, in kg m-2 yr-1.

    Melt is exactly linear in ``K``, so every member of the parameter ensemble
    is this field times its ``K``.  Working at ``K = 1`` keeps the whole
    replication in 2-D fields instead of 120 copies of each.
    """
    draft = static['draft']
    floating = static['floating']

    tf = xr.load_dataset(tf_file)[tf_var]
    so = xr.load_dataset(so_file)[so_var]

    tf_draft = tf.sel(z=draft, method='nearest').where(floating)
    so_draft = so.sel(z=draft, method='nearest').where(floating)

    return local_quadratic_melt(1.0, tf_draft, so_draft, static['slope'])


def _scaled(aggregate_unit, k_values):
    """
    Turn a K = 1 aggregate into the (p1, p2, ...) ensemble aggregate.

    Exact, because the melt -- and therefore every linear aggregate of it -- is
    proportional to K.
    """
    k = xr.DataArray(k_values, dims=['p1'], coords={'p1': k_values})
    out = (aggregate_unit * k).expand_dims({'p2': np.ones(1)})
    return out.transpose('p1', 'p2', ...)


def build_terms(data_root=DEFAULT_DATA_ROOT, k_values=K_VALUES, verbose=True):
    """
    Build the four objective-function terms on the ISMIP 8 km grid.

    Returns
    -------
    dict
        ``t{n}_model``, ``t{n}_obs_mean``, ``t{n}_obs_sigma`` for n in 1..4,
        plus ``weights`` and ``static``.
    """
    param, _ = _paths(data_root)
    static = load_static(data_root)
    area, mask = static['area'], static['mask']
    basins, bfrn = static['basins'], static['bfrn']

    if verbose:
        print(
            f'mean draft slope: {static["slope"]:.6f} rad '
            f'(sin = {np.sin(static["slope"]):.6f})'
        )

    # ---- J1 and J2: present-day climatology -----------------------------
    clim_dir = os.path.join(
        data_root, 'obs', 'ocean', 'climatology', 'zhou_annual_06_nov'
    )
    pd_unit = unit_melt(
        os.path.join(
            clim_dir,
            'tf',
            'v3',
            'tf_AIS_obs_ocean_climatology_zhou_annual_06_nov_v3_1972-2024.nc',
        ),
        os.path.join(
            clim_dir,
            'so',
            'v4',
            'so_AIS_obs_ocean_climatology_zhou_annual_06_nov_v4_1972-2024.nc',
        ),
        static,
    )
    if verbose:
        print('built present-day ensemble')

    t1_unit = integrate_by_group(
        pd_unit, area, mask, basins, CELL_DIMS, group_dim='basins'
    )
    t1_model = _scaled(t1_unit, k_values)

    t2_unit = integrate_by_group(
        pd_unit,
        area,
        mask,
        bfrn['BFRN_bins'],
        CELL_DIMS,
        group_dim='BFRN_bins',
    )
    t2_model = _scaled(t2_unit, k_values)

    # ---- J3: ocean-model cold/warm pairs ---------------------------------
    model_dir = os.path.join(param, 'ocean_modelling_data')
    means = {}
    for state in ('cold', 'warm'):
        per_model = []
        for prefix, label in OCEAN_MODELS:
            stem = prefix.format(state=state)
            melt = unit_melt(
                os.path.join(model_dir, f'{stem}TF.nc'),
                os.path.join(model_dir, f'{stem}S.nc'),
                static,
            )
            # regional models only constrain the basins they cover
            if label in WEDDELL_ONLY:
                melt = melt.where(basins == 14)
            elif label in AMUNDSEN_WEDDELL_ONLY:
                melt = melt.where((basins == 9) | (basins == 14))
            agg = average_by_group(
                melt, area, mask, basins, CELL_DIMS, group_dim='basins'
            )
            per_model.append(agg.expand_dims({'model': [label]}))
            if verbose:
                print(f'  {state:4s} {label}')
        means[state] = xr.concat(per_model, dim='model')

    t3_model = _scaled(means['warm'] - means['cold'], k_values)

    # ---- J4: Amundsen observations ---------------------------------------
    obs_dir = os.path.join(param, 'ocean_observations_data')
    per_year = []
    for year in OBS_YEARS:
        melt = unit_melt(
            os.path.join(obs_dir, f'Obs_{year}_TF.nc'),
            os.path.join(obs_dir, f'Obs_{year}_S.nc'),
            static,
        )
        agg = integrate_by_group(
            melt,
            area,
            mask,
            static['region_label'],
            CELL_DIMS,
            group_dim='region',
        )
        per_year.append(agg.expand_dims({'year': [year]}))
    t4_unit = xr.concat(per_year, dim='year')
    # drop the '' label used for cells outside PIG and Dotson
    t4_unit = t4_unit.sel(region=[r for r in t4_unit.region.values if r])
    t4_model = _scaled(t4_unit, k_values)

    t4_obs = static['t4_obs']
    t4_model = t4_model.where(t4_obs.melt_rate.notnull())
    t4_model = t4_model.reindex_like(t4_obs.melt_rate)
    if verbose:
        print('built observational ensemble')

    # ---- targets and weights ---------------------------------------------
    melt_imbie = static['melt_imbie']
    t1_obs_mean = xr.DataArray(
        melt_imbie['BMR (Gt/yr)'].values.astype(float),
        dims=['basin'],
        coords={'basin': np.arange(len(melt_imbie))},
    )
    t1_obs_sigma = xr.DataArray(
        melt_imbie['BMR uncert (Gt/yr)'].values.astype(float),
        dims=['basin'],
        coords={'basin': np.arange(len(melt_imbie))},
    )

    t3_obs_mean = (
        static['warm_target'].melt_rate - static['cold_target'].melt_rate
    )
    t3_obs_sigma = np.sqrt(
        static['warm_target'].melt_rate_uncert ** 2
        + static['cold_target'].melt_rate_uncert ** 2
    )

    t1_weights = xr.DataArray(
        np.ones(t1_model.sizes['basins']),
        dims=['basins'],
        coords={'basins': t1_model.basins.values},
    )
    t2_weights = xr.DataArray(
        (bfrn['BFRN_medians'] / bfrn['BFRN_median']).values,
        dims=['BFRN_bins'],
        coords={'BFRN_bins': static['buttressing_target'].BFRN_bins.values},
    )
    t3_weights = xr.DataArray(
        np.ones((t3_model.sizes['model'], t3_model.sizes['basins'])),
        dims=['model', 'basins'],
        coords={
            'model': t3_model.model.values,
            'basins': t3_model.basins.values,
        },
    )
    keep = (
        'mathiot',
        'naughten_ais_1',
        'jourdain_naughten',
        'naughten_naughten',
    )
    t3_weights = t3_weights.where(t3_weights.model.isin(keep), other=0)

    t4_weights = xr.DataArray(
        np.ones(
            (t4_obs.melt_rate.sizes['region'], t4_obs.melt_rate.sizes['year'])
        ),
        dims=['region', 'year'],
        coords={'region': t4_obs.region.values, 'year': t4_obs.year.values},
    )
    t4_weights = t4_weights.where(t4_weights.region == 'pig', other=0)
    t4_weights = t4_weights.where(
        (t4_weights.year == 2009) | (t4_weights.year == 2012), other=0
    )

    return dict(
        t1_model=t1_model,
        t1_obs_mean=t1_obs_mean,
        t1_obs_sigma=t1_obs_sigma,
        t1_weights=t1_weights,
        t2_model=t2_model,
        t2_obs_mean=static['buttressing_target']['melt_mean'],
        t2_obs_sigma=static['buttressing_target']['melt_mean_err'],
        t2_weights=t2_weights,
        t3_model=t3_model,
        t3_obs_mean=t3_obs_mean,
        t3_obs_sigma=t3_obs_sigma,
        t3_weights=t3_weights,
        t4_model=t4_model,
        t4_obs_mean=t4_obs.melt_rate,
        t4_obs_sigma=t4_obs.melt_rate_uncert,
        t4_weights=t4_weights,
        static=static,
    )


def run(
    data_root=DEFAULT_DATA_ROOT, sample_size=100000, seed=None, verbose=True
):
    """
    Build the terms, run the optimisation, and report the percentiles.

    Returns
    -------
    dict
        ``min_p1``, plus the 5th/50th/95th percentiles and the mode.
    """
    toolbox = load_upstream_toolbox()
    terms = build_terms(data_root=data_root, verbose=verbose)

    if seed is not None:
        np.random.seed(seed)

    if verbose:
        print(f'running objective function with {sample_size} samples ...')
    min_p1, min_p2 = toolbox.calculate_objective_function(
        sample_size,
        RESO,
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
    edges = np.append(K_VALUES[0] * 0.5, K_VALUES + 1e-7)
    counts, _ = np.histogram(min_p1, bins=edges)

    result = {
        'min_p1': min_p1,
        'p5': float(np.percentile(min_p1, 5)),
        'median': float(np.median(min_p1)),
        'p95': float(np.percentile(min_p1, 95)),
        'mode': float(K_VALUES[int(np.argmax(counts))]),
    }

    if verbose:
        published = {'p5': 4.75e-5, 'median': 8.5e-5, 'p95': 13.75e-5}
        print()
        print(f'{"":10s} {"ours":>12s} {"published":>12s}')
        for key in ('p5', 'median', 'p95'):
            print(f'{key:10s} {result[key]:12.3e} {published[key]:12.3e}')
        print(f'{"mode":10s} {result["mode"]:12.3e}')
    return result


def main(argv=None):
    """CLI entry point for ``python -m mali_melt_calib replicate``."""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', default=DEFAULT_DATA_ROOT)
    parser.add_argument('--sample-size', type=int, default=100000)
    parser.add_argument('--seed', type=int, default=None)
    args = parser.parse_args(argv)

    run(
        data_root=args.data_root,
        sample_size=args.sample_size,
        seed=args.seed,
    )
    return 0


if __name__ == '__main__':
    sys.exit(main())
