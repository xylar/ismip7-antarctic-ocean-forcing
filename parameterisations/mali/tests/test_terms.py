"""
Tests for the mesh-agnostic objective-function terms.

The central claim being tested is that passing a *uniform* cell area
reproduces the upstream structured-grid toolbox exactly, so that the same code
can serve the ISMIP 8 km grid and MALI's variable-resolution mesh.  The
upstream formulas are reimplemented here directly from
``parameterisations/parameter_selection_toolbox.py`` so the comparison is
against that code's semantics rather than against our own restatement.
"""

import numpy as np
import pandas as pd
import pytest
import xarray as xr
from mali_melt_calib.terms import (
    KG_PER_GT,
    average_by_group,
    calculate_term1,
    calculate_term3,
    integrate_by_group,
    uniform_area,
)

RESO = 8000.0
NY, NX = 12, 10
NBASIN = 4


@pytest.fixture
def grid():
    """A small structured ensemble with basins and a floating mask."""
    rng = np.random.default_rng(20260908)
    y = np.arange(NY, dtype=float) * RESO
    x = np.arange(NX, dtype=float) * RESO
    p1 = np.array([1.0e-5, 2.0e-5, 3.0e-5])
    p2 = np.array([1.0])

    melt = rng.uniform(-50.0, 900.0, size=(len(p1), len(p2), NY, NX))
    ds = xr.Dataset(
        {
            'melt_rate': (
                ('p1', 'p2', 'y', 'x'),
                melt,
            )
        },
        coords={'p1': p1, 'p2': p2, 'y': y, 'x': x},
    )

    basins = xr.DataArray(
        rng.integers(0, NBASIN, size=(NY, NX)).astype(float),
        dims=('y', 'x'),
        coords={'y': y, 'x': x},
        name='basins',
    )
    mask = xr.DataArray(
        rng.random((NY, NX)) > 0.3,
        dims=('y', 'x'),
        coords={'y': y, 'x': x},
        name='mask',
    )
    return ds, basins, mask


def upstream_integrate(ds, mask, basins, reso):
    """Upstream calculate_term1/2/4 aggregation, verbatim in form."""
    cvt = reso**2 / 1e12
    out = (
        ds['melt_rate'].where(mask, np.nan).groupby(basins).sum(skipna=True)
        * cvt
    )
    return out.where(out != 0, np.nan)


def upstream_average(ds, mask, basins):
    """Upstream calculate_term3 aggregation, verbatim in form."""
    out = ds['melt_rate'].where(mask, np.nan).groupby(basins).mean()
    return out.where(out != 0, np.nan)


def test_integrate_matches_upstream_on_uniform_grid(grid):
    ds, basins, mask = grid
    area = uniform_area(ds['melt_rate'], RESO, ('y', 'x'))

    ours = integrate_by_group(
        ds['melt_rate'], area, mask, basins, ('y', 'x'), group_dim='basins'
    )
    theirs = upstream_integrate(ds, mask, basins, RESO)

    assert ours.sizes['basins'] == NBASIN
    xr.testing.assert_allclose(
        ours.transpose('p1', 'p2', 'basins'),
        theirs.transpose('p1', 'p2', 'basins'),
    )


def test_average_matches_upstream_on_uniform_grid(grid):
    ds, basins, mask = grid
    area = uniform_area(ds['melt_rate'], RESO, ('y', 'x'))

    ours = average_by_group(
        ds['melt_rate'], area, mask, basins, ('y', 'x'), group_dim='basins'
    )
    theirs = upstream_average(ds, mask, basins)

    xr.testing.assert_allclose(
        ours.transpose('p1', 'p2', 'basins'),
        theirs.transpose('p1', 'p2', 'basins'),
    )


def test_integral_is_area_weighted_not_cell_counted(grid):
    """A variable-area mesh must not reduce to counting cells."""
    ds, basins, mask = grid
    rng = np.random.default_rng(7)
    area = xr.DataArray(
        rng.uniform(0.25, 4.0, size=(NY, NX)) * RESO**2,
        dims=('y', 'x'),
        coords={'y': ds.y, 'x': ds.x},
    )

    ours = integrate_by_group(
        ds['melt_rate'], area, mask, basins, ('y', 'x'), group_dim='basins'
    )
    uniform = integrate_by_group(
        ds['melt_rate'],
        uniform_area(ds['melt_rate'], RESO, ('y', 'x')),
        mask,
        basins,
        ('y', 'x'),
        group_dim='basins',
    )
    assert not np.allclose(ours.values, uniform.values)

    # check one basin explicitly against a hand-rolled sum
    b = 0
    sel = (basins == b) & mask
    expect = float(
        (ds['melt_rate'].isel(p1=0, p2=0) * area).where(sel).sum() / KG_PER_GT
    )
    got = float(ours.isel(p1=0, p2=0).sel(basins=b))
    assert got == pytest.approx(expect, rel=1e-12)


def test_unstructured_matches_structured_after_flattening(grid):
    """
    The same cells presented as a 1-D 'nCells' mesh must give the same answer.

    This is the property that lets one implementation serve both the ISMIP grid
    and the MALI mesh.
    """
    ds, basins, mask = grid
    area2d = uniform_area(ds['melt_rate'], RESO, ('y', 'x'))

    structured = integrate_by_group(
        ds['melt_rate'], area2d, mask, basins, ('y', 'x'), group_dim='basins'
    )

    ncells = NY * NX
    flat = xr.Dataset(
        {
            'melt_rate': (
                ('p1', 'p2', 'nCells'),
                ds['melt_rate'].values.reshape(
                    ds.sizes['p1'], ds.sizes['p2'], ncells
                ),
            )
        },
        coords={'p1': ds.p1, 'p2': ds.p2},
    )
    area1d = xr.DataArray(np.full(ncells, RESO**2), dims=('nCells',))
    basins1d = xr.DataArray(basins.values.reshape(ncells), dims=('nCells',))
    mask1d = xr.DataArray(mask.values.reshape(ncells), dims=('nCells',))

    unstructured = integrate_by_group(
        flat['melt_rate'],
        area1d,
        mask1d,
        basins1d,
        ('nCells',),
        group_dim='basins',
    )

    np.testing.assert_allclose(
        structured.transpose('p1', 'p2', 'basins').values,
        unstructured.transpose('p1', 'p2', 'basins').values,
    )


def test_empty_group_becomes_nan(grid):
    """A basin with no contributing cells must drop out, not read as zero."""
    ds, basins, mask = grid
    area = uniform_area(ds['melt_rate'], RESO, ('y', 'x'))
    # mask out basin 1 entirely
    mask = mask & (basins != 1)

    ours = integrate_by_group(
        ds['melt_rate'], area, mask, basins, ('y', 'x'), group_dim='basins'
    )
    assert bool(ours.sel(basins=1).isnull().all())
    assert bool(ours.sel(basins=0).notnull().any())


def test_term1_targets_passed_through(grid):
    ds, basins, mask = grid
    area = uniform_area(ds['melt_rate'], RESO, ('y', 'x'))
    melt_obs = pd.DataFrame(
        {
            'BMR (Gt/yr)': [10.0, 20.0, 30.0, 40.0],
            'BMR uncert (Gt/yr)': [1.0, 2.0, 3.0, 4.0],
        }
    )

    model, obs_mean, obs_sigma = calculate_term1(
        ds, area, mask, basins, melt_obs, ('y', 'x')
    )
    assert model.sizes['basins'] == NBASIN
    np.testing.assert_allclose(obs_mean.values, [10.0, 20.0, 30.0, 40.0])
    np.testing.assert_allclose(obs_sigma.values, [1.0, 2.0, 3.0, 4.0])


def test_term3_is_warm_minus_cold(grid):
    ds, basins, mask = grid
    area = uniform_area(ds['melt_rate'], RESO, ('y', 'x'))
    cold = ds
    warm = ds.copy(deep=True)
    warm['melt_rate'] = warm['melt_rate'] * 2.0

    targets = xr.Dataset(
        {
            'melt_rate': ('basins', np.zeros(NBASIN)),
            'melt_rate_uncert': ('basins', np.ones(NBASIN)),
        },
        coords={'basins': np.arange(NBASIN)},
    )

    model, obs_mean, obs_sigma = calculate_term3(
        cold, warm, targets, targets, area, mask, basins, ('y', 'x')
    )

    cold_mean = average_by_group(
        cold['melt_rate'], area, mask, basins, ('y', 'x')
    )
    xr.testing.assert_allclose(model, cold_mean)
    np.testing.assert_allclose(
        obs_sigma.values, np.sqrt(2.0) * np.ones(NBASIN)
    )


def test_stack_cells_rejects_missing_dim(grid):
    ds, _, _ = grid
    area = uniform_area(ds['melt_rate'], RESO, ('y', 'x'))
    with pytest.raises(ValueError, match='not found in dimensions'):
        integrate_by_group(ds['melt_rate'], area, None, area, ('nCells',))
