"""
Tests for the quadratic local melt parameterisation.

The key property is agreement with ``multimelt``, which the protocol's own
worked example uses: if our formula and theirs disagree, any difference in the
calibrated ``K`` would be an artefact of the implementation rather than of the
mesh or the ice-sheet model.  ``multimelt`` is optional, so those tests skip
when it is not installed.
"""

import numpy as np
import pytest
import xarray as xr
from mali_melt_calib.quadratic import (
    F_CORIOLIS,
    ICE_DENSITY,
    MELT_FACTOR,
    SECONDS_PER_YEAR,
    draft_slope,
    local_quadratic_melt,
    mean_slope,
    u_factor,
)


@pytest.fixture
def fields():
    rng = np.random.default_rng(42)
    ny, nx = 9, 7
    coords = {'y': np.arange(ny) * 8000.0, 'x': np.arange(nx) * 8000.0}
    tf = xr.DataArray(
        rng.uniform(-1.0, 4.0, size=(ny, nx)), dims=('y', 'x'), coords=coords
    )
    so = xr.DataArray(
        rng.uniform(33.0, 35.0, size=(ny, nx)), dims=('y', 'x'), coords=coords
    )
    draft = xr.DataArray(
        -rng.uniform(50.0, 1500.0, size=(ny, nx)),
        dims=('y', 'x'),
        coords=coords,
    )
    return tf, so, draft


def test_melt_is_linear_in_k(fields):
    """The property the whole ensemble strategy rests on (findings F1)."""
    tf, so, _ = fields
    slope = 0.005

    m1 = local_quadratic_melt(1.0, tf, so, slope)
    m7 = local_quadratic_melt(7.0, tf, so, slope)

    xr.testing.assert_allclose(m7, 7.0 * m1)


def test_melt_sign_follows_thermal_forcing(fields):
    """Positive TF melts, negative TF refreezes."""
    tf, so, _ = fields
    melt = local_quadratic_melt(1.0e-5, tf, so, 0.005)

    # index rather than .where(), which would leave NaNs that compare False
    warm = melt.values[tf.values > 0.0]
    cold = melt.values[tf.values < 0.0]
    assert warm.size > 0 and cold.size > 0
    assert (warm > 0.0).all()
    assert (cold < 0.0).all()


def test_u_factor_matches_definition():
    salinity = xr.DataArray([34.0])
    expected = (3974.0 / 334000.0) * 0.000786 * (9.81 / (2.0 * F_CORIOLIS))
    np.testing.assert_allclose(u_factor(salinity).values, expected * 34.0)


def test_matches_multimelt(fields):
    """Bit-for-bit agreement with the protocol's reference implementation."""
    mf = pytest.importorskip('multimelt.melt_functions')
    tf, so, _ = fields
    k = 8.5e-5
    slope = 0.0051117

    ours = local_quadratic_melt(k, tf, so, slope)
    theirs = (
        mf.quadratic_mixed_slope(k, MELT_FACTOR, tf, tf, u_factor(so), slope)
        * SECONDS_PER_YEAR
        * ICE_DENSITY
    )
    xr.testing.assert_allclose(ours, theirs)


def test_constants_match_multimelt():
    """Guard against multimelt changing its constants under us."""
    const = pytest.importorskip('multimelt.constants')
    assert MELT_FACTOR == const.melt_factor
    assert SECONDS_PER_YEAR == const.yearinsec
    assert F_CORIOLIS == const.f_coriolis


def test_draft_slope_matches_multimelt(fields):
    """Our slope reproduces multimelt's finite-difference scheme."""
    pf = pytest.importorskip('multimelt.plume_functions')
    _, _, draft = fields
    dx = dy = 8000.0

    ours = draft_slope(draft, dx, dy)

    xslope = pf.check_slope_one_dimension(
        draft, draft.shift(x=1), draft.shift(x=-1), dx
    )
    yslope = pf.check_slope_one_dimension(
        draft, draft.shift(y=1), draft.shift(y=-1), dy
    )
    theirs = np.arctan(np.sqrt(xslope**2 + yslope**2))

    xr.testing.assert_allclose(ours, theirs)


def test_mean_slope_is_a_scalar_over_floating_ice(fields):
    _, _, draft = fields
    floating = xr.DataArray(
        np.arange(draft.size).reshape(draft.shape) % 3 == 0,
        dims=draft.dims,
        coords=draft.coords,
    )
    value = mean_slope(draft, floating, 8000.0, 8000.0)
    assert isinstance(value, float)
    assert 0.0 < value < np.pi / 2


def test_semi_local_variant_uses_the_average(fields):
    """Passing thermal_forcing_avg gives protocol Eq. (2) rather than (1)."""
    tf, so, _ = fields
    tf_avg = tf.mean()

    local = local_quadratic_melt(1.0e-5, tf, so, 0.005)
    semi = local_quadratic_melt(
        1.0e-5, tf, so, 0.005, thermal_forcing_avg=tf_avg
    )
    assert not np.allclose(local.values, semi.values)

    expected = local * abs(tf_avg) / abs(tf)
    xr.testing.assert_allclose(semi, expected)
