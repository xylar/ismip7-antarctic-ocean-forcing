"""
The ISMIP7 "quadratic local" melt parameterisation, Burgard et al. (2022).

Protocol Eq. (1):

.. code-block:: none

    m = K sin(theta) (rho_o/rho_i) (c_o/L_i)^2 beta_S S_loc
        (g/(2|f|)) |TF_loc| TF_loc

This is the reference implementation used both to replicate the published
8 km calibration (see :mod:`mali_melt_calib.replicate`) and as the
specification for the MALI Fortran version on branch
``mali/add-burgard-melt-param``.

It is written to agree bit-for-bit with ``multimelt.melt_functions
.quadratic_mixed_slope``, which is what the protocol's own worked example uses,
so that any difference in the calibrated ``K`` is attributable to the mesh and
the ice-sheet model rather than to the melt formula.

Note on the Coriolis parameter: multimelt uses a *constant* ``f = 1.4e-4``
rather than a latitude-varying one.  On the MALI mesh we could use ``latCell``
instead, but that would change ``K``; see ``protocol-and-toolbox-questions.md``
A4.
"""

from __future__ import annotations

import numpy as np

#: seconds per year, as used by multimelt
SECONDS_PER_YEAR = 31556926.080000002

#: (rho_sw * c_pw) / (rho_i * L_i), K^-1 -- multimelt's ``melt_factor``
MELT_FACTOR = 0.013338444158574889

#: specific heat capacity of seawater, J kg^-1 K^-1
C_PO = 3974.0

#: latent heat of fusion of ice, J kg^-1
L_I = 334000.0

#: haline contraction coefficient, from Lazeroms et al.
BETA_S = 0.000786

#: gravitational acceleration, m s^-2
GRAVITY = 9.81

#: Coriolis parameter, s^-1 (constant, as in multimelt)
F_CORIOLIS = 0.00014

#: ice density used by the protocol's worked example, kg m^-3
ICE_DENSITY = 918.0


def u_factor(salinity):
    """
    The velocity-scale factor of Jenkins et al. (2018).

    ``(c_o / L_i) * beta_S * g / (2 |f|) * S_loc``

    Parameters
    ----------
    salinity : xarray.DataArray or float
        Practical salinity at the ice draft.

    Returns
    -------
    Same type as ``salinity``.
    """
    return (
        (C_PO / L_I) * BETA_S * (GRAVITY / (2.0 * abs(F_CORIOLIS))) * salinity
    )


def local_quadratic_melt(
    k, thermal_forcing, salinity, slope, thermal_forcing_avg=None
):
    """
    Melt rate from the quadratic local parameterisation, in kg m-2 yr-1.

    Parameters
    ----------
    k : float or xarray.DataArray
        The calibration parameter ``K``.
    thermal_forcing : xarray.DataArray
        Local thermal forcing at the ice draft, in K.
    salinity : xarray.DataArray
        Practical salinity at the ice draft.
    slope : float or xarray.DataArray
        Ice-draft slope angle in radians (positive).  A scalar gives the
        "constant Antarctic-mean slope" variant; a field gives the
        slope-dependent one.
    thermal_forcing_avg : xarray.DataArray, optional
        Thermal forcing to use in the ``|TF|`` factor.  Defaults to
        ``thermal_forcing``, giving the *local* form; pass a shelf- or
        basin-average for the *semi-local* form of protocol Eq. (2).

    Returns
    -------
    xarray.DataArray
        Melt rate in kg m-2 yr-1, positive for melting.

    Notes
    -----
    Melt is exactly linear in ``k``, which is what allows the calibration to
    use one model run per ocean state rather than one per (state, parameter)
    pair; see ``survey-findings.md`` F1.
    """
    if thermal_forcing_avg is None:
        thermal_forcing_avg = thermal_forcing

    melt_m_per_s = (
        k
        * MELT_FACTOR
        * u_factor(salinity)
        * thermal_forcing
        * abs(thermal_forcing_avg)
        * np.sin(slope)
    )
    return melt_m_per_s * SECONDS_PER_YEAR * ICE_DENSITY


def draft_slope(draft, dx, dy, x_dim='x', y_dim='y'):
    """
    Ice-draft slope angle on a structured grid, in radians.

    Reproduces the centred-difference scheme of
    ``multimelt.plume_functions.check_slope_one_dimension``, including its
    one-sided fallbacks at NaN neighbours and its substitution of zero where
    the slope cannot be computed at all.

    Parameters
    ----------
    draft : xarray.DataArray
        Ice draft (negative below sea level), on a regular grid.
    dx, dy : float
        Grid spacing in m.
    x_dim, y_dim : str, optional
        Names of the horizontal dimensions.

    Returns
    -------
    xarray.DataArray
        Slope angle in radians.
    """
    slope_x = _one_dimensional_slope(draft, x_dim, dx)
    slope_y = _one_dimensional_slope(draft, y_dim, dy)
    return np.arctan(np.sqrt(slope_x**2 + slope_y**2))


def _one_dimensional_slope(draft, dim, spacing):
    """One-sided-tolerant centred difference, as in multimelt."""
    shifted_minus = draft.shift({dim: -1})
    shifted_plus = draft.shift({dim: 1})

    both = (shifted_minus - shifted_plus) / np.sqrt((2.0 * spacing) ** 2)
    right = (draft - shifted_plus) / np.sqrt(spacing**2)
    left = (shifted_minus - draft) / np.sqrt(spacing**2)

    slope = both.combine_first(right).combine_first(left)
    return slope.where(np.isfinite(slope), 0.0)


def mean_slope(draft, floating, dx, dy, x_dim='x', y_dim='y'):
    """
    Antarctic-mean ice-draft slope angle over floating ice, in radians.

    This is the "constant slope" used by the protocol's worked example.  Note
    that it is a mean of *angles*, not of ``sin(theta)``, and that its value is
    resolution dependent; both points are raised in
    ``protocol-and-toolbox-questions.md`` A4.

    Parameters
    ----------
    draft : xarray.DataArray
        Ice draft on a regular grid.
    floating : xarray.DataArray
        Boolean mask, True on floating ice.
    dx, dy : float
        Grid spacing in m.
    x_dim, y_dim : str, optional
        Names of the horizontal dimensions.

    Returns
    -------
    float
        Mean slope angle in radians.
    """
    slope = draft_slope(draft, dx, dy, x_dim=x_dim, y_dim=y_dim)
    return float(slope.where(floating).mean())
