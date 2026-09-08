"""
Check MALI's melt field against the Python reference implementation.

Two independent things can go wrong in the Fortran: the vertical interpolation
of thermal forcing and salinity to the ice draft, and the melt expression
itself.  Comparing MALI's melt against
:func:`mali_melt_calib.quadratic.local_quadratic_melt` evaluated on **MALI's
own** ``TFdraft`` and ``Sdraft`` isolates the second from the first, so a
disagreement points at one or the other rather than at both together.

Run against a directory produced by :mod:`mali_melt_calib.ensemble`::

    python -m mali_melt_calib.verify <run_dir>
"""

from __future__ import annotations

import os
import sys

import numpy as np
import xarray as xr

from mali_melt_calib.quadratic import MALI, local_quadratic_melt

#: seconds per year, matching MALI's li_constants scyr.  Both sides of the
#: comparison must use the same year, or they disagree by 0.066%.
MALI_SCYR = MALI.seconds_per_year


def read_namelist_floats(path, keys):
    """Pull a few real-valued options out of a run's namelist."""
    values = {}
    with open(path) as handle:
        for line in handle:
            if '=' not in line:
                continue
            key, value = line.split('=', 1)
            key = key.strip()
            if key in keys:
                values[key] = value.strip().strip('\'"')
    return values


def melt_from_mali(ds):
    """
    Melt rate in kg m-2 yr-1, positive for melting, from a MALI output file.

    MALI stores ``floatingBasalMassBal`` in kg m-2 s-1 with melting negative,
    which is the opposite sign convention to the calibration toolbox.
    """
    bmb = ds['floatingBasalMassBal']
    if 'Time' in bmb.dims:
        bmb = bmb.isel(Time=-1)
    return -bmb * MALI_SCYR


def verify_run(run_dir, rtol=1.0e-10):
    """
    Compare MALI's melt with the reference formula on the same drafts.

    Parameters
    ----------
    run_dir : str
        A melt-diagnostic run directory containing ``output_melt.nc`` and
        ``namelist.landice``.
    rtol : float, optional
        Relative tolerance for the comparison.  The two should agree to
        round-off, since the reference is evaluated with MALI's own constants
        on MALI's own drafts.

    Returns
    -------
    dict
        Summary statistics.
    """
    out_file = os.path.join(run_dir, 'output_melt.nc')
    ds = xr.open_dataset(out_file, decode_times=False, decode_timedelta=False)

    params = read_namelist_floats(
        os.path.join(run_dir, 'namelist.landice'),
        {
            'config_ismip7_melt_K',
            'config_ismip7_melt_sin_slope',
            'config_ismip7_melt_coriolis',
            'config_ismip7_melt_semi_local',
        },
    )
    melt_k = float(params['config_ismip7_melt_K'])
    sin_slope = float(params['config_ismip7_melt_sin_slope'])
    coriolis = float(params['config_ismip7_melt_coriolis'])
    semi_local = params['config_ismip7_melt_semi_local'].lower() == '.true.'

    def last(name):
        da = ds[name]
        return da.isel(Time=-1) if 'Time' in da.dims else da

    tf = last('ismip6shelfMelt_TFdraft')
    so = last('ismip7shelfMelt_Sdraft')
    mali_melt = melt_from_mali(ds)

    # melting cells only; elsewhere MALI writes zero by construction
    active = (tf != 0.0) | (so != 0.0)

    constants = MALI
    if coriolis != constants.coriolis:
        constants = type(constants)(
            **{**constants.__dict__, 'coriolis': coriolis}
        )

    if semi_local:
        raise NotImplementedError(
            'the semi-local form needs the basin means to be recomputed here; '
            'verify the local form first'
        )

    # arcsin so that local_quadratic_melt's own sin() recovers the value MALI
    # was given; at these angles the two are equal to one part in 1e5 anyway
    reference = local_quadratic_melt(
        melt_k, tf, so, np.arcsin(sin_slope), constants=constants
    )

    diff = (mali_melt - reference).where(active)
    denom = abs(reference).where(active)
    rel = (abs(diff) / denom.where(denom > 0.0)).values
    rel = rel[np.isfinite(rel)]

    result = {
        'n_active': int(active.sum()),
        'max_abs_diff': float(abs(diff).max()),
        'max_rel_diff': float(rel.max()) if rel.size else float('nan'),
        'mali_melt_range': (
            float(mali_melt.where(active).min()),
            float(mali_melt.where(active).max()),
        ),
        'passed': bool(rel.size and rel.max() < rtol),
        'K': melt_k,
        'sin_slope': sin_slope,
        'semi_local': semi_local,
    }
    return result


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        print(__doc__)
        return 2
    run_dir = argv[0]

    result = verify_run(run_dir)
    print(f'run: {run_dir}')
    print(
        f'  K = {result["K"]:.4e}, sin(slope) = {result["sin_slope"]:.6f}, '
        f'semi-local = {result["semi_local"]}'
    )
    print(f'  melting cells:        {result["n_active"]}')
    lo, hi = result['mali_melt_range']
    print(f'  MALI melt kg/m2/yr:   {lo:.2f} .. {hi:.2f}')
    print(f'  max abs difference:   {result["max_abs_diff"]:.3e}')
    print(f'  max rel difference:   {result["max_rel_diff"]:.3e}')
    print(f'  {"PASS" if result["passed"] else "FAIL"}')
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
