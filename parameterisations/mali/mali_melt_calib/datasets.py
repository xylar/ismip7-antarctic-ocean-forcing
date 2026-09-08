"""
Registry of the ISMIP7 ocean states used in the calibration.

One source of truth for *which* ocean states feed *which* objective-function
term, shared by the 8 km replication (:mod:`mali_melt_calib.replicate`) and the
MALI-mesh workflow (:mod:`mali_melt_calib.forcing`), so the two cannot drift
apart.

The sets follow protocol Table 2 and the focus group's 31 July 2026 update:

* ``minimal`` -- the states marked X in Table 2: present-day climatology, two
  circum-Antarctic cold/warm pairs, and PIG in 2009 and 2012.  7 states.
* ``recommended`` -- adds the two regional cold/warm pairs marked +, which the
  focus group strongly suggests.  11 states.  **This is the default.**
* ``all`` -- every state distributed for ISMIP7.  26 states.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional, Tuple

#: default root of the ISMIP7 AIS datasets on LCRC
DEFAULT_DATA_ROOT = '/lcrc/group/e3sm/ac.xylar/ismip7/forcing-data/data/AIS'

#: ISMIP grid resolution of the calibration datasets, in km
ISMIP_RESOLUTION_KM = 8

#: ocean-model datasets: file prefix, toolbox label, basins constrained.
#: ``basins`` is None where the simulation is circum-Antarctic; otherwise it
#: lists the ISMIP7 (0-based) basins the regional domain actually covers, and
#: melt outside them must be discarded.
OCEAN_MODELS = [
    ('Mathiot_NEMO_{state}_v3_', 'mathiot', None),
    ('Timmermann_FESOM_{state}_v3_', 'timmermann', (14,)),
    ('Naughten_FESOM_ACCESS_{state}_v2_', 'naughten_ais_1', None),
    ('Naughten_FESOM_MMM_{state}_v2_', 'naughten_ais_2', None),
    ('Jourdain-Naughten_NEMO-MITgcm_{state}_', 'jourdain_naughten', (9, 14)),
    ('Naughten_MITamu-MITwed_{state}_', 'naughten_naughten', (9, 14)),
]

#: all years of Amundsen near-ice-shelf ocean observations
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

#: minimal J3 set (Table 2, marked X)
MINIMAL_MODELS = ('mathiot', 'naughten_ais_1')

#: J3 set the focus group recommends (adds those marked +)
RECOMMENDED_MODELS = (
    'mathiot',
    'naughten_ais_1',
    'jourdain_naughten',
    'naughten_naughten',
)

#: J4 years the protocol suggests as a minimum, a cold and a warm PIG state
RECOMMENDED_OBS_YEARS = (2009, 2012)


@dataclass(frozen=True)
class OceanState:
    """
    One ocean state to force the melt module with.

    Attributes
    ----------
    name : str
        Unique, filesystem-safe identifier, e.g. ``mathiot_cold``.
    kind : str
        ``'climatology'``, ``'model'`` or ``'obs'``.
    tf_file, so_file : str
        Paths to the thermal-forcing and salinity files on the ISMIP grid.
    term : str
        Which objective-function term this state feeds.
    label : str, optional
        Toolbox model label, for ``'model'`` states.
    state : str, optional
        ``'cold'`` or ``'warm'``, for ``'model'`` states.
    year : int, optional
        Observation year, for ``'obs'`` states.
    basins : tuple of int, optional
        ISMIP7 (0-based) basins the state constrains; None means all.
    """

    name: str
    kind: str
    tf_file: str
    so_file: str
    term: str
    label: Optional[str] = None
    state: Optional[str] = None
    year: Optional[int] = None
    basins: Optional[Tuple[int, ...]] = None


def climatology_files(data_root=DEFAULT_DATA_ROOT):
    """Paths to the Zhou et al. present-day climatology TF and salinity."""
    clim = os.path.join(
        data_root, 'obs', 'ocean', 'climatology', 'zhou_annual_06_nov'
    )
    stem = 'AIS_obs_ocean_climatology_zhou_annual_06_nov'
    return (
        os.path.join(clim, 'tf', 'v3', f'tf_{stem}_v3_1972-2024.nc'),
        os.path.join(clim, 'so', 'v4', f'so_{stem}_v4_1972-2024.nc'),
    )


def ocean_states(data_root=DEFAULT_DATA_ROOT, subset='recommended'):
    """
    The ocean states to run the melt module for.

    Parameters
    ----------
    data_root : str, optional
        Root of the ISMIP7 AIS datasets.
    subset : {'minimal', 'recommended', 'all'}, optional
        Which set of states to include; see the module docstring.

    Returns
    -------
    list of OceanState
    """
    if subset not in ('minimal', 'recommended', 'all'):
        raise ValueError(
            f"subset must be 'minimal', 'recommended' or 'all', got '{subset}'"
        )

    param = os.path.join(data_root, 'parameterisations', 'ocean')
    model_dir = os.path.join(param, 'ocean_modelling_data')
    obs_dir = os.path.join(param, 'ocean_observations_data')

    if subset == 'minimal':
        keep_models = MINIMAL_MODELS
        keep_years = RECOMMENDED_OBS_YEARS
    elif subset == 'recommended':
        keep_models = RECOMMENDED_MODELS
        keep_years = RECOMMENDED_OBS_YEARS
    else:
        keep_models = tuple(label for _, label, _ in OCEAN_MODELS)
        keep_years = tuple(OBS_YEARS)

    tf_file, so_file = climatology_files(data_root)
    states = [
        OceanState(
            name='climatology',
            kind='climatology',
            tf_file=tf_file,
            so_file=so_file,
            term='J1,J2',
        )
    ]

    for prefix, label, basins in OCEAN_MODELS:
        if label not in keep_models:
            continue
        for state in ('cold', 'warm'):
            stem = prefix.format(state=state)
            states.append(
                OceanState(
                    name=f'{label}_{state}',
                    kind='model',
                    tf_file=os.path.join(model_dir, f'{stem}TF.nc'),
                    so_file=os.path.join(model_dir, f'{stem}S.nc'),
                    term='J3',
                    label=label,
                    state=state,
                    basins=basins,
                )
            )

    for year in OBS_YEARS:
        if year not in keep_years:
            continue
        states.append(
            OceanState(
                name=f'obs_{year}',
                kind='obs',
                tf_file=os.path.join(obs_dir, f'Obs_{year}_TF.nc'),
                so_file=os.path.join(obs_dir, f'Obs_{year}_S.nc'),
                term='J4',
                year=year,
                basins=(9,),
            )
        )

    return states


def missing_files(states):
    """Return the (state, path) pairs whose input files are not present."""
    missing = []
    for state in states:
        for path in (state.tf_file, state.so_file):
            if not os.path.exists(path):
                missing.append((state.name, path))
    return missing
