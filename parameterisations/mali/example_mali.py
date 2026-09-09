#!/usr/bin/env python
"""
Worked example: calibrating MALI's sub-shelf melt module for ISMIP7.

Unlike the other examples in this directory, the melt rates here are **not**
computed independently in Python.  They come from MALI itself, on MALI's own
variable-resolution mesh, which is what ``parameterisations/README.md`` asks
for::

    In ISMIP7/meltMIP, the melt rate should be calculated through the
    ice-sheet model code and grid instead.

Two parameterisations are compared, both of which MALI implements:

* **local** -- protocol Eq. (1), the Burgard et al. (2022) local quadratic,
  ``config_basal_mass_bal_float = 'ismip7'``, calibrating ``K``
* **semi-local** -- protocol Eq. (2).  With the salinity held constant this is
  *algebraically identical* to the ISMIP6 non-local method, so MALI provides it
  as ``config_basal_mass_bal_float = 'ismip6'``, calibrating ``gamma0``

Both are exactly linear in their parameter for a fixed geometry, so the
parameter sweep is a rescaling of **one run per ocean state** rather than one
run per (state, parameter) pair.

Usage
-----
::

    python example_mali.py --ensemble-dir <dir> --masks <file>

``--ensemble-dir`` holds one subdirectory per ocean state for each
parameterisation (``ismip7/<state>/`` and ``ismip6/<state>/``), each
containing the ``output_melt.nc`` of a single-timestep MALI run.  See
``README.md`` for how those are produced.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

from mali_melt_calib.calibrate import (
    load_mali_static,
    load_targets,
    run_optimisation,
    toolbox_terms,
    unit_aggregates,
)
from mali_melt_calib.datasets import DEFAULT_DATA_ROOT

#: the parameter grid of the published quadratic example
K_VALUES = np.arange(0.25e-5, 3.025e-4, 0.25e-5)

#: gamma0 grid covering the same relative range, so that the two
#: parameterisations are explored with equal resolution.  The mapping is the
#: gamma0 equivalent to the published K = 8.5e-5 at a reference salinity.
GAMMA0_PER_K = 11519.0 / 8.5e-5

#: how the published quadratic example weights the terms: four ocean models for
#: J3, and PIG alone in 2009 and 2012 for J4
PUBLISHED = {}

#: everything the protocol distributes
FULL = dict(t3_models=None, t4_regions=None, t4_years=None)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--ensemble-dir', required=True,
        help='directory holding ismip7/<state>/ and ismip6/<state>/ runs',
    )
    parser.add_argument(
        '--masks', required=True,
        help='ISMIP7 masks remapped onto the MALI mesh',
    )
    parser.add_argument('--data-root', default=DEFAULT_DATA_ROOT)
    parser.add_argument('--sample-size', type=int, default=100000)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument(
        '--weighting', choices=['published', 'full', 'both'], default='both',
        help='term weighting; "published" reproduces the weighting of the '
             'quadratic example, "full" uses every dataset distributed',
    )
    args = parser.parse_args(argv)

    forms = {
        'local (K)': (
            os.path.join(args.ensemble_dir, 'ismip7'), 8.5e-5, K_VALUES,
        ),
        'semi-local (gamma0)': (
            os.path.join(args.ensemble_dir, 'ismip6'), 11519.0,
            K_VALUES * GAMMA0_PER_K,
        ),
    }
    for name, (path, _, _) in forms.items():
        if not os.path.isdir(path):
            parser.error(f'{name}: no run directory at {path}')

    static = load_mali_static(
        args.masks, os.path.join(args.ensemble_dir, 'ismip7', 'climatology')
    )
    targets = load_targets(args.data_root)

    weightings = (
        [('published', PUBLISHED), ('full', FULL)]
        if args.weighting == 'both'
        else [(args.weighting,
               PUBLISHED if args.weighting == 'published' else FULL)]
    )

    units = {
        name: unit_aggregates(path, param, static)
        for name, (path, param, _) in forms.items()
    }
    print(f'ocean models: {list(units["local (K)"]["t3"].model.values)}')
    print(f'observation years: '
          f'{[int(y) for y in units["local (K)"]["t4"].year.values]}\n')

    for label, kwargs in weightings:
        print(f'=== {label} weighting ===')
        print(f'{"":22}{"5th":>12}{"50th":>12}{"95th":>12}{"mode":>12}')
        for name, (_, _, grid) in forms.items():
            terms = toolbox_terms(units[name], targets, grid, **kwargs)
            result = run_optimisation(
                terms, grid, sample_size=args.sample_size, seed=args.seed
            )
            print(f'{name:22}' + ''.join(
                f'{result[k]:12.4g}'
                for k in ('p5', 'median', 'p95', 'mode')
            ))
        print()

    print('For reference, the published 8 km calibration of the local form '
          'gives\n  K = 4.75e-5 / 8.5e-5 / 1.375e-4')
    return 0


if __name__ == '__main__':
    sys.exit(main())
