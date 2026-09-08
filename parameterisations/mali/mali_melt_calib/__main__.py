"""
Command-line entry point: ``python -m mali_melt_calib <step>``.

The steps mirror the workflow in ``mali-melt-calibration-plan.md`` §3.2.  Steps
that are not written yet exit with a clear message naming the plan phase they
belong to, rather than failing obscurely.
"""

from __future__ import annotations

import argparse
import sys

STEPS = {
    'inputs': ('phase 1', 'check/stage the ISMIP7 and MALI input datasets'),
    'mesh': (
        'phase 3',
        'build the MALI-mesh file with basins, BFRN bins and masks',
    ),
    'forcing': (
        'phase 4',
        'remap ISMIP 8 km TF and salinity onto the MALI mesh',
    ),
    'ensemble': (
        'phase 6',
        'generate and submit the single-timestep MALI runs',
    ),
    'collect': (
        'phase 6',
        'assemble MALI output into toolbox ensemble datasets',
    ),
    'calibrate': ('phase 7', 'run the objective-function optimisation'),
    'report': ('phase 7', 'plots and chosen parameter values'),
}


def build_parser():
    parser = argparse.ArgumentParser(
        prog='mali_melt_calib',
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest='step', metavar='STEP')
    for name, (phase, help_text) in STEPS.items():
        sub.add_parser(name, help=f'({phase}) {help_text}')
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.step is None:
        parser.print_help()
        return 2

    phase, help_text = STEPS[args.step]
    print(
        f"'{args.step}' ({help_text}) is not implemented yet; it is "
        f'{phase} of the plan. See mali-melt-calibration-plan.md §6.',
        file=sys.stderr,
    )
    return 1


if __name__ == '__main__':
    sys.exit(main())
