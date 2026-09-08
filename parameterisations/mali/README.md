# MALI

Calibration of the [MALI](https://github.com/MALI-Dev/E3SM) (MPAS-Albany Land Ice) sub-shelf
melt module following the ISMIP7 Antarctic ice–ocean protocol, Sect. 4.2.

This is the MALI counterpart to the quadratic, PICO and LADDIE examples in the parent
directory. It uses the same parameter-selection toolbox
(`../parameter_selection_toolbox.py`) and the same ISMIP7 calibration datasets, but computes
melt rates with MALI itself, on MALI's unstructured, variable-resolution Antarctic mesh,
rather than with a standalone parameterisation on the ISMIP 8 km grid.

## Contents

* [`mali-melt-calibration-plan.md`](mali-melt-calibration-plan.md) — the plan: what the
  calibration involves, the approach, the decisions taken and still open, and the phased
  schedule. **Start here.**
* [`survey-findings.md`](survey-findings.md) — the evidence behind the plan: what exists
  already in MALI, Compass, MPAS-Tools and the ISMIP7 datasets, plus the results verified
  during the survey. Read this when you want to check a conclusion rather than take it on
  trust.
* [`protocol-and-toolbox-questions.md`](protocol-and-toolbox-questions.md) — running log of
  points where the protocol manuscript or this repository's toolbox was ambiguous or
  inconsistent, split into manuscript feedback for the focus group and concrete repository
  changes to propose.

Tools will be added here as they are written; the proposed layout is in §4.4 of the plan.

## Status

Draft plan under review. Several decisions still need input from the MALI developers and the
ISMIP7 AIS Ocean Focus Group — see §5 of the plan.
