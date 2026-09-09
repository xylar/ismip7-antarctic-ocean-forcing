# MALI example: calibration on an unstructured, variable-resolution mesh

A worked example of the §4.2 calibration for **MALI** (MPAS-Albany Land Ice), contributed as a
fourth example alongside the quadratic, PICO and LADDIE ones.

It differs from those in two ways, and both are the point of contributing it.

**The melt rates come from the ice-sheet model.** The top-level `README.md` says:

> Note that here we calculate melt rates independently. In ISMIP7/meltMIP, the melt rate should
> be calculated through the ice-sheet model code and grid instead.

That is what this example does. Every melt field is produced by MALI's own Fortran melt module,
on MALI's own mesh, from MALI's own geometry — not by a Python reimplementation.

**The mesh is unstructured and variable-resolution.** The toolbox's term functions assume a
structured grid: a scalar `reso` for the cell area, `groupby` over a 2-D field, and
`stack(grid=('x','y'))`. The Antarctic MALI mesh used here refines from 20 km inland to 4 km at
the grounding line, so cell area is a *field*, not a constant, and an area-weighted mean is not
the same as a plain mean. `terms.py` here is a mesh-agnostic rewrite that takes a cell-area
array and a generic cell dimension, and reduces to the existing behaviour on a uniform grid.

## Both parameterisations

MALI implements both quadratic forms the protocol offers, and the example calibrates each:

| | protocol | MALI option | parameter |
|---|---|---|---|
| **local** | Eq. (1), Burgard et al. (2022) | `config_basal_mass_bal_float = 'ismip7'` | `K` |
| **semi-local** | Eq. (2) | `config_basal_mass_bal_float = 'ismip6'` | `gamma0` |

The second row deserves a note, because it surprised us. Written out with the basin temperature
correction in place, the two are

```
ISMIP6 non-local:    melt = C6 * (TF_loc + dT) * |<TF> + dT|
Burgard semi-local:  melt = C7 * <S> * (TF_loc + dT) * |<TF> + dT|
```

— the same functional form, differing only in how the constant is decomposed. **With the
salinity held constant, Eq. (2) *is* the ISMIP6 non-local parameterisation**, and every `K` has
an exactly equivalent `gamma0`. With a basin-mean salinity they differ only by the per-basin
spread of `<S>`, measured at ±0.8% on this mesh, which a fitted `dT_b` largely absorbs. So MALI
implements Eq. (2) by keeping its existing ISMIP6 code path rather than adding a second one.
The **local** form is the one that genuinely differs.

## Results

Calibrated on the Antarctic 4–20 km MALI mesh, 28 ocean states, 100,000 samples:

| weighting | | 5th | 50th | 95th |
|---|---|---:|---:|---:|
| as published | local, `K` | 5.00e-5 | 9.25e-5 | 1.45e-4 |
| as published | semi-local, `gamma0` (m/yr) | 7792 | 1.73e4 | 2.91e4 |
| everything distributed | local, `K` | 5.00e-5 | 7.75e-5 | 1.15e-4 |
| everything distributed | semi-local, `gamma0` (m/yr) | 7453 | 1.22e4 | 2.10e4 |

"As published" reproduces the weighting of `parameter_selection_quadratic_example.ipynb` — four
ocean models for J3, PIG alone in 2009 and 2012 for J4. "Everything distributed" uses all seven
ocean models and all thirteen observation years for both shelves.

**The local result sits on the published distribution.** Against the published 8 km numbers of
4.75e-5 / 8.5e-5 / 1.375e-4, the 5th and 50th percentiles agree to within a single grid step of
the parameter sweep. An independent Fortran implementation, on a variable-resolution unstructured
mesh, with masks remapped from the ISMIP grid and melt aggregated on the model's own floating
cells, lands where the protocol's own Python does on its own grid. We take that as evidence the
mesh-agnostic generalisation is faithful rather than merely self-consistent.

Two observations that may interest the focus group:

* Widening from the published weighting to everything distributed **tightens** both
  distributions and lowers the medians — more data, more constraint.
* Both distributions are **bimodal**: a narrow spike at the low end and a broad hump above it.
  That is the random term weighting acting on terms that disagree — J1 and J3 prefer a smaller
  parameter than J2 and J4 — rather than anything about MALI.

## What is here

| file | what it is |
|---|---|
| `example_mali.py` | the example: calibrates both forms and prints the table above |
| `mali_melt_calib/terms.py` | **mesh-agnostic** `calculate_term1..4`, area-weighted |
| `mali_melt_calib/quadratic.py` | protocol Eq. (1) as a reference implementation, with the constants collected in one place |
| `mali_melt_calib/calibrate.py` | aggregation on a MALI mesh and the toolbox driver |
| `mali_melt_calib/datasets.py` | registry of the ocean states and which term each feeds |
| `mali_melt_calib/replicate.py` | reproduces the published 8 km calibration, as a check |
| `tests/` | 16 tests, including that the mesh-agnostic terms reproduce the structured-grid answers on a uniform-area mesh |

## Running it

```bash
python example_mali.py --ensemble-dir <dir> --masks <masks-on-mali-mesh.nc>
```

`--ensemble-dir` holds `ismip7/<state>/` and `ismip6/<state>/`, each with the `output_melt.nc`
of a single-timestep MALI run for that ocean state.

Producing those runs needs MALI and the remapped forcing, and is being packaged as a
[Compass](https://github.com/MPAS-Dev/compass) test group so that it is reproducible rather than
bespoke. Two things about it are worth knowing even if you never run MALI:

* Melt is **exactly** proportional to `K` (and to `gamma0`) for a fixed geometry and fixed
  `dT_b`, to machine precision. So the parameter ensemble needs **one run per ocean state**, not
  one per (state, parameter) pair — 28 runs rather than about 1300. Any ice-sheet model whose
  melt module is linear in its calibration parameter can do the same.
* `dT_b` is zero throughout, following §4.2.1's second option, in which the corrections are
  calibrated *after* the parameter rather than before.

## Validation

`replicate.py` runs the published 8 km calibration through this code path and reproduces all
three published percentiles exactly, which is what gives us confidence that the generalised term
functions are equivalent on a structured grid.

MALI's melt field was checked against `quadratic.py` evaluated on MALI's own draft-level thermal
forcing and salinity, agreeing to 7e-16; and MALI's vertical interpolation of the forcing to the
ice draft was checked independently on all four of its code paths, agreeing to machine precision.

## Feedback

Working through the protocol and the toolbox produced a list of points where the text or the code
could be clearer, or where we think there may be a defect. Those are collected separately and
will be passed to the focus group rather than buried here.
