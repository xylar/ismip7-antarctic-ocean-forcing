# Handoff: porting the MALI ISMIP7 melt calibration into Compass

For the agent implementing a new **ISMIP7 calibration test group** in Compass.

The prototype in this directory answered the science questions and produced calibrated
parameters; it was never meant to be the production workflow. This document says what exists,
what to reuse, what to rebuild, and — most usefully — which mistakes have already been made so
they need not be made twice.

Read [`mali-melt-calibration-plan.md`](mali-melt-calibration-plan.md) for the plan and
[`survey-findings.md`](survey-findings.md) for the evidence. Findings are cited below as F1,
F6, … and each is worth reading before touching the corresponding area.

---

## 1. Where everything is

| What | Path | Branch | Pushed to |
|---|---|---|---|
| **Compass worktree (yours)** | `/home/ac.xylar/compass/add-ismip7-calibration` | `add-ismip7-calibration` | not yet |
| Compass main checkout | `/home/ac.xylar/compass/main` | `main` | — |
| MALI with the Burgard melt module | `.../MALI-melt-calibration/mali/add-burgard-melt-param` | `mali/add-burgard-melt-param` | `xylar/E3SM` |
| Mask tool | `~/mpas_work/MPAS-Tools/add-ismip7-mali-masks` | `add-ismip7-mali-masks` | `xylar/MPAS-Tools` |
| Calibration prototype + docs | `.../MALI-melt-calibration/add-mali-melt-calibration` | `add-mali-melt-calibration` | `xylar/ismip7-antarctic-ocean-forcing` (temporary; see below) |
| Ensemble output (28 states x 2 forms) | `.../MALI-melt-calibration/work/ensemble/` | — | — |
| Remapped forcing (28 states) | `.../MALI-melt-calibration/work/forcing/` | — | — |
| ISMIP7 source datasets | `/lcrc/group/e3sm/ac.xylar/ismip7/forcing-data/data/AIS/` | — | — |

The worktree is branched from `origin/main` at `b60c7fb26`.

`add-mali-melt-calibration` is a **temporary** branch: it records how this work was done and
carries the feedback document for the AIS Ocean Forcing focus group, and will not be merged
upstream.  The mesh-agnostic analysis (`terms.py`, `quadratic.py`, `calibrate.py`) will be
offered to `ismip7-antarctic-ocean-forcing` separately, from a new focused branch, as the MALI
worked example.  So do not assume any of that code is upstream when you reach for it.

**Environment.** Xylar is deploying a Compass pixi environment and load script that includes
Albany. Until then, MALI can be built without Albany for melt-only work
(`config_velocity_solver = 'none'`), which is all the calibration needs.

---

## 2. First task: find the shared framework

Before writing the calibration test group, work out what `ismip7_forcing`, `ismip7_run` and the
new calibration group should share, and lift it into a framework module. The aim is to use
what exists and to change other repositories only where necessary.

Concrete starting points, all in `compass/landice/`:

* **`tests/ismip7_forcing/create_mapfile.py`** — `build_mapping_file()`, built on
  `mpas_tools.scrip.from_mpas` plus `ESMF_RegridWeightGen`. **This may make the MPAS-Tools mask
  tool unnecessary** — see §4.
* **`tests/ismip7_forcing/fracture/remap_utils.py`** — `extrapolate_source()`,
  `open_rename_and_trim()`, `add_xtime_and_write()`. Remapping-adjacent helpers living in a
  sub-test-case; obvious framework candidates.
* **`tests/ismip7_forcing/ice_sheet_params.py`** — shared parameters, already factored out.
* **`iceshelf_melt.py`** — has `calc_mean_TF(geometry_file, forcing_file)`, the basin-mean
  thermal forcing. Directly relevant to the ISMIP6/semi-local form.
* **`extrapolate.py`**, **`util.py`**, **`mesh.py`** — check before writing anything new.
* **`tests/ismip7_forcing/ismip7_forcing.cfg`** — the config-option conventions to follow
  (`base_path_ismip7`, `mali_mesh_file`, `mali_mesh_name`, ESMF task count, …).

Report what you find before refactoring: consolidation that touches `ismip7_forcing` or
`ismip7_run` affects Trevor's and others' work and should be agreed first.

---

## 3. MALI build

Use **`mali/add-burgard-melt-param`** rather than the Compass `MALI-Dev` submodule.

Convenient fact: the submodule currently points at `ee0f74cc35`, which is *exactly* the
merge-base of that branch. Pointing it at the branch is a clean fast-forward of two commits, not
a divergence.

The branch adds `config_basal_mass_bal_float = 'ismip7'`: the Burgard **local** quadratic,
protocol Eq. (1). The semi-local form was deliberately **not** implemented — with constant
salinity it is algebraically identical to the ISMIP6 non-local method MALI already has (F7).
Use `'ismip6'` for that behaviour. Both paths exist, so the port does not need the
local-vs-semi-local question settled before it starts (§8).

---

## 4. Does the mask tool survive?

`interpolate_ismip7_masks_to_mali.py` (MPAS-Tools) remaps the ISMIP7 masks onto a MALI mesh with
pyremap and writes `ismip7BasinNumber`, `ismip6shelfMelt_basin`, `ismip7BFRNBin`,
`ismip7FloatingMask`, `ismip7ShelfRegion` and a zero `ismip6shelfMelt_deltaT`, plus a
bidirectional cross-check against `regionCellMasks`.

**Decide whether it is still needed.** Compass's `build_mapping_file()` may cover the remapping,
in which case the tool reduces to the mask-derivation logic and could live in the test group. Do
not keep it out of momentum; equally, do not lose the cross-check, which is what caught the
basin-numbering trap.

---

## 5. Sourcing `parameter_selection_toolbox.py`

The objective function comes from `parameter_selection_toolbox.py` in
`ismip7-antarctic-ocean-forcing`. There is no conda-forge package, and the module would not be
in one anyway. Options, to be decided and documented:

1. **Git submodule** of `ismip7-antarctic-ocean-forcing` in Compass — exact provenance, but adds
   a submodule for one file.
2. **Vendored copy with a recorded tag** — the upstream repository is supposed to tag each
   toolbox update, so record the tag and a checksum, and add a test that the copy still matches.
3. **Reimplement** the objective in Compass against the documented equations — the most work and
   the greatest risk of silently diverging from the published numbers.

The prototype imports it by path (`replicate.load_upstream_toolbox()`), which is fine for a
prototype and not for Compass. Whatever you choose, keep the ability to reproduce the published
8 km percentiles exactly (§9) — that is the check that the toolbox is being driven correctly.

---

## 6. What to build

**Masks** on the MALI mesh (§4).

**Forcing.** Remap ISMIP7 thermal forcing **and salinity** onto the mesh.
`ismip7_forcing/ocean_thermal/process_thermal_forcing.py` currently handles **thermal forcing
only — there is no salinity path**. The Burgard local form needs `so`; ISMIP7 distributes it for
every model and scenario. If MALI ends up using the semi-local/ISMIP6 form, salinity is not
needed, so sequence this after the §8 decision if you want to avoid possibly wasted work.

**Ensemble.** One single-timestep MALI run per ocean state. **Not** one per parameter value:
melt is *exactly* proportional to `K` (F1, measured to machine precision in F9), so the
parameter sweep is a scaling of one run. That is 28 runs, not ~1300. Do not "improve" this away.

**Aggregation and objective.** Aggregate melt to basins, BFRN bins and PIG/Dotson regions, then
drive the toolbox. `terms.py`, `quadratic.py` and `calibrate.py` here are mesh-agnostic and are
the intended reusable pieces; `forcing.py`, `ensemble.py` and `verify.py` are prototypes for run
setup and are what the test group replaces.

**The `dT_b` fit**, protocol §4.2.1 option 2 — after parameter selection, not before. The
toolbox's `optimise_deltaT()` does it: per basin, grid-search `dT` in [-2, 2] minimising
|basin-integrated melt − IMBIE observed|. Never implemented here, deliberately.

**Validation** using Compass's baseline-comparison framework (§9).

**Docs and registration.** Add
`docs/developers_guide/landice/test_groups/ismip7_calibration.rst` alongside the existing
`ismip7_forcing.rst` / `ismip7_run.rst`, and register the group in `compass/landice/__init__.py`.

---

## 7. Traps already hit

Each of these cost real time here.

**Basin numbering is off by one.** ISMIP7 `basinNumber` is 0-based (0-15; basin 9 = Eastern
Amundsen, 14 = Ronne-Filchner); MALI `ismip6shelfMelt_basin` is 1-based (1-16). Using one where
the other is expected mis-assigns every basin and still produces plausible numbers (F2).

**Aggregation must be area-weighted.** On a 4-20 km mesh a plain mean is wrong. Integrals use
`melt * areaCell`; means are area-weighted.

**The diagnostic run evolves the geometry.** Melt is applied over the timestep and thins the ice
by up to 1.67 m, so `lowerSurface`/`thickness` in the output are *post*-step while `TFdraft` was
computed *pre*-step. Any check pairing output geometry with output melt is wrong by up to a metre
of draft; reconstruct the initial draft from the mesh file instead (F8).

**MPAS reads dimensions from the input stream only.** `config_nISMIP6OceanLayers` must be set
from the forcing file, or the 3-D fields are allocated against a zero-length dimension and the
run dies trying to allocate hundreds of GB.

**Forcing files must not contain a variable whose name matches a dimension** (e.g. a variable
`nISMIP6OceanLayers`), and `Time` must be written **UNLIMITED**. Either breaks MPAS's reader in
ways that are painful to diagnose.

**The partition file must be named `graph.info.part.<ntasks>`**, matching
`config_block_decomp_file_prefix`. A symlink called `graph.info` is silently not found.

**Seconds per year.** MALI's `scyr` is 31,536,000 (365-day); multimelt uses 31,556,926.08
(365.2422-day). The 0.066% difference will show up in any comparison against the reference
implementation.

**The published J4 weighting uses 2 of 18 available observations** — PIG only, 2009 and 2012 —
and excludes Dotson entirely. As weighted, J4 cannot discriminate between melt modules; with
Dotson it can (F14). Getting this wrong reverses which parameterisation appears better; it did
here, and the error survived one round of review.

**The objective normalises each term by its own median** over the parameter ensemble, so a
minimised `I` is **not comparable between melt modules** — only between parameter values within
one module (feedback item A11).

**Solve the optimal scale, do not scan it.** For a weighted L1 misfit the minimiser is a
weighted median of `t_i/u_i` with weights `w_i|u_i|`. A scale *grid* silently clips and produces
spurious "improvements" (F16, method note).

**Activate the environment.** `pyremap`, `mpas_tools` and friends shell out to `ncks`,
`ncremap`, `ESMF_RegridWeightGen`. Calling python by absolute path leaves those off `PATH` and
they fail with a bare `FileNotFoundError`.

---

## 8. Still open

**Local or semi-local?** Evidence points slightly to local, not decisively:

* J2 (buttressing bins) favours local, 96.8% of target draws; J4 favours local, 98.7% (F14)
* J3 (warm−cold sensitivity) favours semi-local at its own optimum, but **reverses** at the
  calibrated parameter, where local is better by 17.3% (F15)
* local extrapolates *more conservatively* to the extreme warm anchor (x26.8 vs x33.6) (F15)
* J1 does not discriminate either way

Both forms are available in MALI, so the port does not need this settled. It does determine
whether salinity forcing is required (§6).

**The initial condition** is a database pointer; it changes the numbers only slightly.

**A negative result worth not repeating.** Generalising to `melt ~ K|TF_loc|^p|<TF>_b|^q` gives
no robust improvement: only J3 cross-validates, at (0.40, 1.80), beating semi-local by 0-8%; J2
and J4 have no interior optimum and their fits do not generalise (F16). No exogenous predictor
available — TF, stratification, water column, bed depth — explains *where* the forms fail. Shelf
area does (r = +0.85) but is endogenous and unusable. The missing physics is cavity overturning,
which no power law of TF can stand in for.

---

## 9. Numbers the port should reproduce

Use these as baselines; each was verified here.

| Check | Expected |
|---|---|
| MALI melt vs the Python reference, on MALI's own `TFdraft`/`Sdraft` | agreement to **7e-16** (F6) |
| MALI's vertical interpolation, all four code paths | machine precision (F8) |
| Melt linearity in `K` | total melt / `K` identical to all digits (F9) |
| Replication of the published 8 km calibration | **4.75e-5 / 8.5e-5 / 1.375e-4** exactly |
| MALI mesh, published weighting, 28 states | K = 5.00e-5 / 9.25e-5 / 1.45e-4 |
| **MALI mesh, full weighting, 28 states** | **K = 5.00e-5 / 7.75e-5 / 1.15e-4** (F12, F14) |
| ISMIP6 non-local on the MALI mesh, full weighting | gamma0 = 7453 / 12200 / 21000 |

Figures in the workspace root: `mali_melt_parameter_distribution.png`,
`mali_melt_term_comparison.png`, `mali_melt_pq_landscape.png`, `melt_local_vs_nonlocal.png`.
