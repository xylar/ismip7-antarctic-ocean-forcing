# Plan: calibrating the MALI sub-shelf melt parameterization for ISMIP7

**Branch:** `ismip7-antarctic-ocean-forcing` @ `add-mali-melt-calibration`
**Status:** draft for review. Nothing blocks starting (see §6); the open questions in §5.2
matter for later phases.
**Date:** 2026-09-08

**Companion documents**

* [`survey-findings.md`](survey-findings.md) — what was inspected and measured to reach the
  conclusions here. The plan states decisions; the findings document holds the evidence.
* [`protocol-and-toolbox-questions.md`](protocol-and-toolbox-questions.md) — running log of
  points of confusion, for feedback to the AIS Ocean Forcing focus group.

*Path conventions:* bare relative paths are relative to this repository's root; paths
beginning with `/` are absolute on LCRC.

---

## 1. Goal

Follow the ISMIP7 AIS ice–ocean protocol (Reese et al., §4.2) to calibrate the free
parameter(s) of MALI's sub-shelf melt module, and produce **at least three parameter values**
(5th / 50th / 95th percentile of the optimal-parameter distribution) for the ISMIP7 MALI
projections.

Secondary goals:

* Full provenance — every input dataset downloaded or derived by scripted, re-runnable steps.
* Tools delivered as an installable **Python package with a CLI**, not notebooks.
* No hard dependency on Compass; use its `landice/ismip7_forcing` and `landice/ismip7_run`
  test groups as reference implementations.

---

## 2. What the calibration involves

A two-stage procedure. Stage 1 is the MALI-specific work; stage 2 is cheap and already
implemented upstream.

### Stage 1 — build the melt ensemble

Run MALI with fixed, present-day geometry for a **single time step**, once per ocean forcing
state. The states, and the objective-function term each feeds (protocol Table 2):

| Group | Datasets | Term | Count |
|---|---|---|---|
| Present-day climatology | Zhou et al. observational climatology | J1, J2 | 1 |
| Ocean models, circum-Antarctic | Mathiot_NEMO, Naughten_FESOM_ACCESS (**minimum**), Naughten_FESOM_MMM, Timmermann_FESOM | J3 | 8 (cold+warm) |
| Ocean models, regional | Jourdain-Naughten_NEMO-MITgcm, Naughten_MITamu-MITwed (**recommended**), Haid_FESOM | J3 | 6 (cold+warm) |
| Amundsen observations | PIG/Dotson 1994–2020; **minimum PIG 2009 + 2012** | J4 | up to 13 |

Minimal set 7 ocean states; **recommended set 11** (§5.1, Q6); everything 26.

Output: melt rate in kg m⁻² yr⁻¹, positive = melting, on the MALI mesh.

**We need only one MALI run per ocean state, not per (state, parameter) pair.** Melt is
exactly linear in `gamma0` (and in `K`), so the full 120-value parameter ensemble follows by
scaling — ~15 runs instead of ~1,300. Derivation and caveats in
[`survey-findings.md`](survey-findings.md) F1; the linearity is to be confirmed numerically
before we rely on it.

### Stage 2 — parameter optimisation

`parameterisations/parameter_selection_toolbox.py` aggregates modelled melt into four terms —
**J1** basin-integrated melt, **J2** melt by buttressing bin, **J3** warm−cold basin-mean melt
difference, **J4** PIG/Dotson integrated melt per observation year — then draws 100,000
samples of the term weights and target uncertainties, minimising
`I = Σ a_i · J_i / median(J_i)` each time. The resulting distribution of optimal parameter
values gives the 5th / 50th / 95th percentiles.

The toolbox assumes a uniform structured grid, so the four terms need area-weighted,
`nCells`-based equivalents for MALI's mesh (findings 1.7). That is the main new code.

---

## 3. Approach

### 3.1 Geometry: calibrate on the un-relaxed mesh

The protocol's two §4.2 bullets are sometimes read as being in tension. They are not: bullet 1
constrains the **geometry**, bullet 2 constrains the **code, grid and parameters**. Both can be
satisfied at once, and should be.

For MALI that means using the **un-relaxed** geometry. A 10-year relaxation does not move the
geometry far, but the drift is model-specific, and reducing exactly that model-to-model
variation is the point of asking everyone to calibrate against a common present-day geometry.
This costs nothing: the E3SM inputdata meshes are pre-relaxation and already essentially
observed (findings F3 — thickness median −12.4 m vs Bedmap3, volume within 1.2%). **No
geometry replacement and no new mesh are needed.**

The relaxed initial condition is worth one extra run as a *sensitivity*, to confirm
quantitatively that the shift is small.

Caveat regardless of geometry: J1, J2 and J4 are *integrated* melt, so any ice-shelf area
mismatch enters the calibrated parameter directly; J3 (basin-*mean*) is much less sensitive.
This is intended by the protocol, but we should report MALI's per-basin shelf area against the
observed area alongside the results.

### 3.2 Workflow

```
 (0) environment          pixi env for the calibration tools
 (1) inputs               link/verify local ISMIP7 data; confirm the ISMIP7 grid definition
 (2) mesh preparation     build a MALI-mesh file (un-relaxed, present-day geometry) with:
                            geometry (thickness, bedTopography, lowerSurface)
                            ismip6shelfMelt_basin   (ISMIP7 IMBIE2 numbering, 16 basins)
                            BFRN bins               (10 bins, remapped from ISMIP 8 km)
                            floating mask + areaCell
                            PIG/Dotson region mask for J4
 (3) forcing remap        11 (or 26) ISMIP 8 km TF fields -> MALI mesh, 30 z-levels
                          (bilinear, mirroring compass ismip7_forcing/ocean_thermal)
 (4) MALI ensemble        one 1-timestep, velocity-solver-off run per ocean state
                          (+ 3 gamma0 values on one state to verify linearity)
 (5) collect              assemble pd / cold / warm / obs ensembles in the toolbox layout
                          (p1, p2, [model|year], nCells), melt kg m^-2 yr^-1 positive = melt
                          (MALI: melt = -floatingBasalMassBal * s_yr)
 (6) calibrate            unstructured-mesh J1..J4 + objective function, 100,000 samples
 (7) report               distribution plot, per-term diagnostics (protocol Fig. 5/7),
                          chosen 5th/50th/95th values, per-basin shelf area, optional deltaT
```

Steps 5–7 are minutes of CPU; step 4 is small; step 3 is the only one needing real throughput
(ESMF weight generation, once per mesh).

### 3.3 MALI run configuration for step 4

```
config_velocity_solver          = 'none'      ! melt does not depend on velocity
config_basal_mass_bal_float     = 'ismip6'    ! or the new Burgard option, per Q3
config_thermal_solver           = 'none'
config_thermal_calculate_bmb    = .false.
config_dt                       = <one short step>
config_run_duration             = <one step>
config_ocean_data_extrapolation = ?           ! Q7
```

Streams: `input` = prepared mesh file; an initial-only TF stream with
`ismip6shelfMelt_3dThermalForcing` + `ismip6shelfMelt_zOcean`; an `ismip7_params` stream with
`gamma0`/`deltaT`/`basin`; an output stream at the end of the step with
`floatingBasalMassBal`, `ismip6shelfMelt_TFdraft`, `cellMask`, `areaCell`, `thickness`,
`lowerSurface`, `bedTopography`.

A diagnostic-only run will not work — MALI computes melt inside the timestep, not in the
initial diagnostic solve (findings 1.1).

### 3.4 Software design

A Python package + CLI in this directory. Mask generation goes upstream to MPAS-Tools instead,
so the generated files carry provenance (§5.1, Q4; §5.2, Q11).

```
parameterisations/mali/
  mali-melt-calibration-plan.md      # this document
  survey-findings.md                 # evidence
  protocol-and-toolbox-questions.md  # focus-group feedback log
  pyproject.toml / pixi.toml
  mali_melt_calib/
    config/                          # .cfg defaults and per-run overrides
    inputs.py                        # dataset registry, fetch, checksums, provenance
    mesh.py                          # thin driver over the MPAS-Tools mask tool
    forcing.py                       # 8 km TF -> MALI mesh remap
    ensemble.py                      # generate + submit the MALI runs
    collect.py                       # MALI output -> toolbox ensemble datasets
    terms.py                         # unstructured-mesh J1..J4  <-- new, upstreamable
    calibrate.py                     # objective function driver
    report.py                        # plots + chosen parameter values
    __main__.py                      # `mali-melt-calib <step>` CLI
  tests/
```

`terms.py` takes a cell-area array rather than assuming `reso**2`, so the same functions serve
structured grids and unstructured meshes and can be offered upstream.

### 3.5 Feedback to the focus group

We are the first unstructured-mesh model through this protocol, so we will keep hitting points
where the manuscript or toolbox is ambiguous. These are logged as they arise in
[`protocol-and-toolbox-questions.md`](protocol-and-toolbox-questions.md) — Part A for the
manuscript (to pass to Ronja), Part B for concrete repository changes (to fold into our
upstream PR). Entries record what we checked and concluded, so the focus group can tell a real
problem from a documentation gap.

---

## 4. What we still need to acquire

Almost nothing — the calibration datasets, the MALI mesh, the masks and the Globus tooling are
all already on LCRC (findings 1.4–1.6).

1. **`ISMIP7_8km-60m_AIS_grid_ocean_v3.nc`** (46 MB, Globus) — *probably not required.* ESMF
   weight generation needs only x/y coordinates, which
   `grid_and_mapping.create_ismip7_grid_file('AIS', 8, ...)` writes directly, and
   `i7aof.grid.ismip` is a third route. Worth fetching once to confirm all three agree.
2. **The relaxed initial condition** `relaxed_10yrs_4km.nc` (NERSC) — only for the §3.1
   sensitivity run, not the primary path.
3. **Nothing new** for the Burgard form if Q3 goes ahead: salinity is already present for
   every ocean state (`*_S.nc`), and latitude and ice-draft slope are derivable on the mesh.

---

## 5. Decisions

### 5.1 Settled

| # | Decision |
|---|---|
| **Q4** | **Basin numbering.** ISMIP7 is 0-based, MALI 1-based; MALI basin 10 = ISMIP7 basin 9 (findings F2). Derive the basin field by remapping the ISMIP7 mask so numbering is right by construction, and keep the cross-tabulation as a regression test. The remapping tool goes in MPAS-Tools/Compass for provenance, not in this package. |
| **Q5** | **Basin ΔT.** Calibrate with ΔT = 0, then optionally fit ΔT_b to J1 afterwards. Matches the production parameter file (ΔT = 0 everywhere) and preserves the linearity in findings F1. |
| **Q6** | **Dataset selection.** J3 with Mathiot_NEMO + Naughten_FESOM_ACCESS (all basins) and Jourdain-Naughten + Naughten_MITamu-MITwed (basins 9 and 14 only); J4 with PIG 2009 and 2012 only. |
| **Q8** | **Machine.** Chrysalis — the mesh does not need moving (findings 1.4). |
| **Q9** | **Tool location.** `parameterisations/mali/`, this directory. |
| **Q10** | **Unstructured meshes.** Contribute the MALI code upstream as another worked example, including the mesh-agnostic J1–J4 implementation. |
| **Q12** | **Branch name.** `add-mali-melt-calibration`, matching this repository's branch style. |
| — | **Geometry.** Calibrate on the un-relaxed mesh; see §3.1. |

### 5.2 Open

| # | Question | Blocks | For |
|---|---|---|---|
| **Q3** | Implement the Burgard et al. (2022) quadratic in MALI? | Phase 5 | Matt / Trevor |
| **Q11** | Where does the ISMIP7 mask tool live in MPAS-Tools? | Phase 3 upstream | Matt / Trevor |
| **Q7** | Is MALI's ocean extrapolation needed after remapping? | Phase 4 | Matt / Trevor |
| **Q1** | Which MALI branch is authoritative for ISMIP7? | Phase 5 | Matt / Trevor |
| **Q13** | Does a plain (non-Albany) MALI build suffice? | Phase 5 | Matt / Trevor |
| **Q2** | Which initial-condition vintage? | Phase 3 re-run | Matt / Trevor |
| **A3, A4** | Manuscript items still unresolved | — | Focus group |

**Q3 — Burgard et al. (2022) quadratic.** ISMIP7 explicitly recommends it; MALI has only the
ISMIP6 form (findings 1.1). Options: (a) calibrate `gamma0` in the existing form — defensible,
no code change, but not the recommended formulation; (b) add an option implementing Eqs.
(1)/(2) with parameter `K`; (c) (b) plus reporting the ISMIP6-form `gamma0` for comparison.
Xylar leans towards (b) or (c) and is checking with colleagues. **This is the highest-priority
answer** — it shapes the ensemble. If it goes ahead: who writes it; which variants (local vs.
semi-local, constant vs. local slope); a new `config_basal_mass_bal_float` value or a
sub-option; the 3-D **salinity** input stream Eq. (1) needs, which compass's ocean step does
not currently remap; and any date by which the melt module must be frozen.

**Q11 — MPAS-Tools mask tool.** Proposal: a new `landice/mesh_tools_li/`
`interpolate_ismip7_masks_to_mali.py` writing the ISMIP7-numbered basin field, BFRN bins,
floating mask and PIG/Dotson mask onto a MALI mesh, reusing `grid_and_mapping.py` with
`build_mapping_file` generalised so the ISMIP grid can be the source. Should that helper move
out of `output_processing_li/`, given an input-preparation tool would now call it? And should
`tune_ismip6_melt_deltat.py` gain ISMIP7's three-product J1 target (findings 1.3), or should
we do the ΔT fit here and leave it alone?

**Q7 — Ocean extrapolation.** ISMIP7 TF is already extrapolated into cavities on the 8 km
grid. After remapping to the MALI mesh, is `config_ocean_data_extrapolation` needed, or is the
remapped field already gap-free? What does the production setup do?

**Q1 — MALI branch.** `MALI-Dev/E3SM` @ `develop` has the ISMIP6 quadratic and recent
ISMIP7-adjacent merges. Is there a `matthewhoffman/*` or `trhille/*` branch with ISMIP7 melt
changes not yet on `develop`?

**Q13 — Albany.** With `config_velocity_solver = 'none'` the runs should not need Albany.
Confirm, and point us at a current plain-MALI build on Chrysalis if one exists.

**Q2 — IC vintage.** All five vintages are geometrically identical (findings F3), so this is
about matching the projections: which vintage do they use, do they add a relaxation step, and
is a BedMachine v3 / Bedmap3-based successor mesh in progress?

---

## 6. Phased work plan

**Nothing blocks starting.** Phases 0–2 need no answers; phase 3 can start on
`ais_4to20km.20250625.nc` and switch vintage later as a config change.

| Phase | Work | Blocked by | Output |
|---|---|---|---|
| **0. Feedback log** | append to `protocol-and-toolbox-questions.md` throughout; do not defer to the end | — | notes for Ronja; items for the PR |
| **1. Scaffold** | pixi env; package skeleton + CLI; config system; dataset registry with checksums | — | `mali-melt-calib inputs` runs green |
| **2. Terms on unstructured meshes** | area-weighted `calculate_term1..4`; tests reproducing the structured-grid answers on a uniform-area mesh; end-to-end replication of the published quadratic numbers (median K = 8.5e-5, 5th = 4.75e-5, 95th = 13.75e-5) as a regression test | — | `terms.py` + tests |
| **3. Mesh preparation** | `interpolate_ismip7_masks_to_mali.py` in MPAS-Tools reusing `grid_and_mapping.py`; assemble the mesh file from the un-relaxed vintage, asserting the F2 basin mapping | Q11 (review), Q2 (vintage only) | MPAS-Tools PR + `mali-melt-calib mesh` |
| **4. Forcing remap** | 8 km → MALI-mesh remap of 11 (then 26) TF fields, plus `so` if Q3 → (b)/(c) | Q7 | `mali-melt-calib forcing` |
| **5. MALI ensemble** | run directories, namelists/streams, job scripts; 3-value linearity check; production runs | **Q3**, Q1, Q13 | melt fields |
| **6. Calibration + report** | assemble ensembles; 100,000-sample optimisation; protocol Fig. 5/7 equivalents; per-basin shelf area; relaxed-IC sensitivity; optional ΔT | — | parameter values + plots |
| **7. Upstream PR** | the MALI example, the mesh-agnostic terms, and the Part B toolbox fixes | — | PR to `ismip7-antarctic-ocean-forcing` |

**Start with phase 2.** It needs no MALI runs and no open answers, and reproducing the
published quadratic numbers through our own code path validates the whole of stage 2 while
giving us the regression test the unstructured generalisation must not break. It should also
settle feedback item A3 empirically: the published percentiles came out of the toolbox as it
is, so matching them confirms which sampling the published numbers actually used.

Write phase 4 to remap a *list* of variables rather than hardcoding `tf`, so that a late "yes"
on Q3 adds salinity by configuration rather than by rework.

---

## 7. Risks

| Risk | Impact | Mitigation |
|---|---|---|
| Basin numbering used inconsistently | Silently wrong J1/J3/J4 — worst failure mode here, since results still look plausible | Confirmed real (findings F2). Derive basins by remapping the ISMIP7 mask; assert the cross-tabulation in the tool |
| Q3 answered late | Phase 5 rework, possibly phase 4 | Design phase 4 for a variable list; the toolbox is parameter-agnostic so phases 2 and 6 are unaffected either way |
| Linearity assumption wrong | Ensemble under-sampled | Explicit 3-value numerical check before relying on it |
| TF gaps after remapping | Spurious zero/invalid melt near grounding lines | Compare MALI-mesh `TFdraft` against 8 km TF at draft; check `config_invalid_value_TF`; enable MALI extrapolation if needed |
| Ice-shelf area mismatch biases J1/J2/J4 | Systematically shifted parameter distribution | Intended by the protocol; bounded by calibrating on un-relaxed present-day geometry and reported as per-basin modelled vs. observed area |
| Wrong IC vintage | Calibration against a geometry the projections do not use | Vintages are geometrically identical (findings F3); keep the IC a config option and re-run phases 3–6, which are cheap |
| Focus-group data revisions (a `_v4` appears) | Rework | Dataset registry pins versions and checksums; a refresh is a config change |
