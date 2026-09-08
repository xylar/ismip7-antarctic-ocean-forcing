# Plan: calibrating the MALI sub-shelf melt parameterization for ISMIP7

**Branch:** `ismip7-antarctic-ocean-forcing` @ `add-mali-melt-calibration`
**Status:** active. All phases are unblocked — the questions that remain are being carried
as working assumptions (§5.2) so that coding and testing can proceed through phase 5.
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
config_basal_mass_bal_float     = 'ismip6'    ! and the new Burgard option, per Q3
config_thermal_solver           = 'none'
config_thermal_calculate_bmb    = .false.
config_dt                       = <one short step>
config_run_duration             = <one step>
config_ocean_data_extrapolation = .false.     ! Q7 assumption: TF is pre-extrapolated
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

### 3.5 Possible direction change: porting the MALI-specific parts to Compass

Noted 2026-09-08, **nothing to do now.** The rest of the MALI ISMIP7 work — `ismip7_forcing`
and `ismip7_run` — lives in Compass, so a self-contained package here is the odd one out.
Porting the MALI-specific orchestration into Compass would be more intuitive for the other
MALI developers, who already know where to look for it.

If that happens, the split would be roughly:

* **moves to Compass** — the pieces that mirror what is already there: preparing the mesh
  file, remapping the calibration forcing (a variant of `ismip7_forcing/ocean_thermal`), and
  generating and running the single-timestep ensemble (a variant of `ismip7_run`).
* **stays here** — `terms.py`, which is mesh-agnostic and belongs with the toolbox it
  generalises; `quadratic.py`, the reference implementation; and the worked example and its
  documentation.

That last point is why this is a possible move rather than a mistake to correct: the
`parameterisations/mali/` example serves ice-sheet modellers who are not MALI developers and
will not have Compass, sitting alongside the quadratic, PICO and LADDIE examples. It earns
its place regardless of where the MALI plumbing ends up.

Deferring the decision costs little — the modules are already separated along roughly that
line, and `datasets.py` keeps the ocean-state definitions in one place either way.

### 3.6 Feedback to the focus group

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

### 5.2 Working assumptions — proceeding on a best guess, to be confirmed

Xylar's direction (2026-09-08): rather than wait, code and test through phase 5 on our best
guess for each. Each row names what we assume and what it costs if the answer differs.

| # | Assumption | If wrong |
|---|---|---|
| **Q1** | `MALI-Dev/E3SM` @ `develop` is authoritative. Work on branch `mali/add-burgard-melt-param` | Rebase onto the right branch |
| **Q2** | Calibrate on `ais_4to20km.20250625.nc`, the newest vintage | Re-run phases 3–6; all vintages are geometrically identical (findings F3), so the risk is near zero |
| **Q3** | **Implement it.** Add the Burgard et al. (2022) *local* quadratic (Eq. 1) with a constant Antarctic-mean slope as a new `config_basal_mass_bal_float` option; calibrate `K` and also report the ISMIP6-form `gamma0` | Wasted MALI development, but the calibration pipeline is unchanged — the toolbox is parameter-agnostic, so we fall back to calibrating `gamma0` |
| **Q7** | `config_ocean_data_extrapolation = .false.`. **Supported by evidence:** after bilinear remap to the MALI mesh the thermal forcing is 100% finite over all 385,379 cells, so there are no gaps for MALI to fill | Re-run phase 4 with extrapolation on; cheap |
| **Q11** | pyremap does the remapping (§5.4); the driver script goes in `MPAS-Tools/landice`, branch `add-ismip7-mali-masks` | Move the script to Compass |
| **Q13** | Build MALI **without** Albany, since `config_velocity_solver = 'none'` | Use the existing Albany build (below) |
| **A3** | ~~Assumption~~ **confirmed**: continuous `U(0,1)` weights reproduce all three published percentiles exactly (phase 2). The manuscript text is what is out of date | n/a |
| **A4** | Constant slope, computed on our own mesh per the notebook's advice that it "should match the geometry". On Bedmap3 at 8 km this gives sinθ = 0.0051117, matching the paper's 0.005. Report the value explicitly | Re-run with the published constant; `K` is not comparable across slope conventions |

Verified while setting these up: Chrysalis + gnu + openmpi is an Albany-supported combination
(`compass/deploy/albany_supported.txt`), and `/lcrc/soft/climate/compass/chrysalis/spack/`
already holds both `dev_compass_1_2_0_gnu_openmpi` and `..._gnu_openmpi_albany`. So the Q13
fallback costs nothing if we need `FO` later.

**Still to ask, in priority order:** Q3 (do the MALI developers agree?), Q11 (script name
and home), Q7, Q2, Q13, Q1. Plus A4 and A7 to the focus group; A3 is now answered by our
replication. All are recorded in
[`protocol-and-toolbox-questions.md`](protocol-and-toolbox-questions.md).

### 5.3 Q3 in detail — what we are implementing

Protocol Eq. (1), the "quadratic local" form:

```
m = K sinθ (ρ_o/ρ_i) (c_o/L_i)² β_S S_loc (g / (2|f|)) |TF_loc| TF_loc
```

Scope for the first pass:

* **local** (Eq. 1), not semi-local (Eq. 2) — the protocol accepts either, and local matches
  the worked example in this repository.
* **constant** Antarctic-mean slope, not the locally-varying one. Slope enters as a scalar so
  the linearity in findings F1 is preserved.
* a new value of `config_basal_mass_bal_float` rather than a sub-option of `'ismip6'`, so the
  existing ISMIP6 path is untouched and both can be run for comparison.
* new registry fields for `K` and the constant slope; `S_loc` needs a **3-D salinity input
  stream** alongside the existing TF one, which is the largest piece of new plumbing —
  compass's ocean step currently remaps only `tf`.
* `f` **constant** at 1.4e-4, matching the reference implementation, rather than from
  `latCell`. A latitude-varying `f` would shift `K` away from the published value; raised as
  feedback A7.

**A local slope was tried in MALI before and rejected.** `MALI-Dev/E3SM` PR #48 added a
per-cell `shelfBaseSlope` field to the ISMIP6 quadratic and was closed unmerged: it melted
holes in shelves as TF rose, through a feedback between basal slope and melt rate. Our
constant scalar slope cannot reproduce that — it is mathematically a rescaling of `K` — so
the finding supports the scoping above rather than contradicting it. See findings §1.8. It
does mean the locally-varying slope variant of Eq. (1) should not be added without
revisiting that PR.

**Algorithm from the protocol, constants from MALI.** The formula is protocol Eq. (1), but
`cp_seawater`, `latent_heat_ice`, `gravity` and the ice and ocean densities are MALI's own,
from `li_constants` and the namelist — not the reference implementation's. They differ
slightly (that implementation uses `L_i = 334` kJ/kg against MALI's 335, and `g = 9.81`
against 9.80616), so at the same `K` melt differs by a few tenths of a percent. That is far
inside the parametric spread the calibration explores — the 5th-to-95th range spans a factor
of about three — and `K` absorbs any constant offset regardless. Agreement with the reference
implementation is therefore a **sanity check, not a requirement**.

The one new physical constant, the haline contraction coefficient `beta_S = 7.86e-4` PSU⁻¹,
is added to `li_constants` where MALI's other physical constants live.

`K` and the slope defaults are the published 8 km values, flagged **provisional** in the
Registry: both were obtained with the reference implementation's constants on a different
grid, and `K` absorbs the slope, the Coriolis parameter and the constants, so production runs
should use a value calibrated with MALI on the mesh being run — which is what the rest of
this plan produces.

### 5.4 Q11 in detail — the mask tool

Use **pyremap**, not a hand-rolled ESMF wrapper. Generalising
`grid_and_mapping.build_mapping_file` is dropped: it hardcodes MALI as the source
(`-s mali_scripfile -d ismip7_scripfile`) and pyremap removes the need for it.

```python
remapper = Remapper(ntasks=..., map_filename=..., method='neareststod')
remapper.src_from_proj(ismip_grid_file, 'ismip8km', proj_str=EPSG3031)
remapper.dst_from_mpas(mali_mesh_file, 'ais4to20km')
remapper.build_map(logger=logger)
masks_on_mali = remapper.remap_numpy(masks)
```

pyremap is already an `i7aof` dependency, so this adds nothing new. Verified that all four
ISMIP7 mask files carry real projected coordinates (−3040 km…+3040 km at 8 km, 761 points),
which is what `src_from_proj` needs. Nearest-neighbour (`neareststod`) throughout, since
every field is categorical.

The script should also emit the cross-tabulation against the existing ISMIP6
`regionCellMasks`, so the off-by-one in findings F2 is asserted rather than rediscovered.

Portability nit to raise in the PR: `build_mapping_file` decides whether to use `srun` from
`hostname.startswith('nid')`, which is Cori/Perlmutter-specific and would not fire on
Chrysalis.

---

## 6. Phased work plan

**Nothing is blocked.** The remaining questions are carried as working assumptions (§5.2),
each cheap to revise. Repositories are already set up:

```
add-mali-melt-calibration/                    ismip7-antarctic-ocean-forcing
mali/add-burgard-melt-param/                  E3SM (MALI-Dev)
~/mpas_work/MPAS-Tools/add-ismip7-mali-masks  MPAS-Tools
```

| Phase | Work | Blocked by | Output |
|---|---|---|---|
| **0. Feedback log** | append to `protocol-and-toolbox-questions.md` throughout; do not defer to the end | — | notes for Ronja; items for the PR |
| **1. Scaffold** ✅ | pixi env; package skeleton + CLI | — | `pixi.toml`, `mali_melt_calib/` |
| **2. Terms on unstructured meshes** ✅ | area-weighted `calculate_term1..4`; tests reproducing the structured-grid answers on a uniform-area mesh; **replication of the published quadratic numbers reproduces all three percentiles exactly** | — | `terms.py`, `quadratic.py`, `replicate.py`, 16 tests |
| **3. Mesh preparation** ✅ | `interpolate_ismip7_masks_to_mali.py` in `MPAS-Tools/landice/mesh_tools_li`, using pyremap; basins, BFRN bins, floating mask and PIG/Dotson regions on the mesh, with a bidirectional cross-check against `regionCellMasks` | — | MPAS-Tools `f7407bc6`; `work/mesh/ais_4to20km_ismip7_masks.nc` |
| **4. Forcing remap** ✅ | `forcing.py` remaps TF **and salinity** for the 11 recommended states (26 available) via pyremap, in compass's field names and dimension order; `datasets.py` is the shared ocean-state registry | — | `work/forcing/ocean_forcing_*.nc` |
| **5. MALI melt module** | implement the Burgard local quadratic in MALI (§5.3) on `mali/add-burgard-melt-param`; build without Albany; verify against `quadratic.py` on the same drafts | — | MALI branch |
| **6. MALI ensemble** | run directories, namelists/streams, job scripts; 3-value linearity check; production runs | — | melt fields |
| **7. Calibration + report** | assemble ensembles; 100,000-sample optimisation; protocol Fig. 5/7 equivalents; per-basin shelf area; relaxed-IC sensitivity; optional ΔT | — | parameter values + plots |
| **8. Upstream PRs** | the MALI example and mesh-agnostic terms here; the mask tool to MPAS-Tools; the melt module to MALI-Dev | — | three PRs |

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
| Q3 turns out to be unwanted | The MALI melt-module work is wasted | Contained to its own branch; the calibration pipeline is unchanged, since the toolbox is parameter-agnostic and we fall back to calibrating `gamma0` |
| Linearity assumption wrong | Ensemble under-sampled | Explicit 3-value numerical check before relying on it |
| TF gaps after remapping | Spurious zero/invalid melt near grounding lines | Compare MALI-mesh `TFdraft` against 8 km TF at draft; check `config_invalid_value_TF`; enable MALI extrapolation if needed |
| Ice-shelf area mismatch biases J1/J2/J4 | Systematically shifted parameter distribution | Intended by the protocol; bounded by calibrating on un-relaxed present-day geometry and reported as per-basin modelled vs. observed area |
| Wrong IC vintage | Calibration against a geometry the projections do not use | Vintages are geometrically identical (findings F3); keep the IC a config option and re-run phases 3–6, which are cheap |
| Focus-group data revisions (a `_v4` appears) | Rework | Dataset registry pins versions and checksums; a refresh is a config change |
