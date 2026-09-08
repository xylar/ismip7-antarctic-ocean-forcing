# Plan: calibrating the MALI sub-shelf melt parameterization for ISMIP7

**Branch:** `ismip7-antarctic-ocean-forcing` @ `add-mali-melt-calibration`
**Status:** draft for review — several decisions still need input from Matt Hoffman /
Trevor Hillebrand / the AIS Ocean Focus Group (see §5).
**Date:** 2026-09-08

*Path conventions used below:* bare relative paths such as
`parameterisations/parameter_selection_toolbox.py` are relative to the root of this
repository. Paths beginning with `/` are absolute on LCRC. MALI source files are given
relative to an `E3SM` checkout (`MALI-Dev/E3SM` @ `develop`); the working clone for this
effort is `/lcrc/group/e3sm/ac.xylar/ismip7/MALI-melt-calibration/MALI-Dev`, alongside which
the survey below was carried out.

---

## 1. Goal

Follow the ISMIP7 AIS ice–ocean protocol (Reese et al., `ISMIP7_AIS_ice_ocean_protocol.pdf`,
Sect. 4.2) to calibrate the free parameter(s) of MALI's sub-shelf melt module, and produce
**at least three parameter values** (5th / 50th / 95th percentile of the optimal-parameter
distribution) for use in the ISMIP7 MALI projections.

Secondary goals:

* Full provenance: every input dataset is downloaded/derived by scripted, re-runnable steps.
* Tools committed to the branch as an installable **Python package with a CLI**, not notebooks.
* Avoid a hard dependency on Compass; use its `landice/ismip7_forcing` and
  `landice/ismip7_run` test groups as reference implementations.

---

## 2. How the calibration actually works

The protocol (§4.2) is a two-stage procedure. Understanding the split matters, because
stage 1 is the expensive/MALI-specific part and stage 2 is cheap and already written.

### Stage 1 — build the melt ensemble (this is the work)

Run MALI with **fixed, present-day geometry** for **a single time step**, once for every
combination of:

* **ocean forcing state** — one of ~26 thermal-forcing (TF) fields provided by ISMIP7, and
* **melt parameter value** — a sampling of the module's free parameter.

The ocean states, and which objective-function term each feeds (protocol Table 2):

| Group | Datasets | Term | Count |
|---|---|---|---|
| Present-day climatology | Zhou et al. observational climatology | J1, J2 | 1 |
| Ocean models, circum-Antarctic | Mathiot_NEMO, Naughten_FESOM_ACCESS (**minimum**), Naughten_FESOM_MMM, Timmermann_FESOM | J3 | 8 (cold+warm) |
| Ocean models, regional | Jourdain-Naughten_NEMO-MITgcm, Naughten_MITamu-MITwed (**strongly suggested**), Haid_FESOM | J3 | 6 (cold+warm) |
| Amundsen observations | PIG/Dotson years 1994–2020; **minimum = PIG 2009 + 2012** | J4 | up to 13 |

Minimal set = 1 + 4 + 2 = **7 ocean states**; recommended set (incl. the two regional
datasets the focus group added on 31 July 2026) = **11 ocean states**; everything = 26.

Output of stage 1: melt rate in **kg m⁻² yr⁻¹** (positive = melting) on the MALI mesh, for
each (ocean state, parameter value).

### Stage 2 — parameter optimisation (already implemented upstream)

`parameterisations/parameter_selection_toolbox.py` in this repo:

1. Aggregates modelled melt into four targets:
   * **J1** basin-integrated melt (Gt yr⁻¹) over the 16 IMBIE2 basins, vs. Paolo/Davison/Adusumilli.
   * **J2** melt integrated over 10 buttressing (BFRN) bins, weighted by median BFRN.
   * **J3** warm−cold basin-mean melt difference (kg m⁻² yr⁻¹) vs. the ocean models.
   * **J4** PIG/Dotson integrated melt (Gt yr⁻¹) per observation year.
2. Draws 100,000 samples of the term weights `a_i ~ U(0,1)` and of the target uncertainties
   (normal, per-basin/bin/model/year), minimising
   `I = Σ a_i · J_i / median(J_i)` each time.
3. Returns 100,000 optimal parameter values → take the **5th / 50th / 95th percentiles**
   (and the mode, for exploration).

### Key simplification: the ensemble over parameter values is (exactly) free

MALI's ISMIP6 quadratic is

```
BMB = -gamma0 · (ρ_sw c_p / (ρ_i L_i))² / s_yr · ρ_i · (TF_draft + ΔT_b) · |<TF>_b + ΔT_b| + offset
```

(`mpas_li_iceshelf_melt.F:1369-1374`). For a fixed geometry and a fixed ΔT field, `TF_draft`,
`<TF>_b`, `cellMask` and `connectedOceanMask` are all independent of `gamma0` — so **melt is
exactly linear in `gamma0`**. The same is true of the Burgard et al. (2022) form, which is
linear in `K`.

**Consequence:** we need only **one MALI run per ocean state** (≈11 runs), and the full
120-value parameter ensemble is obtained by scaling. I will *verify* this numerically by
running 3 values of `gamma0` for one ocean state and confirming exact linearity, rather than
assuming it. This turns a ~1,300-run campaign into ~15 runs.

This linearity breaks if we introduce basin-wide ΔT *inside* the optimisation. The protocol
offers two orderings (§4.2.1); **option (2) — compute ΔT after the parameter optimisation, as
in ISMIP6 — preserves linearity** and is what I propose (see §5, Q5).

---

## 3. What already exists (survey results)

### 3.1 MALI code

* The working clone is `MALI-Dev/E3SM` on `develop`, clean, with `matthewhoffman`, `trhille`
  and `xylar` remotes configured.
* `config_basal_mass_bal_float = 'ismip6'` →
  `iceshelf_melt_ismip6` in
  `components/mpas-albany-landice/src/mode_forward/mpas_li_iceshelf_melt.F:1167`.
  This is the **ISMIP6 *non-local* quadratic** (Jourdain et al. 2020), with:
  * `ismip6shelfMelt_gamma0` — scalar, m yr⁻¹
  * `ismip6shelfMelt_deltaT` — per-cell (nCells), K
  * `ismip6shelfMelt_basin` — per-cell integer basin index
  * `ismip6shelfMelt_offset` — per-cell, kg m⁻² s⁻¹
  * `ismip6shelfMelt_3dThermalForcing` (nISMIP6OceanLayers × nCells) + `ismip6shelfMelt_zOcean`
  * `ismip6shelfMelt_TFdraft` — diagnostic, TF interpolated to the ice draft
* **It is not the Burgard et al. (2022) update** that ISMIP7 recommends. The ISMIP7 cheat
  sheet is explicit: *"The only recommendation is to avoid linear formulations and to update
  the ISMIP6 formulation following Burgard et al. (2022)."* Protocol Eqs. (1)–(2) add
  `sinθ`, `β_S·S_loc`, `g/(2|f|)` and change `(ρ_o/ρ_i)` from squared to first power.
  → **This is the single biggest open item (§5, Q3).**
* `config_velocity_solver = 'none'` exists → melt-only runs are cheap.
* `li_basal_melt_floating_ice` is called *inside* the timestep
  (`mpas_li_time_integration_fe_rk.F:288`), **not** in the initial diagnostic solve, so a
  `config_stop_time == config_start_time` diagnostic-only run would produce no melt. We must
  take one short real timestep (or add a small MALI hook — not proposed).

### 3.2 Compass (`/home/ac.xylar/compass/main`, branch `main`)

Active and directly relevant:

* `compass/landice/tests/ismip7_forcing/ocean_thermal/process_thermal_forcing.py` — remaps 3D
  TF from the ISMIP 8 km × 30-level grid onto a MALI mesh, writing
  `ismip6shelfMelt_3dThermalForcing`, `ismip6shelfMelt_zOcean`, `ismip6shelfMelt_zBndsOcean`.
  **This is exactly the remapping we need**, just driven by calibration TF files instead of
  CMIP scenario files.
* `compass/landice/tests/ismip7_run/ismip7_ais/` — the ISMIP7 AIS production setup:
  `namelist.landice`, `streams.landice.template`, `set_up_experiment.py`. It tells us the
  intended production configuration:
  * `config_velocity_solver = 'FO'`, `config_basal_mass_bal_float = 'ismip6'`, `config_dt = 1 month`
  * `init_cond_path = .../ISMIP6-2300/initial_conditions/AIS_4to20km_20230105/relaxation_0TGmelt_10yr/relaxed_10yrs_4km.nc`
  * `melt_params_path = .../basin_and_coeff_DeltaT_quadratic_non_local_gamma14500.nc`
  * an `ismip7_params` input stream carrying `gamma0`, `deltaT`, `basin`
* Those specific compass paths are on NERSC CFS, but **the mesh and equivalent inputs are
  already on LCRC** in the E3SM inputdata tree,
  `/lcrc/group/e3sm/{data,public_html}/inputdata/glc/mpasli/mpas.ais4to20km/`:
  * `ais_4to20km.{20230105,20240708,20241224,20250411,20250625}.nc` — initial-condition
    vintages, `nCells = 385,379`, carrying `thickness`, `bedTopography`, `areaCell`,
    `xCell`/`yCell`/`latCell` (no `lowerSurface`/`cellMask` — those are derived at init)
  * `ais_4to20km.20230105.scrip.nc` — SCRIP description, ready for ESMF weight generation
  * `ais_4to20km_region_mask.20230105.nc` — `regionCellMasks(nCells, 16)`, named
    `ISMIP6 Basin A-Ap` … `K-A` (the standard ISMIP6 16-basin partition)
  * `ais_4to20km_tf_params.20250724.nc` — `ismip6shelfMelt_gamma0 = 14500`,
    **`ismip6shelfMelt_deltaT = 0` everywhere**, `ismip6shelfMelt_basin` in 0…16
  * `mpasli.graph.info.240507.part.{64…3840}` — graph partitions

  So phases 3–5 are **not** blocked on a NERSC transfer. What remains open is *which*
  vintage is the ISMIP7 initial condition, and whether the 10-year relaxation applied in the
  compass config is baked into any of these files or applied separately (Q2).

### 3.3 ISMIP7 data — almost everything is already on LCRC

Root: `/lcrc/group/e3sm/ac.xylar/ismip7/forcing-data/data/AIS/`

| Subtree | Contents | Status |
|---|---|---|
| `parameterisations/ocean/ocean_modelling_data/` | 59 files: `{Model}_{cold,warm}_{T,S,TF,m}.nc` + `melt_{cold,warm}_target_term3_v2.nc` | complete |
| `parameterisations/ocean/ocean_observations_data/` | 41 files: `Obs_{year}_{T,S,TF}.nc` + `melt_observations_target_term4.nc` | complete |
| `parameterisations/ocean/imbie2/` | `basin_numbers_ismip{1,2,4,8,16,32,64}km_v2.nc` | complete |
| `parameterisations/ocean/bfrns/` | `BFRN_ismip{...}km_v2.nc` | complete |
| `parameterisations/ocean/floatingmasks/` | `floatingmask_ismip{...}km.nc` | complete |
| `parameterisations/ocean/meltobs/` | `Melt_Paolo_Davison_Adusumilli_imbie2.csv`, `melt_target_term2_v3.nc` | complete |
| `parameterisations/ocean/shelfmask/` | `shelf_mask_ismip8km.nc` (PIG id 110, Dotson id 97) | complete |
| `obs/ocean/climatology/zhou_annual_06_nov/{tf,so,thetao}/` | present-day climatology | complete |
| `obs/ocean/topography/{bedmap3,BedMachineAntarctica-v3}/v3/` | 8 km topo: `bed`, `draft`, `surface`, `thickness`, `*_frac` masks | complete |
| `obs/ocean/IMBIE-basins/v3/` | basin polygons/masks | complete |
| `grid/ocean/ISMIP7/8km-60m/v3/` | `ISMIP7_8km-60m_AIS_grid_ocean_v3.nc` (46 MB) | **missing** |

All TF/S fields are on the ISMIP 8 km polar-stereographic grid, `761 × 761 × 30 levels`.
Diffed against the Globus manifest (`ismip7-source-coop/manifests/ISMIP7_AIS.json.gz`,
109,834 files / 3.5 TB): the only gaps in the subtrees we need are the one grid file above
and `parameterisations/fracture/ant_surf_roughness_*.nc` (not needed for melt calibration).

### 3.4 Globus access is already solved

`/lcrc/group/e3sm/ac.xylar/ismip7/forcing-data/ismip7-source-coop/` is a pixi project with
`globus-cli` and working scripts (`transfer_to_lcrc.py`, `globus_survey.py`,
`refresh_manifest.py`). Endpoint: **`ccc9bbd2-4091-4e35-addd-eeb639cf5332`** (`GHub_upload`),
tree root `/ISMIP7/`. I will reuse this rather than build new Globus tooling; the calibration
package will just declare the manifest paths it needs and shell out to the existing transfer
script (or `globus transfer`), with checksums recorded.

### 3.5 The toolbox assumes a structured grid — this is our main code gap

`parameter_selection_toolbox.py` hard-codes uniform-cell-area, 2-D-raster assumptions:

* `cvt_m = reso**2 / 1e12` — constant cell area (`calculate_term1`, `2`, `4`)
* `.groupby(basins_m).sum()` over a 2-D `(y, x)` field
* `obs_ensemble.stack(grid=('x','y'))` in `calculate_term4`

MALI's mesh is **unstructured and variable-resolution (4–20 km)**, so all four term
calculations need area-weighted, `nCells`-based equivalents. This is a contained, testable
piece of work and a genuinely useful upstream contribution (the toolbox README already says
"only regular grids are supported at the moment").

---

## 4. Proposed approach

### 4.1 The geometry question — recommendation: do both

The protocol pulls in two directions:

* §4.2 bullet 1: *"modellers are asked to use a recent present-day geometry, e.g., BedMap3 or
  BedMachine3, for the calibration rather than an ice-sheet model initial state with
  substantially different grounding line position, ice thickness, or ice shelf extent."*
* §4.2 bullet 2: *"Ideally, the melt module is calibrated and evaluated in the same setting as
  it is used in the core ISMIP7 ice-sheet simulations, i.e. with the same code and on the same
  grid."*

The MALI ISMIP7 initial condition (`relaxed_10yrs_4km.nc`) is a 10-year relaxation, so it is
*not* the present-day geometry. Because the melt runs are cheap (§2), I propose we do not
choose:

* **Config A (primary): MALI ISMIP7 AIS mesh + ISMIP7 initial condition.** Same mesh, same
  code, same masks as the projections. Satisfies bullet 2.
* **Config B (sensitivity): same mesh, geometry replaced by Bedmap3 interpolated onto it.**
  `thickness`, `bedTopography` (and hence `lowerSurface`, `cellMask`) taken directly from the
  8 km Bedmap3 product already on disk, no relaxation. Satisfies bullet 1.

The difference between the two calibrated parameter distributions *is* the geometry-induced
bias, and reporting it is far more useful than picking one and hoping. **A new mesh is almost
certainly not needed** — we reuse the production mesh and swap the geometry fields.

Caveat to watch: J1 and J4 are *integrated* melt (Gt yr⁻¹), so a mismatch in ice-shelf area
between MALI and observations biases them directly. J3 (basin-*mean* kg m⁻² yr⁻¹) is much less
sensitive. Config B largely removes this; Config A will show it.

### 4.2 Workflow

```
 (0) environment          pixi env for the calibration tools
 (1) inputs               link/verify local ISMIP7 data; globus-fetch the 8km grid file;
                          obtain the MALI ISMIP7 mesh + initial condition
 (2) mesh preparation     for each of Config A / B, build a MALI-mesh file carrying:
                            geometry (thickness, bedTopography, lowerSurface)
                            ismip6shelfMelt_basin   (IMBIE2, 16 basins)
                            BFRN bins               (10 bins, remapped from ISMIP 8 km)
                            floating mask + areaCell
                          + the PIG/Dotson region mask for J4
 (3) forcing remap        11 (or 26) ISMIP 8 km TF fields -> MALI mesh, 30 z-levels
                          (bilinear, mirroring compass ismip7_forcing/ocean_thermal)
 (4) MALI ensemble        one 1-timestep, velocity-solver-off MALI run per ocean state
                          (+ 3 gamma0 values on one state to verify linearity)
 (5) collect              assemble pd_ensemble / cold_ensemble / warm_ensemble /
                          obs_ensemble in the toolbox's (p1, p2, [model|year], nCells) layout,
                          melt in kg m^-2 yr^-1, positive = melting
                          (MALI: melt = -floatingBasalMassBal * s_yr)
 (6) calibrate            unstructured-mesh J1..J4 + objective function, 100,000 samples
 (7) report               distribution plot, per-term diagnostic plots (protocol Fig. 5/7),
                          chosen 5th/50th/95th values, optional basin deltaT
```

Steps 5–7 are minutes of CPU; step 4 is small; step 3 is the only step needing real
throughput (ESMF weight generation, once per mesh).

### 4.3 MALI run configuration for step 4

```
config_velocity_solver        = 'none'      ! melt does not depend on velocity
config_basal_mass_bal_float   = 'ismip6'    ! (or the new Burgard option, per Q3)
config_thermal_solver         = 'none'
config_thermal_calculate_bmb  = .false.
config_dt                     = <one short step>
config_run_duration           = <one step>
config_ocean_data_extrapolation = ?         ! Q7 — TF is pre-extrapolated on the 8km grid
```

Streams: `input` = prepared mesh file; an initial-only TF stream carrying
`ismip6shelfMelt_3dThermalForcing` + `ismip6shelfMelt_zOcean`; an `ismip7_params` stream with
`gamma0`/`deltaT`/`basin`; an output stream at the end of the step with
`floatingBasalMassBal`, `ismip6shelfMelt_TFdraft`, `cellMask`, `areaCell`, `thickness`,
`lowerSurface`, `bedTopography`.

### 4.4 Software design

A Python package + CLI, no notebooks. Proposed location inside the branch: **this directory**,
`parameterisations/mali/` (sibling of the PICO/LADDIE/quadratic examples), or `i7aof/mali/` if
we prefer it inside the installed package — **see §5, Q9**.

```
parameterisations/mali/
  mali-melt-calibration-plan.md    # this document
  pyproject.toml / pixi.toml       # pixi env (user preference)
  mali_melt_calib/
    config/                        # .cfg defaults, per-config (A/B) overrides
    inputs.py                      # dataset registry, globus fetch, checksums, provenance
    mesh.py                        # basins/BFRN/floating mask/region mask -> MALI mesh
    forcing.py                     # 8km TF -> MALI mesh remap (pyremap/ESMF)
    ensemble.py                    # generate + submit the MALI runs
    collect.py                     # MALI output -> toolbox ensemble datasets
    terms.py                       # unstructured-mesh J1..J4  <-- new, upstreamable
    calibrate.py                   # objective function driver (reuses upstream toolbox)
    report.py                      # plots + chosen parameter values
    __main__.py                    # `mali-melt-calib <step>` CLI
  tests/
```

`terms.py` should be written so the same functions serve both structured and unstructured
input (area array supplied explicitly rather than assumed), so it can be offered upstream to
replace the `cvt_m = reso**2` assumption.

### 4.5 Feedback to the AIS Ocean Forcing focus group

Being the first unstructured-mesh model through this protocol, we will keep hitting points
where the manuscript or the toolbox is ambiguous, inconsistent, or silent. These are worth
capturing as they arise rather than reconstructing later, and they are useful to the focus
group in two different ways. Both live in
[`protocol-and-toolbox-questions.md`](protocol-and-toolbox-questions.md), started alongside
this plan:

* **Part A — manuscript.** Points of confusion in the ISMIP7 AIS ice–ocean protocol paper:
  ambiguities, internal inconsistencies, and places where guidance is missing for an ISM
  whose grid or initialisation differs from the 8 km structured assumption. Collected as a
  standalone markdown document to pass back to Ronja later. It is *feedback*, not a change
  request — some items will turn out to be our misreading, and that is itself useful signal
  about where the text could be clearer.
* **Part B — repository.** Points of confusion in `parameterisations/` — the toolbox, its
  README, and the worked examples — where the fix is a concrete change. These are proposed
  as part of our PR adding the MALI example, so each item should say what the change would
  be, not just what confused us.

Every entry should record what we checked and what we concluded, so the focus group can see
whether it is a real problem or a documentation gap. Where we found a genuine
manuscript/code discrepancy, both sides are cited.

Both parts are seeded with what the initial survey already turned up (five manuscript items,
five repository items), including one confirmed manuscript/code discrepancy in how the J3/J4
inclusion pre-factors are sampled.

---

## 5. Decisions

### 5.0 Resolved (2026-09-08)

**Q5 — Basin ΔT: confirmed.** Calibrate `gamma0` (or `K`) with ΔT = 0, then optionally fit
ΔT_b to J1 afterwards. This is the expected ISMIP7 workflow for MALI, is consistent with the
production `ais_4to20km_tf_params.20250724.nc` (ΔT = 0 everywhere), and preserves the exact
linearity in §2.

**Q6 — Dataset selection: confirmed.** J3 with Mathiot_NEMO + Naughten_FESOM_ACCESS (all
basins) and Jourdain-Naughten + Naughten_MITamu-MITwed (basins 9 and 14 only); J4 with PIG
2009 and 2012 only.

**Q9 — Tool location: confirmed** `parameterisations/mali/` — this directory.

**Q12 — Branch name: resolved.** Development moved to `add-mali-melt-calibration`, matching
this repository's branch-naming style (`fix-date-line`, `update-deps`, `replace-fortran`) and
the convention that each worktree directory is named for its branch. The earlier
`mali-melt-calibration-notebooks` branch remains at the initial plan commit and is no longer
used.

**Q10 — Unstructured meshes: confirmed.** The MALI code will be contributed upstream as
another worked example alongside the quadratic, PICO and LADDIE ones, including the
area-weighted, mesh-agnostic J1–J4 implementation.

**Q4 — Basin index mapping: resolved by inspection; still needs a provenance-carrying tool.**
Verified directly against the data:

* ISMIP7 `basin_numbers_ismip8km_v2.nc` is **0-based, 0…15**, covering every grid cell (no
  "outside" value). Confirmed geographically: PIG and Dotson both fall in basin **9**
  (Eastern Amundsen) and basin **14** is Ronne-Filchner, matching the protocol's usage.
* MALI `ismip6shelfMelt_basin` is **1-based, 1…16** (plus a single stray cell labelled 0),
  equal to `regionCellMasks` column index + 1.
* Cross-tabulating the two over the 296,538 ice cells of `ais_4to20km.20250625.nc`:
  **MALI = ISMIP7 + 1 for 100.0% of ice cells**, with per-basin purity ≥ 98.8% for 15 of 16
  basins. The exception is MALI 14 / ISMIP7 13 (Antarctic Peninsula, `Ipp-J`) at 85.2%,
  where the ISMIP6-2022 and IMBIE2-v3 boundaries genuinely differ.

  So MALI basin 10 = ISMIP7 basin 9 (Eastern Amundsen) and MALI basin 15 = ISMIP7 basin 14
  (Ronne-Filchner). **Using MALI's index where an ISMIP7 basin number is expected is off by
  one** — precisely the silent-failure mode flagged in §8.

**Decision:** derive the basin field on the MALI mesh by remapping the ISMIP7 mask, so the
numbering is correct by construction, and keep the cross-tabulation above as a regression
test. Per Xylar, that remapping tool belongs in **MPAS-Tools or Compass**, not in this
package, so that the generated basin file carries proper provenance; this package will call
it. Which of the two, and its interface, is a new question for the MALI/Compass developers
(Q11).

**Q8 — Machine: effectively settled.** The mesh does not need moving (§3.2), so the
calibration runs on Chrysalis. Remaining sub-question folded into Q11 below: since the
calibration runs set `config_velocity_solver = 'none'`, they should not need Albany at all —
only `'L1L2'`, `'FO'` and `'Stokes'` require external dycores — so a plain MALI build should
suffice. To be confirmed at build time.

### 5.1 Still open — for Matt Hoffman / Trevor Hillebrand (MALI + Compass)

**Q1 — Which MALI branch is authoritative for ISMIP7?** `MALI-Dev/E3SM` @ `develop` has the ISMIP6
quadratic and recent ISMIP7-adjacent merges (`ismip-flux-time-avg`, fracture work). Is there a
`matthewhoffman/*` or `trhille/*` branch with ISMIP7 melt changes not yet on `develop`?

**Q2 — Which initial-condition vintage?** The `ais_4to20km` mesh and five IC vintages
(20230105 → 20250625) are already on LCRC in the E3SM inputdata tree, along with the SCRIP
file, graph partitions, region mask, and a `tf_params` file. Compass instead points at a
NERSC file, `AIS_4to20km_20230105/relaxation_0TGmelt_10yr/relaxed_10yrs_4km.nc` — the
20230105 mesh after a 10-year relaxation. Which vintage is the ISMIP7 AIS initial condition,
and is the relaxation part of it? Should I reproduce the relaxation on LCRC or fetch the
relaxed file? Is a BedMachine v3 / Bedmap3-based successor mesh in progress?

**Q3 — Do we implement the Burgard et al. (2022) quadratic in MALI?** ISMIP7 explicitly
recommends it; MALI currently has only the ISMIP6 form. Options:
  (a) calibrate `gamma0` in the existing ISMIP6 non-local form (no code change; the toolbox is
      parameter-agnostic, so this is fully defensible, just not the recommended formulation);
  (b) add an option to `config_basal_mass_bal_float` implementing Eqs. (1)/(2) with
      parameter `K`, local vs. semi-local, constant vs. local slope;
  (c) (b) but calibrate `K` and *also* report the ISMIP6-form `gamma0` for comparison.

Xylar's leaning is that the development should be part of this plan — i.e. (b) or (c) — but
this is **not yet decided**; he is checking with colleagues. This decision gates the ensemble
design, so it is the highest-priority answer needed.

If it does go ahead, the follow-on questions are:
  * Who writes the MALI code — us on a branch, or Matt/Trevor?
  * Which variants to support: local (Eq. 1) vs. semi-local (Eq. 2); constant Antarctic-mean
    slope vs. local slope updated each step. The protocol accepts any; supporting local +
    constant slope is the minimum, and matches the worked example in this repo.
  * Naming: a new `config_basal_mass_bal_float = 'ismip7'` option, or a sub-option of
    `'ismip6'`? New registry fields are needed for `K` and the slope.
  * Salinity: Eq. (1) needs `S_loc` at the draft, so MALI needs a 3-D salinity input stream
    alongside the existing TF one. The ISMIP7 datasets provide `*_S.nc` for every ocean state,
    but the compass `ismip7_forcing` ocean step currently remaps only `tf`.
  * Timeline relative to the ISMIP7 projections — is there a date by which the melt module
    must be frozen?

**Q11 — Where should the ISMIP7 basin-mask tool live, and does a plain MALI build suffice?**
Two follow-ons from §5.0:
  * The tool that rasterises/remaps the ISMIP7 IMBIE2 basin mask (and the BFRN bins, floating
    mask and PIG/Dotson region mask) onto a MALI mesh should live in **MPAS-Tools or
    Compass** for provenance. Which, and is there an existing entry point to extend —
    `mpas_tools.landice`? a step in `compass/landice/tests/ismip7_forcing/`? The ISMIP6 basin
    mask on this mesh was made with a bespoke `compass run custom` in 2022 (per the file's
    history attribute), so there may be no reusable path yet.
  * Do the calibration runs need an Albany-linked MALI build? With
    `config_velocity_solver = 'none'` they should not, which would make the ensemble much
    easier to build and run on Chrysalis. Is that right, and is there a current plain-MALI
    build on Chrysalis to use?

**Q7 — Ocean data extrapolation.** ISMIP7 TF is already extrapolated into cavities on the 8 km
grid. After bilinear remap to the MALI mesh, do we need MALI's
`config_ocean_data_extrapolation`, or does compass's remap already produce a
gap-free `ismip6shelfMelt_3dThermalForcing`? What does the production ISMIP7 setup do?

### 5.2 Still open — for the ISMIP7 AIS Ocean Focus Group

Collected as we go in
[`protocol-and-toolbox-questions.md`](protocol-and-toolbox-questions.md) rather than here —
see §4.5. That document has two parts: manuscript points to pass back to Ronja, and toolbox
points to fold into our upstream PR.

---

## 6. Datasets: have / need

### Already on LCRC — no action

Everything under `/lcrc/group/e3sm/ac.xylar/ismip7/forcing-data/data/AIS/parameterisations/`
and `.../obs/ocean/` (§3.3). This covers all TF/S forcing, all J1–J4 targets, the IMBIE2 basin
masks, BFRN bins, floating masks, the PIG/Dotson shelf mask, and both Bedmap3 and
BedMachine v3 topography on the ISMIP 8 km grid.

### To acquire

1. **`ISMIP7_8km-60m_AIS_grid_ocean_v3.nc`** (46 MB) — Globus,
   `/ISMIP7/AIS/grid/ocean/ISMIP7/8km-60m/v3/`. Needed as the source-grid description for
   ESMF weight generation. (Alternative: generate it with `i7aof.grid.ismip` — check
   equivalence.)
2. **The relaxed MALI initial condition** — `relaxed_10yrs_4km.nc`, currently NERSC-only,
   *if* the ISMIP7 IC is the relaxed state rather than one of the LCRC inputdata vintages.
   Blocked on Q2. The mesh itself, five IC vintages, the SCRIP file, graph partitions, region
   mask and `tf_params` are all already on LCRC (§3.2) — nothing to fetch there.
3. *(Config B only)* nothing new — Bedmap3 on the 8 km grid is already local.
4. *(Only if Q3 → option b)* nothing new; the Burgard form needs salinity (`so`) at the draft,
   which is already present for every ocean state (`*_S.nc`), plus latitude for the Coriolis
   parameter and the ice-draft slope, both derivable on the MALI mesh.

### Reference documents

Not committed here (third-party PDFs); local copies live in the work area at
`/lcrc/group/e3sm/ac.xylar/ismip7/MALI-melt-calibration/`:

* `ISMIP7_AIS_ice_ocean_protocol.pdf` — Reese et al., the ISMIP7 AIS ice–ocean protocol.
  Sect. 4.2 is the calibration protocol implemented here.
* `reference-docs/ISMIP7_Protocol_Cheat_Sheet.pdf` — 15 pp., Google Drive
  `1pRTauVlAWHr80I5qA_wrfCmwAzvM0tvq`. Source of the "avoid linear formulations, update the
  ISMIP6 formulation following Burgard et al. (2022)" recommendation and the 31 July 2026
  focus-group update on the new regional ocean datasets.
* `reference-docs/ISMIP7_Globus_CLI_Workflow.pdf` — 7 pp., Google Drive
  `14hAtVqtPu1upCiEPAE0vA_MaDJB9AZ5r`. Endpoint IDs and transfer recipes.

The two Drive documents are linked from <https://www.ismip.org/research/ismip7>. Their `/view`
URLs are not directly fetchable; use
`curl -L "https://drive.google.com/uc?export=download&id=<ID>"`.

---

## 7. Phased work plan

Phases 1–3 are unblocked and can start immediately. Phase 5 is the one that depends on Q3.

| Phase | Work | Blocked by | Output |
|---|---|---|---|
| **0. Feedback log** | start `protocol-and-toolbox-questions.md`; append throughout, do not defer to the end | — | manuscript notes for Ronja; toolbox items for the PR |
| **1. Scaffold** | pixi env; package skeleton + CLI; config system; dataset registry with checksums; fetch the 8 km grid file via Globus | — | branch commit; `mali-melt-calib inputs` runs green |
| **2. Terms on unstructured meshes** | area-weighted `calculate_term1..4`; unit tests that reproduce the structured-grid answers when given a uniform-area mesh; end-to-end replication of the published quadratic-example numbers (median K = 8.5e-5, 5th = 4.75e-5, 95th = 13.75e-5) on the 8 km grid as a regression test | — | `terms.py` + tests; upstreamable PR |
| **3. Mesh preparation** | remap IMBIE2 basins, BFRN bins, floating mask, PIG/Dotson mask onto the MALI mesh; build Config A and Config B geometry files. Basin mapping verified (§5.0); the remap tool itself lands in MPAS-Tools/Compass | Q11 for the tool's home; Q2 only for the final IC choice — can start now on `ais_4to20km.20250625.nc` | `mali-melt-calib mesh` + an upstream tool |
| **4. Forcing remap** | 8 km → MALI-mesh remap of 11 (then 26) TF fields (and `so` too if Q3 → (b)/(c)); reuse compass `process_thermal_forcing` logic standalone; SCRIP file already on LCRC | Q7 | `mali-melt-calib forcing` |
| **5. MALI ensemble** | run-directory generation, namelists/streams, job scripts; 3-value linearity check; production runs | **Q3**, Q1, Q11 | `mali-melt-calib ensemble` + melt fields |
| **6. Calibration + report** | assemble ensembles, run the 100,000-sample optimisation, produce protocol Fig. 5/7 equivalents for Config A and B; optional per-basin ΔT (Q5: after optimisation) | — | parameter values + plots + write-up |
| **7. Upstream PR** | contribute the MALI example, the mesh-agnostic terms, and the Part B toolbox fixes | — | PR to `ismip7-antarctic-ocean-forcing` |

**Suggested immediate next step:** Phase 2. Reproducing the protocol's published quadratic
numbers on the 8 km structured grid, through our own code path, validates the whole stage-2
chain and gives us a regression test that the unstructured generalisation must not break. It
needs no MALI runs and no answers to any of the open questions.

Phase 3 can proceed in parallel now that the mesh is confirmed on LCRC — starting on
`ais_4to20km.20250625.nc` and switching IC vintage later is a config change, not rework.

---

## 8. Risks

| Risk | Impact | Mitigation |
|---|---|---|
| MALI lacks the Burgard (2022) formulation (Q3) | Calibration is of `gamma0` in the ISMIP6 form, not the ISMIP7-recommended form | Decide early; the toolbox is parameter-agnostic, so option (a) is a valid fallback and the pipeline is unchanged |
| Wrong IC vintage chosen (Q2) | Calibration is done against a geometry the projections do not use | Mesh + 5 IC vintages are on LCRC, so this is a re-run, not a re-acquisition; keep the IC a config option and re-run phases 3–6 (cheap) once Q2 is settled |
| MALI/ISMIP7 basin numbering disagrees | Silently wrong J1/J3/J4 — worst failure mode here, because it produces plausible numbers | **Confirmed real:** MALI is 1-based, ISMIP7 0-based, so MALI basin 10 = ISMIP7 basin 9 (§5.0). Derive basins by remapping the ISMIP7 mask and keep the cross-tabulation as a regression test |
| Ice-shelf area mismatch biases J1/J4 | Systematically shifted parameter distribution | Config B (present-day Bedmap3 geometry) quantifies it; report both |
| TF gaps after remap to the MALI mesh | Spurious zero/invalid melt near grounding lines | Compare MALI-mesh `TFdraft` against 8 km TF at draft; check `config_invalid_value_TF` handling; use MALI's extrapolation if needed |
| Linearity assumption wrong | Ensemble under-sampled | Explicit 3-value numerical check before relying on it |
| Focus-group data revisions (a `_v4` appears) | Rework | Dataset registry pins versions + checksums; a refresh is a config change, not a code change |

---

## 9. Reference paths

In this repository:

```
Upstream toolbox     parameterisations/parameter_selection_toolbox.py
Worked example       parameterisations/parameter_selection_quadratic_example.ipynb
Toolbox data notes   parameterisations/README.md
This plan            parameterisations/mali/mali-melt-calibration-plan.md
```

In `MALI-Dev/E3SM` @ `develop`:

```
MALI melt code       components/mpas-albany-landice/src/mode_forward/mpas_li_iceshelf_melt.F:1167
MALI registry        components/mpas-albany-landice/src/Registry.xml:381 (config_basal_mass_bal_float)
```

In `MPAS-Dev/compass` @ `main`:

```
Compass forcing      compass/landice/tests/ismip7_forcing/
Compass AIS run      compass/landice/tests/ismip7_run/ismip7_ais/
```

On LCRC:

```
ISMIP7 data (local)  /lcrc/group/e3sm/ac.xylar/ismip7/forcing-data/data/AIS/
MALI mesh + IC       /lcrc/group/e3sm/public_html/inputdata/glc/mpasli/mpas.ais4to20km/
                       ais_4to20km.20250625.nc            (newest IC vintage, nCells=385379)
                       ais_4to20km.20230105.scrip.nc      (SCRIP, for ESMF weights)
                       ais_4to20km_region_mask.20230105.nc (16 ISMIP6 basins)
                       ais_4to20km_tf_params.20250724.nc  (gamma0=14500, deltaT=0)
                       mpasli.graph.info.240507.part.*    (graph partitions)
Globus tooling       /lcrc/group/e3sm/ac.xylar/ismip7/forcing-data/ismip7-source-coop/
Existing conda env   /gpfs/fs1/home/ac.xylar/chrysalis/miniforge3/envs/ismip7_dev
Work area            /lcrc/group/e3sm/ac.xylar/ismip7/MALI-melt-calibration/
                       ISMIP7_AIS_ice_ocean_protocol.pdf
                       reference-docs/                   (cheat sheet, Globus manual)
                       MALI-Dev/                         (E3SM clone, develop)
                       add-mali-melt-calibration/        (this repo, worktree)
```

Globus endpoint: `ccc9bbd2-4091-4e35-addd-eeb639cf5332` (`GHub_upload`), tree root `/ISMIP7/`.
