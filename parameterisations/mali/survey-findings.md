# Survey findings: MALI, Compass, MPAS-Tools and the ISMIP7 datasets

Evidence behind the decisions in
[`mali-melt-calibration-plan.md`](mali-melt-calibration-plan.md). The plan states
conclusions; this document records what was actually inspected, measured or computed to
reach them, so the conclusions can be checked rather than taken on trust.

**Date of survey:** 2026-09-08.

*Path conventions:* bare relative paths such as
`parameterisations/parameter_selection_toolbox.py` are relative to the root of this
repository. Paths beginning with `/` are absolute on LCRC. Source files in other projects are
given relative to that project's checkout, named in the heading.

---

## Part 1 — What already exists

### 1.1 MALI (`MALI-Dev/E3SM` @ `develop`)

Working clone: `/lcrc/group/e3sm/ac.xylar/ismip7/MALI-melt-calibration/MALI-Dev`, clean, with
`matthewhoffman`, `trhille` and `xylar` remotes configured.

`config_basal_mass_bal_float = 'ismip6'` selects `iceshelf_melt_ismip6` in
`components/mpas-albany-landice/src/mode_forward/mpas_li_iceshelf_melt.F:1167`. This is the
**ISMIP6 non-local quadratic** (Jourdain et al. 2020), driven by:

| Field | Shape | Units |
|---|---|---|
| `ismip6shelfMelt_gamma0` | scalar | m yr⁻¹ |
| `ismip6shelfMelt_deltaT` | nCells | K |
| `ismip6shelfMelt_basin` | nCells | integer index |
| `ismip6shelfMelt_offset` | nCells | kg m⁻² s⁻¹ |
| `ismip6shelfMelt_3dThermalForcing` | nISMIP6OceanLayers × nCells | °C |
| `ismip6shelfMelt_zOcean` | nISMIP6OceanLayers | m |
| `ismip6shelfMelt_TFdraft` | nCells (diagnostic) | °C |

**It is not the Burgard et al. (2022) update** that ISMIP7 recommends. The ISMIP7 cheat sheet
is explicit: *"The only recommendation is to avoid linear formulations and to update the
ISMIP6 formulation following Burgard et al. (2022)."* Protocol Eqs. (1)–(2) add `sinθ`,
`β_S·S_loc` and `g/(2|f|)`, and change `(ρ_o/ρ_i)` from squared to first power.

Two further observations relevant to how we run it:

* `config_velocity_solver = 'none'` exists, and only `'L1L2'`, `'FO'` and `'Stokes'` require
  external dycores — so melt-only runs should not need Albany.
* `li_basal_melt_floating_ice` is called *inside* the timestep
  (`mpas_li_time_integration_fe_rk.F:288`), **not** in the initial diagnostic solve. A
  `config_stop_time == config_start_time` diagnostic-only run would therefore produce no
  melt; we must take one short real timestep.

### 1.2 Compass (`MPAS-Dev/compass` @ `main`, clone `/home/ac.xylar/compass/main`)

* `compass/landice/tests/ismip7_forcing/ocean_thermal/process_thermal_forcing.py` remaps 3-D
  TF from the ISMIP 8 km × 30-level grid onto a MALI mesh, writing
  `ismip6shelfMelt_3dThermalForcing`, `ismip6shelfMelt_zOcean` and
  `ismip6shelfMelt_zBndsOcean`. This is the remapping the calibration needs, just driven by
  calibration TF files rather than CMIP scenario files.
* `compass/landice/tests/ismip7_run/ismip7_ais/` is the ISMIP7 AIS production setup and shows
  the intended configuration: `config_velocity_solver = 'FO'`,
  `config_basal_mass_bal_float = 'ismip6'`, `config_dt` = 1 month, an `ismip7_params` input
  stream carrying `gamma0`/`deltaT`/`basin`, an initial condition at
  `.../AIS_4to20km_20230105/relaxation_0TGmelt_10yr/relaxed_10yrs_4km.nc` and melt parameters
  at `.../basin_and_coeff_DeltaT_quadratic_non_local_gamma14500.nc`.

Those two paths are on NERSC CFS, but equivalents are already on LCRC — see 1.4.

### 1.3 MPAS-Tools (`MPAS-Dev/MPAS-Tools` @ `master`, clone `/home/ac.xylar/mpas_work/MPAS-Tools/master`, `be1954cf`)

Three pieces are directly relevant. Together they mean the MALI-side work is mostly
*extending* MPAS-Tools rather than writing from scratch.

**`landice/mesh_tools_li/tune_ismip6_melt_deltat.py`** — the tool that produced the production
parameter file. On the MALI mesh it already implements:

* `TFdraft` by interpolating 3-D TF to `lowerSurface`, and area-weighted `meanTF` per region
* the MALI ISMIP6 quadratic: `melt = gamma0 * coef * (TFdraft + dT) * |meanTFcell + dT|`
* **area-weighted basin aggregation to Gt yr⁻¹**:
  `(melt[ind] * areaCell[ind]).sum() * rhoi / 1.0e12` — exactly the J1 aggregation that 1.5
  says the toolbox lacks, already written for an unstructured mesh
* a per-region `dT` sweep against an observational target, writing
  `basin_and_coeff_DeltaT_quadratic_non_local_gamma<gamma0>.nc` — the filename the compass
  ISMIP7 AIS config points at, with the same default `gamma0 = 14500`

*Gap:* its target is Paolo et al. alone, spatially averaged over floating cells per region.
ISMIP7's J1 target is the three-product mean (Paolo/Davison/Adusumilli) with the expanded
uncertainty of protocol §A8.

**`landice/output_processing_li/ismip7_postprocessing/grid_and_mapping.py`** — ISMIP7-aware
remapping infrastructure:

* `create_ismip7_grid_file(icesheet, res_km, output_file)` writes the ISMIP7 standard grid
  (AIS EPSG:3031, −3040 km…+3040 km) as x/y coordinate variables — all a SCRIP/ESMF
  source-grid description needs.
* `build_mapping_file(...)` builds SCRIP files via `mpas_tools.scrip.from_mpas.scrip_from_mpas`
  and calls `ESMF_RegridWeightGen`. Wired MALI → ISMIP for post-processing; the calibration
  needs ISMIP → MALI, i.e. source and destination swapped.
* `check_ismip7_grid_file(...)` validates extents.

**`mpas_tools.landice.interpolate.interpolate_to_mpasli_grid`** — structured grid → MALI mesh
with bilinear, barycentric or ESMF-weight methods. Expects CISM or MPAS input conventions, so
ISMIP7 rasters would need adapting.

**What does not exist:** nothing produces the *ISMIP7* IMBIE2 basin field, the BFRN bins or
the PIG/Dotson shelf mask on a MALI mesh. The existing basin field derives from ISMIP6
`regionCellMasks` (built with `geometric_features` in 2022), not from
`basin_numbers_ismip8km_v2.nc`.

### 1.4 The MALI mesh and its inputs are on LCRC

`/lcrc/group/e3sm/{data,public_html}/inputdata/glc/mpasli/mpas.ais4to20km/`:

| File | Contents |
|---|---|
| `ais_4to20km.{20230105,20240708,20241224,20250411,20250625}.nc` | IC vintages, `nCells = 385,379`, with `thickness`, `bedTopography`, `areaCell`, `xCell`/`yCell`/`latCell` (no `lowerSurface`/`cellMask` — derived at init) |
| `ais_4to20km.20230105.scrip.nc` | SCRIP description, ready for ESMF weight generation |
| `ais_4to20km_region_mask.20230105.nc` | `regionCellMasks(nCells, 16)`, named `ISMIP6 Basin A-Ap` … `K-A` |
| `ais_4to20km_tf_params.20250724.nc` | `gamma0 = 14500`, **`deltaT = 0` everywhere**, `basin` in 0…16 |
| `mpasli.graph.info.240507.part.{64…3840}` | graph partitions |

The `public_html` copy additionally has the 20230105 IC and the SCRIP file; the `data` copy
does not. Because compass points at NERSC paths instead, it is easy to conclude these are
unavailable on LCRC — they are not.

### 1.5 ISMIP7 calibration data is almost entirely on LCRC

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
| `obs/ocean/topography/{bedmap3,BedMachineAntarctica-v3}/v3/` | 8 km topography and masks | complete |
| `obs/ocean/IMBIE-basins/v3/` | basin polygons/masks | complete |
| `grid/ocean/ISMIP7/8km-60m/v3/` | `ISMIP7_8km-60m_AIS_grid_ocean_v3.nc` (46 MB) | **absent locally** |

All TF/S fields are on the ISMIP 8 km polar-stereographic grid, 761 × 761 × 30 levels.

Diffed against the Globus manifest `ismip7-source-coop/manifests/ISMIP7_AIS.json.gz`
(109,834 files, 3.5 TB): in the subtrees the calibration needs, the only gaps are that one
grid file and `parameterisations/fracture/ant_surf_roughness_*.nc`, which melt calibration
does not use.

### 1.6 Globus access is already solved

`/lcrc/group/e3sm/ac.xylar/ismip7/forcing-data/ismip7-source-coop/` is a pixi project with
`globus-cli` and working scripts (`transfer_to_lcrc.py`, `globus_survey.py`,
`refresh_manifest.py`), plus gzipped manifests of the whole remote tree under `manifests/`.

Endpoint **`ccc9bbd2-4091-4e35-addd-eeb639cf5332`** (`GHub_upload`), tree root `/ISMIP7/`.

### 1.7 The parameter-selection toolbox assumes a structured grid

`parameterisations/parameter_selection_toolbox.py` hard-codes uniform-cell-area, 2-D-raster
assumptions:

* `cvt_m = reso**2 / 1e12` — a single scalar cell area (`calculate_term1`, `2`, `4`)
* `.groupby(basins_m).sum()` over a 2-D `(y, x)` field
* `obs_ensemble.stack(grid=('x','y'))` in `calculate_term4`

MALI's mesh is unstructured and variable-resolution (4–20 km), so all four terms need
area-weighted, `nCells`-based equivalents. The toolbox README already notes "only regular
grids are supported at the moment". 1.3 gives a working reference implementation of the
area-weighted aggregation to follow.

### 1.8 Prior art: a slope-dependent melt option was tried in MALI and rejected

`MALI-Dev/E3SM` PR #48, *"Add option to include shelf base slope in ISMIP6 melt param."*
(opened 2022-09-12, closed 2022-11-10 without merging). Flagged by Matt Hoffman.

It added `config_ismip6shelfMelt_use_slope` to the existing ISMIP6 non-local quadratic,
multiplying the melt by a **per-cell slope field**:

```
<var name="shelfBaseSlope" type="real" dimensions="nCells Time" ...>
call calc_shelf_base_slope(geometryPool, meshPool)      ! every melt call
floatingBasalMassBal(iCell) = -coef * shelfBaseSlope(iCell) * (TFdraft + deltaT) * abs(...)
```

recomputed each timestep from `lowerSurface` on edges, then smoothed and capped.

Matt's closing comment:

> further testing revealed the slope-dependent form had undesirable properties - it has a
> tendency to melt holes in ice shelves when TF increases and also exhibited evidence of
> undesirable feedbacks between shelf base slope and melt rate. With further theory
> development, it may be possible to implement a variation to this, but for now, we can
> consider this a dead end.

**Why this does not apply to `config_basal_mass_bal_float = 'ismip7'` as implemented.** Both
symptoms follow from the slope being a *geometry-derived field*: melt thins the shelf locally,
which steepens the local basal slope, which increases melt. Our implementation takes the slope
from a **scalar namelist option**, `config_ismip7_melt_sin_slope`, which enters `coef`
alongside `K`; melt therefore depends only on the product `K · sinθ` and the slope cannot
respond to geometry at all. The feedback has no pathway.

This is independent support for the scoping decision in plan §5.3 to use a constant
Antarctic-mean slope rather than the locally varying variant, which protocol Eq. (1) also
permits. **Do not add the local-slope variant without revisiting this PR.**

**What remains open.** Our change relative to ISMIP6 is a different one — non-local to local,
i.e. `TF_local × |TF_basin-mean|` to `TF_local × |TF_local|`. The basin mean acts as a spatial
damper, so removing it sharpens the melt pattern where TF is locally high. Whether Matt's
concern attaches to the slope specifically or to anything that sharpens the pattern is worth
asking him, since he has seen the failure mode. If the latter, protocol Eq. (2) — the
*semi-local* form, which keeps a shelf- or basin-mean factor while adopting the Burgard
constants — is an explicitly permitted fallback that would retain the damping of MALI's
current behaviour.

---

## Part 2 — Verified results

These were computed during the survey rather than read off. Each states the check so it can
be repeated.

### F1. Melt is exactly linear in `gamma0` (and in `K`)

MALI's ISMIP6 quadratic (`mpas_li_iceshelf_melt.F:1369-1374`) is

```
BMB = -gamma0 · (ρ_sw c_p / (ρ_i L_i))² / s_yr · ρ_i · (TF_draft + ΔT_b) · |<TF>_b + ΔT_b| + offset
```

For a fixed geometry and a fixed ΔT field, `TF_draft`, `<TF>_b`, `cellMask` and
`connectedOceanMask` are all independent of `gamma0`. So melt is exactly — not approximately —
linear in `gamma0`. The Burgard et al. (2022) form is likewise linear in `K`.

**Consequence:** one MALI run per *ocean state* suffices (≈11 for the recommended dataset
set); the full 120-value parameter ensemble follows by scaling. That turns a ~1,300-run
campaign into ~15 runs.

**Caveat:** this holds only while ΔT is fixed. Fitting ΔT *inside* the optimisation would
break it — which is one reason the ΔT-after-optimisation ordering matters.

**Still to verify numerically:** run 3 values of `gamma0` for one ocean state and confirm
exact proportionality before relying on it.

### F2. ISMIP7 basins are 0-based; MALI's are 1-based

The single most dangerous discrepancy found, because it produces plausible-looking wrong
answers rather than an error.

* ISMIP7 `basin_numbers_ismip8km_v2.nc` is **0-based, 0…15**, covering every grid cell with
  no "outside" value. Confirmed geographically: sampling the basin field under
  `shelf_mask_ismip8km.nc`, **PIG and Dotson both fall in basin 9** (Eastern Amundsen), and
  basin 14's centroid is Ronne-Filchner — matching the protocol's references to
  "b = 9 (Eastern Amundsen)" and "b = 14 (Ronne-Filchner)".
* MALI `ismip6shelfMelt_basin` is **1-based, 1…16**, plus a single stray cell labelled 0. It
  equals the `regionCellMasks` column index + 1.
* Cross-tabulating the two over the 296,538 ice cells of `ais_4to20km.20250625.nc` (ISMIP7
  field sampled nearest-neighbour at MALI cell centres): **MALI = ISMIP7 + 1 for 100.0% of
  ice cells**, with per-basin purity ≥ 98.8% for 15 of 16 basins. The exception is
  MALI 14 / ISMIP7 13 (Antarctic Peninsula, `Ipp-J`) at 85.2%, where the ISMIP6-2022 and
  IMBIE2-v3 boundaries genuinely differ.

So **MALI basin 10 = ISMIP7 basin 9** and **MALI basin 15 = ISMIP7 basin 14**. Using one
index where the other is expected silently mis-assigns every basin.

Corroborated in the source: `tune_ismip6_melt_deltat.py:345` builds the field as
`regionCells[regionCellMasks[:, reg] == 1] = reg + 1`.

### F3. The LCRC inputdata meshes are pre-relaxation and essentially present-day

Sampling `ais_4to20km.20250625.nc` against Bedmap3 on the ISMIP 8 km grid, over the 295,819
ice cells where both are defined:

| Field | MALI − Bedmap3 |
|---|---|
| thickness | median −12.4 m, MAD 35.9 m, p5 −236 m, p95 +137 m |
| bed | median −13.8 m, MAD 41.0 m |
| total ice volume | 26.11 vs 26.44 ×10⁶ km³ (−1.2 %) |

Residuals of this size are what one expects from BedMachine (MALI's source, via compass
`landice/tests/antarctica/mesh.py` → `add_bedmachine_thk_to_ais_gridded_data`) versus Bedmap3,
plus 8 km sampling of a 4–20 km mesh.

All five IC vintages have **identical ice volume and ice-cell count**; they differ only in
ancillary field counts (60–85 variables). So vintage choice is not a geometry choice.

**Consequence:** the inputdata geometry already satisfies the protocol's present-day-geometry
requirement, so no geometry replacement and no new mesh is needed.

### F4. The published 8 km calibration reproduces exactly through our code path

Running the full ISMIP7 calibration with our area-weighted terms
(:mod:`mali_melt_calib.terms`) driving the upstream objective function, on the ISMIP 8 km
grid with 100,000 samples:

| percentile | ours | published (§4.3.1, Fig. 5) |
|---|---|---|
| 5th | 4.750e-5 | 4.75e-5 |
| 50th | 8.500e-5 | 8.5e-5 |
| 95th | 1.375e-4 | 13.75e-5 |

Reproduce with ``python -m mali_melt_calib.replicate``.

Three things follow:

* Calibration stage 2 is validated end to end before any MALI run exists, and this is the
  regression test the unstructured generalisation must not break.
* The melt formula in ``quadratic.py`` agrees with the one behind the published numbers, so
  a later difference in ``K`` on the MALI mesh is attributable to the mesh and the model
  rather than to the parameterisation.
* Feedback item A3 is answered: the toolbox samples the J3/J4 inclusion pre-factors from a
  continuous ``U(0,1)`` while protocol §4.2.3 describes Bernoulli ``{0,1}`` inclusion, and
  since continuous weighting reproduces the published percentiles, that is what was used.

Incidental measurement: the mean ice-draft slope over floating ice, computed from Bedmap3 on
the 8 km grid with multimelt's finite-difference scheme, is **sin θ = 0.0051117** — matching
the protocol's stated ``sinθ ≈ 0.005``, and confirming that the paper's value is the
geometry-derived mean rather than Burgard et al.'s 2.9e-3.

### F5. The remapped thermal forcing has no gaps on the MALI mesh

Bilinearly remapping the ISMIP7 present-day climatology from the 8 km grid onto
`ais_4to20km.20250625.nc` gives thermal forcing that is **finite at all 385,379 cells**, over
all 30 ocean layers.

This is expected — the ISMIP7 fields are already extrapolated into ice-shelf cavities and
under grounded ice on the 8 km grid (protocol §3.2) — but it is worth having measured, because
it is the evidence behind assuming `config_ocean_data_extrapolation = .false.` (plan §5.2,
Q7). MALI has nothing left to extrapolate.

The check to repeat if the mesh or the remapping method changes:

```python
tf = ds['ismip6shelfMelt_3dThermalForcing'].values
assert np.isfinite(tf).all()
```

### F6. MALI's melt matches the reference formula to round-off

The ISMIP7 local quadratic was run in MALI on the 4-20 km mesh, forced by the present-day
climatology, one timestep with the velocity solver off, at K = 8.5e-5:

| | |
|---|---|
| melting cells | 109,412 |
| MALI melt | 1.81 .. 32,153 kg m⁻² yr⁻¹ |
| max relative difference from the reference | **7.0e-16** |

The comparison evaluates :func:`mali_melt_calib.quadratic.local_quadratic_melt` on **MALI's
own** ``TFdraft`` and ``Sdraft``, which isolates the melt expression from the vertical
interpolation.  Reproduce with ``python -m mali_melt_calib.verify <run_dir>``.

**What it does not test:** the vertical interpolation of thermal forcing and salinity to the
ice draft.  Checking that needs MALI's ``TFdraft`` compared against a Python interpolation of
the 3-D field, which is the natural next step.

**A real inconsistency it caught.** The first comparison disagreed by 6.63e-4 — too large for
round-off, too small for a structural error.  It is entirely the year length: MALI's ``scyr``
is a 365-day year (31,536,000 s, matching its noleap calendar) while the protocol's reference
implementation uses 365.2422 days (31,556,926.08 s).  The ratio is 1.000663562, against an
observed 6.631e-4.  The melt rate in kg m⁻² s⁻¹ is unambiguous; only the conversion to a
per-year rate differs.  ``seconds_per_year`` is now part of the constants set so each
implementation uses its own year.  Negligible physically, but it would have sat in the melt
fields as unexplained noise.

### F7. Semi-local Burgard is degenerate with MALI's existing ISMIP6 method

Implemented, then removed. The reasoning is worth keeping so it is not re-added.

Write both out with the temperature correction in place:

```
ISMIP6 non-local:    melt = C6 * (TF_loc + dT) * |<TF> + dT|
Burgard semi-local:  melt = C7 * <S> * (TF_loc + dT) * |<TF> + dT|

C6 = gamma0 * (rho_sw c_o /(rho_i L_i))^2 * rho_i / s_yr
C7 = K * sin(theta) * (rho_o/rho_i)(c_o/L_i)^2 * beta_S * g/(2|f|) * rho_i
```

Same functional form. Only the decomposition of the constant differs, so **with salinity
held constant the two are algebraically identical** and every `K` has an exactly equivalent
`gamma0`: at the published K = 8.5e-5 that is gamma0 = 11,519 m yr⁻¹, against MALI's
production 14,500.

With a *basin-mean* salinity they are not exactly identical, but the difference is the
per-basin spread of `<S>`, measured on the mesh at **±0.8%** (34.11 to 34.67 across the 16
basins). Since `dT_b` is fitted per basin, that freedom absorbs most of it.

So providing semi-local would mean a second code path, plus a global reduction, computing
what `iceshelf_melt_ismip6` already computes. It was removed on that basis;
`config_basal_mass_bal_float = 'ismip6'` gives that behaviour.

**Removing it does not cost us the local-vs-semi-local comparison.** Because the two are
degenerate, calibrating `gamma0` through MALI's existing `ismip6` code path *is* the
semi-local calibration, and the two can be compared on a fair footing: `gamma0` + 16 `dT_b`
against `K` + 16 `dT_b`, the same number of free parameters, so the minimised objective is a
legitimate model-selection test.  F10 suggests J2 is where the discrimination will live,
since J1 is largely absorbable by the free parameter and `dT_b`.

**The local form is not degenerate.** Replacing `|<TF>|` with `|TF_loc|` changes the spatial
pattern of melt in a way no per-basin constant can reproduce — which is both why it is worth
having and why it carries the risk Matt raised in PR #48 (§1.8).

A corollary worth recording for the salinity question: since `beta_S` multiplies `S` rather
than being a function of it (TEOS-10 gives only 0.04% variation in `beta_S` across
34.2-34.8), melt is *linear* in salinity, so the ±0.8% maps one-to-one with no amplification.
The neglected pressure dependence of `beta_S` is larger — 2.8% from the surface to 1800 dbar,
using the surface value — and unlike a uniform offset it correlates with draft depth, so `K`
does not absorb it.

### F8. MALI's vertical interpolation to the draft is exact

F6 checked the melt *expression* by feeding the Python reference MALI's own ``TFdraft`` and
``Sdraft``, which deliberately says nothing about how those were produced.  Checking the
interpolation independently -- a plain ``np.interp`` in depth, written from the protocol
rather than transliterated from the Fortran -- gives, on the 8 km climatology over 109,412
melting cells:

| code path | cells | max │dTF│ | max │dS│ |
|---|---:|---:|---:|
| interior (linear between layer centres) | 102,603 | 4.4e-16 | 7.1e-15 |
| above the shallowest centre | 508 | 0 | 0 |
| below the deepest centre | 219 | 0 | 0 |
| layer below the draft is beneath the bed | 6,082 | 0 | 0 |

All four paths agree to round-off, so both halves of the melt calculation are now verified.

**The trap this exposed.** The first attempt at this check disagreed by up to 7.9e-3 K, and
the cause was not the interpolation: **the diagnostic run evolves the geometry**.  Melt is
applied over the one-day step, thinning the ice by up to **1.67 m**, so the ``lowerSurface``
and ``thickness`` written to ``output_melt.nc`` are *post*-step while ``TFdraft`` was computed
*pre*-step.  Comparing against the output geometry therefore compares two different drafts.
Reconstructing the initial draft from ``mesh.nc`` (``-rho_i/rho_sw * H``, or the bed where
grounded, with the run's 910/1028) collapsed the discrepancy from 7.9e-3 to 4.4e-16.

This does not affect the melt field, which is computed from the initial geometry, and the
calibration terms use melt, ``areaCell`` and the masks, none of which drift.  But **any
diagnostic that pairs output geometry with output melt is inconsistent by up to a metre of
draft**, which is worth knowing before the Compass port writes its own analysis.

### F9. Linearity in K, measured in MALI rather than argued

F1 derives the linearity from the algebra.  This measures it: three MALI runs on the Zhou
climatology at K = 4.0e-5, 8.5e-5 and 1.3e-4, everything else identical.

| K | total melt (Gt/yr) | max │relative deviation from exact scaling│ |
|---|---:|---:|
| 4.0e-5 | 713.81 | 6.4e-16 |
| 8.5e-5 | 1516.84 | 0 (reference) |
| 1.3e-4 | 2319.87 | 6.4e-16 |

Total melt divided by K is identical to every printed digit (1.78451398566586e7).  Melt is
exactly proportional to K at machine precision, so the ensemble needs **one run per ocean
state**, not one per (state, K) pair -- about 11 runs rather than 1,300.  Anything built on
this workflow can rely on it.

### F10. The local form redistributes melt toward buttressing-relevant ice

The open question from PR #48 (§1.8) is whether the Burgard *local* form gives a melt
pattern MALI can live with.  Two runs on the same Zhou climatology, same mesh, same basin
mask, dT = 0: ISMIP7 local at K = 8.5e-5, and ISMIP6 non-local.  The ISMIP6 field is then
rescaled to the same total melt, so only the **pattern** is compared.

Totals: local 1516.8 Gt/yr; non-local at gamma0 = 11,519 gives 1131.5.  The gamma0 that
matches the local total is **15,443 m/yr**, close to MALI's production 14,500 -- so the
local form at the protocol's published K is in a sensible melt regime, not an extreme one.

At equal total melt the patterns still differ substantially:

* correlation 0.92 over the 109,412 melting cells (unweighted)
* area-weighted RMS difference 1105 kg/m²/yr, **110% of the mean melt** of 1005
* peak local melt 32,153 vs 20,340 kg/m²/yr -- the local form is **58% peakier**

The redistribution is systematic in exactly the direction that matters dynamically.  By BFRN
bin -- the J2 target, bin 9 being the most buttressing-relevant and bin 0 passive ice:

| BFRN bin | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| local / non-local | 0.81 | 0.82 | 0.72 | 0.85 | 0.81 | 0.97 | 0.93 | 1.04 | **1.22** | **1.27** |

The local form moves melt off passive ice and onto buttressing-relevant ice: **+27% in bin 9,
+22% in bin 8, −19% to −28% in bins 0-4.**  Per basin the ratios span 0.86 to 1.54.

**Why this matters for the choice.** The per-basin differences are largely absorbable by a
fitted dT_b, since that is a per-basin constant.  The BFRN redistribution is **not** -- it is
a *within*-basin change in where the melt sits, and no per-basin constant can undo it.  So
switching to the local form is not melt-neutral for ice dynamics even after recalibration:
it puts more melt where buttressing responds most.  Spatially (see
``work/tests/melt_local_vs_nonlocal.png``) the difference concentrates at shelf margins and
calving fronts; Ross and Ronne interiors barely change.

This does not settle the choice, but it makes the skepticism concrete and quantified rather
than qualitative.  Caveats: one ocean state, dT = 0 throughout, and the present-day mesh
geometry.

### F11. Calibrated on the MALI mesh, the local form fits three of the four terms better

The first end-to-end evaluation of the objective function against melt that **MALI itself**
produced: 22 single-timestep runs (11 ocean states x local and semi-local), aggregated on the
MALI mesh, compared against the protocol's own targets.

For each term, the weighted mean absolute error is minimised over the form's own free
parameter (targets at their means, so the comparison is deterministic).  Both forms have the
same degrees of freedom, so lower is better:

| term | what it constrains | local | semi-local | winner |
|---|---|---:|---:|---|
| J1 | melt per IMBIE2 basin | **32.43** | 32.96 | local, by 1.6% |
| J2 | melt per BFRN buttressing bin | **132.20** | 167.85 | local, by 21.2% |
| J3 | warm-cold basin-mean melt | 2855.21 | **2611.84** | semi-local, by 8.5% |
| J4 | PIG/Dotson melt per obs year | **10.32** | 13.22 | local, by 21.9% |

The local form wins the two terms it was expected to -- J2, the buttressing-weighted one that
F10 identified as where the forms genuinely differ, and J4, the Amundsen observations -- and
ties on J1.  It loses J3, the warm-minus-cold sensitivity.

**This survives the masking question.** J1, J2 and J4 are integrals and so are sensitive to
shelf extent, while J3 is an area-weighted mean and is not; since the two forms distribute
melt differently, MALI's extent bias is not purely common-mode.  Repeating everything on the
intersection of MALI's floating cells with the observed ISMIP7 mask changes no winner and
barely moves the margins:

| term | J1 | J2 | J3 | J4 |
|---|---|---|---|---|
| MALI floating | local 1.6% | local 21.2% | semi-local 8.5% | local 21.9% |
| MALI ∩ observed | local 2.5% | local 20.3% | semi-local 6.7% | local 17.5% |

**dT_b = 0 is correct here, not a shortcut.** The protocol selects the melt parameter with
dT_b = 0 in every basin, and calibrates dT_b only *afterwards*, as a separate downstream step.
The toolbox is built that way: `optimise_deltaT`, `select_optimal_deltaT` and
`select_subensemble_using_optimal_deltaT` are separate functions, none of them called from
`calculate_objective_function`, and `p2` is a singleton in the published replication.  So this
comparison should **not** be "improved" by fitting dT_b first -- doing so would depart from the
protocol.

A consequence worth noting: **J1 does not discriminate between the forms either way.**  At
dT_b = 0 it is a genuine constraint but comes out a near-tie (32.43 vs 32.96, 1.6%), and once
dT_b is fitted downstream it collapses towards zero for any parameterisation, since a
per-basin constant is exactly the freedom needed to match a per-basin melt integral.  The
discriminating power lives in J2 and J4 against J3.

**Caveats.** The best-fit parameter differs by term -- K = 4.0e-5 (J1), 1.14e-4 (J2), 4.5e-5
(J3), 1.08e-4 (J4) -- which is precisely why the protocol samples random term weights to
produce a *distribution* rather than one value; the published median 8.5e-5 sits inside that
spread.  This uses the mean targets rather than sampling them.

---

## Reference paths

In this repository:

```
Upstream toolbox     parameterisations/parameter_selection_toolbox.py
Worked example       parameterisations/parameter_selection_quadratic_example.ipynb
Toolbox data notes   parameterisations/README.md
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
Mesh generation      compass/landice/tests/antarctica/mesh.py
```

In `MPAS-Dev/MPAS-Tools` @ `master`:

```
deltaT tuning        landice/mesh_tools_li/tune_ismip6_melt_deltat.py
ISMIP7 grid/mapping  landice/output_processing_li/ismip7_postprocessing/grid_and_mapping.py
structured -> MALI   mpas_tools/landice/interpolate.py  (interpolate_to_mpasli_grid)
```

On LCRC:

```
ISMIP7 data          /lcrc/group/e3sm/ac.xylar/ismip7/forcing-data/data/AIS/
MALI mesh + IC       /lcrc/group/e3sm/public_html/inputdata/glc/mpasli/mpas.ais4to20km/
Globus tooling       /lcrc/group/e3sm/ac.xylar/ismip7/forcing-data/ismip7-source-coop/
Existing conda env   /gpfs/fs1/home/ac.xylar/chrysalis/miniforge3/envs/ismip7_dev
Work area            /lcrc/group/e3sm/ac.xylar/ismip7/MALI-melt-calibration/
                       ISMIP7_AIS_ice_ocean_protocol.pdf
                       reference-docs/                   (cheat sheet, Globus manual)
                       MALI-Dev/                         (E3SM clone, develop)
                       add-mali-melt-calibration/        (this repo, worktree)
```

## Reference documents

Not committed here (third-party PDFs); local copies are in the work area above:

* `ISMIP7_AIS_ice_ocean_protocol.pdf` — Reese et al., the ISMIP7 AIS ice–ocean protocol.
  §4.2 is the calibration protocol implemented here.
* `reference-docs/ISMIP7_Protocol_Cheat_Sheet.pdf` — 15 pp., Google Drive
  `1pRTauVlAWHr80I5qA_wrfCmwAzvM0tvq`. Source of the "avoid linear formulations, update the
  ISMIP6 formulation following Burgard et al. (2022)" recommendation and the 31 July 2026
  focus-group update on the new regional ocean datasets.
* `reference-docs/ISMIP7_Globus_CLI_Workflow.pdf` — 7 pp., Google Drive
  `14hAtVqtPu1upCiEPAE0vA_MaDJB9AZ5r`. Endpoint IDs and transfer recipes.

Both Drive documents are linked from <https://www.ismip.org/research/ismip7>. Their `/view`
URLs are not directly fetchable; use
`curl -L "https://drive.google.com/uc?export=download&id=<ID>"`.
