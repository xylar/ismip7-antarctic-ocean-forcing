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
