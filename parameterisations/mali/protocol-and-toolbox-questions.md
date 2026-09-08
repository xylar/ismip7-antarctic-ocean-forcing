# Points of confusion: ISMIP7 AIS ice–ocean protocol and parameter-selection toolbox

Collected while implementing the calibration for MALI — the first unstructured-mesh,
variable-resolution ice-sheet model to go through this protocol. See
[`mali-melt-calibration-plan.md`](mali-melt-calibration-plan.md) §4.5 for why this document
exists.

**Part A** is feedback on the manuscript, to pass back to Ronja. **Part B** is concrete
changes to this repository, to propose as part of the PR adding the MALI example.

These are points where *we* were unsure. Some will turn out to be our misreading rather than
a defect — that is still useful, since it marks where the text could be more explicit. Each
entry records what we checked and what we concluded.

Status key: **open** (needs a focus-group answer) · **proposed** (we have a concrete
suggestion) · **resolved** (answered; kept for the record).

---

## Part A — manuscript

### A1. Present-day geometry vs. "same setting as the core simulations" — proposed

*Not a contradiction — a clarification.* Our first reading treated these two §4.2 bullets as
being in tension; they are not, and the text could make that harder to misread.

* bullet 1: *"modellers are asked to use a recent present-day geometry, e.g., BedMap3 or
  BedMachine3, for the calibration rather than an ice-sheet model initial state with
  substantially different grounding line position, ice thickness, or ice shelf extent"*
* bullet 2: *"Ideally, the melt module is calibrated and evaluated in the same setting as it
  is used in the core ISMIP7 ice-sheet simulations, i.e. with the same code and on the same
  grid"*

Bullet 1 constrains the **geometry**; bullet 2 constrains the **code, grid and parameters**.
Those are separable, so a model can and should satisfy both at once. The practical
implication — which the text leaves implicit — is that a model whose initial state has been
relaxed or spun up should calibrate on its **pre-relaxation, present-day** geometry, and a
model with no present-day inversion at all (most of them; MALI is unusual here) must
*substitute* a present-day geometry into its initial state for the calibration rather than use
its own initial geometry.

Suggestion: say explicitly that the geometry requirement is independent of the grid/code
requirement, and add a sentence on what groups without a present-day initialisation should do.
That second point is likely to affect more participants than it does us.

### A2. Integrated vs. averaged targets and ice-shelf area bias — proposed

J1, J2 and J4 compare *integrated* melt (Gt yr⁻¹); J3 compares *basin-averaged* melt
(kg m⁻² yr⁻¹). An ISM whose ice-shelf area differs from the observed area therefore carries a
systematic bias in J1/J2/J4 that J3 does not share, and the bias enters directly into the
calibrated parameter.

**Confirmed intended** (Xylar, focus group): the modelled shelf area *is* part of what is
being calibrated against, and modellers should not normalise it away. The paper could say so —
as written, a reader hitting the J1/J2/J4-vs-J3 asymmetry has to guess whether it is
deliberate.

Suggestion: a sentence in §4.2.2 noting that the integrated terms deliberately hold the model
accountable for its ice-shelf extent, while J3 is deliberately area-insensitive so that melt
*sensitivity* is compared independently of geometry errors. This also explains why bullet 1 of
§4.2 (A1) matters as much as it does: geometry error enters the answer through these terms.

### A3. §4.2.3 says the J3/J4 inclusion pre-factors are sampled from {0,1}; the code samples U(0,1) — open

**Confirmed discrepancy between the manuscript and the toolbox.**

Manuscript §4.2.3, for J3: *"for each summand of J3 we randomly sample a pre-factor from
{0,1} that decides whether it is included"*, and the same wording for J4. That describes
Bernoulli inclusion — each dataset/basin term is either in or out.

`parameter_selection_toolbox.py` instead draws a *continuous* weight:

```text
t3_weights_samples = xr.DataArray(
    np.random.uniform(0, 1, size=(len(...model), len(...basins), sample_size)), ...)
...
t4_weights_samples = xr.DataArray(
    np.random.uniform(0, 1, size=(len(...region), len(...year), sample_size)), ...)
```

`np.random.uniform(0, 1)` is continuous on [0, 1); it never yields the two-point distribution
the text describes. The two are not equivalent — continuous weighting always includes every
term at a random strength, whereas Bernoulli inclusion drops terms entirely, which produces a
broader (and differently shaped) parameter distribution.

Which is intended? If the code is right, §4.2.3 should say "sample a weight uniformly from
[0,1]"; if the text is right, the code should use `np.random.randint(0, 2, ...)`. Either way
this affects the published 5th/95th percentiles, so it is worth settling before more groups
calibrate.

### A4. Which constant slope, and measured how — open

Three different values/definitions of the constant slope appear:

* §4.3.1: *"This example assumes an average constant slope around Antarctica (sinθ ≈ 0.005)"*
* `parameter_selection_quadratic_example.ipynb`: `slope = local_alpha_isf.mean()`, i.e. the
  mean of `arctan(sqrt(xslope² + yslope²))` over floating points — a mean of *angles*, not of
  `sinθ`
* the same notebook's comment: *"or np.arcsin(2.9e-3) which is the value used in Burgard et
  al. 2022 but ideally it should match the geometry"*

Three questions: (i) is the quantity to be averaged `sinθ` or `θ`? (ii) which of 0.005,
2.9e-3, or the geometry-derived mean should modellers use? (iii) since the slope is computed
by finite differences on the ice draft, its mean is resolution-dependent — a 4 km mesh gives a
different mean slope than the 8 km grid, and `K` absorbs the difference. Should the constant
slope be *fixed* across models for comparability, or recomputed per model grid?

For an unstructured mesh there is a fourth: the notebook's `shift(x=±1)` slope has no direct
analogue, so we will use MALI's own gradient operator — which is a different numerical
estimate again.

### A5. Basin index base: text sums 1…16, data files are 0…15 — proposed

Eq. (4) writes the J1 sum as `Σ_{b=1}^{16}` and Eq. (5) as `Σ_{i=1}^{10}`, but the
distributed masks are 0-based: `basin_numbers_ismip8km_v2.nc` has `basinNumber` ∈ {0,…,15}
covering every cell, and `BFRN_ismip8km_v2.nc` has `BFRN_bins` ∈ {0,…,9} (plus NaN).

We verified that the *labels* are consistent despite the notation: PIG and Dotson both fall
in basin 9 and basin 14 is Ronne-Filchner, matching the protocol's references to
"b = 9 (Eastern Amundsen)" and "b = 14 (Ronne-Filchner)". So this is editorial, not
substantive — but the mismatch cost us a real check, and §A7's "we generate 9 bins (1-9), and
add afterwards a bin 0" is the only place the 0-based convention is stated.

Suggestion: write the sums as `Σ_{b=0}^{15}` and `Σ_{i=0}^{9}`, or add a sentence in §A2
stating the convention.

### A6. Eqs. (4)–(7) are sums; the toolbox computes means — proposed

The four J terms are written as weighted sums over basins/bins/datasets/years. The toolbox's
`mae()` applies `.mean(dims)`.

We believe this is harmless: each term is normalised by `median(J_i)` over the parameter
ensemble, and the count of contributing elements is fixed across parameter values, so a
constant factor cancels. But it is a silent difference between the equations as written and
the reference implementation, and a reader checking the code against the paper will trip over
it. Suggestion: note in §4.2.2 that the implementation uses the mean, and that the
normalisation makes this equivalent.

---

## Part B — this repository

### B1. `parameterisations/README.md` names the wrong variable for the cold/warm ensembles — proposed

The README says the present-day ensemble carries `melt_rate`:

> Melt rates are saved in variable called "melt_rate" given in kg/m2/a.

but for the cold/warm ensembles:

> They contain corresponding melt rates modelled with your melt module (variable called
> "melt_rates").

The toolbox uses the singular `melt_rate` throughout — `cold_ensemble['melt_rate']`,
`warm_ensemble['melt_rate']` in `calculate_term3`. **Change:** fix the README to `melt_rate`.

### B2. The toolbox assumes a uniform structured grid — proposed

`cvt_m = reso**2 / 1e12` (a single scalar cell area) appears in `calculate_term1`,
`calculate_term2` and `calculate_term4`; `calculate_term4` also does
`stack(grid=('x','y'))`; all four terms `groupby` a 2-D field. The README acknowledges this
("only regular grids are supported at the moment") but the constraint is not visible from the
function signatures.

**Change:** accept a cell-area array instead of a scalar resolution, and operate on a generic
cell dimension so the same functions serve structured grids (area = `reso**2` broadcast) and
unstructured meshes (area = `areaCell`). This is what we are contributing; see plan §4.4.

### B3. `calculate_term1` takes an unused `nBasins` argument — proposed

```text
def calculate_term1(pd_ensemble, mask_m, basins_m, nBasins, cvt_m, MeltData):
```

`nBasins` is never referenced in the body; the basin count is taken from the groupby.
**Change:** drop the parameter (or use it to validate `basins_m`).

### B4. Example notebooks hard-code a personal data path — proposed

`parameter_selection_quadratic_example.ipynb`:

```text
data_path="/media/NAS2/ISMIP7_public/AIS/parameterisations/ocean"
```

with the intended form commented out above it. The Globus-relative layout is then rebuilt with
`os.path.join(data_path, "../..", 'obs', 'ocean', ...)`, which is hard to follow.

**Change:** read the root from an environment variable (say `ISMIP7_DATA`) with the commented
placeholder as the documented default, and define the `obs/` paths from the same root rather
than by walking up two levels.

### B5. Which dataset versions the examples expect is implicit — proposed

Filenames carry versions, and in `ocean_modelling_data` both `v2` and `v3` of several
datasets sit side by side (`Mathiot_NEMO_cold_v2_*` and `_v3_*`, `Timmermann_FESOM_cold_v2_*`
and `_v3_*`). The quadratic example's `model_names` list picks a specific mix — `v3` for
Mathiot and Timmermann, `v2` for the two Naughten FESOM runs, unversioned for the two new
regional datasets — with nothing stating that this is the current combination rather than an
artefact of when the notebook was last run.

(The climatology's `so` `v4` / `tf` `v3` mix is *not* an example of this: only one version of
each is distributed, so there is no choice to make there.)

**Change:** add a short "dataset versions used" table to `parameterisations/README.md`, or a
manifest file the examples read, so a modeller can tell whether they are reproducing the
published numbers. This matters because the 31 July 2026 focus-group update notes the median
`K` shifted when new datasets were added.

---

## Log

| Date | Change |
|---|---|
| 2026-09-08 | Created; seeded with A1–A6 and B1–B5 from the initial survey. |
| 2026-09-08 | A1 reframed: the geometry and grid/code requirements are separable, not in conflict; the ask is a clarification plus guidance for groups without a present-day initialisation. A2 confirmed intended; the ask is a sentence explaining why. Both moved open → proposed. |
