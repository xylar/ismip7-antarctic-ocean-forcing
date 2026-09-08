# Points of confusion: ISMIP7 AIS ice–ocean protocol and parameter-selection toolbox

Collected while implementing the calibration for MALI — the first unstructured-mesh,
variable-resolution ice-sheet model to go through this protocol. See
[`mali-melt-calibration-plan.md`](mali-melt-calibration-plan.md) §3.6 for why this document
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
Those are separable, so a model can and should satisfy both at once.

The practical implications, which the text leaves implicit:

* A model whose initial state has been relaxed or spun up should preferably calibrate on its
  **pre-relaxation** geometry. Not because a short relaxation moves the geometry far — the
  differences are subtle — but because the drift is model-specific, and reducing exactly that
  model-to-model variation is the point of asking everyone to use a common present-day
  geometry.
* A model with no present-day inversion at all (most of them; MALI is unusual here) has no
  present-day geometry to fall back on, and must *substitute* one into its initial state for
  the calibration rather than use its own initial geometry.

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

### A3. §4.2.3 says the J3/J4 inclusion pre-factors are sampled from {0,1}; the code samples U(0,1) — proposed

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

**Settled by replication.** Running the full 8 km calibration through our own code path with
the toolbox's continuous `U(0,1)` weighting reproduces all three published percentiles
exactly:

===========  ===========  ===========
percentile   ours         published
===========  ===========  ===========
5th          4.750e-5     4.75e-5
50th         8.500e-5     8.5e-5
95th         1.375e-4     13.75e-5
===========  ===========  ===========

So the published numbers were produced with continuous weighting, and it is the **manuscript
text that is out of date**, not the code.

**Suggestion:** change §4.2.3 to say that a weight is drawn uniformly from [0,1] for each
summand of J3 and J4, rather than a pre-factor from {0,1}. If Bernoulli inclusion was the
intent, then the code and the published percentiles both need revisiting, which is worth
knowing sooner rather than later.

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

**Partly answered by measurement.** Computing the mean over floating cells of the Bedmap3
draft on the ISMIP 8 km grid, with multimelt's own finite-difference scheme, gives

    mean angle = 0.0051117 rad,  sin(mean angle) = 0.0051117

which matches §4.3.1's stated `sinθ ≈ 0.005`. So the paper's 0.005 *is* the
geometry-derived mean, not Burgard et al.'s 2.9e-3 — those are two different quantities and
the text reads as though they were alternatives.

That also disposes of question (i): at θ ≈ 0.005 rad, `sinθ` and `θ` agree to one part in
10^5, so averaging angles rather than sines makes no practical difference. What remains is
(ii) which value modellers should adopt, and (iii) the resolution dependence, which is the
part that actually bites: `K` absorbs the slope, so a model computing its own mean slope on a
finer mesh gets a systematically different `K` and the two are not comparable.

**Our working assumption:** compute the constant slope on our own mesh, per the notebook's
advice that it "should match the geometry", and report the value alongside `K`.

**A second, separable point: the `sin θ` notation itself is confusing when θ is constant.**
Written as `K sin θ ...`, it reads as a geometric quantity that varies over the ice shelf. In
the constant-slope variant it is nothing of the kind — it is a fixed dimensionless
coefficient multiplying `K`, and only the product `K · sin θ` affects the melt. Two readers
here (an ice-sheet modeller and an ice-ocean modeller) independently found this a stumbling
block, which suggests it is worth a sentence rather than being left implicit.

**Suggestion:** state explicitly whether the constant slope is meant to be fixed across
models (making `K` comparable) or recomputed per model grid (making it not), and if the
latter, ask modellers to report the slope they used. And in the constant-slope case, say
plainly that `sin θ` is then a fixed coefficient absorbed into `K` — or quote the value of
`sin θ` directly rather than as the sine of an angle.

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

### A7. The Coriolis parameter is a constant in the reference implementation — minor

Protocol Eq. (1) contains `g / (2|f|)` with `f` the Coriolis parameter, but
`multimelt.constants` defines a single `f_coriolis = 1.4e-4` and the worked example uses it
everywhere. The manuscript does not say whether a constant or a latitude-varying `f` is
intended.

Quantifying it before making anything of it: `f = 2Ω sin(φ)` gives

| latitude | \|f\| (s⁻¹) |
|---|---|
| 63°S (Peninsula tip) | 1.300e-4 |
| 74°S | 1.402e-4 |
| 85°S (southern Ross/Ronne) | 1.453e-4 |

So `|f|` varies by about **12%** across the latitudes where ice shelves actually sit — not
the large factor one might assume from the full pole-to-equator range — and multimelt's
1.4e-4 corresponds to 73.7°S, near the middle of that span. Since melt goes as `1/|f|`, using
the constant introduces at most a ~6% spatial modulation either side of the mean, and `K`
absorbs the mean.

**This is therefore a documentation point, not a correctness one.** We use the constant
`f = 1.4e-4` to stay consistent with the published calibration.

**Suggestion:** state in §4.1.1 that `f` is taken as constant, and give the value, so
implementers do not each decide independently.

### A8. The observational datasets are not climatology-filled, but the regional model ones are — open

The ISMIP7 datasets treat partial spatial coverage two different ways, and only one is
documented:

| dataset | finite fraction on the 8 km grid |
|---|---|
| `Jourdain-Naughten_NEMO-MITgcm_cold_TF.nc` | 1.0000 |
| `Timmermann_FESOM_cold_v3_TF.nc` | 1.0000 |
| `Obs_2009_TF.nc` | 0.0603 |

§A9 says the regional ocean-model domains have "the remaining basins filled with the ISMIP7
ocean climatology", and the distributed files are indeed complete. §A10 says the
observational profiles are applied "uniformly to the entire Amundsen basin", and the
distributed files are undefined outside it — about 94% missing.

An ice-sheet model needs valid forcing wherever it has ice, so every group using the J4
datasets has to invent a fill, and different choices will give different melt outside the
Amundsen (which does not enter J4, but does have to not crash the model). We fill from the
present-day climatology before remapping, matching how the model datasets were prepared.

**Suggestion:** either distribute the observational files climatology-filled, as the model
files already are, or state in §A10 what modellers should do outside the covered basin. The
former would be more consistent and would remove a step everyone has to reinvent.

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
unstructured meshes (area = `areaCell`). This is what we are contributing; see plan §3.4.

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

### B4b. The worked example uses two different ice densities — proposed

`multimelt.constants` sets `rho_i = 917.0`, which is what its `melt_factor` is built from,
but `parameter_selection_quadratic_example.ipynb` defines its own `ice_density = 918` and
uses that to convert the melt rate from m of ice to kg m⁻² yr⁻¹.

Ice density appears twice in protocol Eq. (1) — once in the `(rho_o / rho_i)` factor and once
in the conversion to a mass flux — and the two should cancel exactly, leaving melt independent
of the value chosen. Using 917 in one place and 918 in the other leaves a spurious factor of
918/917, about 0.1%.

That is far too small to matter for the calibration, and it is folded into the published `K`
in any case. But it took some care to disentangle when writing an implementation intended to
agree with the reference, and it will do so again for the next group.

**Change:** use a single ice density, taken from `multimelt.constants`, in both places.

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
| 2026-09-08 | Added B4b: the worked example uses rho_i = 917 in the melt factor and 918 in the mass conversion, leaving a spurious 0.1% factor where the two should cancel. |
| 2026-09-08 | Added A8: the observational J4 datasets are distributed unfilled (6% coverage) while the regional model datasets are climatology-filled. |
| 2026-09-08 | Extended A4: the sin(theta) notation is confusing when theta is constant, since it is then a fixed coefficient absorbed into K rather than a geometric quantity. |
| 2026-09-08 | A3 settled by replication: continuous U(0,1) reproduces all three published percentiles exactly, so the manuscript text is what needs updating. |
| 2026-09-08 | A4 partly answered by measurement: the mean draft slope on Bedmap3 at 8 km is sin(theta) = 0.0051117, matching the paper's 0.005. Added A7 on the constant Coriolis parameter. |
| 2026-09-08 | A7 corrected: |f| varies ~12% across ice-shelf latitudes (1.300e-4 at 63S to 1.453e-4 at 85S), not a factor of two, and multimelt's 1.4e-4 sits at 73.7S. Downgraded from open to minor. |
| 2026-09-08 | Working assumptions recorded for A3 and A4 so implementation can proceed; both still need a focus-group answer. |
| 2026-09-08 | A1 reframed: the geometry and grid/code requirements are separable, not in conflict; the ask is a clarification plus guidance for groups without a present-day initialisation. A2 confirmed intended; the ask is a sentence explaining why. Both moved open → proposed. |
