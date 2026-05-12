# User guide

This package supports three kinds of users. Pick the lane that matches what you want:

| Lane | You want… | Start here |
|---|---|---|
| 1 | The published quarterly distributional data — straight off the shelf | [§ 1](#1-just-the-data) |
| 2 | Custom moments (quantiles, copula density, counterfactual factor values) computed from the published coefficients | [§ 2](#2-compute-custom-moments-from-coefficients) |
| 3 | To re-run the model with your own settings (different tag, factor count, muted waves, etc.) | [§ 3](#3-run-the-model-yourself) |

---

## 1. Just the data

Every CSV in [`data/synthetic/`](../data/synthetic/) is quarterly **1962-Q3 → 2024-Q1**, joint over consumption, income, and wealth. See [`data/synthetic/README.md`](../data/synthetic/README.md) for column-by-column documentation. The two most useful files for a quick start:

```python
import pandas as pd

# 43 latent factors F_t — most compact summary of the joint dynamics
F = pd.read_csv("data/synthetic/smoothed_factors.csv", index_col="time")

# Reconstructed PSID coefficients (in-sample, HP-trend conventions)
psid = pd.read_csv("data/synthetic/PSID_coefficients_normal.csv", index_col="time")
```

**Trend convention.** Files named `_normal` re-add the HP-filter trend (date-anchored, use inside 1962–2024). Files named `_average` re-add the time-averaged trend (extrapolation-friendly).

**Functional data vs coefficients.** `*_functional_data*.csv` are *inputs* to the state-space model (per-dataset $\hat\phi^j_t$, NaN where the survey wasn't run). `*_coefficients_*.csv` are *outputs* (Kalman-smoother reconstructions, dense across all 248 quarters).

---

## 2. Compute custom moments from coefficients

The published coefficient files store Legendre polynomial weights. To turn them into moments you can plot (decile means, copula densities at arbitrary points, etc.), use the small standalone helpers shipped with the repo:

* Python: [`code/python/reconstruct.py`](../code/python/reconstruct.py) — needs only `numpy` + `pandas`
* Julia: [`code/julia/reconstruct.jl`](../code/julia/reconstruct.jl) — needs only `CSV.jl` + `DataFrames.jl`

Both implement the same orthonormal Legendre basis $Q_o(u) = \sqrt{2o+1}\,P_o(2u-1)$ and agree to machine precision.

### Marginal quantile function

```python
from reconstruct import Reconstruction
r = Reconstruction("data/synthetic/PSID_coefficients_normal.csv")
r.quantile_at("2008-Q3", "consum", [0.1, 0.5, 0.9])
# → array([0.320, 0.734, 1.386])    (relative to per-HH mean)
```

```julia
using .DistributionalReconstruction
r = Reconstruction("data/synthetic/PSID_coefficients_normal.csv")
quantile_at(r, "2008-Q3", :consum, [0.1, 0.5, 0.9])
```

### Copula density

```python
r.copula_density_at("2008-Q3", 0.5, 0.5, 0.5)        # density at the joint median
r.copula_density_grid("2008-Q3", n=30)               # 30³ grid for plotting
r.copula_pmf_grid("2008-Q3")                         # 10³ probability masses (matches the published ciw_*)
```

### Factor → coefficient mapping

If you want to ask *"what would the joint distribution look like if factor 3 were one standard deviation higher?"*, there are two ways to get a coefficient row from a factor vector. **Prefer the projection mode**, which uses the same linear map the model's Kalman smoother uses internally and is therefore exact.

#### Mode A — exact reconstruction via the model's PCA loadings (Recommended)

The state-space estimator already computes a per-dataset projection matrix $\Gamma_j$. We export it (plus the per-coefficient means, stds, and trends) at the end of every model run into a `projection/` subfolder. Once those CSVs are in `data/synthetic/projection/`, do:

```python
from reconstruct import FactorMap
fm = FactorMap.from_projection(
    "data/synthetic/projection/",
    dataset="PSID",          # or "SCF", "CEX"
    trend_kind="normal",     # "normal" = HP-trend, date-anchored; "average" = time-mean trend
)
print(fm.summary())          # FactorMap (projection): dataset=PSID, K=8, n_coefs=1730, …

import numpy as np
F = np.zeros(fm.n_factors); F[0] = 1.0    # counterfactual: factor 1 at +1 SD
fm.quantile_at(F, "consum", [0.1, 0.5, 0.9], date="2008-Q3")
fm.copula_density_at(F, 0.5, 0.5, 0.5, date="2008-Q3")
```

The `date` argument selects the HP-trend row for `trend_kind="normal"`; with `trend_kind="average"` the trend is a constant per coefficient and `date` is ignored.

Julia mirrors the API:

```julia
pm = load_projection_map("data/synthetic/projection/", "PSID"; trend_kind="normal")
println(DistributionalReconstruction.summary(pm))

F = zeros(pm.n_factors); F[1] = 1.0
quantile_at(pm, F, :consum, [0.1, 0.5, 0.9]; date="2008-Q3")
copula_density_at(pm, F, 0.5, 0.5, 0.5; date="2008-Q3")
```

Internally the helper computes the same identity that lives inside `Reconstruction.jl`:

$$\text{coef}_t = \text{stds} \odot (\Gamma \cdot F^{\text{dist}}_t) + \text{means} + \text{trend}_t$$

where $F^{\text{dist}}_t$ is the distributional portion of the smoothed factors (the first $\text{factor\_count}\times 4$ entries of `smoothed_factors.csv`). For the published 8-factor model this is the first 32 columns.

#### Mode B — OLS approximation (fallback, no projection files needed)

If you don't have a `projection/` folder (e.g., you're working only from the public coefficient + smoothed factor CSVs), you can fit a quick OLS map. One regression per coefficient, rows with any NaN dropped:

```python
fm = FactorMap(
    "data/synthetic/smoothed_factors.csv",
    "data/synthetic/PSID_coefficients_normal.csv",
    n_factors=8,                  # 8 distributional factors; use 43 for all
)
print(fm.summary())
# FactorMap (OLS): K=8, T_used=100 (dropped 147 NaN rows of 247), R² median=0.222, …
```

R² is moderate at K=8 (the OLS is a linear approximation to what PCA + Kalman do). Bumping to `n_factors=43` uses the full file and drives the median R² to ~0.70. The OLS mode does **not** require a `date` argument because the means/stds are estimated from the data directly — there's no trend to re-add.

The two modes have the same `.predict`, `.quantile_at`, `.copula_density_at` API; pick whichever input data you have.

---

## 3. Run the model yourself

The entry points are:

```bash
bash master.sh data        # Stage 1: clean raw survey data (~1 h)
bash master.sh estimate    # Stage 2: MCMC estimation (~48–72 h)
bash master.sh results     # Stage 3: post-estimation, figures, tables (~3 h)
```

The single file you actually edit when configuring a new run is [`code/julia/Structures.jl`](../code/julia/Structures.jl) — specifically the [`ModelOptions`](../code/julia/Structures.jl) struct. Its docstring lists every field; the seven you'll touch in practice are:

| Field | What it controls | Typical change |
|---|---|---|
| `tag` | Result-folder suffix | `" my_experiment"` (leading space) — outputs go to `7_Results/<m_label><tag>/...` |
| `measures` | Which household variables to model jointly | Always `["consum", "income", "wealth"]` in the paper |
| `number_of_dfs` | # of distributional factors retained from PCA | Default `8` (~99% of business-cycle variation) |
| `lags` / `agg_lags` | AR depth of states and aggregate factors | Defaults `1` / `4` |
| `blind_to` | Mute a measure in a specific survey | `Dict("CEX" => ["wealth"])` etc. |
| `data_to_mute` | Drop survey waves entirely (validation experiments) | `Dict("SCF" => muted_quarters_between(QuarterlyDate(2004,1), QuarterlyDate(2009,4)))` |
| `estimator` | `SeriesEstimator` polynomial order / output grid | `SeriesEstimator(grid_pcf=14, grid_cop=14, …)` for higher-order |

The block of comments immediately below `ModelOptions` lists tag conventions used in the paper's validation runs (e.g., `" excluding housing cycle wealth"`, `" every 4 years"`, `" 6 factors"`), with the matching `data_to_mute` recipes. Use them as templates.

### Output layout for a fresh tag

When you bump `tag` to something new, the pipeline writes to a brand-new tree (the `mkpath` calls scattered through the post-estimation code ensure every subfolder is created on first write):

```
7_Results/
└── consum_and_income_and_wealth<tag>/
    ├── from_mcmc/
    │   ├── data/                    ← *_coefficients_*.csv, smoothed_factors.csv, *_functional_data*.csv
    │   │   └── projection/          ← PCA loadings Γ, means, stds, trends (one set per dataset)
    │   ├── plots/                   ← quantile / copula plots + correlation tables
    │   └── bayesian_convergence/
    └── other_results/                ← counterfactuals, raw data, likelihood diagnostics
```

The `projection/` subfolder is what feeds `FactorMap.from_projection(...)` — copy it (or the files inside) over to `data/synthetic/projection/` when you want users to use the exact reconstruction.

### Sanity checks before kicking off

* `model_options.tag` is new (won't collide with old results)
* If you set `compare_to_other_est = true`, the baseline run with the comparison tag must already exist
* `data_to_mute` keys match an actual dataset name in `obs_data.files`

When the run finishes, drop the newly generated `*_coefficients_*.csv`, `*_functional_data*.csv`, and `smoothed_factors.csv` into `data/synthetic/` (overwriting the previous publication) and re-run the validation in [`code/python/validate_reconstruction.py`](../code/python/validate_reconstruction.py) to confirm coherence between the coefficient and functional-data outputs.
