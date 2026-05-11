# Synthetic Distributional Data

This folder contains **ready-to-use, high-frequency** estimates of the joint distribution of consumption, income, and wealth for U.S. households from **1962-Q3 to 2024-Q1**. They are the main output of the *Distributional Dynamics* method.

Citation: please cite Bayer, Calderon, and Kuhn (2026). See [`CITATION.cff`](../../CITATION.cff).

## Files

| File | Rows × Cols | What it is |
|---|---|---|
| [`smoothed_factors.csv`](smoothed_factors.csv) | 248 × 44 | **Latent factors $\hat F_t$** (43 factors + time). These are the states of the Kalman smoother and the most compact summary of distributional + aggregate dynamics. |
| [`PSID_functional_data.csv`](PSID_functional_data.csv) | 248 × 1091 | **Raw $\hat\phi^j_t$ from PSID**, with trend retained (levels). Coefficients of the Legendre-polynomial expansion of the joint distribution. NaN in quarters where PSID has no observation. |
| [`PSID_functional_data_detrended.csv`](PSID_functional_data_detrended.csv) | 248 × 1091 | Same as above but **cyclical component only** (HP-trend removed). This is what enters the state-space observation equation. |
| [`SCF_functional_data.csv`](SCF_functional_data.csv) | 248 × 1091 | Raw $\hat\phi^j_t$ from SCF, with trend. NaN in non-survey quarters. |
| [`SCF_functional_data_detrended.csv`](SCF_functional_data_detrended.csv) | 248 × 1091 | Cyclical component of SCF $\hat\phi^j_t$. |
| [`PSID_coefficients_normal.csv`](PSID_coefficients_normal.csv) | 248 × 1731 | **Model-implied PSID coefficients with HP trend re-added.** Use this for the in-sample period (1962-Q3 to 2024-Q1) — trend is date-anchored. |
| [`PSID_coefficients_average.csv`](PSID_coefficients_average.csv) | 248 × 1731 | Same model output, but the **time-averaged** trend is added back. Use this when you want to extrapolate or compare beyond the HP-anchored sample. |
| [`SCF_coefficients_normal.csv`](SCF_coefficients_normal.csv) | 248 × ~1.7 K | SCF reconstruction with HP trend (in-sample). |
| [`SCF_coefficients_average.csv`](SCF_coefficients_average.csv) | 248 × ~1.7 K | SCF reconstruction with averaged trend (extrapolation-friendly). |
| [`CEX_coefficients_normal.csv`](CEX_coefficients_normal.csv) | 248 × ~1.7 K | CEX reconstruction with HP trend (in-sample). |
| [`CEX_coefficients_average.csv`](CEX_coefficients_average.csv) | 248 × ~1.7 K | CEX reconstruction with averaged trend (extrapolation-friendly). |

All time series are quarterly. Time index format: `1962-Q3 … 2024-Q1`.

## Two trend conventions: `_normal` vs `_average`

The state-space model is estimated on the **cyclical** component of the data ($\boldsymbol{\hat\phi}_t$ minus a trend). To reconstruct levels for use, the trend has to be added back. We provide two conventions:

- **`_normal`** — adds back the **HP-filter trend**, which is *date-anchored*. This is the right object when you want the reconstruction at observed survey dates (1962-Q3 through 2024-Q1).
- **`_average`** — adds back the **time-average of the HP trend**, which is a single constant per coefficient. This is the right object when you want to push estimates **outside the sample** (forecasts / nowcasts at dates where the HP trend is not well identified).

If you're not sure which to use: use **`_normal`** inside 1962–2024, **`_average`** otherwise.

## Functional data vs coefficients

- **`*_functional_data*.csv`** are *inputs* to the state-space model — the raw, dataset-specific estimates of the Legendre coefficients $\hat\phi^j_t$, with NaN where the dataset wasn't observed.
- **`*_coefficients_*.csv`** are *outputs* — the model-implied coefficients $\hat\phi_t = \Gamma \hat F_t$ reconstructed via the Kalman smoother, dense across all 248 quarters.

The factor representation ($\hat F_t$ in `smoothed_factors.csv`) is the most compact summary; the coefficient files are the same information mapped back into coefficient space.

## Column naming

Functional-data and coefficient columns are indexed by polynomial order on each margin. Example: `ciw_321.0` is the coefficient on the Legendre polynomial of order $(o_c, o_y, o_w) = (3, 2, 1)$ in (consumption, income, wealth).

The exact ordering follows `code/julia/CreateTimeSeries.jl`. To map a coefficient time series back to a distribution, see `code/julia/Reconstruction.jl`.

## Quick start (Python)

```python
import pandas as pd

# Just want the latent state-space factors?
F = pd.read_csv("smoothed_factors.csv", index_col="time")

# Want the reconstructed PSID distribution (in-sample, with HP trend)?
psid = pd.read_csv("PSID_coefficients_normal.csv", index_col="time")
```

### Compute marginal quantiles and copula density

The `*_coefficients_*.csv` files store Legendre polynomial coefficients. To
evaluate $\Xi^{-1}_{m,t}(u)$ or $dC_t(u_c, u_y, u_w)$ at arbitrary quantile
points, use the standalone helper at [`code/python/reconstruct.py`](../../code/python/reconstruct.py):

```python
from reconstruct import Reconstruction

r = Reconstruction("PSID_coefficients_normal.csv")

# Marginal quantile function for consumption in 2008-Q3 at deciles
r.quantile_at("2008-Q3", "consum", [0.1, 0.5, 0.9])
# -> array([0.320, 0.734, 1.386])  (relative to mean)

# Copula density at the joint median
r.copula_density_at("2008-Q3", 0.5, 0.5, 0.5)
# -> 1.754

# Full 30 x 30 x 30 density grid (handy for heatmaps / iso-surfaces)
grid = r.copula_density_grid("2008-Q3", n=30)   # shape (30, 30, 30)
```

Or as a CLI:

```bash
python code/python/reconstruct.py data/synthetic/PSID_coefficients_normal.csv \
    quantile --date 2008-Q3 --measure consum --u 0.1,0.5,0.9

python code/python/reconstruct.py data/synthetic/PSID_coefficients_normal.csv \
    copula --date 2008-Q3 --u 0.5,0.5,0.5
```

The basis is the orthonormal Legendre family on $[0,1]$: $Q_o(u) = \sqrt{2o+1}\,P_o(2u-1)$. Column layout: $x_1$..$x_{1694}$ are the mutable copula coefficients $\kappa$ (the leading "L-shape" of 34 immutable entries is reconstructed by the script), and $x_{1695}$..$x_{1730}$ are the marginal quantile coefficients $\xi$ — 12 polynomial orders per measure, in alphabetical order: consum, income, wealth.

## Provenance

Each file is generated by the replication pipeline in this repository. See the [main README](../../README.md) for the full pipeline and the paper for methodology.
