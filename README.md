# Replication Package for: "Distributional Dynamics"

> Christian Bayer, Luis Calderon, Moritz Kuhn

## 📥 Synthetic data downloads

**If you just want the ready-to-use, high-frequency distributional series — start here.** All files are quarterly from **1962-Q3 to 2024-Q1**, joint over consumption, income, and wealth.

| Want | File | Size |
|------|------|------|
| The 43 latent factors $\hat F_t$ that summarize all dynamics | [`data/synthetic/smoothed_factors.csv`](data/synthetic/smoothed_factors.csv) | 209 KB |
| Raw $\hat\phi^j_t$ inputs (PSID, with trend / cyclical) | [`PSID_functional_data.csv`](data/synthetic/PSID_functional_data.csv) · [`_detrended`](data/synthetic/PSID_functional_data_detrended.csv) | 2.9 MB · 5.5 MB |
| Raw $\hat\phi^j_t$ inputs (SCF, with trend / cyclical) | [`SCF_functional_data.csv`](data/synthetic/SCF_functional_data.csv) · [`_detrended`](data/synthetic/SCF_functional_data_detrended.csv) | 1.7 MB · 1.7 MB |
| **Reconstructed PSID coefficients** (in-sample, HP trend) | [`PSID_coefficients_normal.csv`](data/synthetic/PSID_coefficients_normal.csv) | 4.5 MB |
| **Reconstructed PSID coefficients** (extrapolation, averaged trend) | [`PSID_coefficients_average.csv`](data/synthetic/PSID_coefficients_average.csv) | 8.3 MB |
| Reconstructed SCF coefficients (HP / averaged trend) | [`_normal`](data/synthetic/SCF_coefficients_normal.csv) · [`_average`](data/synthetic/SCF_coefficients_average.csv) | 2.2 MB each |
| Reconstructed CEX coefficients (HP / averaged trend) | [`_normal`](data/synthetic/CEX_coefficients_normal.csv) · [`_average`](data/synthetic/CEX_coefficients_average.csv) | 2.0 MB · 2.2 MB |

**Trend convention:** `_normal` files add back the **HP-filter trend** (date-anchored — use inside 1962–2024). `_average` files add back the **time-average** of the HP trend (use for extrapolation / dates outside the sample).

**Functional data vs coefficients:** `*_functional_data*.csv` are the *inputs* to the state-space model (per-dataset estimates with NaN where the survey wasn't run). `*_coefficients_*.csv` are the *outputs* of the Kalman smoother — dense across all 248 quarters.

See [`data/synthetic/README.md`](data/synthetic/README.md) for column naming, file formats, and a quick-start snippet.

---

## Overview

This package provides all code and instructions needed to reproduce the results in "Distributional Dynamics." The code estimates a state-space model that combines multiple household surveys (PSID, SCF, CEX, CPS, SIPP) to reconstruct the joint distribution of income, wealth, and consumption at quarterly frequency from 1962 to 2024.

The replication package is organized into three stages:

| Stage | Command | Description | Approx. Runtime |
|-------|---------|-------------|-----------------|
| 1: Data | `bash master.sh data` | Clean raw survey data | ~1 hour |
| 2: Estimation | `bash master.sh estimate` | MCMC estimation (4 chains) | ~48-72 hours |
| 3: Results | `bash master.sh results` | Post-estimation analysis, figures, tables | ~3 hours |

Pre-computed posterior estimates are provided in `output/estimates/`, allowing replicators to skip Stage 2 and go directly from Stage 1 to Stage 3.

## Data Availability Statement

All data used in this paper are publicly available. Some datasets require free registration. No restricted-access or confidential data are used.

| Dataset | Provider | Access | Years | URL |
|---------|----------|--------|-------|-----|
| Panel Study of Income Dynamics (PSID) | University of Michigan | Free registration | 1968-2021 | https://psidonline.isr.umich.edu/ |
| Survey of Consumer Finances (SCF) | Federal Reserve Board | Public download | 1962-2022 | https://www.federalreserve.gov/econres/scfindex.htm |
| Consumer Expenditure Survey (CEX) | Bureau of Labor Statistics | Public download | 1984-2023 | https://www.bls.gov/cex/ |
| Current Population Survey (CPS) | Census Bureau / IPUMS | Free registration | 1964-2023 | https://cps.ipums.org/ |
| Survey of Income and Program Participation (SIPP) | Census Bureau | Public download | 1984-2022 | https://www.census.gov/programs-surveys/sipp.html |
| World Inequality Database (WID) | WID.world | Public download | 1962-2023 | https://wid.world/ |
| Distributional Financial Accounts (DFA) | Federal Reserve Board | Public download | 1989-2023 | https://www.federalreserve.gov/releases/z1/dataviz/dfa/ |
| FRED Economic Data | Federal Reserve Bank of St. Louis | Public download | Various | https://fred.stlouisfed.org/ |
| Forbes 400 | Forbes / Fernholz & Haslberger (2023) | Public | 1985-2024 | https://www.forbes.com/forbes-400/ |

**Forbes 400 details:** The SCF+ is augmented with the Forbes 400. For the years 2021 to 2024, we use data directly from Forbes. For the years 1985 to 2020, we use the per capita dataset of Fernholz and Haslberger (2023), whose observations originate from families of the Forbes 400. This provides more complete coverage (nearly all 400 observations each year) than the Forbes website alone, which has incomplete records from the 1990s through the late 2000s.

**Included in this package:** Synthetic (model-generated) microdata are provided in `data/synthetic/`. Pre-computed posterior parameter vectors are provided in `output/estimates/`.

**Not included:** Raw survey microdata must be downloaded by the replicator from the sources above. See `data/raw/DOWNLOAD_INSTRUCTIONS.md` for detailed variable lists and extraction instructions.

## Computational Requirements

### Software

| Software | Version | Purpose |
|----------|---------|---------|
| Julia | >= 1.10 | Model estimation and all quantitative analysis |
| Stata | >= 17 (MP recommended) | Survey data cleaning and local projections |
| Python | >= 3.10 | Data preprocessing (stationarity, seasonal adjustment) |
| R | >= 4.3 | Hermite series estimation, copula estimation, X-12 |

### Julia Packages

All Julia dependencies are pinned in `code/julia/env/Project.toml` and `code/julia/env/Manifest.toml`. Key packages include: Optimization, AdvancedHMC, Distributions, Interpolations, KernelDensity, JLD2, BlackBoxOptim, JuMP, HiGHS, ChebyshevApprox, RCall.

### Python Packages

pandas, numpy, statsmodels, scipy.

### R Packages

seasonal, x12.

### Hardware

- **Minimum:** 32 GB RAM, 8 CPU cores
- **Recommended:** 64 GB RAM, 16+ cores (estimation benefits from parallelization)
- **Storage:** ~10 GB for data and results

Estimation (Stage 2) was performed on a Linux computing cluster. Post-estimation (Stage 3) can run on a standard workstation.

### Runtime

| Stage | Hardware | Approximate Time |
|-------|----------|-----------------|
| 1: Data | Standard workstation | ~1 hour |
| 2: Estimation | 16-core server, 64 GB RAM | ~48-72 hours |
| 3: Results | Standard workstation | ~3 hours |

Using pre-computed estimates (skipping Stage 2), the total runtime is approximately 4 hours.

## Package Contents

```
├── README.md                     This file
├── LICENSE
├── master.sh                     Master replication script
│
├── code/
│   ├── julia/
│   │   ├── DistributionalDynamics.jl   Main Julia entry point (loads all modules)
│   │   ├── run_estimation.jl           Stage 2: MCMC estimation
│   │   ├── run_postestimation.jl       Stages 3-5: results, figures, tables
│   │   ├── config.jl                   Path configuration
│   │   ├── Structures.jl              Data structures and model options
│   │   ├── DataConstructor.jl         Data loading and preprocessing
│   │   ├── ModelPrep.jl               Functional data preparation and PCA
│   │   ├── Model.jl                   State-space model specification
│   │   ├── MCMC.jl                    MCMC sampling and optimization
│   │   ├── SelectPrior.jl             Minnesota prior construction
│   │   ├── DIMES.jl                   DIME sampler integration
│   │   ├── DIMESampler.jl             DIME sampler implementation
│   │   ├── SeparationStrategy.jl      Distributional/aggregate factor separation
│   │   ├── BlackBoxSSM.jl             Black-box optimization
│   │   ├── EM.jl                      EM algorithm
│   │   ├── RobustProjection.jl        Robust factor projection
│   │   ├── HyperparameterOptimization.jl  Hyperparameter tuning
│   │   ├── Reconstruction.jl          Synthetic distribution reconstruction
│   │   ├── CreateTimeSeries.jl        Time series export (quantiles, shares, levels)
│   │   ├── IntervalEstimation.jl      Confidence interval estimation
│   │   ├── ForecastSSM.jl             Out-of-sample forecasting
│   │   ├── Validation.jl              External validation (DFA, WID, SCF, CPS, ACS)
│   │   ├── Correlations.jl            Correlation analysis
│   │   ├── CorrelationTables.jl       Correlation table export
│   │   ├── FEVD.jl                    Forecast error variance decomposition
│   │   ├── HistoricalDecomposition.jl Historical decomposition
│   │   ├── CounterfactualRuns.jl      Counterfactual exercises
│   │   ├── CyclicalityOfConsumption2.jl  Recession consumption dynamics
│   │   ├── AllPlots.jl                Proof-of-concept and main figures
│   │   ├── CEXFunctions.jl            CEX data utilities
│   │   ├── MDD.jl                     Marginal data density
│   │   ├── SupportingFunctions.jl     Utility functions
│   │   └── env/                       Julia environment (Project.toml, Manifest.toml)
│   ├── stata/
│   │   ├── 00_master_data.do          Master Stata script (Stage 1)
│   │   ├── data_cleaning.do           Main survey data cleaning (SCF, PSID, CEX, ACS, CPS)
│   │   ├── clean_SCF_2022.do          SCF 2022 wave cleaning
│   │   ├── clean_SIPP_panels.do       SIPP panel cleaning
│   │   ├── connect_SIPP_panels.do     SIPP panel concatenation
│   │   ├── sipp_panel_constructor.do  SIPP quarterly aggregation
│   │   ├── master_sipp_file.do        SIPP utilities
│   │   ├── insert_interview_dates.do  PSID interview date formatting
│   │   ├── import_aggregates.do       FRED macro series import
│   │   ├── x12series.do               X-12 seasonal adjustment
│   │   ├── process_WIDq.do            WID quarterly data processing
│   │   ├── prep_micro_data.do         Micro data for local projections
│   │   ├── prep_micro_3D.do           3D micro data preparation
│   │   ├── prep_macro_3D.do           3D correlation measures
│   │   ├── prep_functional_data.do    Functional data preparation
│   │   ├── prep_actualCEX.do          CEX data with income imputation
│   │   ├── ginis_and_consumption_trends.do  Gini coefficients and trends
│   │   ├── other_results.do           Aggregate correlations
│   │   └── SIPP/                      SIPP panel do-files by wave
│   ├── python/
│   │   ├── convert_monthly_to_quarterly.py  Frequency conversion
│   │   ├── make_stationary.py         Stationarity transformations
│   │   ├── make_stationary_x12.py     Post-X12 stationarity tests
│   │   ├── generateForbes400.py       Forbes 400 wealth data
│   │   ├── clean_brake_data.py        Geographic data cleaning
│   │   ├── MEILC_MEGC.py             Multivariate Gini coefficient
│   │   └── multidim_inequality.py     Multidimensional inequality measures
│   └── R/
│       ├── HermiteSeriesEstimator.R   Hermite series density estimation
│       ├── NonParametricCopula.R      Non-parametric copula estimation
│       └── X12_script.R              X-12 seasonal adjustment
│
├── data/
│   ├── raw/                           Raw survey data (user-provided)
│   │   └── DOWNLOAD_INSTRUCTIONS.md
│   ├── processed/                     Cleaned data (generated by Stage 1)
│   ├── synthetic/                     ⭐ READY-TO-USE quarterly distributional series
│   │   ├── README.md                  File-by-file documentation
│   │   ├── smoothed_factors.csv       Latent factors F_t
│   │   ├── *_functional_data*.csv     Inputs: raw phi-hat per dataset
│   │   └── *_coefficients_*.csv       Outputs: Kalman-smoother reconstructions
│   └── aggregates/                    Public macro data
│
├── output/
│   ├── figures/                       Generated figures
│   ├── tables/                        Generated LaTeX tables
│   └── estimates/                     Pre-computed parameter vectors
│
├── slides/                            Conference / seminar slide decks
│   └── brazil/                        Brazil 2026 (full + short 20-min version)
│
├── doc/                               Methodological notes
│   └── linearization_notes.tex
│
└── CITATION.cff                       Machine-readable citation
```

## Instructions to Replicators

### Quick Start (using pre-computed estimates)

1. Install Julia >= 1.10, Stata >= 17, Python >= 3.10, R >= 4.3.
2. Download raw survey data following `data/raw/DOWNLOAD_INSTRUCTIONS.md`.
3. From the project root directory, run:
   ```bash
   bash master.sh data      # Stage 1: ~1 hour
   bash master.sh results   # Stages 3-5: ~3 hours
   ```

### Full Replication (including estimation)

1. Complete Steps 1-2 above.
2. From the project root directory, run:
   ```bash
   bash master.sh all       # All stages: ~50-75 hours
   ```
   Or run stages individually:
   ```bash
   bash master.sh data      # Stage 1: clean data
   bash master.sh estimate  # Stage 2: MCMC estimation
   bash master.sh results   # Stage 3: post-estimation
   ```

### Configuration

All paths are configured in `code/julia/config.jl`. The script auto-detects the project root directory. No manual path editing should be necessary.

### Random Seed

The MCMC sampler uses Julia's default random number generator. For exact replication of the estimation, set `Random.seed!(12345)` at the top of `run_estimation.jl`. Post-estimation results are deterministic given the parameter estimates.

## List of Tables and Figures

### Main Text

| Exhibit | Source Code |
|---------|-------------|
| Figure 1 | TikZ (in main.tex) |
| Figure 2 | `SeriesEstimators.jl` |
| Figure 3 | TikZ (in main.tex) |
| Figure 4 | `AllPlots.jl`: `gen_proof_of_concept_*()` |
| Figure 5 | `AllPlots.jl`: `gen_proof_of_concept_*()` |
| Figure 6 | `ForecastSSM.jl` |
| Figure 7 | `CreateTimeSeries.jl`: `generate_specific_plots()` |
| Figure 8 | Separate HANK replication package (not included) |
| Figure 9 | Stata (TBD) |
| Figure 10 | `OrderAnalysis.jl` |
| Figure 11 | `OrderAnalysis.jl` |
| Figure 12 | `OrderAnalysis.jl` |
| Figure 13 | `DIMESampler.jl` |
| Figure 14 | `CreateTimeSeries.jl` |
| Figure 15 | `Validation.jl` |
| Table 1 | `.tex` (manual) |
| Table 2 | TBD |
| Table 3 | `FEVD.jl` |
| Table 4 | `FEVD.jl` |
| Table 5 | `MDD.jl` |
| Table 6 | TBD |
| Table 7 | TBD |
| Table 8 | TBD |
| Table 9 | TBD |

## References

Bayer, Christian, Luis Calderon, and Moritz Kuhn. "Distributional Dynamics." CEPR Discussion Paper 19829, 2026.

Board of Governors of the Federal Reserve System. *Survey of Consumer Finances*. Various years (1962-2022). https://www.federalreserve.gov/econres/scfindex.htm.

Board of Governors of the Federal Reserve System. *Distributional Financial Accounts*. https://www.federalreserve.gov/releases/z1/dataviz/dfa/.

Federal Reserve Bank of St. Louis. *FRED Economic Data*. https://fred.stlouisfed.org/.

Panel Study of Income Dynamics. *Public Use Dataset*. Produced and distributed by the Survey Research Center, Institute for Social Research, University of Michigan, Ann Arbor, MI. Various years (1968-2021). https://psidonline.isr.umich.edu/.

Fernholz, Ricardo T. and Haslberger, Christoph. "Rising Concentration and Group Lending: The Forbes 400 in General Equilibrium." *Journal of Economic Dynamics and Control*, 2023.

Ruggles, Steven, Sarah Flood, Ronald Goeken, et al. *IPUMS CPS: Version 11.0* [Current Population Survey]. Minneapolis, MN: IPUMS. https://cps.ipums.org/.

U.S. Bureau of Labor Statistics. *Consumer Expenditure Surveys*. Various years (1984-2023). https://www.bls.gov/cex/.

U.S. Census Bureau. *Survey of Income and Program Participation*. Various panels (1984-2008). https://www.census.gov/programs-surveys/sipp.html.

World Inequality Database. *WID.world*. https://wid.world/.
