# Replication Package for: "Distributional Dynamics"

> Christian Bayer, Luis Calderon, Moritz Kuhn

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
| Survey of Income and Program Participation (SIPP) | Census Bureau | Public download | 1984-2008 | https://www.census.gov/programs-surveys/sipp.html |
| American Community Survey (ACS) | Census Bureau / IPUMS | Free registration | 2005-2023 | https://usa.ipums.org/ |
| World Inequality Database (WID) | WID.world | Public download | 1962-2023 | https://wid.world/ |
| Distributional Financial Accounts (DFA) | Federal Reserve Board | Public download | 1989-2023 | https://www.federalreserve.gov/releases/z1/dataviz/dfa/ |
| FRED Economic Data | Federal Reserve Bank of St. Louis | Public download | Various | https://fred.stlouisfed.org/ |
| Forbes 400 | Forbes | Public | 2021-2024 | https://www.forbes.com/forbes-400/ |

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
│   │   ├── plot_HANK.jl               HANK model comparison plots
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
│   ├── synthetic/                     Model-generated microdata (provided)
│   └── aggregates/                    Public macro data
│
└── output/
    ├── figures/                       Generated figures
    ├── tables/                        Generated LaTeX tables
    └── estimates/                     Pre-computed parameter vectors
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

| Exhibit | Description | Source Code | Output Path |
|---------|-------------|-------------|-------------|
| Figure 1 | Estimation workflow | TikZ (in main.tex) | N/A (diagram) |
| Figure 2 | Legendre coefficients | `AllPlots.jl`: `gen_proof_of_concept_copulas()` | `Plots/proof_of_concept/Legendre_coefs_*.pdf` |
| Figure 3 | Copula observational weights | `Correlations.jl` | `Plots/proof_of_concept/copula_weight_*.pdf` |
| Figure 4 | Geometric representation | TikZ (in main.tex) | N/A (diagram) |
| Figure 5 | Quantile function comparison | `AllPlots.jl`: `gen_proof_of_concept_figure()` | `Plots/proof_of_concept/*_quantiles_*.pdf` |
| Figure 6 | KL divergence | `AllPlots.jl`: `gen_proof_of_concept_copulas()` | `Plots/proof_of_concept/KL_divergence_*.pdf` |
| Figure 7 | Out-of-sample predictability (baseline) | `ForecastSSM.jl`: `perform_forecast()` | `7_Results/.../forecasts/extensive/*.pdf` |
| Figure 8 | Out-of-sample predictability (re-estimated) | `CounterfactualRuns.jl` | `7_Results/.../quantiles_levels/*.pdf` |
| Figure 9 | HANK model comparison | `plot_HANK.jl`: `generate_quantiles_shares_levels_HANK()` | `Plots/HANK/*.pdf` |
| Figure 10 | Consumption dynamics across recessions | `CyclicalityOfConsumption2.jl`: `generate_relative_to_peak_plots()` | `Plots/consumption_plots/*.jpg` |
| Figure 11 | Consumption over business cycles | `ginis_and_consumption_trends.do` | `Plots/IRFs_and_Trends/*.pdf` |
| Table 1 | Data coverage | Manual | `Tables/DataCoverage.tex` |
| Table 2 | HANK average correlations | `CorrelationTables.jl` | Inline in main.tex |
| Table 3 | FEVD on distributional factors | `FEVD.jl`: `fevd_dist_vs_agg()` | Inline in main.tex |
| Table 4 | FEVD of group means | `FEVD.jl`: `fevd_dist_vs_agg()` | Inline in main.tex |

### Appendix

| Exhibit | Description | Source Code |
|---------|-------------|-------------|
| Appendix A | Local linear approximation | Analytical (no code) |
| Appendix B | Marginal data density | `MDD.jl` |
| Appendix C | Minnesota prior | `SelectPrior.jl` |
| Appendix D | Bayesian convergence | `MCMC.jl`: `run_diagnostics()` |
| Appendix E | Reconstructing estimates | `Reconstruction.jl`, `CreateTimeSeries.jl` |
| Appendix F | Hyperparameter validation | `HyperparameterOptimization.jl` |
| Appendix G | External validation | `Validation.jl`: `compare_to_external_sources()` |
| Appendix H | Data description | `data_cleaning.do` |
| Appendix I | Correlations | `Correlations.jl`, `CorrelationTables.jl` |

## References

Bayer, Christian, Luis Calderon, and Moritz Kuhn. "Distributional Dynamics." *Econometrica*, forthcoming.

Board of Governors of the Federal Reserve System. *Survey of Consumer Finances*. Various years (1962-2022). https://www.federalreserve.gov/econres/scfindex.htm.

Board of Governors of the Federal Reserve System. *Distributional Financial Accounts*. https://www.federalreserve.gov/releases/z1/dataviz/dfa/.

Federal Reserve Bank of St. Louis. *FRED Economic Data*. https://fred.stlouisfed.org/.

Panel Study of Income Dynamics. *Public Use Dataset*. Produced and distributed by the Survey Research Center, Institute for Social Research, University of Michigan, Ann Arbor, MI. Various years (1968-2021). https://psidonline.isr.umich.edu/.

Ruggles, Steven, Sarah Flood, Matthew Sobek, et al. *IPUMS USA: Version 14.0* [American Community Survey]. Minneapolis, MN: IPUMS. https://usa.ipums.org/.

Ruggles, Steven, Sarah Flood, Ronald Goeken, et al. *IPUMS CPS: Version 11.0* [Current Population Survey]. Minneapolis, MN: IPUMS. https://cps.ipums.org/.

U.S. Bureau of Labor Statistics. *Consumer Expenditure Surveys*. Various years (1984-2023). https://www.bls.gov/cex/.

U.S. Census Bureau. *Survey of Income and Program Participation*. Various panels (1984-2008). https://www.census.gov/programs-surveys/sipp.html.

World Inequality Database. *WID.world*. https://wid.world/.
