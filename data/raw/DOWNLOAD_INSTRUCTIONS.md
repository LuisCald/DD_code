# Data Download Instructions

This file describes how to obtain the raw survey data needed for Stage 1 (data cleaning). All datasets are publicly available.

## 1. Panel Study of Income Dynamics (PSID)

- **URL:** https://psidonline.isr.umich.edu/
- **Registration:** Free account required
- **Years:** 1968-2021 (biennial after 1997)
- **Variables needed:** family income, net worth, housing wealth, financial assets, debt, family size, age, education, employment status, interview date
- **Instructions:** Use the PSID Data Center to create a custom extract. The extraction script `J341380.do` (in `1_Data/PSID/`) documents the exact variable list.
- **Output format:** Stata .do extraction file

## 2. Survey of Consumer Finances (SCF)

- **URL:** https://www.federalreserve.gov/econres/scfindex.htm
- **Registration:** None required
- **Years:** 1962-2022 (triennial)
- **Files needed:**
  - Summary extract public data (`.dta` format) for each wave
  - Replicate weight files (`p{year}_rw1.dta`)
  - Historical SCF (`HSCF_2019.dta`) for pre-1989 waves
- **Variables needed:** income, net worth, financial assets, housing, debt, family characteristics

## 3. Consumer Expenditure Survey (CEX)

- **URL:** https://www.bls.gov/cex/pumd_data.htm
- **Registration:** None required
- **Years:** 1984-2023 (quarterly)
- **Files needed:** Interview survey public-use microdata (PUMD)
- **Variables needed:** total expenditures, food, housing, transportation, healthcare, entertainment, income
- **Note:** Pre-processed CEX file (`CEX_processed.csv`) is expected in `2_Data_processing/`

## 4. Current Population Survey (CPS)

- **URL:** https://cps.ipums.org/
- **Registration:** Free IPUMS account required
- **Years:** 1964-2023 (annual March supplement)
- **Extract name:** `cps_00005`
- **Variables needed:** total family income, earnings, family size, age, education, state

## 5. American Community Survey (ACS)

- **URL:** https://usa.ipums.org/
- **Registration:** Free IPUMS account required
- **Years:** 2005-2023
- **Extract name:** `usa_00009`
- **Variables needed:** household income, family size, age, education, state

## 6. Survey of Income and Program Participation (SIPP)

- **URL:** https://www.census.gov/programs-surveys/sipp.html
- **Registration:** None required
- **Panels:** 1984, 1985, 1986, 1987, 1988, 1990, 1991, 1992, 1993, 1996, 2001, 2004, 2008
- **Variables needed:** monthly income, assets, liabilities, demographics
- **Note:** Panel-specific do-files in `code/stata/SIPP/` handle each wave

## 7. World Inequality Database (WID)

- **URL:** https://wid.world/data/
- **Registration:** None required
- **Country:** United States
- **Variables needed:** Pre-tax national income shares, net personal wealth shares (quarterly)
- **Files:** Download as CSV, save as `income Data WID.csv` and `wealth Data WID.csv`

## 8. FRED Economic Data

- **URL:** https://fred.stlouisfed.org/
- **Registration:** None required
- **Series needed:**
  - `CPIAUCSL` (CPI-U, all items)
  - S&P 500 dividend yield
  - University of Michigan household expectations
  - GDP and components
- **Files:** Download as CSV to `2_Data_processing/`

## 9. Distributional Financial Accounts (DFA)

- **URL:** https://www.federalreserve.gov/releases/z1/dataviz/dfa/
- **Registration:** None required
- **Files:** Download wealth distribution by percentile group, save as `DFA.xlsx`
- **Used for:** External validation only (Appendix G)

## Directory Structure After Download

Place downloaded files in the following locations relative to the project root:

```
1_Data/
├── PSID/
│   └── J341380.do          (PSID extraction script)
├── SCF+/
│   ├── HSCF_2019.dta       (Historical SCF)
│   ├── SCF_2022.dta        (2022 wave)
│   ├── SCF_2022_wealth.dta
│   ├── p19_rw1.dta         (2019 replicate weights)
│   └── p22_rw1.dta         (2022 replicate weights)
├── CEX/
│   └── [CEX PUMD files]
├── SIPP/
│   └── [SIPP panel .dta files by year]
├── shocks/
│   └── all_shocks.dta
└── [CPS and ACS IPUMS extracts]
```
