* ─────────────────────────────────────────────────────────────
* 00_master_data.do — Stage 1: Data cleaning and preparation
* ─────────────────────────────────────────────────────────────
*
* Usage:
*   stata -b do code/stata/00_master_data.do
*
* This master file runs all Stata data processing scripts in
* the correct order. It assumes raw survey data (PSID, SCF,
* CEX, CPS, ACS, SIPP) are in 1_Data/ and outputs processed
* files to 2_Data_processing/.
*
* Prerequisites:
*   - Stata 17+
*   - Raw data downloaded per DOWNLOAD_INSTRUCTIONS.md
* ─────────────────────────────────────────────────────────────

clear all
set more off
set maxvar 32000

* ── Set paths ─────────────────────────────────────────────────
* All paths relative to Distributional_Dynamics/ (project root)
global init_path = "`c(pwd)'"

display "Project root: $init_path"
display "Starting data pipeline at `c(current_date)' `c(current_time)'"

* ── Step 1: Clean SCF 2022 wave ───────────────────────────────
display _n "──── Step 1: Cleaning SCF 2022 ────"
do "${init_path}/5_Code/code/stata/clean_SCF_2022.do"

* ── Step 2: SIPP panel construction ───────────────────────────
display _n "──── Step 2: SIPP panel construction ────"
do "${init_path}/5_Code/code/stata/master_sipp_file.do"
do "${init_path}/5_Code/code/stata/clean_SIPP_panels.do"
do "${init_path}/5_Code/code/stata/sipp_panel_constructor.do"
do "${init_path}/5_Code/code/stata/connect_SIPP_panels.do"

* ── Step 3: Main data cleaning (SCF, PSID, CEX, ACS, CPS) ────
display _n "──── Step 3: Main data cleaning ────"
do "${init_path}/5_Code/code/stata/data_cleaning.do"

* ── Step 4: Import aggregate macro series ─────────────────────
display _n "──── Step 4: Importing aggregates ────"
do "${init_path}/5_Code/code/stata/import_aggregates.do"

* ── Step 5: Seasonal adjustment ───────────────────────────────
display _n "──── Step 5: Seasonal adjustment (X-12) ────"
do "${init_path}/5_Code/code/stata/x12series.do"

* ── Step 6: Process WID quarterly data ────────────────────────
display _n "──── Step 6: Processing WID data ────"
do "${init_path}/5_Code/code/stata/process_WIDq.do"

* ── Step 7: Other results and aggregates ──────────────────────
display _n "──── Step 7: Other results ────"
do "${init_path}/5_Code/code/stata/other_results.do"

display _n "═══ Stage 1 (Stata) finished ═══"
