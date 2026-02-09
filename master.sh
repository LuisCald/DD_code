#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# master.sh — Replication script for Bayer, Calderon, Kuhn
#              "Distributional Dynamics"
# ─────────────────────────────────────────────────────────────
#
# Usage:
#   bash master.sh [stage]
#
# Stages:
#   data       Stage 1: Clean raw survey data (Stata + Python + R)
#   estimate   Stage 2: MCMC estimation (Julia, ~48-72 hrs)
#   results    Stages 3-5: Post-estimation analysis, figures, tables
#   all        Run everything end-to-end (default)
#
# Requirements:
#   - Julia >= 1.10  (with project environment in code/julia/env/)
#   - Stata >= 17
#   - Python >= 3.10 (with pandas, statsmodels, numpy)
#   - R >= 4.3       (with seasonal, x12)
#
# Notes:
#   - Run from the project root (parent of 5_Code/)
#   - Pre-computed estimates in 7_Results/ allow skipping Stage 2
# ─────────────────────────────────────────────────────────────

set -euo pipefail

# ── Resolve project root ──────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# If script is in 5_Code/, project root is one level up
if [[ "$(basename "$SCRIPT_DIR")" == "5_Code" ]]; then
    PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
else
    PROJECT_ROOT="$SCRIPT_DIR"
fi

CODE_DIR="$PROJECT_ROOT/5_Code"
JULIA_DIR="$CODE_DIR/code/julia"
STATA_DIR="$CODE_DIR/code/stata"
PYTHON_DIR="$CODE_DIR/code/python"
R_DIR="$CODE_DIR/code/R"

echo "════════════════════════════════════════════════════════════"
echo "  Distributional Dynamics — Replication Package"
echo "  Project root: $PROJECT_ROOT"
echo "  Started:      $(date)"
echo "════════════════════════════════════════════════════════════"

STAGE="${1:-all}"

# ── Helper ────────────────────────────────────────────────────
run_stage() {
    echo ""
    echo "────────────────────────────────────────────────────────"
    echo "  $1"
    echo "────────────────────────────────────────────────────────"
}

# ══════════════════════════════════════════════════════════════
# Stage 1: Data cleaning (Stata + Python + R)
# ══════════════════════════════════════════════════════════════
run_data() {
    run_stage "Stage 1: Data cleaning and preparation"

    # 1a. Python preprocessing
    echo "[1a] Python: converting monthly series to quarterly..."
    python3 "$PYTHON_DIR/convert_monthly_to_quarterly.py"

    echo "[1a] Python: generating Forbes 400 data..."
    python3 "$PYTHON_DIR/generateForbes400.py"

    echo "[1a] Python: cleaning Brake data..."
    python3 "$PYTHON_DIR/clean_brake_data.py"

    # 1b. Stata data cleaning
    echo "[1b] Stata: running master data pipeline..."
    cd "$PROJECT_ROOT"
    stata-mp -b do "$STATA_DIR/00_master_data.do"

    # 1c. Python stationarity transformations
    echo "[1c] Python: stationarity transformations..."
    python3 "$PYTHON_DIR/make_stationary.py"
    python3 "$PYTHON_DIR/make_stationary_x12.py"

    # 1d. R scripts
    echo "[1d] R: Hermite series estimation..."
    Rscript "$R_DIR/HermiteSeriesEstimator.R"

    echo "[1d] R: non-parametric copula estimation..."
    Rscript "$R_DIR/NonParametricCopula.R"

    echo "[1d] R: X-12 seasonal adjustment..."
    Rscript "$R_DIR/X12_script.R"

    echo "Stage 1 complete."
}

# ══════════════════════════════════════════════════════════════
# Stage 2: Estimation (Julia MCMC)
# ══════════════════════════════════════════════════════════════
run_estimate() {
    run_stage "Stage 2: MCMC estimation (~48-72 hours)"

    cd "$PROJECT_ROOT"
    julia --project="$JULIA_DIR/env" "$JULIA_DIR/run_estimation.jl"

    echo "Stage 2 complete."
}

# ══════════════════════════════════════════════════════════════
# Stages 3-5: Post-estimation (results, figures, tables)
# ══════════════════════════════════════════════════════════════
run_results() {
    run_stage "Stages 3-5: Post-estimation analysis, figures, and tables"

    cd "$PROJECT_ROOT"
    julia --project="$JULIA_DIR/env" "$JULIA_DIR/run_postestimation.jl"

    # Post-estimation Stata analysis (local projections, Ginis)
    echo "[Post-estimation] Stata: preparing micro data for local projections..."
    stata-mp -b do "$STATA_DIR/prep_micro_data.do"
    stata-mp -b do "$STATA_DIR/prep_micro_3D.do"
    stata-mp -b do "$STATA_DIR/prep_macro_3D.do"
    stata-mp -b do "$STATA_DIR/prep_functional_data.do"
    stata-mp -b do "$STATA_DIR/prep_actualCEX.do"
    stata-mp -b do "$STATA_DIR/ginis_and_consumption_trends.do"

    # Multidimensional inequality (Python)
    echo "[Post-estimation] Python: multidimensional inequality..."
    python3 "$PYTHON_DIR/multidim_inequality.py"
    python3 "$PYTHON_DIR/MEILC_MEGC.py"

    echo "Stages 3-5 complete."
}

# ══════════════════════════════════════════════════════════════
# Dispatch
# ══════════════════════════════════════════════════════════════
case "$STAGE" in
    data)
        run_data
        ;;
    estimate)
        run_estimate
        ;;
    results)
        run_results
        ;;
    all)
        run_data
        run_estimate
        run_results
        ;;
    *)
        echo "Usage: bash master.sh [data|estimate|results|all]"
        echo ""
        echo "  data       Clean raw survey data (Stata + Python + R)"
        echo "  estimate   Run MCMC estimation (Julia, ~48-72 hrs)"
        echo "  results    Post-estimation: analysis, figures, tables"
        echo "  all        Everything end-to-end (default)"
        exit 1
        ;;
esac

echo ""
echo "════════════════════════════════════════════════════════════"
echo "  Finished: $(date)"
echo "════════════════════════════════════════════════════════════"
