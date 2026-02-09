# ─────────────────────────────────────────────────────────────
# run_estimation.jl — Stage 2: MCMC estimation
# ─────────────────────────────────────────────────────────────
#
# Usage (from repo root):
#   julia --project=code/julia/env code/julia/run_estimation.jl
#
# This script:
#   1. Loads all packages and source files via DistributionalDynamics.jl
#   2. Prepares functional data from survey microdata
#   3. Runs black-box optimization to find posterior mode
#   4. Runs MCMC (4 chains) starting from the posterior mode
#   5. Saves parameter vectors and chains to 7_Results/
#
# Runtime: ~48-72 hours depending on hardware.
# ─────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "DistributionalDynamics.jl"))

@info "═══ Stage 2: MCMC Estimation ═══"
@info "BASE_PATH = $BASE_PATH"

# ── Step 1: Prepare data ──────────────────────────────────────
@info "Preparing functional data..."
const func_data, time_params, model_elements = estimation_prep(obs_data, model_options)
@info "Data preparation complete."

# ── Step 2: Run optimization + MCMC ───────────────────────────
@info "Starting optimization and MCMC..."
const smoother_results, dv, all_chains, diagnos, label = SSM_optimize(
    model_elements, model_options, mcmc_options, diagnostics_options,
    obs_data, func_data, time_params
)
@info "Estimation complete. Label: $label"
@info "═══ Stage 2 finished ═══"
