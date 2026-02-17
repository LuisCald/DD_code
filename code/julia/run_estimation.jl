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
#   3. Sets up parameters and priors
#   4. Runs black-box optimization + DIME sampler (MCMC)
#   5. Saves parameter vectors and chains to 7_Results/
#
# ─────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "DistributionalDynamics.jl"))

@info "═══ Stage 2: MCMC Estimation ═══"
@info "BASE_PATH = $BASE_PATH"

# ── Step 1: Prepare data ──────────────────────────────────────
@info "Preparing functional data..."
const func_data, time_params, model_elements = estimation_prep(obs_data, model_options)
@info "Data preparation complete."

# ── Step 2: Set parameters and priors ─────────────────────────
@info "Setting parameters and priors..."
const (param_vector, param_sizes, priors, meas_ind, Σ_ids) = set_params(model_elements, time_params, model_options)
@info "Parameters set. Dimension: $(length(param_vector))"

# ── Step 3: Run DIME sampler (black-box opt + MCMC) ───────────
niter = 400
@info "Starting DIME sampler with $niter iterations..."
par_final = run_DIME_sampler(model_elements, niter, param_vector, param_sizes, priors, meas_ind, Σ_ids, model_options)
@info "Estimation complete."

# ── Verify ────────────────────────────────────────────────────
hyperpriors = priors[end-5:end]
logV, alarm = likeli(model_elements, par_final, param_sizes, hyperpriors, Σ_ids, model_options)
@info "Log-likelihood at estimated parameters: $logV"

@info "═══ Stage 2 finished ═══"
