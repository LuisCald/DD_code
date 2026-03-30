# ─────────────────────────────────────────────────────────────
# Just intervals
include(joinpath(@__DIR__, "DistributionalDynamics.jl"))

@info "═══ Stage 2: MCMC Estimation ═══"
@info "BASE_PATH = $BASE_PATH"

# ── Step 1: Prepare data ──────────────────────────────────────
@info "Preparing functional data..."
const func_data, time_params, model_elements = estimation_prep(obs_data, model_options)
@info "Data preparation complete."
