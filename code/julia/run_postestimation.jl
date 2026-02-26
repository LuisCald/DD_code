# ─────────────────────────────────────────────────────────────
# run_postestimation.jl — Stages 3-5: results, figures, tables
# ─────────────────────────────────────────────────────────────
#
# Usage (from repo root):
#   julia --project=code/julia/env code/julia/run_postestimation.jl
#
# This script loads pre-computed parameter estimates and generates
# all post-estimation results, figures, and tables for the paper:
#   - Reconstruction of synthetic distributions
#   - Validation against external sources (DFA, WID, SCF, CPS, ACS)
#   - Forecast error variance decomposition (FEVD)
#   - Historical decomposition
#   - Counterfactual exercises
#   - Correlation tables
#   - Cyclicality of consumption analysis
#   - Out-of-sample forecasts
#
# Requires: completed estimation (Stage 2) with saved parameter
# vectors in 7_Results/.
# ─────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "DistributionalDynamics.jl"))

@info "═══ Stages 3-5: Post-estimation ═══"
@info "BASE_PATH = $BASE_PATH"

# ── Unpack options ────────────────────────────────────────────
@unpack measures, data_cutoffs, tag, case = model_options
label = "3D_$case"
m_label = measures_folder(measures)
kind_of_plots = :mcmc
opttag = "from_mcmc"

# ── Step 1: Prepare data ─────────────────────────────────────
@info "Preparing functional data..."
const func_data, time_params, model_elements = estimation_prep(obs_data, model_options)
const (param_vector, param_sizes, priors, meas_ind, Σ_ids) = set_params(model_elements, time_params, model_options)
hyperpriors = priors[end-5:end]

SSM(par) = -likeli(model_elements, par, param_sizes, hyperpriors, Σ_ids, model_options)[1]
SSM(par, _) = -likeli(model_elements, par, param_sizes, hyperpriors, Σ_ids, model_options)[1]

@unpack tmin, tmax = time_params
@unpack gdp_series = obs_data
@unpack data_sources = func_data
user_t = (deepcopy(tmin), deepcopy(tmax))

# ── Step 2: Load estimated parameters ────────────────────────
@info "Loading estimated parameters..."
# par_final = get_param_vector(measures, kind_of_plots, label, data_cutoffs, tag)
logV, alarm = likeli(model_elements, par_final, param_sizes, hyperpriors, Σ_ids, model_options)
@info "Log-likelihood at loaded parameters: $logV"

# ── Step 3: Reconstruct synthetic distributions ──────────────
@info "Reconstructing synthetic distributions..."
dv, _ = reconstruct_data(
    par_final, param_sizes, hyperpriors, meas_ind, Σ_ids,
    model_elements, obs_data, model_options, time_params, data_sources
)

# ── Step 4: Export time series, validate, generate figures ───
@info "Exporting time series and generating validation plots..."
include("plot_HANK.jl")
within_stat_dict = Dict()
for (c, k) in enumerate(keys(dv))
    if occursin("HANK", tag) && k == "consensus"
        continue
    end

    within_stat_dict[k] = Dict()
    for ty in ["normal"]
        within_stat_dict[k][ty], dv[k][ty] = export_functional_data(
            dv[k][ty], ty, k, kind_of_plots, obs_data, func_data,
            time_params, user_t, model_options, false, true
        )

        if c == length(keys(dv)) && ty == "normal"
            init_path = BASE_PATH
            dict_path = init_path * "/7_Results/$m_label" * "$tag" * "/$opttag/plots/"
            compare_to_data(dv, ty, func_data, obs_data, user_t, time_params, model_options, kind_of_plots, label)
            export_combined_stat_dict_to_latex(within_stat_dict, [measures..., "copula"], dict_path, label)
            if !occursin(" HANK", tag)
                compare_to_external_sources(dv, ty, func_data, obs_data, user_t, time_params, model_options, kind_of_plots, label)
            end
        end
    end
end
@info "Time series export and validation complete."

# ── Step 5: Correlation tables ────────────────────────────────
@info "Generating correlation tables..."
type = "from_mcmc"
export_table_to_tex_with_strings(measures, :from_mcmc)
generate_correlations_table_for_external_comparisons("SCF", measures, tag, type, "cycle")
@info "Correlation tables complete."

# ── Step 6: FEVD ─────────────────────────────────────────────
@info "Computing forecast error variance decomposition..."
A, B_mat, C, D, Ω_var, Ω_corr, Σ = matrisize(par_final[1:end-6], param_sizes)
r_dist = param_sizes[1][1]
q_agg = param_sizes[2][2]

Ω_var[diagind(Ω_var)] = log.(exp.(Ω_var[diagind(Ω_var)]) .+ 1)
mat_Ω_corr = Matrix(Ω_corr)
Ω = Ω_var * mat_Ω_corr * Ω_var'

r, q = size(B_mat)
nₛ = 4r + q
Tval = eltype(A)

Φ = zeros(Tval, nₛ, nₛ)
AI = Matrix{Tval}(I, r, r)
@views begin
    Φ[1:r, 1:r] .= A
    Φ[1:r, 4r+1:4r+q] .= B_mat
    Φ[r+1:2r, 1:r] .= AI
    Φ[2r+1:3r, r+1:2r] .= AI
    Φ[3r+1:4r, 2r+1:3r] .= AI
    Φ[4r+1:end, 1:r] .= C
    Φ[4r+1:end, 4r+1:end] .= D
end

ids_sub = vcat(1:r, 4r+1:4r+q)
Φ_sub = Φ[ids_sub, ids_sub]

tbl = fevd_dist_vs_agg(model_elements, model_options, obs_data, Φ_sub, Ω;
    r, q, horizon=20, factor_names=nothing, shock_order=:agg_first
)
@info "FEVD complete."

# ── Step 7: Historical decomposition ─────────────────────────
@info "Computing historical decomposition..."
smoother_res, _, _ = likeli(model_elements, par_final, param_sizes, hyperpriors, Σ_ids, model_options; smooth=true)
@unpack x_smoothed, dε_smoothed = smoother_res

hd = historical_decomp_factors_blockchol(Φ, r, q, dε_smoothed, x_smoothed, Ω;
    splitting=:group, shock_order=:agg_first
)
@info "Historical decomposition complete. Max reconstruction error: $(hd.maxerr)"

# ── Step 8: Cyclicality of consumption ────────────────────────
@info "Generating cyclicality of consumption analysis..."
try
    generate_relative_to_peak_plots()
    @info "Cyclicality analysis complete."
catch e
    @warn "Cyclicality analysis skipped: $e"
end

# ── Step 9: Out-of-sample forecasts ──────────────────────────
@info "Running out-of-sample forecasts..."
how_much = 1
perform_forecast("SCF", par_final, param_sizes, priors, meas_ind, Σ_ids, how_much,
    obs_data, model_options, model_elements, time_params, user_t, func_data,
    kind_of_plots, ["iterative"]
)

@unpack time_dict, year_vec = time_params
periods_to_remove = muted_quarters_between(QuarterlyDate(1984, 1), QuarterlyDate(2021, 4))
filtering_criteria = Dict("periods_to_remove" => periods_to_remove)
perform_forecast("CEX", par_final, param_sizes, hyperpriors, meas_ind, Σ_ids,
    filtering_criteria, obs_data, model_options, model_elements, time_params,
    user_t, func_data, kind_of_plots, ["extensive", "data_only"]
)
@info "Forecasts complete."

@info "═══ Stages 3-5 finished ═══"
