# To just generate the plots 
# cd("/Users/lc/Dropbox/Distributional_Dynamics/5_Code")
include("DistributionalDynamics.jl")

@unpack measures, data_cutoffs, tag = model_options
println(tag)
label = "3D_A non-diag"
m_label = measures_folder(measures)

const func_data, time_params, model_elements = estimation_prep(obs_data, model_options);
const (param_vector, param_sizes, priors, meas_ind, Σ_ids) = set_params(model_elements, time_params, model_options)
println("finished setting priors")
SSM(par) = -likeli(model_elements, par, param_sizes, priors, meas_ind, Σ_ids, model_options)[1]
SSM(par, _) = -likeli(model_elements, par, param_sizes, priors, meas_ind, Σ_ids, model_options)[1]

# Generate plots 
kind_of_plots = :mcmc
par_final = run_black_box_opt(SSM, param_vector, param_sizes, priors, measures)
println("finished optimization")
# par_final     = get_param_vector(measures, kind_of_plots, label, data_cutoffs, tag) #TODO: adapt code to actual number of chains


niter = param_sizes[1][1] * 100
par_final = run_DIME_sampler(model_elements, niter, param_vector, param_sizes, priors, meas_ind, Σ_ids, model_options)
println("finished optimization")

# min_param_vectors now contains the parameter vector for each chain that minimizes the log probability
# store_optim_estimate(par_final, label, m_label, data_cutoffs, tag)

# Perform Variance Decomposition
smoother_res, logV, alarm = likeli(model_elements, par_final, param_sizes, priors, meas_ind, Σ_ids, model_options; smooth=true)
@unpack x_smoothed, x_filtered = smoother_res               # F̂_t   (nF × T)
@unpack u = model_elements

X_choice = x_smoothed
# --- 1.  VAR(1) residuals and their covariance -----------------------------
# e = X_choice[:, 2:end] .- L * X_choice[:, 1:end-1]     # (r+q) × (T-1)
A, B, D, Ω, _ = matrisize(par_final, param_sizes)
Ω[diagind(Ω)] = log.(exp.(Ω[diagind(Ω)]) .+ 1)
Ω_full = Ω # Diagonal(diag(cov(e; dims = 2)))                               # (r+q) × (r+q)

nx = size(A, 1)                         # = r  (number of factors)
nq = size(D, 1)                         # = q  (idiosyncratic AR processes)

Ω_x = Ω_full[1:nx, 1:nx]     # Var(η^x_t)
Ω_ε = Ω_full[nx+1:end, nx+1:end]     # Var(η^ε_t)
# (cross blocks are zero in most DSGE/FA set-ups; grab them too if not)

# --- 2.  Reduced-form factor VAR -------------------------------------------
C = zeros(nq, nx)                       # still zero under your spec
A_star = A + B * ((I - D) \ C)          # == A here because C == 0
IA = inv(I - A_star)                   # (I - A*)^{-1}

# --- 3.  “Multiplier” from ε-innovations into x_t ---------------------------
G = B * ((I - D) \ I(nq))                 # = B (I‒D)^{-1}

# --- 4.  Unconditional variance of factors and shares ----------------------
var_F = IA * Ω_x * IA' + IA * G * Ω_ε * G' * IA'

part_f = tr(IA * Ω_x * IA')       # contribution of η^x shocks
part_y = tr(IA * G * Ω_ε * G' * IA')    # contribution of η^ε shocks
tot_var = part_f + part_y

share_f = part_f / tot_var
share_y = part_y / tot_var

# New procedure for variance decomposition
@unpack u = model_elements
@unpack case = model_options
smoother_res, logV, alarm = likeli(model_elements, par_final, param_sizes, priors, meas_ind, Σ_ids, model_options; smooth=true)
@unpack x_smoothed, x_filtered = smoother_res               # F̂_t   (nF × T)
@unpack u = model_elements
X_choice = x_filtered

for i in axes(x_filtered, 1)
    Plots.plot(x_filtered[i, :], label="F_$i", xlabel="Time", ylabel="Value", title="Filtered State Variable $i")
    Plots.plot!(x_smoothed[i, :], label="F̂_$i", linestyle=:dash)
    if i > 8
        Plots.plot!(u[i-8, :], label="Y_$i", linestyle=:dot)
    end
    Plots.savefig("filtered_state_variable_$i.pdf")
end


# off_diagonal_indices = setdiff(1:19^2, diagind(A))
# θ_cop[off_diagonal_indices] .= 0  
n = param_sizes[1][1]
b = param_sizes[2][2]
A = reshape(par_final[1:n*n], (n, n))
B = reshape(par_final[n*n+1:n*n+n*b], (n, b))
Ω = Diagonal(par_final[n*n+n*b+1:n*n+n*b+n])
Σ = Diagonal(par_final[end-param_sizes[4][1]+1:end])




@unpack tmin, tmax = time_params
@unpack gdp_series = obs_data
@unpack data_sources = func_data
# user_t     = (Dict("year" => 1990, "quarter" => 1), Dict("year" => 2021, "quarter" => 4))
user_t = (deepcopy(tmin), deepcopy(tmax))
# init_path                               = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
# all_chains = jldopen(init_path * "/2D_income_and_wealth_A non-diag.jld2", "r")
# @unpack par_mcmc, parameter_chain, bounds = all_chains["all_chains"]

opttag = "from_mcmc"
# opttag = "from_optimization"

#TODO: decrease maximum of process eror, mess with aggs 
logV, alarm = likeli(model_elements, par_final, param_sizes, priors, meas_ind, Σ_ids, model_options)

local dv
if tag == "Γ estimated"
    A_new, B_new, Ω_new, Δ_new, G_new, likeli_vec, Δ_log = run_EM_algorithm(param_vector, param_sizes, meas_ind, Σ_ids, model_elements, model_options)
    dv, _ = reconstruct_data_short(A_new, B_new, Ω_new, Δ_new, G_new, likeli_vec, Δ_log, model_elements, obs_data, model_options, time_params, data_sources; reconstruction_to_show=false, dε_smoothed=false)
else
    dv, _ = reconstruct_data(par_final, param_sizes, priors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources)
end

within_stat_dict = Dict()
include("CreateTimeSeries.jl")
include("Validation.jl")
for (c, k) in enumerate(keys(dv))
    within_stat_dict[k] = Dict()
    for ty in ["normal", "average"]
        within_stat_dict[k][ty], dv[k][ty] = export_functional_data(dv[k][ty], ty, k, kind_of_plots, obs_data, func_data, time_params, user_t, model_options, false, true)

        if c == length(keys(dv)) && ty == "normal"
            init_path = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
            dict_path = init_path * "/7_Results/$m_label" * "$tag" * "/$opttag/plots/"
            compare_to_data(dv, ty, func_data, obs_data, user_t, time_params, model_options, kind_of_plots, label)
            export_combined_stat_dict_to_latex(within_stat_dict, [measures..., "copula"], dict_path, label) # [measures..., "copula"]
            compare_to_external_sources(dv, ty, func_data, obs_data, user_t, time_params, model_options, kind_of_plots, label)
        end
    end
end
println("done!")

# Generate Correlations Table 
export_table_to_tex_with_strings(measures, type)

# Generate parameter vectors based on the MCMC chains
# param_mat = vcat(par_final, DIME_chains[end, rand(1:50, 4), :])
# while size(param_mat, 2) < 50
#     param_mat_int = copy(par_final)
#     while size(param_mat_int, 2) < 5
#         pv        = min_param_vectors[rand(1:end), :]
#         if all(isnan.(pv))
#             nothing
#         else
#             param_mat_int = hcat(param_mat_int, pv)
#         end
#     end
#     # Take average 
#     param_mat = hcat(param_mat, mean(param_mat_int, dims=2))
# end


# Generate microdata implicates
# include("CreateTimeSeries.jl")
# generate_microdata_implicates(200, "SCF", param_sizes, priors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources, tag)
# generate_microdata_implicates(200, "CEX", param_sizes, priors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources, tag)
# generate_microdata_implicates(draws, "PSID", param_sizes, priors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources, tag)





# # Generate Consumption plots 
# generate_relative_to_peak_plots()


# Performing the different forecasts 
# how_much = 1 # removes X observation(s) from SCF
# perform_forecast("SCF", par_final, param_sizes, priors, meas_ind, Σ_ids, how_much, obs_data, model_options, model_elements, time_params, user_t, func_data, kind_of_plots, ["iterative"])


# filtering_criteria = Dict("periods" => 20)
filtering_criteria = Dict("dates" => ("2004-Q4", "2009-Q4"))
perform_forecast("SCF", par_final, param_sizes, priors, meas_ind, Σ_ids, filtering_criteria, obs_data, model_options, model_elements, time_params, user_t, func_data, kind_of_plots, ["extensive", "all_data"])

@unpack time_dict = time_params
nd = length(collect(keys(time_dict[1])))
periods_to_remove = setdiff(collect(1:nd), collect(2:4:nd))
cex_years = sort(collect(keys(time_dict[1])))[periods_to_remove]
dvec = []
for i in cex_years
    push!(dvec, QuarterlyDate(i, 4))
end
filtering_criteria = Dict("periods_to_remove" => dvec)

perform_forecast("CEX", par_final, param_sizes, priors, meas_ind, Σ_ids, filtering_criteria, obs_data, model_options, model_elements, time_params, user_t, func_data, kind_of_plots, ["extensive", "data_only"])


