# To just generate the plots 
# cd("/Users/lc/Dropbox/Distributional_Dynamics/5_Code")
include("DistributionalDynamics.jl")

@unpack measures, data_cutoffs, tag = model_options
println(tag)
label = "3D_A non-diag"
m_label = measures_folder(measures)

const func_data, time_params, model_elements = estimation_prep(obs_data, model_options);
println("finished prepping data")
const (param_vector, param_sizes, priors, meas_ind, Σ_ids) = set_params(model_elements, time_params, model_options)
println("finished setting priors")
hyperpriors = priors[end-5:end] #TODO: hardcoded
SSM(par) = -likeli(model_elements, par, param_sizes, hyperpriors, Σ_ids, model_options)[1]
SSM(par, _) = -likeli(model_elements, par, param_sizes, hyperpriors, Σ_ids, model_options)[1]

# Generate plots 
kind_of_plots = :mcmc
par_final = get_param_vector(measures, kind_of_plots, label, data_cutoffs, tag) #TODO: adapt code to actual number of chains

@unpack tmin, tmax = time_params
@unpack gdp_series = obs_data
@unpack data_sources = func_data
user_t = (deepcopy(tmin), deepcopy(tmax))
opttag = "from_mcmc"
logV, alarm = likeli(model_elements, par_final, param_sizes, hyperpriors, Σ_ids, model_options)

local dv
if tag == "Γ estimated"
    A_new, B_new, Ω_new, Δ_new, G_new, likeli_vec, Δ_log = run_EM_algorithm(param_vector, param_sizes, meas_ind, Σ_ids, model_elements, model_options)
    dv, _ = reconstruct_data_short(A_new, B_new, Ω_new, Δ_new, G_new, likeli_vec, Δ_log, model_elements, obs_data, model_options, time_params, data_sources; reconstruction_to_show=false, dε_smoothed=false)
else
    dv, _ = reconstruct_data(par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources)
end

within_stat_dict = Dict()
include("Reconstruction.jl")
include("CreateTimeSeries.jl")
include("Validation.jl")
include("ModelPrep.jl")
include("plot_HANK.jl")
for (c, k) in enumerate(keys(dv))
    if occursin("HANK", tag) && k == "consensus"
        continue
    end

    within_stat_dict[k] = Dict()
    for ty in ["normal"]
        # for ty in ["normal", "average"]
        within_stat_dict[k][ty], dv[k][ty] = export_functional_data(dv[k][ty], ty, k, kind_of_plots, obs_data, func_data, time_params, user_t, model_options, false, true)

        if c == length(keys(dv)) && ty == "normal"
            init_path = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
            dict_path = init_path * "/7_Results/$m_label" * "$tag" * "/$opttag/plots/"
            @info "Comparing to data"
            compare_to_data(dv, ty, func_data, obs_data, user_t, time_params, model_options, kind_of_plots, label)
            export_combined_stat_dict_to_latex(within_stat_dict, [measures..., "copula"], dict_path, label) # [measures..., "copula"]
            if !occursin(" HANK", tag)
                compare_to_external_sources(dv, ty, func_data, obs_data, user_t, time_params, model_options, kind_of_plots, label)
            end
        end
    end
end
println("done!")

# /home/luisc/Distributional_Dynamics/7_Results/illiqd_and_income_and_liquid HANK full 1/from_mcmc/data
@unpack func_dict = func_data
economy_number = split(tag, " ")[end]
jldsave("/home/luisc/Distributional_Dynamics/7_Results/illiqd_and_income_and_liquid$(tag)/from_mcmc/data/HANK_full_$(economy_number).jld2"; data=func_dict)

@load "/home/luisc/Distributional_Dynamics/7_Results/illiqd_and_income_and_liquid HANK full 1/from_mcmc/plots/correlations_HANK full 1.jld2" corr_dict

# how_much = 1 # removes X observation(s) from SCF
# perform_forecast("SCF", par_final, param_sizes,  hyperpriors, meas_ind, Σ_ids, how_much, obs_data, model_options, model_elements, time_params, user_t, func_data, kind_of_plots, ["iterative"])

@unpack time_dict = time_params
@unpack year_vec = time_params
periods_to_remove = muted_quarters_between(QuarterlyDate(1984, 1), QuarterlyDate(2021, 4))
cex_quarters = QuarterlyDate(1984, 1):Quarter(1):QuarterlyDate(2021, 4)
dvec = filter(x -> x ∉ periods_to_remove, cex_quarters)

filtering_criteria = Dict("periods_to_remove" => periods_to_remove)
perform_forecast("CEX", par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, filtering_criteria, obs_data, model_options, model_elements, time_params, user_t, func_data, kind_of_plots, ["extensive", "data_only"])
