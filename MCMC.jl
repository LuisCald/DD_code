# TODO: Resources - http://www.stat.columbia.edu/~gelman/research/published/taumain.pdf
function SSM_optimize(model_elements, model_options, mcmc_options, diagnostics_options, obs_data, func_data, time_params)
    """Searches for the paramater vector that maximizes the likelihood of our SSM model."""
    # Define initial draw and prior for parameters
    @unpack estimator, estimation_object, measures, case, equivalized, bottom_coded, data_cutoffs, tag = model_options
    @unpack grid_cop, grid_pcf = estimator
    @unpack tmin, tmax = time_params
    @unpack gdp_series, agg_data, df_vec = obs_data
    @unpack data_sources = func_data

    param_vector, param_sizes, priors, meas_ind, Σ_ids = set_params(model_elements, time_params, model_options)

    println(param_sizes)
    dimension = length(measures)
    equiv = equivalized == true ? "eq" : ""
    botcod = !isempty(bottom_coded) ? "bc" : ""
    label = "$dimension" * "D" * "_$case" # * "_$estimation_object" * "_$equiv" * "$botcod"
    m_label = measures_folder(measures)
    ndim = length(param_vector)

    # Define bounds 
    # lb, ub       = construct_bounds(param_vector, param_sizes, case)

    SSM(par) = -likeli(model_elements, par, param_sizes, priors, meas_ind, Σ_ids, model_options)[1]
    SSM(par, _) = -likeli(model_elements, par, param_sizes, priors, meas_ind, Σ_ids, model_options)[1]

    # function SSM(model_elements, param_sizes, priors, meas_ind, model_options)
    #     # Convert matrix into a sequence of vectors 
    #     # vp = [par[:, i] for i in axes(par, 2)]

    #     function _SSM(par)   
    #         return likeli(model_elements, par, param_sizes, priors, meas_ind, model_options)[1]
    #     end
    #     # _SSM.(vp)
    # end

    # logpdf = SSM(model_elements, param_sizes, priors, meas_ind, model_options)

    # # Run MCMC using Gregor's paper 
    # LogProbParallel(x) = pmap(logpdf, eachslice(x, dims=2))  # running an ensemble in parallel (ensemble of chains)

    # nchain             = ndim * 5 # a sane default
    # init_chains        = draw_from_prior(param_sizes, priors, nchain)
    # niter              = 3000

    # # off you go sampling
    # chains, lprobs, propdist = RunDIME(LogProbParallel, init_chains, niter, progress=true, aimh_prob=0.15)
    # return chains, lprobs, propdist


    # Run optimization 
    # par_final      = get_param_vector(measures, :mcmc, label, data_cutoffs, tag) 
    par_final = run_black_box_opt(SSM, param_vector, param_sizes, priors, measures)
    store_optim_estimate(par_final, label, m_label, data_cutoffs, tag)

    # par_final      = run_optimizer(SSM, param_vector, label, m_label, dimension, grid, case)
    posterior_like = SSM(par_final)
    @info("After optimization, the posterior likelihood is $posterior_like")

    # Generate Results from optimization 
    # user_t     = (deepcopy(tmin), deepcopy(tmax))


    # Run MCMC 
    @info("Collecting estimates of unobserved factors at posterior mode ...")
    smoother_results, _ = likeli(model_elements, par_final, param_sizes, priors, meas_ind, Σ_ids, model_options, true)

    @info("Entering MCMC step for $tag...")
    all_chains, mcmc_acceptance_rate = run_MCMC(par_final, model_elements, model_options, mcmc_options, param_sizes, priors, meas_ind, Σ_ids, label, m_label, tag)
    init_path = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
    jldsave(init_path * "/3D_" * m_label * "_$tag.jld2"; all_chains)
    @unpack par_mcmc, parameter_chain, bounds = all_chains  # These are parameter bounds 
    # mcmc_like           = SSM(par_mcmc)

    # likelihoods and diagnostics 
    @info("After optimization, the posterior likelihood is $posterior_like")
    # @info("After MCMCs, the likelihood at the posterior mean is $mcmc_like")

    # Generate results from MCMC
    # dv, _                 = reconstruct_data(par_mcmc, param_sizes, priors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources)
    # @info("Generating reconstruction bounds")
    # local posterior_bounds_dict
    # try
    #     posterior_bounds_dict = generate_reconstruction_bounds(parameter_chain, par_mcmc, param_sizes, priors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources, mcmc_options)
    #     for (c, k) in enumerate(keys(dv))
    #         dv[k] = export_functional_data(dv[k], k, :mcmc, obs_data, func_data, time_params, user_t, model_options, posterior_bounds_dict[k], true)        
    #         if c == length(keys(dv))
    #             compare_to_data(dv, func_data, obs_data, user_t, time_params, model_options, :mcmc, label)
    #             compare_to_external_sources(dv, func_data, obs_data, user_t, time_params, model_options, :mcmc, label, posterior_bounds_dict)
    #         end
    #     end
    # catch ee
    #     @warn("Could not generate reconstruction bounds")
    #     for (c, k) in enumerate(keys(dv))
    #         dv[k] = export_functional_data(dv[k], k, :mcmc, obs_data, func_data, time_params, user_t, model_options, false, true)        
    #         if c == length(keys(dv))
    #             compare_to_data(dv, func_data, user_t, time_params, model_options, :mcmc, label)
    #             compare_to_external_sources(dv, func_data, user_t, time_params, model_options, :mcmc, label)

    #         end
    #     end
    # end

    # for (c, k) in enumerate(keys(dv))
    #     dv[k] = export_functional_data(dv[k], k, :mcmc, obs_data, func_data, time_params, user_t, model_options, posterior_bounds_dict[k])        
    #     if c == length(keys(dv))
    #         compare_to_external_sources(dv, func_data, user_t, time_params, model_options, :mcmc, label, posterior_bounds_dict)
    #         compare_to_data(dv, func_data, user_t, time_params, model_options, :mcmc, label)

    #     end
    # end

    diagnos = run_diagnostics(parameter_chain, mcmc_acceptance_rate, diagnostics_options)
    analyze_mcmc(diagnos)

    return smoother_results, dv, all_chains, diagnos, label
end

function generate_reconstruction_bounds(parameter_chain, par_mcmc, param_sizes, priors, meas_ind, Σ_ids, model_elements, obs_data, estimator, model_options, time_params, data_sources, mcmc_options)
    @unpack chains = mcmc_options
    @unpack measures, number_of_dfs, reconstruction_to_show = model_options
    @unpack grid_pcf, grid_cop = estimator
    @unpack gdp_series = obs_data
    @unpack tmin, tmax = time_params

    # Defines size of sample 
    param_length = size(parameter_chain, 2)
    posterior_dist = vcat([parameter_chain[:, :, i] for i in 1:chains]...)
    n_draws = 2000

    # For percentile functions adjustment 
    gdp_series[!, "date"] = QuarterlyDate.(gdp_series[!, "time"])
    correction_names = [meas * "_per_hh" for meas in measures]
    filter!(row -> row.date >= QuarterlyDate(tmin["year"], tmin["quarter"]), gdp_series)
    filter!(row -> row.date <= QuarterlyDate(tmax["year"], tmax["quarter"]), gdp_series)
    select_series = select(gdp_series, correction_names)

    # Getting additional params 
    X, _ = reconstruct_data(par_mcmc, param_sizes, priors, meas_ind, Σ_ids, model_elements, obs_data, estimator, model_options, time_params, data_sources, reconstruction_to_show=reconstruction_to_show)  # Dictionary of split data 
    M = length(measures)
    d = grid^(M) + grid * M # X["consensus"] is in split data format 
    T = size(X["consensus"][1])[end]
    cop_size = tuple(vcat([grid for i in 1:M], [T])...)

    # Containers 
    recon = zeros(d, T, n_draws)  # container of reconstructions 
    cop = zeros(grid^M, T)
    pcfs = zeros(grid * M, T)
    θ̄s = zeros(1, param_length)
    bounds = Dict()
    ds = Dict() # data samples 

    all_reconstructions = vcat(data_sources, ["consensus"])

    for s in all_reconstructions
        bounds[s] = Dict()
        # Initialize arrays to store the 10th and 90th percentiles
        bounds[s]["lb"] = zeros(d, T)
        bounds[s]["ub"] = zeros(d, T)

        # Main container, data samples 
        ds[s] = Vector{AbstractArray}(undef, n_draws)
    end

    @info("For the point confidence intervals, we have $(n_draws) samples")
    # Generating the mean of each sample 
    for draw in 1:n_draws
        psample = view(posterior_dist, sample(1:size(posterior_dist, 1), 3 * param_length, replace=true, ordered=true), :)
        θ̄s .= mean(psample, dims=1)
        reconstructed, _ = reconstruct_data(θ̄s, param_sizes, priors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources, reconstruction_to_show="SCF")

        for data_source in all_reconstructions
            cop .= reshape(reconstructed[data_source][1], (grid^M, T))
            pcfs .= reconstructed[data_source][2]
            split_pcfs = [pcfs[I, :] for I in Iterators.partition(axes(pcfs, 1), grid)]
            for t in 1:T
                pcfs[:, t] .= vcat([split_pcfs[m][:, t] .* select_series[t, correction_names[m]] for m in eachindex(split_pcfs)]...)
                # pcfs[:, t]        .= pcfs[:, t] .* gdp_series[t, :real_gdp_pc]
            end

            ds[data_source][draw] = vcat(cop, pcfs)
        end
    end

    # Generate the data samples, should be 240 by T by length(samples) for 2D case 
    for data_source in all_reconstructions
        recon .= reshape(collect(Iterators.flatten(ds[data_source])), (d, T, n_draws))

        # Calculate the lower and upper bound for each element. This returns a lower bound dataset and an upper bound dataset
        for i in axes(recon, 1)
            for j in axes(recon, 2)
                bounds[data_source]["lb"][i, j] = quantile(recon[i, j, :], 0.025)[1, 1]
                bounds[data_source]["ub"][i, j] = quantile(recon[i, j, :], 0.975)[1, 1]
            end
        end
    end

    # Reformat Bounds to have copulas and percentile functions separate
    posterior_bounds_dict = Dict()
    for data_source in all_reconstructions
        posterior_bounds_dict[data_source] = Dict()
        for bound in collect(keys(bounds[data_source]))
            pcfs = bounds[data_source][bound][(grid*grid+1):end, :]
            copulas = reshape(bounds[data_source][bound][1:(grid*grid), :], cop_size)  # TODO: needs to be generalized to 3D case 
            levels, shares = generate_shares_levels(pcfs, model_options, gdp_series)
            posterior_bounds_dict[data_source][bound] = create_time_series_dictionary([copulas, pcfs, levels, shares], grid, measures)
        end
    end

    return posterior_bounds_dict
end


function draw_from_prior(param_sizes, priors, nchain)
    nf = param_sizes[1][1]  # number of factors
    ny = param_sizes[2][2]  # number of controls
    nΣ = param_sizes[6][1]

    # Draw from prior for state equation + shocks
    A_B_C_D = rand(priors[1], nchain)  # a matrix 
    Ωf = hcat([rand(priors[2], nf) for _ in 1:nchain]...)
    Ωy = hcat([rand(priors[3], ny) for _ in 1:nchain]...)
    Σ = hcat([rand.(priors[4:3+nΣ]) for _ in 1:nchain]...)
    θ = vcat(A_B_C_D, Ωf, Ωy, Σ)
    return θ
end

# function generate_all_plots(dv, type, func_data, gdp_series, time_params, user_time_params, label, model_options, data_bounds=false)
#     @unpack estimation_object, measures, grid = model_options
#     @unpack tmin, tmax = time_params  
#     D                  = length(measures)


#     if estimation_object == "levels and percentile functions"
#         # Organize data into their topologies 
#         for k in collect(keys(dv))
#             levels, quantiles, shares = break_univariate_data(dv[k], gdp_series, time_params, model_options)
#             # Break data into levels and quantiles 
#             # split_data = [vcat([dv[k][i] for i in 1:D]...), vcat([dv[k][i] for i in D+1:2*D]...)]
#             # g, T       = size(split_data[1])  # dimensions of one measure

#             # # Objects of interest 
#             # levels     = [split_data[1][I,:] for I in Iterators.partition(Base.OneTo(D*grid), grid)]  # split by measure 
#             # quantiles  = [split_data[2][I,:] for I in Iterators.partition(Base.OneTo(D*grid), grid)] 
#             # shares     = [zeros(grid, T) for _ in 1:D]

#             # corrected_levels, _, _ = break_univariate_data(dv[k], gdp_series, time_params, model_options)  # N x T x M 

#             # # This will be computed in the order of the measures 
#             # for i in eachindex(shares)
#             #     for t in 1:T
#             #         # First transform into  
#             #         shares[i][:, t] = corrected_levels[i][:,t] ./ sum(corrected_levels[i][:,t])
#             #     end
#             # end
#             dv[k]                     = create_time_series_dictionary([levels, quantiles, shares], model_options)
#         end

#         if data_bounds != false 
#             for k in collect(keys(data_bounds))
#                 for b in ["lb", "ub"]
#                     levels, quantiles, shares = break_univariate_data(data_bounds[k][b], gdp_series, time_params, model_options)
#                     data_bounds[k][b]         = create_time_series_dictionary([levels, quantiles, shares], model_options)
#                 end
#             end
#             generate_specific_plots(dv["consensus"], func_data, time_params, user_time_params, model_options, type, data_bounds["consensus"])
#         else
#             generate_specific_plots(dv["consensus"], func_data, time_params, user_time_params, model_options, type)
#         end
#         # Plot the reconstructions, also doing so against the WID 
#         compare_to_WID(dv, func_data, user_time_params, time_params, model_options, type, label, data_bounds)


#     elseif estimation_object == "copulas and percentile functions"
#         # No need to break here and no comparison to the WID for obvious reasons 
#         for k in collect(keys(dv))
#             dv[k]                   = create_time_series_dictionary(dv[k], model_options)
#         end

#         if data_bounds != false 
#             for k in collect(keys(data_bounds))
#                 for b in ["lb", "ub"]
#                     data_bounds[k][b]            = create_time_series_dictionary(data_bounds[k][b], model_options)
#                 end
#             end
#             generate_specific_plots(dv["consensus"], func_data, time_params, user_time_params, model_options, type, data_bounds["consensus"])
#         else
#             generate_specific_plots(dv["consensus"], func_data, time_params, user_time_params, model_options, type)
#         end

#     end
# end


function measures_folder(measures)
    label = ""
    for (m, meas) in enumerate(sort(measures))
        if m < length(measures)
            label = label * meas * "_and_"
        else
            label = label * meas
        end
    end
    return label
end

function data_tag(sources)
    label = ""

    new_sources = deepcopy(sources)

    # Change 'SCF2016' in sources to 'SCF'
    for (m, source) in enumerate(new_sources)
        if source == "SCF2016"
            new_sources[m] = "SCF"
        end
    end

    for (m, source) in enumerate(sort(new_sources))
        if m < length(new_sources)
            label = label * source * "_and_"
        else
            label = label * source
        end
    end
    return label
end


function run_optimizer(SSM, param_vector, label, m_label, dimension, grid, case)
    local par_final
    if grid < 10 || dimension <= 2 || case == "diag"
        par_final = run_small_optimizer(SSM, param_vector, label, m_label)

    elseif grid >= 10 || dimension >= 3
        par_final = run_large_optimizer(SSM, param_vector, dimension, label, m_label)

    end
    return par_final
end

function get_initializer()
    init_path = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
    other_path = init_path * "/7_Results/income_and_wealth/from_optimization/solution2D_diag_shares and levels_eqbc.csv"
    param_vector = vec(Matrix(CSV.read(other_path, DataFrame, header=0)))
    return param_vector
end

function run_small_optimizer(SSM, param_vector, label, m_label)
    OptOpt = Optim.Options(
        x_tol=0.0,
        f_tol=0.0,
        g_abstol=0.0,
        g_reltol=0.0,
        iterations=300000,
        store_trace=true,
        show_trace=true,
        show_every=100000
    )
    # Optimize 
    @info("Entering optimization step for $label ...")
    @info("Optimizing via NelderMead ...")
    opti = optimize(SSM, param_vector, NelderMead(), OptOpt)
    par_final = Optim.minimizer(opti)


    # @info("Entering optimization step for $label ...")
    # @info("Shortly Optimizing via Nelder Mead to get a nice initial point ...")
    # # OptOpt           = Optim.Options(
    #                             x_tol = 0.0,
    #                             f_tol = 0.0,
    #                             g_abstol = 0.0, 
    #                             g_reltol = 0.0,
    #                             iterations = 100, 
    #                             store_trace = true,
    #                             show_trace = true,
    #                             show_every = 100 
    #                             )
    # opti             = optimize(SSM, param_vector, LBFGS(), autodiff = :forward, OptOpt)
    # f_par_final      = Optim.minimizer(opti)


    # @info("Optimizing via reverse auto-differentiation for speedy AD, LBFGS ...")
    # OptOpt       = Optim.Options(
    #                                 x_tol = 0.0,
    #                                 f_tol = 0.0,
    #                                 g_abstol = 0.0, 
    #                                 g_reltol = 0.0,
    #                                 iterations = 2000, # 2000
    #                                 store_trace = true,
    #                                 show_trace = true,
    #                                 show_every = 100
    #                                 )  
    # # # Using ReverseDiff (sensitive to initial values)
    # tape     = ReverseDiff.GradientTape(SSM, param_vector)
    # g!(G, x) = ReverseDiff.gradient!(G, tape, x)

    # opti             = optimize(SSM, g!, param_vector, LBFGS(linesearch=LineSearches.BackTracking()), OptOpt)  # doesn't work because log() of a negative number 
    # par_final        = Optim.minimizer(opti)

    # # Nelder-Mead 
    # OptOpt       = Optim.Options(
    #                             x_tol = 0.0,
    #                             f_tol = 0.0,
    #                             g_abstol = 0.0, 
    #                             g_reltol = 0.0,
    #                             iterations = 30000, 
    #                             store_trace = true,
    #                             show_trace = true,
    #                             show_every = 1000 
    #                             )

    # # Optimize 
    # @info("Optimizing via NelderMead to finish ...")
    # opti            = optimize(SSM, r_par_final, NelderMead(), OptOpt)
    # par_final       = Optim.minimizer(opti)

    # @info("Optimizing via reverse auto-differentiation for speedy AD, LBFGS ...")
    # OptOpt       = Optim.Options(
    #                                 x_tol = 0.0,
    #                                 f_tol = 0.0,
    #                                 g_abstol = 0.0, 
    #                                 g_reltol = 0.0,
    #                                 iterations = 200, 
    #                                 store_trace = true,
    #                                 show_trace = true,
    #                                 show_every = 100
    #                                 )  
    # # # Using ReverseDiff (sensitive to initial values)
    # tape     = ReverseDiff.GradientTape(SSM, par_final)
    # g!(G, x) = ReverseDiff.gradient!(G, tape, x)
    # # g!(x) = ReverseDiff.gradient(SSM, x)
    # opti             = optimize(SSM, g!, par_final, LBFGS(linesearch=LineSearches.BackTracking()), OptOpt)  # doesn't work because log() of a negative number 
    # par_final      = Optim.minimizer(opti)
    # posterior_like  = SSM(par_final)


    # if posterior_like > 1e7
    #     # Nelder-Mead 
    #     OptOpt       = Optim.Options(
    #                         x_tol = 0.0,
    #                         f_tol = 0.0,
    #                         g_abstol = 0.0, 
    #                         g_reltol = 0.0,
    #                         iterations = 10000, # 10000
    #                         store_trace = true,
    #                         show_trace = true,
    #                         show_every = 1000 
    #                         )
    #     # Optimize 
    #     @info("Optimizing via NelderMead again ...")
    #     opti         = optimize(SSM, nd_par_final, NelderMead(), OptOpt)
    #     par_final    = Optim.minimizer(opti)

    # # Reverse-Differentiation 
    # @info("Optimizing via reverse auto-differentiation, LBFGS ...")
    # OptOpt       = Optim.Options(
    #                                 x_tol = 0.0,
    #                                 f_tol = 0.0,
    #                                 g_abstol = 0.0, 
    #                                 g_reltol = 0.0,
    #                                 iterations = 500, 
    #                                 store_trace = true,
    #                                 show_trace = true,
    #                                 show_every = 100
    #                                 )  
    # # # Using ReverseDiff (sensitive to initial values)
    # tape     = ReverseDiff.GradientTape(SSM, par_final)
    # g!(G, x) = ReverseDiff.gradient!(G, tape, x)
    # # g!(x) = ReverseDiff.gradient(SSM, x)
    # opti           = optimize(SSM, g!, par_final, LBFGS(linesearch=LineSearches.BackTracking()), OptOpt)  # doesn't work because log() of a negative number 
    # par_final      = Optim.minimizer(opti)
    #     posterior_like = SSM(par_final)
    # end
    # f    = OptimizationFunction(SSM, Optimization.AutoReverseDiff())
    # prob = OptimizationProblem(f, par_final)  
    # sol  = solve(prob, LBFGS(linesearch=LineSearches.BackTracking()))
    # par_final    = sol.u

    # opti         = optimize(SSM, param_vector, LBFGS(), autodiff = :forward, OptOpt)
    store_optim_estimate(par_final, label, m_label)
    return par_final
end

function run_large_optimizer(SSM, param_vector, dimension, label, m_label)
    @info("Entering optimization step for $label ...")
    @info("Shortly Optimizing via Nelder-Mead to get a nice initial point ...")

    # local OptOpt
    # if dimension == 2
    OptOpt = Optim.Options(
        x_tol=0.0,
        f_tol=0.0,
        g_abstol=0.0,
        g_reltol=0.0,
        iterations=300000, #300000, 
        store_trace=true,
        show_trace=true,
        show_every=10000
    )
    # elseif dimension == 3
    #     OptOpt           = Optim.Options(
    #                                 x_tol = 0.0,
    #                                 f_tol = 0.0,
    #                                 g_abstol = 0.0, 
    #                                 g_reltol = 0.0,
    #                                 iterations = 300000, #300000, 
    #                                 store_trace = true,
    #                                 show_trace = true,
    #                                 show_every = 100000
    #                                 )
    # end

    # Optimize 
    opti = optimize(SSM, param_vector, NelderMead(), OptOpt)
    # nd_par_final    = Optim.minimizer(opti)
    par_final = Optim.minimizer(opti)
    posterior_like = SSM(par_final)

    # Optimize 
    # @info("Try to optimize via reverse auto-differentiation, LBFGS ...")
    # OptOpt       = Optim.Options(
    #                                 x_tol = 0.0,
    #                                 f_tol = 0.0,
    #                                 g_abstol = 0.0, 
    #                                 g_reltol = 0.0,
    #                                 iterations = 2000, #2000
    #                                 store_trace = true,
    #                                 show_trace = true,
    #                                 show_every = 100
    #                                 )  
    # # Using ReverseDiff (sensitive to initial values)
    # # local par_final, posterior_like
    # # try 
    #     tape          = ReverseDiff.GradientTape(SSM, param_vector)
    #     # compiled_tape = ReverseDiff.compile(tape)
    #     g!(G, x)      = ReverseDiff.gradient!(G, tape, x)

    #     # g!(x) = ReverseDiff.gradient(SSM, x)
    #     opti             = optimize(SSM, g!, param_vector, LBFGS(linesearch=LineSearches.BackTracking()), OptOpt)  # doesn't work because log() of a negative number 
    #     par_final        = Optim.minimizer(opti)
    #     posterior_like   = SSM(par_final)

    #     # OptOpt       = Optim.Options(
    #     #     x_tol = 0.0,
    #     #     f_tol = 0.0,
    #     #     g_abstol = 0.0, 
    #     #     g_reltol = 0.0,
    #     #     iterations = 50000, 
    #     #     store_trace = true,
    #     #     show_trace = true,
    #     #     show_every = 1000
    #     #     )  

    #     # opti             = optimize(SSM, par_final, NelderMead(), OptOpt)
    #     # par_final        = Optim.minimizer(opti)
    #     # posterior_like   = SSM(par_final)
    # catch e 
    #     println(e)
    #     @info("Trying again Reverse differentiation")
    #     try 
    #         tape     = ReverseDiff.GradientTape(SSM, param_vector)
    #         g!(G, x) = ReverseDiff.gradient!(G, tape, x)
    #         # tape          = ReverseDiff.GradientTape(SSM, param_vector)
    #         # compiled_tape = ReverseDiff.compile(tape)
    #         # g!(G, x)      = ReverseDiff.gradient!(G, compiled_tape, x)

    #         # g!(x) = ReverseDiff.gradient(SSM, x)
    #         opti             = optimize(SSM, g!, param_vector, LBFGS(), OptOpt)  # doesn't work because log() of a negative number 
    #         par_final        = Optim.minimizer(opti)
    #         posterior_like   = SSM(par_final)
    #     catch ee
    #         @info("Trying Forward differentiation")
    #         try 
    #             # tape     = ReverseDiff.GradientTape(SSM, param_vector)
    #             # g!(G, x) = ReverseDiff.gradient!(G, tape, x)
    #             # tape          = ReverseDiff.GradientTape(SSM, param_vector)
    #             # compiled_tape = ReverseDiff.compile(tape)
    #             # g!(G, x)      = ReverseDiff.gradient!(G, compiled_tape, x)

    #             # g!(x) = ReverseDiff.gradient(SSM, x)
    #             opti             = optimize(SSM, param_vector, LBFGS(), autodiff = :forward, OptOpt)  # doesn't work because log() of a negative number 
    #             par_final        = Optim.minimizer(opti)
    #             posterior_like   = SSM(par_final)
    #         catch eee
    #             println(eee)
    #             @info("Optimizing via NelderMead ...")
    #             OptOpt       = Optim.Options(
    #                 x_tol = 0.0,
    #                 f_tol = 0.0,
    #                 g_abstol = 0.0, 
    #                 g_reltol = 0.0,
    #                 iterations = 30000, 
    #                 store_trace = true,
    #                 show_trace = true,
    #                 show_every = 1000
    #                 )  
    #             local opti
    #             try 
    #                 @info("Trying again NelderMead with par_final")
    #                 opti             = optimize(SSM, par_final, NelderMead(), OptOpt)
    #             catch e
    #                 @info("Trying again NelderMead with paramvector")
    #                 opti             = optimize(SSM, param_vector, NelderMead(), OptOpt)
    #             end 
    #             par_final        = Optim.minimizer(opti)
    #             posterior_like   = SSM(par_final)
    #         end
    #     end
    # end
    #     # Trying again with Nelder-Mead 
    #     # @info("Optimization via reverse auto-differentiation, LBFGS failed. Trying NelderMead again ...")
    #     # # Optimize
    #     # OptOpt       = Optim.Options(
    #     #     x_tol = 0.0,
    #     #     f_tol = 0.0,
    #     #     g_abstol = 0.0, 
    #     #     g_reltol = 0.0,
    #     #     iterations = 30000, 
    #     #     store_trace = true,
    #     #     show_trace = true,
    #     #     show_every = 1000
    #     #     )  

    #     # # opti             = optimize(SSM, nd_par_final, LBFGS(), autodiff = :forward, OptOpt)
    #     # opti             = optimize(SSM, nd_par_final, NelderMead(), OptOpt)
    #     # par_final        = Optim.minimizer(opti)
    #     # posterior_like   = SSM(par_final)
    # # end
    # @info("The posterior likelihood is     $posterior_like     ")

    # if posterior_like > 1e7
    #     @info("Optimizing via NelderMead again ...")
    #     OptOpt       = Optim.Options(
    #         x_tol = 0.0,
    #         f_tol = 0.0,
    #         g_abstol = 0.0, 
    #         g_reltol = 0.0,
    #         iterations = 50000, 
    #         store_trace = true,
    #         show_trace = true,
    #         show_every = 1000
    #         )  

    #     opti         = optimize(SSM, param_vector, NelderMead(), OptOpt)
    #     par_final    = Optim.minimizer(opti)
    # end

    store_optim_estimate(par_final, label, m_label)
    return par_final
end

function construct_bounds(param_vector, param_sizes, case)
    local lb, ub
    if case == "diag"
        A_length = param_sizes[1][1]
        B_length = param_sizes[2][1] * param_sizes[2][2]
        Ω_length = param_sizes[3][1]
        Σ_length = param_sizes[4][1]
        lb = vcat(repeat([-Inf], A_length + B_length), repeat([0], Ω_length + Σ_length))
        ub = repeat([Inf], length(param_vector))

    elseif case == "A non-diag"
        A_length = param_sizes[1][1] * param_sizes[1][2]
        B_length = param_sizes[2][1] * param_sizes[2][2]
        Ω_length = param_sizes[3][1]
        Σ_length = param_sizes[4][1]
        lb = vcat(repeat([-Inf], A_length + B_length), repeat([0], Ω_length + Σ_length))
        ub = repeat([Inf], length(param_vector))
    elseif case == "A, Σ non-diag"
        A_length = param_sizes[1][1] * param_sizes[1][2]
        B_length = param_sizes[2][1] * param_sizes[2][2]
        Ω_length = param_sizes[3][1]
        Σ_length = param_sizes[4][1] * param_sizes[4][2]
        lb = vcat(repeat([-Inf], A_length + B_length), repeat([0], Ω_length + Σ_length))
        ub = repeat([Inf], length(param_vector))
    end
    return lb, ub
end

function get_hessian(label, m_label, tag)
    init_path = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
    inv_hessian = Matrix(CSV.read(init_path * "/7_Results/$m_label" * "$tag" * "/from_mcmc/hessians/hessian_$label.csv", header=false, DataFrame))
    return inv_hessian
end


function run_MCMC(par_final, model_elements, model_options, mcmc_options, param_sizes, priors, meas_ind, Σ_ids, label, m_label, tag)
    # Compute hessian and draws for MCMC 
    # println(Dates.format(now(), "HH:MM"))
    # inv_hessian            = get_hessian(label, m_label, tag) 
    println(Dates.format(now(), "HH:MM"))
    inv_hessian = construct_hessian(par_final, model_elements, param_sizes, priors, meas_ind, Σ_ids, model_options, mcmc_options)
    store_hessian(inv_hessian, label, m_label, tag)
    println(Dates.format(now(), "HH:MM"))
    initial_draws, sym_hes = multi_chain_init(model_elements, par_final, param_sizes, priors, meas_ind, Σ_ids, model_options, mcmc_options, inv_hessian)

    # Perform MCMC 
    # smoother_results, _                   = likeli(model_elements, par_final, param_sizes, priors, meas_ind, model_options, true) 
    chains, mcmc_acceptance_rate = sampler(model_elements, initial_draws, sym_hes, priors, model_options, mcmc_options, param_sizes, meas_ind, Σ_ids)
    println(Dates.format(now(), "HH:MM"))
    @unpack par_mcmc = chains

    # Store results and diagnostics 
    @info("Finished sampling. Storing posterior means ...")
    store_posterior_mean(par_mcmc, label, m_label, tag)

    return chains, mcmc_acceptance_rate
end

# elseif case == "A non-diag"   #&& estimation_object == "copulas and percentile functions"
#     @info("Entering constrained gradient-based optimization step ...") 
#     # Define State-Space Model function 

#     A_length = param_sizes[1][1] * param_sizes[1][2]
#     B_length = param_sizes[2][1] * param_sizes[2][2]
#     lb       = vcat(repeat([-Inf], A_length + B_length), repeat([0], param_sizes[3][1] + param_sizes[4][1]))
#     ub       = repeat([Inf], length(param_vector))

#     # Specifying some options 
#     OptOpt       = Optim.Options(
#                                 x_tol = 1.0e-6,
#                                 f_tol = 1.0e-6, 
#                                 iterations = 150000, 
#                                 store_trace = true,
#                                 show_trace = true,
#                                 show_every = 100
#                                 )
#     # opti         = optimize(SSM, lb, ub, param_vector, Fminbox(LBFGS()), autodiff = :forward, OptOpt)
#     opti         = optimize(SSM, param_vector,  NelderMead(), OptOpt)
#     # opti         = optimize(SSM, param_vector, SimulatedAnnealing(), OptOpt)

#     par_final    = Optim.minimizer(opti)
#     DelimitedFiles.writedlm("gradient_solution_$label.csv",  par_final, ',')  # save for future runs 
# end
# Posterior likelihood at posterior mode, par_final 

# For large parameter spaces, reverse diff is better 
# tape     = ReverseDiff.GradientTape(SSM, repeat([.01], length(param_vector)))
# g!(G, x) = ReverseDiff.gradient!(G, tape, x)            
# opti     = optimize(SSM, g!, lb, ub, param_vector, Fminbox(LBFGS()), OptOpt)
# f    = OptimizationFunction(SSM, Optimization.AutoReverseDiff())
# prob = OptimizationProblem(f, param_vector, lb=lb, ub=ub)  # the claim is that it doesnt work with constraints yet, but that's not the issue. 
# # prob = OptimizationProblem(f, param_vector)
# sol  = solve(prob, BFGS())
# par_final    = sol.u


function construct_hessian(par_final, model_elements, param_sizes, priors, meas_ind, Σ_ids, model_options, mcmc_options)
    @unpack compute_hessian = mcmc_options

    local inv_hessian
    SSM(par) = -likeli(model_elements, par, param_sizes, priors, meas_ind, Σ_ids, model_options)[1]
    SSM(par, _) = -likeli(model_elements, par, param_sizes, priors, meas_ind, Σ_ids, model_options)[1]

    if compute_hessian == true
        @info("Computing Hessian ...")
        try
            @info("Trying first with Reverse Differentiation ...")
            tape = ReverseDiff.HessianTape(SSM, par_final)
            inv_hessian = inv(ReverseDiff.hessian(SSM, tape, par_final))
        catch e
            try
                @info("Trying now with other derivative methods ...")
                tape = ReverseDiff.compile(ReverseDiff.GradientTape(SSM, par_final))
                f_grad = ReverseDiff.gradient!(SSM, tape, par_final)
                inv_hessian = inv(ForwardDiff.jacobian(f_grad, par_final))
            catch ee
                @info("Computing Hessian failed with Reverse Diff. Trying now with Finite Differences ...")
                # ∂²f∂p∂p     = TwiceDifferentiable(pp -> SSM(pp), par_final, autodiff =:forward)
                # # old
                n = length(par_final)
                ∂²f∂p∂p = zeros(Float64, n, n)
                # buffer  = Array(Float64, n)
                # FiniteDiff.hessian!(∂²f∂p∂p, pp -> SMM(pp), par_final, buffer)  
                # inv_hessian = inv(∂²f∂p∂p)

                # new
                FiniteDiff.finite_difference_hessian!(∂²f∂p∂p, pp -> SSM(pp), par_final)
                inv_hessian = inv(∂²f∂p∂p)
                # ∂²f∂p∂p     = TwiceDifferentiable(pp -> SSM(pp), par_final)
                # inv_hessian = inv(Optim.hessian!(∂²f∂p∂p, par_final))  # When the negative log-likelihood is minimized, the negative Hessian is returned. So this is the (-H)⁻¹            
            end
        end
    else
        @info("Hessian taken as the identity matrix ...")
        inv_hessian = Matrix{Float64}(I, length(par_final), length(par_final))
    end
    return inv_hessian
end


function sampler(model_elements, initial_draws, inv_hessian, priors, model_options, mcmc_options, param_sizes, meas_ind, Σ_ids)
    @unpack sampler = mcmc_options
    @info("Beginning to sample from the posterior ... The sampler of choice is $sampler")
    #    if sampler == "MH"
    return rwmh(model_elements, initial_draws, inv_hessian, priors, param_sizes, model_options, mcmc_options, meas_ind, Σ_ids)

    #    elseif sampler == "BAT sampler"
    #         return BAT_rwmh(model_elements, initial_draws, inv_hessian, priors, param_sizes, model_options, mcmc_options, meas_ind)

    #    elseif sampler == "Gibbs"
    #         return gibbs_sampler(model_elements, smoother_results, initial_draws, priors, model_options, param_sizes)
    #     end
end

function store_optim_estimate(params, label, m_label, data_cutoffs, tag)
    end_year = data_cutoffs["end"] != "" ? data_cutoffs["end"][1:4] : "all"
    init_path = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()

    DelimitedFiles.writedlm(init_path * "/7_Results/$m_label" * "$tag" * "/from_mcmc/parameter_vectors/solution" * "$label" * "_$end_year" * ".csv", params, ',')
end

function store_posterior_mean(posterior_mean, label, m_label, tag)
    for j in eachindex(posterior_mean)
        init_path = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
        DelimitedFiles.writedlm(init_path * "/7_Results/$m_label" * "$tag" * "/from_mcmc/parameter_vectors/posterior_mean_" * "$label" * "$j.csv", posterior_mean[j], ',')
    end
end

function store_hessian(hessian, label, m_label, tag)
    init_path = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
    DelimitedFiles.writedlm(init_path * "/7_Results/$m_label" * "$tag" * "/from_mcmc/hessians/hessian_" * "$label.csv", hessian, ',')
end

scale_hessian(inv_hessian, scale) = scale .* inv_hessian
make_symmetric(mat) = (mat .+ mat) ./ 2

# TODO: We can also optimize in parallel, take the different posterior modes as the initializers. (https://www.imprs-tp.mpg.de/71348/Ballnus_al17.pdf)
function multi_chain_init(model_elements, param_vector, param_sizes, priors, meas_ind, Σ_ids, model_options, mcmc_options, inv_hessian=false)
    """Initializes multiple independent runs of the MCMC i.e., chains. 

    We perturb the final set of parameters, check that indeed they are in the support and initialize the MCMC algorithm. Once for each chain. 

    Purpose: to diagnose multimodality. In theory, starting do not matter, but it does matter for speed. Hence, we use this for both methods. 
    """

    @unpack sampler, mhscale, chains = mcmc_options
    @info("Initializing $chains chains")

    # init_scale = mhscale^2
    # scaled_hes = scale_hessian(inv_hessian, init_scale)
    sym_hes = make_symmetric(inv_hessian)
    # if isposdef(sym_hes) == false 
    #     @info("finding nearest spd of hessian.")
    #     sym_hes .= nearest_spd(sym_hes)
    # end

    prop_dist = MvNormal(zeros(length(param_vector)), Matrix(I, size(sym_hes, 1), size(sym_hes, 1)))
    # prop_dist  = MvTDist(length(param_vector) - 1, zeros(length(param_vector)), Matrix(I, size(sym_hes, 1), size(sym_hes, 1)))
    # objects for the struct 
    draws = Vector{Vector{Float64}}(undef, chains)
    draws[1] = copy(param_vector)
    draw = zeros(length(param_vector))
    scale = sym_hes   #Matrix(cholesky(sym_hes).L) .* mhscale^2 

    # For each chain, get a vector of starting values. This vector will consist of 3 blocks 
    init_success = 1
    init_iter = 1
    while init_success < chains
        if init_iter >= 100 && mod(init_iter, 50) == 0
            scale .= scale .* mhscale
        elseif init_iter > 250 && mod(init_iter, 10) == 0
            @info("at $init_iter")
            scale .= scale .* mhscale
        end
        draw .= scale * rand(prop_dist) .+ param_vector
        _, alarm = likeli(model_elements, draw, param_sizes, priors, meas_ind, Σ_ids, model_options)

        # Check if draw in support and/or generated a non-weird value for the likelihood 
        if alarm == false
            init_success += 1
            init_iter += 1
            draws[init_success] = draw
        else
            init_iter += 1
        end
    end
    return draws, sym_hes
end




"""
Runs Gibbs Sampler, sampling the conditional posterior of the blocks (matrices or group of parameters). 
"""
function gibbs_sampler(model_elements, smoother_results, initial_draws, priors, model_options, param_sizes)

    @unpack x_filtered = smoother_results
    @unpack agg_count, factor_count, chains, mcmc_jsd_draws, case, lags, constant, thinning_steps = model_options
    @unpack MV, G = model_elements

    # Dimensions
    T = size(x_filtered, 2)
    n_draws = nsave + nburn
    state_count = agg_count + factor_count

    # Priors 
    A_prior = mean(priors[1])
    V_prior = cov(priors[1])

    Ω_prior = cov(priors[2])
    v₁_prior = prior[2].df

    Σ_prior = cov(priors[3])
    v₂_prior = prior[3].df

    # Measurements 
    measures = vcat(MV...)              # N x T

    # Define containers 
    F_draws = zeros(nsave, state_count, T + 1, chains)
    A_draws = zeros(nsave, length(A_prior), chains)
    Ω_draws = zeros(nsave, length(Ω_prior), chains)
    Σ_draws = zeros(nsave, length(Σ_prior), chains)

    chain_means = Vector{Vector{Float64}}(undef, chains)

    # Initial conditions   #TODO: dimension check for Σ
    A, Ω, Σ = inflate_parameter_space(initial_draws[c], param_sizes, case)
    F = smoothing_distribution(A, Ω, param_sizes, smoother_results, state_count, T, mcmc_jsd_draws, case)
    ŷ, _ = collapse_observational_vector(measures, G, Σ)

    # For each chain, draw from the conditional posteriors in a sequential fashion 
    for c in 1:chains
        j = 1
        for i in range(1; step=1, length=n_draws * thinning_step)
            F .= smoothing_distribution(A, Ω, param_sizes, smoother_results, state_count, T, mcmc_jsd_draws, case)
            Ω .= rand(posterior_of_Ω(F, Ω_prior, v₁_prior))
            Σ .= rand(posterior_of_Σ(F, ŷ, G, Σ_prior, v₂_prior))
            A .= rand(posterior_of_A(F, A_prior, V_prior, ŷ, state_count, lags, constant, T))

            if i <= nburn || i > nburn && mod(i, thinning_step) != 0
                nothing
            elseif i > nburn && mod(i, thinning_step) == 0
                F_draws[:, j, c] = vec(F)
                Ω_draws[:, j, c] = vec(Ω)
                Σ_draws[:, j, c] = vec(Σ)
                A_draws[:, j, c] = vec(A)
                j += 1
            end
        end
    end

    # Store posterior moments. We use Chains() to be able to generate diagnostics
    return GibbsPosteriorResults(A_draws, Ω_draws, Σ_draws)
end

function posterior_of_A(F_draw, A_prior, V_prior, measures, state_count, lags, constant, T)
    X = lagmatrix(F_draw', lags, constant)  # size(F_draw) = (n, T), size(X) = (T - lags, k)
    σ² = var_estimation(F_draw, lags, constant, T)  # Already diagonal 
    V_prior = [Diagonal(diag(V_prior)[(k-1)*state_count+1:k*state_count]) for k in 1:state_count]

    # Posterior Containers 
    V_post = Vector{SparseMatrixCSC{Float64}}(undef, state_count)
    A_post = zeros(Float64, state_count, var_count)

    # To avoid invertibility issues, we operate on each equation, one by one 
    for k in axes(X, 2)
        V_post[k] = inv(inv(V_prior[k]) + inv(σ²[k]) * X' * X) # size(X) = T-1 by k
        conformable_X = X[(!isnan).(measures[k, :]), :]
        conformable_Y = measures[k, :][k, (!isnan).(measures[k, :])]
        A_post[k, :] = V_post[k] * (inv(V_prior[k]) * A_prior[k, :] + inv(σ²[k]) * conformable_X' * conformable_Y)
    end

    if isposdef(V_post) == false
        V_post .= nearest_spd(Matrix(blockdiag(V_post...)))
    end

    return MvNormal(vec(A_post), V_post)
end


function posterior_of_Ω(F_draw, Ω_prior, v₁_prior)
    # Second, define moments of the posterior for the State VCV         
    comp1 = zeros(size(Ω_prior))

    for t in axes(F_draw, 2)
        diff = F_draw[:, t+1] - F_draw[:, t]
        comp1 += diff * diff'
    end
    Ω_bar = Ω_prior + comp1
    v₁ = v₁_prior + T
    return InverseWishart(v₁, Ω_bar)
end


function posterior_of_Σ(F_draw, measures, G, Σ_prior, v₂_prior)
    comp2 = zeros(size(Σ_prior))

    for t in axes(F_draw, 2)
        diff = measures[:, t] - G[t] * F_draw[:, t]
        comp2 += diff * diff'
    end
    Σ_bar = Σ_prior + comp2
    v₂ = v₂_prior + T
    return InverseWishart(v₂, Σ_bar)
end

function rwmh(model_elements, initial_draws, inv_hessian, priors, param_sizes, model_options, mcmc_options, meas_ind, Σ_ids)
    """Runs the random walk metropolis hastings algorithm. 
    http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/ebooks/html/csa/node27.html#SECTION06133000000000000000
    """

    @unpack mhscale, nsave, nburn, chains, adaptive_rwmh, thinning_step = mcmc_options
    init_scales = ones(chains)
    converged = false
    d = length(initial_draws[1])

    # Define MH container of draws 
    # MH             = [zeros(nburn, d) for _ in eachindex(initial_draws)]
    MH = zeros(nsave, d, length(initial_draws))

    scales = [inv_hessian for _ in eachindex(initial_draws)]
    # prop_dist     = MvNormal(zeros(d), Matrix(I, d, d))
    prop_dist = MvTDist(d - 1, zeros(d), Matrix(I, size(inv_hessian, 1), size(inv_hessian, 1)))
    accept_cont = zeros(chains)
    cycles_done = 0

    # Run rwmh cycles until convergence  
    while converged == false
        Threads.@threads for c in 1:chains
            MH[:, :, c], scales[c], accept_cont[c], init_scales[c] = run_rwmh_cycle(MH[:, :, c], prop_dist, scales[c], initial_draws[c], nburn, nsave, model_elements, param_sizes, priors, meas_ind, Σ_ids, model_options, d, cycles_done, init_scales[c], accept_cont[c], c)
        end

        # Check convergence 
        converged = check_convergence(MH)
        cycles_done += 1

        # If converged, run final cycle 
        if converged == true
            @info("All done tuning. Drawing from stationary distribution.")

            # Threads.@threads for c in 1:chains
            #     MH[c], scales[c], accept_cont[c], init_scales[c] = run_rwmh_cycle(MH[c], prop_dist, scales[c], initial_draws[c], nburn, nsave, model_elements, param_sizes, priors, meas_ind, model_options, d, cycles_done, init_scales[c], accept_cont[c], c, true)
            #     @info("Done drawing from stationary distribution.")
            # end

            # Convert MH to array 
            # iters  = size(MH[1], 1)
            # MH     = reshape(reduce(hcat, MH), (iters, d, :))

            # return RWMHPosteriorResults(MH[burned+1:thinning_step:end, :, :]), accept_cont
            return RWMHPosteriorResults(MH[1:thinning_step:end, :, :]), accept_cont
        else
            @info("Finished ($cycles_done) cycle. Starting next cycle.")
        end
    end
end


function run_rwmh_cycle(MH, prop_dist, scale, initial_draw, draws, n_save, model_elements, param_sizes, priors, meas_ind, Σ_ids, model_options, d, cycle, init_scale, accept_rate, c, final=false)
    # Configurations
    n_draws = draws + n_save
    accepted = round(accept_rate * (cycle * n_draws), digits=0)  # every cycle does 'n_draws'
    start = cycle * n_draws + 1                               # the start of the cycle 
    T = (cycle + 1) * n_draws
    initial_draw = cycle >= 1 ? copy(MH[cycle*nsave, :]) : initial_draw  # TODO: for each additional cycle, it gets costly ... 

    if cycle >= 1
        MH = vcat(MH, zeros(n_save, d))
    end

    # Define first draw and run likeli to get old posterior
    old_posterior, alarm = likeli(model_elements, initial_draw, param_sizes, priors, meas_ind, Σ_ids, model_options)
    u = zeros(d)
    xhat_star = zeros(d)

    # Acceptance-rejection algorithm
    local accep
    for i in range(start=start, stop=T; step=1)
        u .= rand(prop_dist)
        xhat_star .= initial_draw .+ scale * u
        new_posterior, alarm = likeli(model_elements, xhat_star, param_sizes, priors, meas_ind, Σ_ids, model_options)
        accprob = min(exp(new_posterior - old_posterior), 1.0)


        if alarm == false && rand() .<= accprob
            # MH[i, :]       = xhat_star
            initial_draw .= xhat_star
            # posterior[a]   = old_posterior 
            old_posterior = new_posterior
            accepted += 1
            # else
            # MH[i, :]       = MH[i-1, :]
            # posterior[a]   = posterior[a-1]
        end
        # a +=1

        if i >= start + draws
            idx = (i - (start + draws)) + (cycle * n_save + 1)
            MH[idx, :] = initial_draw
        end


        if i == T
            @info("-----------------------")
            @info("Chain $c")
            @info("Acceptance Rate: ", accepted / i)
            @info("Number of draws:", i)
            @info("Posterior Likelihood:", -old_posterior)
            @info("-----------------------")
        end

        if mod(i, 100) == 0 && final == false
            scale, init_scale = tune_scale(scale, init_scale, accepted, i)
        elseif final == true
            nothing
        end

        if i == T
            accep = accepted / T
        end
    end
    return MH, scale, accep, init_scale
end

# function adaptive_tuning()
# Matrix(cholesky(inv_hessian).L) .* init_scale
# η(i)          = min(1, d * i^(-2/3))
#             # Compute new scale matrix 
#         # accpt_r = accepted / i
#         # if adaptive_rwmh == true && i > 2 * nburn
#         # if adaptive_rwmh == true && (accpt_r < .20 || accpt_r > .40)
#         #     M      = scale * (I + η(i) .* (accprob - .30) .* (u * u') ./ norm(u)^2 ) * scale'
#         #     M     .= (M .+ M') .* 0.5
#         #     if isposdef(M) == true
#         #         scale .= cholesky(M).L
#         #     else
#         #         M     .= nearest_spd(M)
#         #         scale .= cholesky(M).L
#         #     end
#         # end
# end

function check_convergence(chains)
    """Brooks and Gelman (1997), with corrected PSRF, accounting for sampling variability."""
    # Gelman-Rubin, from Stata
    df = size(chains, 2) - 1
    N = size(chains, 1)
    M = size(chains, 3)

    # Compute variances of chains. Then, compute their mean 
    θ̄ₘ = [vec(mean(chains[:, :, c], dims=1)) for c in 1:M]
    σ²ₘ = [vec(var(chains[:, :, c], dims=1)) for c in 1:M]
    θ̄ = mean(hcat(θ̄ₘ...), dims=2)

    # Between and within variances 
    C = sum([θ̄ₘ[i] .- θ̄ for i in 1:M])
    B = N .* (C .^ 2) ./ (M - 1)
    W = mean(σ²ₘ)

    # Gelman-Rubin 
    Ḃ = ((M + 1) .* (B ./ (M * N)))
    Ẇ = ((N - 1) .* (W ./ N))
    V = Ḃ + Ẇ
    R = sqrt.(((df + 3) / (df + 1)) .* (V ./ W))

    s = 100 * sum(R .> 1.1) / length(R)
    @info("$s" * "%" * " of the parameters have an gelman stat of greater than 1.2.")
    converged = s == 0.0
    return converged
end

function tune_scale(scale, init_scale, accepted, i)
    accpt_r = accepted / i
    if accpt_r > 0.2 && accpt_r <= 0.5
        return scale, init_scale
    else
        scale .= scale ./ init_scale
        if accpt_r <= 0.001
            init_scale *= 0.1
        elseif accpt_r > 0.001 && accpt_r <= 0.05
            init_scale *= 0.5
        elseif accpt_r > 0.05 && accpt_r <= 0.2
            init_scale *= 0.9
        elseif accpt_r > 0.5 && accpt_r <= 0.75
            init_scale *= 1.1
        elseif accpt_r > 0.75 && accpt_r <= 0.90
            init_scale *= 2
        elseif accpt_r > 0.90
            init_scale *= 10
        end
        return scale .* init_scale, init_scale
    end
end








# # Unpack some parameters
# @unpack nsave, nburn, lags, constant, chains, compute_hessian, n_jsd_draws = model_options

# endo_count = size(x_filtered, 2)
# T          = size(x_filtered, 1)

# # Sample from the smoothing distribution. Once for each chain 
# jsd = smoothing_distribution(initial_draws, param_sizes, smoother_results, endo_count, T, n_jsd_draws)


# var_count  = endo_count * lags + exo_count + sum(constant) 
# ν_post     = size(x_trajectory, 1) + endo_count  # T + ν

# # Posterior containers
# a_post = zeros(Float64, endo_count, var_count)
# V_post = Vector{SparseMatrixCSC{Float64}}(undef, endo_count)

# # Define the data objects 
# X     = lagmatrix(x_trajectory, lags, constant)
# SIGMA = var_estimation(x_trajectory, lags, constant)  # Already diagonal 

# # Parameters for posterior draws 
# ntot  = nsave + nburn

# # Containers, pre-allocated, for the draws 
# post_dist = set_up_mcmc(model_options)
# @unpack ALPHA, alpha_draws, OMEGA_draws, SIGMA_draws = post_dist

# # Initialize some matrices 
# V_prior  = Vector{Matrix{Float64}}(undef, endo_count)
# V_prior .= [Diagonal(full_mat[k,:]) for k in 1:endo_count]

# # For each chain, sample from the conditional posterior of each block
# for chain = 1:chains 
#     for i = 1:ntot 
#         # Looping over each equation, k. One per endogenous variable = factor 
#         for k=1:endo_count
#             V_post[k]                 = inv(inv(V_prior[k]) + inv(SIGMA[k]) * X' * X)
#             a_post[k,:]               = V_post[k] * (inv(V_prior[k]) * prior_Γ[k,:] + inv(SIGMA[k]) * X' * x_trajectory[:,k])  
#         end
#         V_post = nearest_spd(Matrix(blockdiag(V_post...))) 

#         alpha  = vec(a_post) + Matrix(cholesky(V_post))' * randn(var_count * endo_count, 1)  # Draws from a normal distribution 
#         ALPHA .= reshape(alpha, var_count, endo_count)  # TODO: should only keep the diagonal elements here 
#         S_post = S_prior + (x_trajectory - X * ALPHA)' * (x_trajectory - X * ALPHA)  
#         S_post = Matrix(Hermitian(S_post))

#         if i > nburn
#             @views alpha_draws[i-nburn, :, chain]  = alpha
#             @views OMEGA_draws[i-nburn, :, chain]  = vec(inv(rand(Wishart(ν_post, nearest_spd(inv(S_post))))))  # Posterior of SIGMA|ALPHA, Data ~ iW(inv(S_post),v_post)
#         end
#         # Convert to diagonal 
#         SIGMA .= diag(FULL_SIG)
#     end
# end