# Hyperparameter selection
# cd("/Users/lc/Dropbox/Distributional_Dynamics/5_Code")
# include("DistributionalDynamics.jl")

# const func_data, time_params, model_elements                   = data_prep(obs_data, model_options);


# Define priors on hyperparameters, using Cauchy distributions for all of them 
# (except for the persistence of the state LOM, which is a Beta distribution)
# Note that the Cauchy distribution has a very heavy tail, which is why it is used here
# The Beta distribution is a distribution on the interval [0,1], which is why it is used here

# no kappa 3
# κ_0   = minn_params[1]  # specifying the prior variance of coefficients that correspond to own lags of endogenous variables
# κ_1   = minn_params[2]  # specifying the size of the prior variance of endogenous variables, which do not correspond to own lags
# κ_2   = minn_params[3]  # specifying the size of the prior variance of non-deterministic exogenous terms
# κ_3   = minn_params[4]  # specifying the size of the prior variance of deterministic terms
# κ_4   = minn_params[5]  # specifying the function exponent of h(lags) = lags^κ_4. "lag decay rate"
# κ_5   = minn_params[6]  # persistence of state LOM

# SSM(par)     = -likeli(model_elements, par, param_sizes, hyperpriors, Σ_ids, model_options)[1]
# SSM(par,_)   = -likeli(model_elements, par, param_sizes, hyperpriors, Σ_ids, model_options)[1]

# hyperpriors = [
#     # Minnesota parameters
#     Normal(.5, 1), # these will be exponentiated to ensure positivity 
#     Normal(.5, 1), 
#     Normal(.5, 1), 
#     truncated(Normal(0.5, 1); lower=0, upper=2), # lag decay rate cannot be some large number and at least 1 = exp(0)
#     Normal(.5, 1),

#     # Error parameters 
#     Normal(.85, 2), # mean of the state VCV
#     Normal(.85, 2), # variance of the state VCV
#     Normal(.85, 2), # mean of the measurement VCV
#     Normal(.85, 2), # variance of the measurement VCV
#     ]



function draw_from_hyperprior(hyperpriors, nchain)
    """Draws from the hyperprior distribution. """

    hyperprior_draws = zeros(length(hyperpriors), nchain)

    for i in eachindex(hyperpriors)
        hyperprior_draws[i, :] = rand(hyperpriors[i], nchain)
    end

    return hyperprior_draws
end

function hyperparameter_optimization(hyperpriors, model_elements, time_p, model_options)
    @unpack prior = model_options
    @unpack case, estimation_object, measures, tag = model_options
    @unpack n_less_than_one, MV = model_elements
    @unpack hyperparameters = prior

    # First, check if the file exists first 
    init_path = BASE_PATH
    m_label = measures_folder(measures)
    hyper_file_name = init_path * "/7_Results/" * m_label * "$tag" * "/from_mcmc/bayesian_convergence/" * "hyperparameter_opt.jld2"
    file_exists = isfile(hyper_file_name)

    if file_exists
        hyper_par_obj = jldopen(hyper_file_name, "r")

        # Setting the hyperparameters
        hyper_to_export = copy(hyper_par_obj["hyper_par_final"])

        hyperparameters[1] = 0.05 #log(exp(hyper_par_final[1])+1) # κ0
        hyperparameters[2] = hyper_to_export[1] #cdf.(Normal(0,1), hyper_par_final[1]) #1 / (1 + exp(-hyper_par_final[1]))
        hyperparameters[6] = hyper_to_export[2] #2 .* (cdf.(Normal(0,1), hyper_to_export[3])) - 1
        hyperparameters[7] = hyper_to_export[3]
        hyper_to_export[4] = log(exp(hyper_to_export[4]) + 1) #−log(1 − (cdf.(Normal(0,1), hyper_to_export[4])))
        hyper_to_export[5] = log(exp(hyper_to_export[5]) + 1) #−log(1 − (cdf.(Normal(0,1), hyper_to_export[5])))
        @pack! prior = hyperparameters


        return hyper_to_export, hyper_par_obj["d_chains"], hyper_par_obj["lprobs"]
    end

    # Save the old params
    old_params = deepcopy(hyperparameters)

    # Define the objective function
    function objective_function(hyperparams)
        # # All hyperparameters
        # κ_0 = minn_params[1]  # specifying the prior variance of coefficients that correspond to own lags of endogenous variables
        # κ_1 = minn_params[2]  # specifying the size of the prior variance of endogenous variables, which do not correspond to own lags
        # κ_2 = minn_params[3]  # specifying the size of the prior variance of non-deterministic exogenous terms
        # κ_3 = minn_params[4]  # specifying the size of the prior variance of deterministic terms
        # κ_4 = minn_params[5]  # specifying the function exponent of h(lags) = lags^κ_4. "lag decay rate"
        # κ_5 = minn_params[6]  # persistence of state LOM
        # κ_6 = minn_params[7]  # persistence of state LOM for aggregates

        # Replace the minnesota params with the hyperpriors 
        hyperparameters[1] = 0.05 # log(exp(hyperparams[1])+1)
        hyperparameters[2] = hyperparams[1]
        hyperparameters[6] = hyperparams[2]
        hyperparameters[7] = hyperparams[3]

        # Estimate the hyperprior
        hyperprior, meas_ind, alarm, logP = get_hyperprior(model_elements, model_options, hyperparams, hyperpriors)

        if alarm >= 1
            return -1e10
        end

        # Estimate the model
        param_vector, param_sizes, Σ_ids, priors = define_parameter_space(model_elements, model_options, hyperprior)

        # Return the objective function
        logL, alarm = likeli(model_elements, param_vector, param_sizes, priors, meas_ind, Σ_ids, model_options) # must be negative for blackbox

        return logP + logL # must be negative for blackbox
    end

    # Run MCMC using Gregor's paper 
    LogProbParallel(x) = pmap(objective_function, eachslice(x, dims=2))  # running an ensemble in parallel (ensemble of chains)

    nchain = length(hyperpriors) * 6 # a sane default
    init_chains = draw_from_hyperprior(hyperpriors, nchain) # n_p by n_chains
    niter = 1000

    # off you go sampling
    DIME_chains, lprobs, _ = RunDIME(LogProbParallel, init_chains, niter, progress=true, aimh_prob=0.1)
    # hyper_par_final = mean(DIME_chains[end, :, :][:, :], dims=1)
    hyper_par_final = find_mode(DIME_chains, lprobs)

    # Save file of DIME_chains, lprobs, and hyper_par_final
    mkpath(init_path * "/7_Results/" * m_label * "$tag" * "/from_mcmc/bayesian_convergence")
    jldsave(init_path * "/7_Results/" * m_label * "$tag" * "/from_mcmc/bayesian_convergence/" * "hyperparameter_opt.jld2"; d_chains=DIME_chains, lprobs=lprobs, hyper_par_final=hyper_par_final)


    log_LProbs = -1 .* log.(-1 .* lprobs)

    Plots.plot(log_LProbs[:, :], color="orange4", lw=2, alpha=0.1, xformatter=:latex, yformatter=:latex, xlabel=L"\textrm{Iterations}", ylabel=L"\textrm{Log-likelihood}", legend=false)
    Plots.plot!(maximum(log_LProbs) * ones(size(log_LProbs, 1)), color="blue3", lw=2)
    Plots.savefig(init_path * "/7_Results/" * m_label * "$tag" * "/from_mcmc/bayesian_convergence/log_probs_hyper.pdf")

    # Setting the hyperparameters
    hyperparameters[1] = 0.05 #log(exp(hyper_par_final[1])+1) # κ0
    hyperparameters[2] = hyper_par_final[1] #cdf.(Normal(0,1), hyper_par_final[1]) #1 / (1 + exp(-hyper_par_final[1]))
    hyperparameters[6] = hyper_par_final[2] #2 .* (cdf.(Normal(0,1), hyper_par_final[3])) - 1 
    hyperparameters[7] = hyper_par_final[3]
    hyper_par_final[4] = log(exp(hyper_par_final[4]) + 1) #−log(1 − (cdf.(Normal(0,1), hyper_par_final[4]))) 
    hyper_par_final[5] = log(exp(hyper_par_final[5]) + 1) #−log(1 − (cdf.(Normal(0,1), hyper_par_final[5])))
    @pack! prior = hyperparameters

    return hyper_par_final, DIME_chains, lprobs
end


function get_hyperprior(model_elements, model_options, hyperparams, hyperpriors)
    """Implements an independent Normal-Cauchy prior. """

    @unpack agg_count, factor_count, pcs, u, MV, proj = model_elements
    @unpack estimator, measures, constant, lags, prior, measurement_error, estimation_object, case, errors_process, pre_multiply = model_options
    @unpack grid_pcf, grid_cop = estimator
    @unpack hyperparameters = prior

    number_of_dfs = length(MV)
    dimension = length(measures)

    # Construct prior hyperparameters  
    prior_mean, A_prior, B_prior, C_prior, D_prior, V_prior = minnesota_prior(hyperparameters, pcs, u, lags, estimator)  # using minnesota to define hyperparameters in normal dist. 
    n_objs = (dimension + 1)
    alarm = false

    local priors
    try
        # println(hyperparams)
        ρ_F = hyperparams[2]
        ρ_Y = hyperparams[3]
        ϕ_F = log(exp(1 - ρ_F^2) - 1)
        ϕ_Y = log(exp(1 - ρ_Y^2) - 1)

        # Remap variance params, since they can be negative
        v_F = log(exp(hyperparams[4]) + 1)
        v_Y = maximum([log(exp(hyperparams[5]) + 1), 0.1]) # Otherwise, it shrinks this extremely and returns basically negative infinity for prioreval()
        corr_Ω_shape = 2 #1 + log(1 + exp(hyperparams[6])) # shape parameter of 2, meaning mild prior towards identity matrix

        priors = [
            MvNormal(prior_mean, V_prior),
            Normal(ϕ_F, v_F),
            Normal(ϕ_Y, v_Y),
            LKJCholesky(factor_count + agg_count, corr_Ω_shape)
        ] # shape parameter of 2, meaning mild prior towards identity matrix

        ϕₘ = log(exp(1) - 1)
        for _ in 1:number_of_dfs
            for _ in 1:n_objs
                push!(priors, Normal(ϕₘ, 2))
            end
        end

        # Add tight priors for the aggregates
        sigma2_star = 1 / 2000                # series used to construct agg factors
        phi_Y_star = log(exp(sigma2_star) - 1) # ≈ sigma2_star
        s_phi_Y = 0.02                      # super-tight
        for _ in 1:agg_count
            push!(priors, Normal(phi_Y_star, s_phi_Y))
        end

    catch ee
        # println(ee)
        alarm = true
        return [], [], alarm, -1e10
    end


    n_param = number_of_dfs * n_objs + agg_count
    Ω_var, Ω_corr, Σ_prior = set_shock_priors(priors, factor_count, agg_count, n_param)
    param_vector = [A_prior, B_prior, C_prior, D_prior, Ω_var, Ω_corr, Σ_prior]

    # Getting indices for measurement error 
    meas_ind = extract_meas_ind(estimator, dimension)  # TODO: no need to worry about this 

    logP, alarm2 = hyper_prioreval(hyperparams, hyperpriors)

    return Prior(priors, param_vector), meas_ind, alarm + alarm2, logP
end