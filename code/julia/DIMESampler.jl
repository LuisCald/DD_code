
# Trying Gregor's method 
function SSM(model_elements, param_sizes, hyperpriors, meas_ind, Σ_ids, model_options)
    # Convert matrix into a sequence of vectors 
    # vp = [par[:, i] for i in axes(par, 2)]

    function _SSM(par)
        return likeli(model_elements, par, param_sizes, hyperpriors, Σ_ids, model_options)[1]
    end
    # _SSM.(vp)
end

function run_DIME_sampler(model_elements, niter, param_vector, param_sizes, priors, meas_ind, Σ_ids, model_options; bbo_opttime=length(param_vector) * 3)
    @unpack measures, data_cutoffs, tag = model_options
    label = "3D_A non-diag"
    m_label = measures_folder(measures)

    # First run tenative optimization
    hyperpriors = priors[end-5:end] #TODO: hardcoded
    BBO(par) = -likeli(model_elements, par, param_sizes, hyperpriors, Σ_ids, model_options)[1]
    par_final = run_black_box_opt(BBO, param_vector, param_sizes, priors, measures; opttime=bbo_opttime)

    # A_new,B_new,Ω_new,Δ_new,G_new, likeli_vec, Δ_log = run_EM_algorithm(param_vector, param_sizes, meas_ind, Σ_ids, model_elements, model_options)
    # @unpack G = model_elements
    # old_G = deepcopy(G)
    # G = deepcopy(G_new)
    # @pack! model_elements = G

    #         # Save the matrices first
    #         init_path = BASE_PATH
    #         file_name = init_path * "/7_Results/$m_label" * "$tag" * "/from_mcmc/parameter_vectors/solution" * "$label" * ".jld2" 

    #         JLD2.save(file_name, "matrices", [A_new,B_new, Ω_new, Δ_new, G_new, likeli_vec, Δ_log])

    #         Plots.scatter(likeli_vec[2:end], msc=:grey, mc=:white, xformatter=:latex, yformatter=:latex, label="", ylabel=L"\textrm{Log\,\,Likelihood}", xlabel=L"\textrm{Iterations}", legend=:best)
    #         Plots.plot!(likeli_vec[2:end], lc=:blue, label=L"\textrm{Convergence\,\,Path}")
    #         Plots.savefig(init_path * "/7_Results/$m_label" * "$tag" * "/other_results/likelihood_convergence" * ".pdf")

    #         Plots.scatter(Δ_log[2:end], msc=:grey, mc=:white, xformatter=:latex, yformatter=:latex, label="", ylabel=L"\textrm{\Delta\,\,\mathcal{l}_{\bar{k}}}", xlabel=L"\textrm{Iterations}", legend=:best)
    #         Plots.plot!(Δ_log[2:end], lc=:orange, label=L"\textrm{Convergence\,\,Path}")
    #         Plots.savefig(init_path * "/7_Results/$m_label" * "$tag" * "/other_results/delta_log" * ".pdf")


    # # Generate new parameter vector  
    # Δ_new_vec = convert_Δ_new_to_vec(Δ_new, model_options, model_elements)
    # par_final = convert.(Float64, vcat(A_new[:], B_new[:], diag(Ω_new), Δ_new_vec))

    log_likeli = SSM(model_elements, param_sizes, hyperpriors, meas_ind, Σ_ids, model_options)

    # Run MCMC using Gregor's paper 
    LogProbParallel(x) = pmap(log_likeli, eachslice(x, dims=2))  # running an ensemble in parallel (ensemble of chains)
    
    hyper_params = par_final[end-5:end] #TODO: hardcoded
    nchain = length(param_vector) * 4 # a sane default
    priors = get_priors(model_elements, model_options, hyper_params)
    push!(priors, hyperpriors...) # add hyperpriors to the end of the priors list
    init_chains = draw_from_prior(param_sizes, priors, nchain)

    # Replace several chains with the best candidate from the optimization
    init_chains[:, 1:10:end, :] .= par_final

    # off you go sampling
    DIME_chains, lprobs, propdist = RunDIME(LogProbParallel, init_chains, niter, progress=true, aimh_prob=0.1)

    # Save just the last 500 iterations  
    last_25p = round(Int, niter * 0.25)
    to_keep = minimum([500, last_25p])
    DIME_chains = DIME_chains[end-to_keep+1:end, :, :] #     chains = Array{Float64,3}(undef, niter, nchain, ndim)

    # Save the chains, the log probabilities, and the proposal distribution
    init_path = BASE_PATH
    jldsave(init_path * "/posterior_draws" * "/" * m_label * "_$tag.jld2"; d_chains=DIME_chains, lprobs=lprobs, propdist=propdist)
    # jldopen(init_path * "/posterior_draws" * "/" * m_label * "_$tag" * ".jld2", "r")

    # Store paramter vector
    par_final = find_mode(DIME_chains, lprobs)
    store_optim_estimate(par_final, label, m_label, data_cutoffs, tag)

    # par_final3 = mean(DIME_chains[end, :, :][:, :], dims=1) # for 'additional factors', its not 3 dimensions ... only kept last iteration  which is fine

    # Generate plot of chains 
    log_LProbs = -1 .* log.(-1 .* lprobs)

    Plots.plot(log_LProbs[:, :], color="orange4", alpha=0.05, xformatter=:latex, yformatter=:latex, xlabel=L"\textrm{Iterations}", ylabel=L"\textrm{Log-likelihood}", legend=false)
    Plots.plot!(maximum(log_LProbs) * ones(size(log_LProbs, 1)), color="blue3", lw=1)
    Plots.savefig(init_path * "/7_Results/" * m_label * "$tag" * "/from_mcmc/bayesian_convergence/log_probs.pdf")

    # MDD 
    try
        run_mdd(tag, model_elements, m_label, param_sizes, priors, meas_ind, Σ_ids, model_options)
    catch e
        @warn "MDD computation failed with error: $e"
    end

    # Plotting hyperparameter marginals
    generate_marginals_for_hyperparameters(DIME_chains, m_label, tag)

    return par_final
end


function generate_marginals_for_hyperparameters(DIME_chains, m_label, tag)
    # Subset to hyperparameters
    DIME_chains = DIME_chains[:, :, end-5:end]

    init_path = BASE_PATH
    for i in axes(DIME_chains, 3)
        # Identify the bottom 5% quantile along the dimension 
        long_chain = vec(DIME_chains[:, :, i])
        qb = quantile(long_chain, 0.05)
        qt = quantile(long_chain, 0.95)
        filt_chain = long_chain[(long_chain.>qb).&(long_chain.<qt)]
        # if exp_vec[i] == 1
        #     filt_chain = exp.(filt_chain)
        # end
        Plots.plot(filt_chain, lc=:black, color=:orange, fa=0.8, xformatter=:latex, yformatter=:latex, lt=:barhist, legend=false, xtickfontsize=14, ytickfontsize=14, guidefontsize=14)
        Plots.savefig(init_path * "/7_Results/" * m_label * "$tag" * "/from_mcmc/bayesian_convergence/" * "hyperparameter_optimization_$i.pdf")
    end
end

function find_mode(DIME_chains, lprobs)
    to_keep = size(DIME_chains, 1) - 1
    lprobs2 = lprobs[end-to_keep+1:end, :]

    n_chains = size(DIME_chains, 2)

    # Step 2: Initialize variables to track the best (maximum) log-probability and the corresponding parameter vector
    best_log_prob = -Inf  # Start with the lowest possible log-probability (negative infinity)
    best_param_vector = nothing  # Placeholder for the best parameter vector

    # Step 3: Iterate over each chain to find the parameter vector with the highest log-probability (least negative)
    for chain_idx in 1:n_chains
        # Get valid log-probabilities for this chain (if no mask is needed, skip mask filtering)
        valid_log_probs = lprobs2[:, chain_idx]

        # Find the index of the maximum log-probability (i.e., the least negative)
        max_log_prob_idx = argmax(valid_log_probs)

        # Compare to the global best log-probability
        if valid_log_probs[max_log_prob_idx] > best_log_prob
            # Update the best log-probability and corresponding parameter vector
            best_log_prob = valid_log_probs[max_log_prob_idx]
            best_param_vector = DIME_chains[max_log_prob_idx, chain_idx, :]
        end
    end

    return best_param_vector
end

# # Procedure to open the JLD2 file and store the parameter vector
# init_path     = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
# postd = jldopen(init_path * "/posterior_draws" * "/" * m_label * "_$tag.jld2")
# # par_final = find_mode(postd["d_chains"], postd["lprobs"])
# par_final = mean(postd["d_chains"][end, :, :][:,:], dims=1) # for 'additional factors', its not 3 dimensions ... only kept last iteration  which is fine
# store_optim_estimate(par_final, label, m_label, data_cutoffs, tag)

# Generate a par_final which is the average of all the posterior draws 


function find_top_10_unique(DIME_chains, lprobs)
    to_keep = size(DIME_chains, 1) - 1
    lprobs2 = lprobs[end-to_keep+1:end, :]

    n_chains = size(DIME_chains, 2)
    n_params = size(DIME_chains, 3)

    # Step 2: Initialize arrays to track the top 10 best unique log-probabilities and parameter vectors
    best_log_probs = fill(-Inf, 10)  # Array for the 10 highest unique log-probabilities
    best_param_vectors = [zeros(n_params) for _ in 1:10]  # Array for storing the corresponding parameter vectors

    # Step 3: Iterate over each chain to find the top 10 unique parameter vectors with the highest log-probabilities
    for chain_idx in 1:n_chains
        # Get valid log-probabilities for this chain
        valid_log_probs = lprobs2[:, chain_idx]

        # For each valid log-probability, check if it is in the top 10 and unique
        for i in 1:length(valid_log_probs)
            # Check if the current parameter vector is unique among the top 10
            new_param_vector = DIME_chains[i, chain_idx, :]
            is_unique = all(x -> !isequal(x, new_param_vector), best_param_vectors)

            if valid_log_probs[i] > minimum(best_log_probs) && is_unique
                # Replace the lowest in the top 10 if the current log-probability is higher and the parameter vector is unique
                min_idx = argmin(best_log_probs)
                best_log_probs[min_idx] = valid_log_probs[i]
                best_param_vectors[min_idx] = new_param_vector
            end
        end
    end

    # Remove any remaining `nothing` entries in case fewer than 10 unique vectors were found
    non_nothing_indices = findall(!isnothing, best_param_vectors)
    best_log_probs = best_log_probs[non_nothing_indices]
    best_param_vectors = best_param_vectors[non_nothing_indices]

    # Sort the top unique entries based on log-probabilities in descending order
    sorted_indices = sortperm(best_log_probs, rev=true)
    best_log_probs = best_log_probs[sorted_indices]
    best_param_vectors = best_param_vectors[sorted_indices]

    return best_param_vectors, best_log_probs
end


# log_LProbs = -1 .* log.(-1 .* hyper_par_obj["lprobs"])

# Plots.plot(log_LProbs[:,:], color="orange4", alpha=.05, xformatter=:latex, yformatter=:latex, xlabel=L"\textrm{Iterations}", ylabel=L"\textrm{Log-likelihood}", legend=false)
# Plots.plot!(maximum(log_LProbs)*ones(size(log_LProbs, 1)), color="blue3", lw=1)
# Plots.savefig("test.pdf")
