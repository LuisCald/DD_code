# # Import models which are comparable 
# model_dict = Dict()
# models = [tag] # " less AF", " more AF", " 6 factors", " additional factors",
# # models = [""]
# k = tag

# init_path     = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
# for model in models 
#     model_dict[model] = jldopen(init_path * "/posterior_draws" * "/" * m_label * "_$model.jld2") 
# end

# function laplace_approximation(model, model_elements, m_label, param_sizes, priors, meas_ind, Σ_ids, model_options)
#     """
#     H is the negative hessian at the posterior mode.
#     """
#     init_path = BASE_PATH
#     model_post = jldopen(init_path * "/posterior_draws" * "/" * m_label * "_$model.jld2")
#     θ_mode = find_mode(model_post["d_chains"], model_post["lprobs"])

#     # Implementation of the Laplace approximation
#     SSM(par) = -likeli(model_elements, par, param_sizes, hyperpriors, Σ_ids, model_options)[1]
#     SSM(par, _) = -likeli(model_elements, par, param_sizes, hyperpriors, Σ_ids, model_options)[1]

#     n = length(θ_mode)
#     ∂²f∂p∂p = zeros(Float64, n, n)

#     # # new
#     # FiniteDiff.finite_difference_hessian!(∂²f∂p∂p, pp -> SSM(pp), θ_mode)
#     # inv_hessian = inv(∂²f∂p∂p) # inverse of the negative hessian

#     # # Store the hessian 
#     # jldsave(init_path * "/Hessians" * "/$m_label" * "_$model.jld2"; ∂²f∂p∂p)
#     chain_sizes = size(model_post["d_chains"])

#     chain = reshape(model_post["d_chains"], chain_sizes[1] * chain_sizes[2], chain_sizes[end])
#     inv_hess = cov(chain, dims=1) # inverse of the negative hessian

#     log_det_inv_hess = generate_log_det_hessian(inv_hess)

#     lprobs_max = likeli(model_elements, θ_mode, param_sizes, priors, meas_ind, Σ_ids, model_options)[1]

#     mdd = 0.5 * n * log(2 * pi) + 0.5 * log_det_inv_hess + lprobs_max
#     println(mdd)

#     jldsave(init_path * "/MDD" * "/" * "laplace_$model.jld2"; mdd)
# end


# function generate_log_det_hessian(inv_hess)
#     DG = Diagonal(ones(size(inv_hess, 1)))
#     for i in [0.00001, 0.0001, 0.001, 0.01, 0.1, 0.2, 0.5]
#         pert_inv_hess = inv_hess .+ i .* DG
#         D = det(pert_inv_hess)
#         println("Determinant of the hessian: ", D)
#         println("Condition number", cond(pert_inv_hess))
#         if !isapprox(D, 0.0; atol=1e-8)
#             return log(det(pert_inv_hess))
#         end
#     end
# end

# function compute_harmonic_mean_for_all_models!(model_dict, models)
#     model_res = Dict()
#     for model in models
#         chain = model_dict[model]["d_chains"] # iterations by chains by parameters
#         it = size(chain, 1)
#         lprobs = model_dict[model]["lprobs"][end-it+1:end, :] # iterations by chains

#         model_res[model] = mdd_harmonic_mean(chain, lprobs)
#     end

#     # jldsave this dictionary
#     jldsave(init_path * "/MDD" * "/" * "mdd.jld2"; model_res)
#     println(model_res)

#     return model_res
# end

function run_mdd(tag, model_elements, m_label, param_sizes, priors, meas_ind, Σ_ids, model_options)
    # Import chains
    init_path = BASE_PATH
    model_mcmc = jldopen(init_path * "/posterior_draws" * "/" * m_label * "_$tag.jld2")

    # Compute the harmonic mean first
    mdd_hm = mdd_harmonic_mean(model_mcmc)

    # Compute the bridge sampler second
    mdd_bs = bridge_sampler(model_mcmc, model_elements, param_sizes, priors, meas_ind, Σ_ids, model_options)

    # Save in a jld2 file
    init_path = BASE_PATH
    mkpath(init_path * "/7_Results/MDD")
    jldsave(init_path * "/7_Results/MDD" * "/" * "mdd_$tag.jld2"; mdd_hm=mdd_hm, mdd_bs=mdd_bs)
end


function mdd_harmonic_mean(model_mcmc; α=0.05, verbose=false, debug=false)

    chains = model_mcmc["d_chains"] # (100, 2272, 568), iterations by chains by parameters
    it = size(chains, 1)
    lprobs = model_mcmc["lprobs"][end-it+1:end, :] # (100, 2272), iterations by chains


    # Remove chains and the respective lprobs if lprobs is less than -100000
    valid_indices = findall(lprobs .> -100000)
    chains = chains[valid_indices, :]
    lprobs = lprobs[valid_indices]

    # Flatten the chain if it has more than 2 dimensions
    if ndims(chains) > 2
        chains = reshape(chains, size(chains)[1] * size(chains)[2], size(chains)[end])
    end

    # Flatten the log probabilities
    lprobs = vec(lprobs)

    # Compute mean and covariance of the chain
    cmean = mean(chains, dims=1) |> vec
    ccov = cov(chains, dims=1)
    cicov = inv(ccov)  # Inverse of covariance matrix

    nsamples = size(chains, 1)
    d = size(chains, 2)  # Dimension of parameter vector

    # Function to evaluate each chunk
    function runner(chunk::Array{Float64,2})
        res = fill(NaN, size(chunk, 1))

        for i in 1:size(chunk, 1)
            drv = chunk[i, :]
            drl = lprobs[i]

            if (drv - cmean)' * cicov * (drv - cmean) < quantile(Chisq(d), 1 - α)
                res[i] = logpdf(MvNormal(cmean, ccov), drv) - drl
            end
        end

        return res
    end

    # Apply runner function to chain
    mls = runner(chains)
    mls = mls[(!isnan).(mls)]
    mls = mls[mls.<1e12.&&mls.!=-Inf]

    # Numerical stability adjustment
    maxllike = maximum(mls)
    imdd = log(mean(exp.(mls .- maxllike))) + maxllike
    println("max likelihood: ", maxllike)
    println("MDD portion: ", imdd - maxllike)

    return -imdd
end



function bridge_sampler(model_mcmc, model_elements, param_sizes, priors, meas_ind, Σ_ids, model_options; cores=1, tol=1e-10, r0=0)
    bridge_dict = Dict()
    samples = model_mcmc["d_chains"]

    # Split the samples into two parts
    nperchain, nchains, _ = size(samples)

    # fit_index = collect(1:round(Int, nperchain / 2))
    iter_index = 1 #collect(round(Int, nperchain / 2) + 1:nperchain)

    # samples_4_fit = vcat([samples[fit_index, i, :] for i in 1:nchains]...)
    samples_4_iter = vcat([samples[nperchain:nperchain, i, :] for i in 1:nchains]...) # making one long chain from last iteration
    lprobs_sub = model_mcmc["lprobs"][nperchain:nperchain, :][:]
    println(size(samples_4_iter))

    # Fit proposal distribution 
    N1 = size(samples_4_iter, 1)
    N2 = N1
    # m = mean(samples_4_fit, dims=1)[:]
    # V = cov(samples_4_fit)

    # Generate samples from the proposal distribution
    proposal_dist = model_mcmc["propdist"] #MvNormal(m, V)
    gen_samples = rand(proposal_dist, N2)'

    # Evaluate proposal distribution for posterior and generated samples
    q12 = Vector{Float64}(undef, N1)
    q22 = Vector{Float64}(undef, N1)
    Threads.@threads for i in 1:N1 # assumption N1=N2
        q12[i] = logpdf(proposal_dist, samples_4_iter[i, :])
        q22[i] = logpdf(proposal_dist, gen_samples[i, :])
    end

    # Evaluate unnormalized posterior for posterior and generated samples
    q11 = copy(lprobs_sub)
    q21 = Vector{Float64}(undef, N1)
    # @Threads.threads for i in axes(samples_4_iter, 1)
    #     push!(q11, likeli(model_elements, gen_samples[i, :], param_sizes, priors, meas_ind, Σ_ids, model_options))
    # end

    hyperpriors = priors[end-5:end]

    Threads.@threads for j in 1:N1
        q21[j] = likeli(model_elements, gen_samples[j, :], param_sizes, hyperpriors, Σ_ids, model_options)[1]
    end

    # Remove weird likelihood values and ensure everything is the same size 
    scale = 100
    condq = q21 .> -100000
    q21_c = q21[condq]
    q22_c = q22[condq]
    q11_c = q11[condq]
    q12_c = q12[condq]

    # Run iterative scheme to estimate the marginal likelihood
    l1 = q11_c .- q12_c # subtraction since in log terms
    l2 = q21_c .- q22_c
    lstar = median(l1)
    s1 = N1 / (N1 + N2)
    s2 = N2 / (N1 + N2)
    r = 0
    tol = 1e-10
    criterion_val = tol + 1
    i = 0

    # while criterion_val > tol
    #     r_old = r
    #     A = (l2 .- lstar) 
    #     B = (l1 .- lstar) 
    #     numerator = sum(exp.(A ./ (s1 .* exp.(A) .+ s2 .* r))) # for numerical stability # exponential to put things back on the same scale
    #     denominator = sum(1 ./(s1 .* exp.(B) .+ s2 .* r))
    #     r = (N1 / N2) * numerator / denominator
    #     criterion_val = abs((r - r_old) / r)
    #     i += 1
    #     println("Iteration: $i -- Log marginal likelihood estimate: $(round((log(r) + lstar) * scale, digits=6))")
    # end

    while criterion_val > tol
        r_old = r
        A = (l2 .- lstar) ./ scale
        B = (l1 .- lstar) ./ scale
        C = scale .* exp.(A)
        D = scale .* exp.(B)
        E = s1 .* D
        F = s2 .* r
        numerator = sum(C ./ (D .+ F)) # for numerical stability # exponential to put things back on the same scale
        denominator = sum(1 ./ (E .+ F))
        r = (N1 / N2) * numerator / denominator
        criterion_val = abs((r - r_old) / r)
        i += 1
        println("Iteration: $i -- Log marginal likelihood estimate: $(round((log(r) + lstar), digits=6))")
    end

    bridge_dict["logml"] = log(r) + lstar


    return bridge_dict
end
