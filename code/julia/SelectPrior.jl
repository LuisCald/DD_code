# TODO: Implement Hierarchical prior: 
# https://watermark.silverchair.com/rest_a_00483.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAz0wggM5BgkqhkiG9w0BBwagggMqMIIDJgIBADCCAx8GCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQM_0w1babP8v82BixJAgEQgIIC8NsTFrwY45o9ASRyIg7rxyg7n2PfwZez0n8YXhoFL2mYEnZyXS0TU6IGPZyOSLQ9Lrx9EnvBQ27f0Kbn2e9k4SbBFDiF5yLtBcC5WXadXYMIwLmWu8QmXYVvLt_kFVNsOsp_p44IUFKiVNbdy4BCqJEAOrQmiN0LF4Ys8KUnoIbgqEOjVhloBUwKdVa7KxU_54hmZoSMD41FavBQTmBeoF7pf7TYK4q8AJZw1ePGFjkIM1rJ47DV42SaPg4Q5uzzuEzJTs3nUOOtlWGjGFi5BqVeKeSNTHh90asADdlwjS3oiMLL5RV0qSV3508qaHwbPfh5Wz9r7ou6ztx7ZnaYmnvH3G5YVrvDC8MF51N9soctjAn2p7aI8tjqR5gfcLSGm2hhDRvQ6kO3656O-wOq47YBHbnWuzJZK_hePsJYqzFZ5gr6s3Lg-iV_ZIxXcbqz99L2_eFLE8Yv1rEFz8lqOBjqBfUA6KzRfyPJgokS3o_3DZVsYJ7kAPq2_ui5viC_wrRs2hsdz5Lc3iYav6ZELVuLNVWgHapZa8TciecqFrTsC35Wpbql7mFeXTe3MSK3iL08ymyQOprBBXj8dspi5NAeYeoOyNw53kl-NjDzeRQxsWTpPA_2UxUHs2H5v7WgQhB1DHDiEHOQt26pXhRiUSmDzkKTrXgy2BhUGj2AftEUelBrZNjd040YP7NjVJego_qEktUfhmtUEmfJRA2Lxd3gp1oAsxEEso_SI2vidoCVn8ZkqZjCqys8cOxGUSEYzK74IL3hcAIk9LJDsQhI-6rCtZNRRXW_qBT6X_JxSM77b3zL7jIGuBFivVt6oMEhrUNnV0mFHiFgIMbbokNbXzvMK17wZ7iiaWDBXRyeDC1rbq3Dg84M1viE7gtbor8xmCJJCrt6bIPGRWO9-xG56fPS-e-aV5RhJhMYixxSHZTVe4Apcrln2YFlnHkdHplcWcpk_A69OGisTbgz0nfCa_MFL1n7le-QoeqiQkNjBVmA 
# https://www.ecb.europa.eu/pub/pdf/scpwps/ecbwp1494.pdf
# estimation of linear systems using gibbs: http://users.isy.liu.se/rt/schon/Publications/WillsSLN2012.pdf
# nice section on priors: https://pdf.sciencedirectassets.com/273539/1-s2.0-C20120135663/1-s2.0-B9780444627315000154/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEIn%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJIMEYCIQC8QHwu9DdxD35Xdqyis0bHey2ZJruVZcYz5NdGEcEDcwIhAJLeo8FV%2FNNYbuKnSdgw%2FAPK0Gv3SUphrNmI%2BrEpzOF9KtsECPL%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEQBRoMMDU5MDAzNTQ2ODY1Igw%2FtcyZEe7drG4k8X8qrwTm9tjEaXAnxFzUngblagO9loXgYlXB0OLk8ZoC%2Fnr08y%2F3r%2B28gTIRyJoGJ8yjsf8C66Rw1oRyImz%2F01pVRx9ZAJEHyRsl3EuxKa2CsrziOrSMfLyJI6HCe3g2CCFBJ3g047MPnKSiLb8zQZap68XvBLvFimy7Q84iCljC4m13q3dE2z8UEowo%2FdfIEEI7L%2BttgrPYQOGxbvJuemLHBiy9hJsf8SVwHkZSgRv4Qp%2BvHypSa5LGpanZs66s2klJ3YRiLGb4mgcQ%2FqW6nCkZbX2PhCNDrB%2BuPjnWRwydtfY5Bw9ttxzCRn4hdR3ejyTtbBqgxvgfnZMwKKN4oN7mcbcYG0bouAhExBL8yXeKTzbLc8dSoiL5n8nsY4tR2KwcF86fQfE0cVLPdPTCtq1KR8wSwa6KDxdRBFObhmtVqkD0GyG8gTocy1Ic6QA%2F46VmIgAU8dOcT4r1TmBRmPZ6fwfpnyiX9V4GqCxdOZwf2CtRwVHFyuoXbqZtBU3Px9%2BfhuYHAcVwIN%2Fa2CR0jxtP1QJZfOCyLcDyVUYSWhSUNIb65lrq050kA1whuetaaNJdnRuk7Jji%2BgRyuz2RwugbSsidFrfVEjz84BAX%2FmuYIz75Pga6Y1a%2B1NdXdHz%2FoLRW%2B%2FZRO%2BMxZKTVHm%2ByiddRVVEUY4NPK8JvSCtLd6qBRwg%2FoJCF1IZmZb9STJlU%2BAKUSA2ejpJOog3k6LDXZNqFDm9XXGWTjvxfokGSsqP1r8tAMJKx9JcGOqgBRoIlLPWVdV06Rg2rKUUm%2FrW0JCt1wkyDiSMKH0dJkR0ylSBvdvh9peA6FrQpcapCp7aLPLxGzT%2BRn6yMRdqFS7bCobarph2xnyRSs%2FftFCoOMNZiL1Eh23pn9QxIVaBP08po8diVY9qaUmCrcUDwu5kL5rDOXql1nTelAz95wXGU0Ij7XfCCuHFrzEosHykCthesaHEY3mWqew6r6V9TPpX7ecuaPbgr&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20220817T174114Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYWPZNFBIK%2F20220817%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=5c76d0c664b1dc46a79c6cc3e8524596cef8e8ee8ebaff822a794d571516a4eb&hash=fb295db63f67d767d4b415d527ad93bd3aa59e4eee475bdf6080c50004111218&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=B9780444627315000154&tid=spdf-4f08ad4e-0dbe-4531-95a8-91160e96218c&sid=8153f7853696a04cb88a53a-8b5057ce5062gxrqb&type=client&ua=53570d525355550d5b&rr=73c4330c5c7fb8d6
# nice section on priors: https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf
# nice stuff on estimation: https://joshuachan.org/papers/large_BVAR.pdf
# TODO: Inverse Gamma? https://stats.stackexchange.com/questions/482824/bayesian-estimation-of-the-variance
# # TODO: the issue is that the wishart is generating a prior likelihood that is positive...
# # Scaled Inverse Wishart: https://www.tandfonline.com/doi/pdf/10.1198/016214508000000724?casa_token=Q3nW55JKV_gAAAAA:XRsJ9h7RgMGAcdZTFQBS-5KQrnYV4-VA_7-S7tnhh_iA47TmJ6FfSiLJ0qwFnEswsX0ItcD-_IoJBRkw
# A bunch of useful formulas: # https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf
# Other links to make a distribution 
# https://github.com/JuliaStats/Distributions.jl/blob/master/src/univariate/discrete/discretenonparametric.jl
# https://github.com/gragusa/BayesianTools.jl/blob/master/src/ProductDistributions/ProductDistributions.jl

# TODO: MAIN: https://arxiv.org/pdf/1408.4050.pdf - has hamiltonian monte carlo stuff 

# An obvious alternative would be to choose independent scaled inverted chi-squared distributions for each of the variances, 
# as this is the commonly used conjugate prior for a variance. In the real application discussed in Section 3.3 we found the 
# log normal prior more appealing, as it was more difficult to deal with the tail behavior of the inverted chi-squared 
# with low degrees of freedom.

function estimate_prior(model_elements::StateSpaceModel, model_options::ModelOptions)
    """Implements an independent Normal-Cauchy prior. """

    # Unloading
    @unpack agg_count, factor_count, pcs, u, MV, proj = model_elements
    @unpack estimator, measures, constant, lags, prior, measurement_error, estimation_object, case, errors_process, pre_multiply = model_options
    @unpack grid_pcf, grid_cop = estimator
    @unpack hyperparameters = prior

    aggregate_rep = hasproperty(model_options, :aggregate_rep) ? model_options.aggregate_rep : :as_states

    number_of_dfs = length(MV)
    dimension = length(measures)

    if aggregate_rep == :as_inputs
        # ─── :as_inputs ─── No C, D; Ω is factor-only; 4 hyperparameters ───
        hyperpriors = [
            Normal(0, 5),   # κ_1  (prior variance scaling)
            Normal(0, 5),   # κ_5  (factor persistence)
            Normal(0, 1),   # σ_ϕ  (variance of factor persistence prior)
            Normal(0, 2)    # ν    (LKJ shape)
        ]

        hyper_params = [0.5, 0.5, 0.5, 3.0]

        # For Minnesota prior
        hyperparameters[1] = 0.05
        hyperparameters[2] = hyper_params[1]
        hyperparameters[6] = hyper_params[2]
        hyper_params[3] = log(exp(hyper_params[3]) + 1)

        # Get initial parameter values (no C, D)
        prior_mean, A_prior, B_prior, _, _, V_prior = minnesota_prior(hyperparameters, pcs, u, lags, estimator; aggregate_rep=:as_inputs)

        # For measurement error (no aggregate entries)
        n_objs = (dimension + 1)
        n_param = n_objs * number_of_dfs

        # ───  Factor persistence hyper-parameter only  ─────────────────────────
        ρ_F = hyperparameters[6]
        ϕ_F = log(exp(1 - ρ_F^2) - 1)

        σ_ϕ = hyper_params[3]   # already softplus-corrected
        corr_Ω_shape = 1 + log(1 + exp(hyper_params[4]))

        # Defining initial priors (no aggregate persistence, factor_count-only LKJ)
        priors = [
            MvNormal(prior_mean, V_prior),           # [1] coefficients of A,B
            Normal(ϕ_F, σ_ϕ),                        # [2] prior for factor persistence
            LKJCholesky(factor_count, corr_Ω_shape),  # [3] prior for factor correlations (factor_count only)
        ]

        # Measurement error priors (no aggregate entries)
        ϕₘ = log(exp(1) - 1)
        for _ in 1:number_of_dfs
            for _ in 1:n_objs
                push!(priors, Normal(ϕₘ, 2))
            end
        end

        # Pushing hyperpriors to the end of the priors list
        push!(priors, hyperpriors...)

        # Get initial parameters: Errors (factor-only Ω)
        Ωf_mode = repeat([StatsBase.mode(priors[2])], factor_count)
        Ω_var = Ωf_mode
        Ω_corr = lower_offdiag(diagm(ones(factor_count)))
        Σ_mode = [StatsBase.mode(priors[3+i]) for i in 1:n_param]
        Σ_prior = Diagonal(Σ_mode)

        # Final parameter vector (6 blocks: A, B, Ω_var, Ω_corr, Σ, H)
        param_vector = [A_prior, B_prior, Ω_var, Ω_corr, Σ_prior, hyper_params]

        # Getting indices for measurement error
        meas_ind = extract_meas_ind(estimator, dimension)

        return Prior(priors, param_vector), meas_ind

    else
        # ─── :as_states (default) ─── Original behavior ───────────────────────

        # Define hyperpriors
        hyperpriors = [
            # Minnesota parameters
            # Uniform(0.2, 0.99), # specifying the size of the prior variance of endogenous variables, which do not correspond to own lags
            # Uniform(-0.99, 0.99), # persistence of state LOM
            # Uniform(-0.99, 0.99), # persistence of state LOM for aggregates
            Normal(0, 5),
            Normal(0, 5),
            Normal(0, 5),
            Normal(0, 1), # variance of the variances - factors
            Normal(0, 1), # variance of the variances - aggregates
            Normal(0, 2) # hyperparameter for the shape of correlations of LKJ prior, ν
        ]

        # Initial hyperparameters
        hyper_params = [0.5, 0.5, 0.5, 0.5, 0.5, 3.0] # initial values for the hyperparameters

        # For Minnesota prior
        hyperparameters[1] = 0.05
        hyperparameters[2] = hyper_params[1]
        hyperparameters[6] = hyper_params[2]
        hyperparameters[7] = hyper_params[3]
        hyper_params[4] = log(exp(hyper_params[4]) + 1)
        hyper_params[5] = log(exp(hyper_params[5]) + 1)

        # Get initial parameter values
        prior_mean, A_prior, B_prior, C_prior, D_prior, V_prior = minnesota_prior(hyperparameters, pcs, u, lags, estimator)  # TODO: assumes aggs only have 1 lag

        # For measurement error
        n_objs = (dimension + 1)
        n_param = n_objs * number_of_dfs + agg_count

        # ───  Extract the two persistence hyper-parameters  ─────────────────────────
        ρ_F = hyperparameters[6]        # κ₅  = persistence target for the  r-factor block
        ρ_Y = hyperparameters[7]        # κ₆  = persistence target for the  q-aggregate block

        #  Same log-transform you already used for ρ_F, now applied to both blocks
        ϕ_F = log(exp(1 - ρ_F^2) - 1)    # ⇒ Normal prior on ϕ_F  implies Beta prior on ρ_F
        ϕ_Y = log(exp(1 - ρ_Y^2) - 1)    # ⇒ identical mapping for the aggregate persistence

        # ───  Prior list  ──────────────────────────────────────────────────────────
        σ_ϕ = hyper_params[4] # These are already corrected
        σ_ϕ_Y = hyper_params[5]

        corr_Ω_shape = 1 + log(1 + exp(hyper_params[6])) # shape parameter of 2, meaning mild prior towards identity matrix

        # Defining initial priors
        priors = [
            MvNormal(prior_mean, V_prior),   # coefficients of A,B,C,D
            Normal(ϕ_F, σ_ϕ),               # prior for factor persistence
            Normal(ϕ_Y, σ_ϕ_Y),               # prior for aggregate persistence
            LKJCholesky(factor_count + agg_count, corr_Ω_shape),  # prior for factor correlations
        ]

        # For measurement equation
        ϕₘ = log(exp(1) - 1) # when softplus applied to it, it should be 1, which is our prior mean on the measures
        for _ in 1:number_of_dfs
            for _ in 1:n_objs
                push!(priors, Normal(ϕₘ, 2)) # was (5,10)
            end
        end

        # Add tight priors for the aggregates
        sigma2_star = 1 / 2000                # series used to construct agg factors
        phi_Y_star = log(exp(sigma2_star) - 1) # ≈ sigma2_star
        s_phi_Y = 0.02                      # super-tight
        for _ in 1:agg_count
            push!(priors, Normal(phi_Y_star, s_phi_Y))
        end

        # Pushing hyperpriors to the end of the priors list
        push!(priors, hyperpriors...)

        # Get initial parameters: Errors
        Ω_var, Ω_corr, Σ_prior = set_shock_priors(priors, factor_count, agg_count, n_param)

        # Final parameter vector
        param_vector = [A_prior, B_prior, C_prior, D_prior, Ω_var, Ω_corr, Σ_prior, hyper_params]

        # Getting indices for measurement error
        meas_ind = extract_meas_ind(estimator, dimension)  # TODO: no need to worry about this

        return Prior(priors, param_vector), meas_ind
    end
end



function define_lower_bounds(boot_noise_processes)
    # Extracting all dataset names
    dataset_names = sort(collect(keys(boot_noise_processes)))

    # Initialize a dictionary to store the scaled values
    lower_bounds = Dict()

    # Iterate for each row (object)
    for dataset in dataset_names
        # Collect means for each object
        # a = [filter(isfinite, boot_noise_processes[dataset][k, :]) for k in axes(boot_noise_processes[dataset], 1)]
        lower_bounds[dataset] = nanminimum(boot_noise_processes[dataset])
    end

    return lower_bounds
end

function scale_measurement_errors(data_vector, dimension, measurement_default_hyperparameter, α=10)
    # Extracting all dataset names
    dataset_names = sort(collect(keys(data_vector)))

    # Initialize a dictionary to store the scaled values
    scale_dict = Dict()

    n_objs = dimension + 1
    r̄ = Dict()

    # Iterate for each row (object)
    for dataset in dataset_names
        # Collect object-mean
        r̄ⱼ = [nanmean(data_vector[dataset][row, :]) for row in 1:n_objs]

        # Take the average across all means
        r̄[dataset] = nanmean(r̄ⱼ)
    end

    min_r̄ = nanminimum([r̄[dataset] for dataset in dataset_names])
    max_r̄ = nanmaximum([r̄[dataset] for dataset in dataset_names])
    ub = measurement_default_hyperparameter * α

    # Calculate the ratio and update the scale_dict
    for dataset in dataset_names
        # println("dataset: ", dataset, " r̄: ", r̄[dataset], " min_r̄: ", min_r̄, " max_r̄: ", max_r̄, " ub: ", ub)
        # if max_r̄ - min_r̄ == 0
        #     scale_dict[dataset] = measurement_default_hyperparameter
        # else
        scale_dict[dataset] = measurement_default_hyperparameter + ((r̄[dataset] - min_r̄) * (ub - measurement_default_hyperparameter) / (max_r̄ - min_r̄))   # e.g., 1 + 0.2 = 1.2 = the factor to multiply the default with later
        # end

        # scale_dict[dataset] = r̄[dataset] / min_r̄ # TODO: unlikely, but this one does not have an upper bound 
    end

    return scale_dict
end

function scale_measurement_errors_all(data_vector, dimension)
    # Extracting all dataset names
    dataset_names = sort(collect(keys(data_vector)))

    # Initialize a dictionary to store the scaled values
    scale_dict = Dict()
    r̄ = Dict()

    n_objs = dimension + 1

    # Iterate for each row (object)
    for dataset in dataset_names
        # Collect object-mean
        scale_dict[dataset] = [nanmean(data_vector[dataset][row, :]) for row in 1:n_objs]
        # println(scale_dict[dataset])
    end


    # Calculate the ratio and update the scale_dict
    for k in 1:n_objs
        # Find minimum of all objects 
        min_r̄ = nanminimum([scale_dict[dataset][k] for dataset in dataset_names])

        for dataset in dataset_names
            scale_dict[dataset][k, 1] = scale_dict[dataset][k] ./ min_r̄
        end
    end

    return scale_dict
end



# Issue is: I need to know where the income quantiles fall in the measurement matrix, per dataset 
# If income is in the first dimension, then whenever the first dimension changes, it gets a new param, thats it. The next thing to know is where in the vector do the income quantiles lie 
# So, early in the code, I can store the locations of the first  # TODO: is this still true? Need to remind myself  

function extract_meas_ind(estimator, dimension)
    """Extracts the indices of the measurements, removing the immutable part."""

    @unpack grid_pcf, grid_cop = estimator

    # Original Matrix 
    cop_size = tuple([grid_cop for i in 1:dimension]...)
    orig = zeros(cop_size...)
    cart_ind = CartesianIndices(size(orig))

    f(c) = sum((==(1)).(c.I)) >= length(c.I) - 1 # immutable 
    to_keep = cart_ind[filter(!f, CartesianIndices(size(orig)))]

    meas_ind = zeros(length(to_keep))

    for j in eachindex(meas_ind)
        meas_ind[j] = to_keep[j].I[1]
    end

    return meas_ind
end


function var_estimation_old(pcs, MV, u, time_p, lags)
    @unpack tmin, tmax, year_vec, tot_periods = time_p

    # Extract the aggregate factors for the corresponding year and q4
    u_df = DataFrame(u', :auto)
    u_df[:, "dates"] = QuarterlyDate(tmin["year"], (tmin["quarter"])):Quarter(1):QuarterlyDate(tmax["year"], tmax["quarter"])
    u_df[:, "year"] = Dates.year.(u_df[:, "dates"])
    u_df[:, "quarter"] = Dates.quarterofyear.(u_df[:, "dates"])

    fully_informed_indices = []
    for j in eachindex(MV)
        for i in axes(MV[j], 2)
            if sum(isnan.(MV[j][:, i])) == 0
                push!(fully_informed_indices, i)
            end
        end
    end
    sort!(fully_informed_indices)

    # Extract the indices from the dates column that correspond to thiese fully informed indices 
    fully_informed_df = u_df[fully_informed_indices, :]


    c_year_v = sort(comp_years)
    u_years = unique(Int.(c_year_v))
    d_years = Int.([k for (k, v) in countmap(c_year_v) if v > 1])  # duplicate years 

    filter!(row -> row.year <= maximum(u_years), u_df)
    filter!(row -> row.quarter == 4, u_df)
    filter!(row -> row.year in u_years, u_df)

    # Duplicate the row if it has a duplicate year 
    agg_mat = Matrix(u_df)
    k = 1
    for i in axes(agg_mat, 1)
        if u_df[k, "year"] ∉ d_years
            nothing
            k += 1
        else
            insert!.(eachcol(u_df), k + 1, agg_mat[i, :])
            k += 2
        end
    end

    # Remove the date columns 
    select!(fully_informed_df, Not(:dates))
    select!(fully_informed_df, Not(:year))
    select!(fully_informed_df, Not(:quarter))

    # Define OLS components 
    T = size(pcs, 2)
    Y = pcs[:, 1+lags:end]'   # one step ahead, T x N 
    C = Matrix(fully_informed_df)[1:end-lags, :]  # T x N
    X = hcat(pcs[:, 1:end-lags]', C) # lag 

    # Export a latex table with the ols results of the low frequency 
    # X   = pcs[:, 1:end-lags]'
    # T = size(con1, 1) 
    # t = collect(1:T)
    # TT = hcat(ones(T), t, 2 .* t.^2 ./ length(t))
    # stats_dict = Dict("factors" => zeros(3, size(Y, 2)), "aggregates" => zeros(3, size(Y, 2)), "trend" => zeros(3, size(Y, 2)))
    β_aggs = zeros(size(Y, 2), size(C, 2))
    SX = hcat(X[:, 1:5], X[:, size(Y, 2)+1:end])
    for i in axes(Y, 2)
        # # on factors 
        # ols                         = lm(X[:, 1:size(Y, 2)], Y[:, i]) 
        # A_diag[i] = coef(ols)[i]
        # stats_dict["factors"][1, i] = r2(ols)
        # stats_dict["factors"][2, i] = mean(coef(ols))
        # stats_dict["factors"][3, i] = mean(abs.(coef(ols)))  

        # on aggregates 
        # ols                            = lm(X[:, size(Y, 2)+1:end], Y[:, i])
        ols = lm(SX, Y[:, i])

        β_aggs[i, :] = coef(ols)[6:end]

        # Analyzing which aggregates are the most important
        # r2ols = lm(X[:, size(Y, 2)+1:end], Y[:, i])
        # println(r2(r2ols))
        # println(coeftable(r2ols))
        # stats_dict["aggregates"][1, i] = r2(ols)
        # stats_dict["aggregates"][2, i] = mean(coef(ols))
        # stats_dict["aggregates"][3, i] = mean(abs.(coef(ols)))  


        # # on trend 
        # ols                       = lm(TT, Y[:, i])
        # stats_dict["trend"][1, i] = r2(ols)
        # stats_dict["trend"][2, i] = mean(coef(ols))
        # stats_dict["trend"][3, i] = mean(abs.(coef(ols)))  
    end

    # Create plots 
    # Plots.plot(axes(Y,2), stats_dict["factors"][1, :], label = L"R^2\, \textrm{of\, factors\, on\, factors}", legend = :topright)
    # Plots.plot!(axes(Y,2), stats_dict["aggregates"][1, :], label = L"R^2\, \textrm{of\, factors\, on\, aggregates}", legend = :topright)
    # Plots.plot!(axes(Y,2), stats_dict["trend"][1, :], label = L"R^2\, \textrm{of\, factors\, on\, trend}", legend = :topright)
    # Plots.savefig("r2_factors.png")


    # Plots.plot(axes(Y,2), stats_dict["factors"][2, :], label = L"β̄\, \textrm{ of\, factors\, on\, factors}", legend = :topright)
    # Plots.plot!(axes(Y,2), stats_dict["aggregates"][2, :], label = L"β̄\, \textrm{ of\, factors\, on\, aggregates}", legend = :topright)
    # Plots.plot!(axes(Y,2), stats_dict["trend"][2, :], label = L"β̄\, \textrm{ of\, factors\, on\, trend}", legend = :topright)
    # Plots.savefig("coefs_factors.png")

    # Plots.plot(axes(Y,2), stats_dict["factors"][3, :], label = L"|β̄|\, \textrm{ of\, factors\, on\, factors}", legend = :topright)
    # Plots.plot!(axes(Y,2), stats_dict["aggregates"][3, :], label = L"|β̄|\, \textrm{ of\, factors\, on\, aggregates}", legend = :topright)
    # Plots.plot!(axes(Y,2), stats_dict["trend"][3, :], label = L"|β̄|\, \textrm{ of\, factors\, on\, trend}", legend = :topright)
    # Plots.savefig("abs_coefs_factors.png")


    # Perform the least squares estimation. 
    # A_OLS    = inv(X' * X) * (X' * Y)
    A_OLS = inv(X' * X) * (X' * Y)
    var_Ω = diag(((Y - X * A_OLS)' * (Y - X * A_OLS)) ./ (T - size(pcs, 1) + 1))


    # Computing the standard deviations of the controls directly 
    # var_Ω[2] = diag(cov(con))

    # println(describe(A_OLS[:]))
    # println(describe(var_Ω[1]))
    # println(describe(var_Ω[2]))


    return var_Ω, β_aggs #vcat(var_Ω...)

end


function var_estimation(MV, pcs, u, time_p, lags, proj)
    @unpack tmin, tmax, year_vec, tot_periods = time_p

    # First estimate the tuning parameters of the measurement error 
    # Σ = extract_measurement_error_params(MV, dimension, grid, errors_process)
    # if errors_process ==  "one per object, per dataset"
    #     Σ = estimate_measurement_VCV(MV, dimension, grid, errors_process, time_p)
    #     return Σ
    # end

    # Linear interpolate points 
    # Find the dataset with the most data 
    data_lengths = [size(MV[j][:, mapslices(col -> all((!isnan).(col)), MV[j], dims=1)[:]])[2] for j in eachindex(MV)]
    id_max = findmax(data_lengths)[2]
    LI_df = similar(MV[id_max])

    # Time index by dataset
    t = collect(1:1:tot_periods)

    # Interpolate this dataset for all rows 
    for i in axes(MV[id_max], 1)
        mask = (!isnan).(MV[id_max][i, :])
        interp_linear = linear_interpolation(t[mask], MV[id_max][i, :][mask], extrapolation_bc=Line())
        LI_df[i, :] .= interp_linear.(t)
    end

    LI_df = convert.(Float64, LI_df) # otherwise, things will be wrong numerically (since it initially returns a simpleratio type)

    # Use the projection matrix from the PCA to project the interpolated series onto the factors
    DF = proj' * LI_df                # distributional factors 
    Y = DF[:, 1+lags:end]
    # Y_l = DF[:, 1:end-lags]

    V = VAR(Y', 4, false)

    # # Aggregates 
    # C   = transpose(u)[1:end-lags, :]
    # X   = vcat(Y_l, C')
    # X   = copy(Y_l)
    # A   = pinv(X * X') * (X * Y')  # given the data structure, I had to transpose everything  (42, 19)
    # Ω   = diag(((Y' - X' * A)' * (Y' - X' * A)) ./ (tot_periods - size(X, 1) + 1))
    # # B   = transpose(A[size(Y, 1)+1:end, :])
    Ω = diag(V.Σ)

    # # To analyze in Stata 
    # some_mat = hcat(Y', C)
    # CSV.write("A.csv", DataFrame(some_mat, :auto))

    return Ω
end

function add_trends(X)
    T = size(X, 1)
    t = collect(1:T)
    X = hcat(ones(T), t, 2 .* t .^ 2 ./ length(t), X)
    # X = hcat(ones(T), t, X)
    return X
end


function minnesota_prior(minn_params, pcs, controls, lags, estimator; aggregate_rep=:as_states)
    """ Returns priors for A and Q. Treatment from https://documentation.sas.com/doc/en/etsug/15.2/etsug_varmax_details30.htm"""

    @unpack grid_pcf, grid_cop = estimator

    n_factors = size(pcs, 1)  # number of factors
    n_aggs = size(controls, 1) # number of aggregates

    # Unload parameters. From: https://cran.r-project.org/web/packages/bvartools/bvartools.pdf pg.59
    κ_0 = minn_params[1]  # specifying the prior variance of coefficients that correspond to own lags of endogenous variables
    κ_1 = minn_params[2]  # specifying the size of the prior variance of endogenous variables, which do not correspond to own lags
    # κ_2 = minn_params[3]  # specifying the size of the prior variance of non-deterministic exogenous terms
    # κ_3 = minn_params[4]  # specifying the size of the prior variance of deterministic terms
    κ_4 = minn_params[5]  # specifying the function exponent of h(lags) = lags^κ_4. "lag decay rate"
    κ_5 = minn_params[6]  # persistence of state LOM
    κ_6 = minn_params[7]  # persistence of state LOM for aggregates

    var_F = VAR_fit(pcs', lags)
    var_Y = VAR_fit(controls', lags)
    var_Ω = vcat(var_F, var_Y)  # variance of the factors and controls (needed for scaling in both modes)

    # Set mean
    A_mean = diagm([ones(n_factors)..., zeros(n_factors * (lags - 1))...]) .* κ_5
    B_mean = zeros(n_factors, n_aggs)

    # Prior variance matrix (Minnesota logic)
    # ----------------   prior variances for endogenous lags   ----------------
    total_state = n_factors + n_aggs
    int_mat = zeros(n_factors, total_state * lags)

    # V_prior for the top half of the matrix (A and B blocks)
    for row in 1:n_factors
        s = 1
        for col in 1:(total_state*lags)
            scaling_factor = var_Ω[row] / var_Ω[col]
            if scaling_factor == 1.0
                int_mat[row, s:s+lags-1] = [(κ_0 / l^κ_4) for l in 1:lags]
            else
                int_mat[row, s:s+lags-1] = [(κ_0 * κ_1 / l^κ_4) * scaling_factor for l in 1:lags]
            end
            s += lags
        end
    end

    if aggregate_rep == :as_inputs
        # ─── :as_inputs ─── No C, D; prior covers only A, B ─────────────────
        prior_mean = vcat(A_mean[:], B_mean[:])
        V_all = vec(int_mat)  # only A, B variance entries
        V_prior = diagm(V_all)

        return prior_mean, A_mean, B_mean, nothing, nothing, V_prior

    else
        # ─── :as_states (default) ─── Original behavior with C, D ───────────
        C_mean = zeros(n_aggs, n_factors)
        D_mean = ones(n_aggs) .* κ_6
        prior_mean = vcat(A_mean[:], B_mean[:], C_mean[:], D_mean[:])

        # V_prior for C
        VC_prior = zeros(n_aggs, n_factors)
        for row in n_factors+1:n_factors+n_aggs
            s = 1
            for col in 1:n_factors*lags
                scaling_factor = var_Ω[row] / var_Ω[col]
                if scaling_factor == 1.0
                    VC_prior[row-n_factors, s:s+lags-1] = [(κ_0 / l^κ_4) for l in 1:lags]
                else
                    VC_prior[row-n_factors, s:s+lags-1] = [(κ_0 * κ_1 / l^κ_4) * scaling_factor for l in 1:lags]
                end
                s += lags
            end
        end

        # V_prior just for D
        # ----------------   prior variances for aggregates   ----------------
        VD_prior = [κ_0 for _ in 1:n_aggs]
        V_all = vcat(vec(int_mat), vec(VC_prior), VD_prior)
        V_prior = diagm(V_all)

        return prior_mean, A_mean, B_mean, C_mean, D_mean, V_prior
    end
end

function hyper_prioreval(par, priors)
    AB_cond = 0

    for i in eachindex(par)
        AB_cond += all(insupport(priors[i], par[i]))
    end


    log_priorval = 0.0
    alarm = false
    if AB_cond == length(par)
        for i in eachindex(par)
            log_priorval += sum(logpdf(priors[i], par[i]))  # very costly unfortunately
        end
    else
        # println("not in support")
        alarm = true
        log_priorval = -1.e9
    end

    return log_priorval, alarm
end


function prioreval(par, priors, param_sizes; aggregate_rep=:as_states)
    """Evaluates the parameters at their prior distribution.

    p(A, Π) = p(A) * p(Π)
    """

    if aggregate_rep == :as_inputs
        # ─── :as_inputs ─── No aggregate persistence prior; priors layout:
        #   [1] MvNormal(AB), [2] Normal(ϕ_F), [3] LKJ(factor_count), [4...] Σ, [end-n_hyper+1:end] hyperpriors
        σ²_Ω = diag(par[2])     # factor-only Ω variances
        σ²_Σ = diag(par[4])

        AB_cond = all(insupport(priors[1], par[1]))
        Ωf_cond = all([insupport(priors[2], σ²_Ω[i]) for i in eachindex(σ²_Ω)])
        Ω_corr_cond = insupport(priors[3], par[3])
        Σ_cond = all([insupport(priors[3+i], σ²_Σ[i]) for i in eachindex(σ²_Σ)])
        hyper_cond = all([insupport(priors[end-length(par[5])+i], par[5][i]) for i in eachindex(par[5])])

        if AB_cond && Ωf_cond && Σ_cond && Ω_corr_cond && hyper_cond
            log_priorval = 0.0
            alarm = false

            log_priorval = sum(logpdf(priors[1], par[1]))
            log_priorval += sum([logpdf(priors[2], σ²_Ω[i]) for i in eachindex(σ²_Ω)])
            log_priorval += logpdf(priors[3], par[3])
            log_priorval += sum([logpdf(priors[3+i], σ²_Σ[i]) for i in eachindex(σ²_Σ)])
            log_priorval += sum([logpdf(priors[end-length(par[5])+i], par[5][i]) for i in eachindex(par[5])])
        else
            alarm = true
            log_priorval = -1.e9
        end

        return log_priorval, alarm

    else
        # ─── :as_states (default) ─── Original behavior ─────────────────────
        nf = param_sizes[1][1]
        σ²_Ω = diag(par[2])
        σ²_Ωf = σ²_Ω[1:nf]
        σ²_Ωy = σ²_Ω[nf+1:end]
        σ²_Σ = diag(par[4])

        ABCD_cond = all(insupport(priors[1], par[1]))
        Ωf_cond = all([insupport(priors[2], σ²_Ωf[i]) for i in eachindex(σ²_Ωf)])
        Ωy_cond = all([insupport(priors[3], σ²_Ωy[i]) for i in eachindex(σ²_Ωy)])
        Ω_corr_cond = insupport(priors[4], par[3])
        Σ_cond = all([insupport(priors[4+i], σ²_Σ[i]) for i in eachindex(σ²_Σ)])

        hyper_cond = all([insupport(priors[end-length(par[5])+i], par[5][i]) for i in eachindex(par[5])])

        if ABCD_cond && Ωf_cond && Ωy_cond && Σ_cond && Ω_corr_cond && hyper_cond
            # split covariance matrices into standard deviations and correlation matrices
            log_priorval = 0.0
            alarm = false

            # Log-prior value
            log_priorval = sum(logpdf(priors[1], par[1]))  # very costly unfortunately
            log_priorval += sum([logpdf(priors[2], σ²_Ωf[i]) for i in eachindex(σ²_Ωf)])
            log_priorval += sum([logpdf(priors[3], σ²_Ωy[i]) for i in eachindex(σ²_Ωy)])
            log_priorval += logpdf(priors[4], par[3])
            log_priorval += sum([logpdf(priors[4+i], σ²_Σ[i]) for i in eachindex(σ²_Σ)])

            # Log-hyperprior value
            log_priorval += sum([logpdf(priors[end-length(par[5])+i], par[5][i]) for i in eachindex(par[5])])


        else
            # println("not in support")
            alarm = true
            log_priorval = -1.e9
        end

        return log_priorval, alarm
    end
end


function estimate_measurement_VCV(MV, dimension, grid, errors_process, time_p)
    Σ = Vector{Matrix}(undef, length(MV))
    Σ_lb = Vector{Matrix}(undef, length(MV))

    @unpack tmin, tmax, year_vec, tot_periods = time_p

    for j in eachindex(MV)
        # Estimate the VCV matrix
        Σ[j] = nancov(MV[j], dims=2) # it computes the VCV over the columns, so, we transpose 

        # Returns zero instead of NaNs for unobserved variances. Replace them 
        replace!(Σ[j], -0.0 => NaN)
    end

    # Return the mean of the variances for each dataset and object // same for variance 
    prior_μ = [zeros(dimension + 1) for _ in eachindex(MV)] # one for each object 
    prior_σ = [zeros(dimension + 1) for _ in eachindex(MV)] # one for each object 

    # Define the different set of rows 
    imm_portion = (grid + (dimension - 1) * (grid - 1))
    cop_rows = 1:grid^dimension-imm_portion
    quantile_rows = [I for I in Iterators.partition(length(cop_rows)+1:length(cop_rows)+grid*dimension, grid)]
    row_sets = [cop_rows, quantile_rows...]


    # Defining the moments of the distribution of the VARIANCE of the measurement error
    for j in eachindex(prior_μ)
        diag_Σ = diag(Σ[j])
        for k in eachindex(prior_μ[j])
            # First, find the rows which are NaN
            prior_μ[j][k] = nanmean(diag_Σ[row_sets[k]])          # mean of the variances, by object  
            prior_σ[j][k] = nan_variance(diag_Σ[row_sets[k]])      # variance of the variances, by object 
        end
    end

    # Now, it may be the case that this dataset is missing some parts of the distribution, so we insert the average of the others 
    # Important for scaling later 
    # TODO: No longer necessary! 
    # obj_vec_μ = []
    # obj_vec_σ = []

    # for k in eachindex(prior_μ[1]) # the "1" doesnt matter 
    #     push!(obj_vec_μ, nanmean([prior_μ[j][k] for j in eachindex(prior_μ)]))
    #     push!(obj_vec_σ, nanmean([prior_σ[j][k] for j in eachindex(prior_σ)]))
    # end

    # for j in eachindex(prior_μ)
    #     for k in eachindex(prior_μ[j])
    #         if isnan(prior_μ[j][k])
    #             prior_μ[j][k] = obj_vec_μ[k]
    #             prior_σ[j][k] = obj_vec_σ[k]
    #         end
    #     end
    # end

    # Now we want to define the lower bound. I need to first linear interpolate the MV[j] for all periods in the model estimation 
    # t                = collect(1:1:tot_periods)

    # # Interpolate this dataset for all rows 
    # for j in eachindex(MV)
    #     LI_df        = similar(MV[j])
    #     for i in axes(MV[j], 1)
    #         mask             = (!isnan).(MV[j][i,:])
    #         if sum(mask) >= 2 
    #             interp_linear    =  linear_interpolation(t[mask], MV[j][i, :][mask], extrapolation_bc=Line())
    #             LI_df[i, :]     .=  interp_linear.(t)
    #         else 
    #             LI_df[i, :]     .=  MV[j][i, :]
    #         end
    #     end
    #     LI_df                =  convert.(Float64, LI_df) # otherwise, things will be wrong numerically (since it initially returns a simpleratio type)

    #     # compute VCV of this object 
    #     Σ_lb[j] = nancov(LI_df, dims=2) # it computes the VCV over the columns, so, we transpose so we get VCV of rows 
    # end

    # # Return the mean of the variances for each dataset and object // same for variance 
    # prior_μ_lb = [zeros(dimension+1) for _ in eachindex(MV)] # one for each object 

    # # Defining the moments of the distribution of the VARIANCE of the measurement error
    # for j in eachindex(prior_μ_lb)
    #     for k in eachindex(prior_μ_lb[j])
    #         # First, find the rows which are NaN
    #         prior_μ_lb[j][k] = nan_mean(diag(Σ_lb[j])[row_sets[k]])          # mean of the variances, by object  
    #     end
    # end

    if errors_process == "one per object, per dataset"
        # For each dataset (j), remove NaNs from vector 
        # for j in eachindex(MV)
        #     mask = .!isnan.(prior_μ[j])
        #     prior_μ[j] = prior_μ[j][mask]
        #     prior_σ[j] = prior_σ[j][mask]
        # end

        # For each dataset (j), replace NaNs with 1000
        # for j in eachindex(MV)
        #     replace!(prior_μ[j], NaN=>1000) # As a flag, but will not be in the estimation 
        #     replace!(prior_σ[j], NaN=>1000)
        # end

        # return [prior_μ, prior_σ, prior_μ_lb]
        return [prior_μ, prior_σ]

    elseif errors_process == "average"
        # For now, return the mean of the means and the mean of the variances 
        prior_μ_means = nan_mean(reduce(hcat, prior_μ))
        prior_σ_means = nan_mean(reduce(hcat, prior_σ))
        prior_μ_lb_means = nan_mean(reduce(hcat, prior_μ_lb))

        return [prior_μ_means, prior_σ_means, prior_μ_lb_means]
    end

end

function extract_measurement_error_params(MV, dimension, grid, errors_process)
    # Variance of each row, for each dataset 
    σ_cont = zeros(size(MV[1], 1), length(MV))

    for j in eachindex(MV)
        for i in axes(MV[j], 1)
            mask = .!isnan.(MV[j][i, :])
            σ_cont[i, j] = var(MV[j][i, mask])
        end

        # # Plot the variance distribution for each dataset 
        # p = barhist(σ_cont[i, :], bins=5)
        # Plots.savefig(p, "variance_distribution_dataset_$j.png")
    end

    # Return the mean of the variances for each dataset and object // same for variance 
    prior_μ = [zeros(dimension + 1) for _ in eachindex(MV)] # one for each object 
    prior_σ = [zeros(dimension + 1) for _ in eachindex(MV)] # one for each object 

    # Define the different set of rows 
    imm_portion = (grid + (dimension - 1) * (grid - 1))
    cop_rows = 1:grid^dimension-imm_portion
    quantile_rows = [I for I in Iterators.partition(length(cop_rows)+1:length(cop_rows)+grid*dimension, grid)]
    row_sets = [cop_rows, quantile_rows...]


    # Defining the moments of the distribution of the VARIANCE of the measurement error
    for j in eachindex(prior_μ)
        for k in eachindex(prior_μ[j])
            # First, find the rows which are NaN
            prior_μ[j][k] = nan_mean(σ_cont[row_sets[k], j])          # mean of the variances, by object  
            prior_σ[j][k] = nan_variance(σ_cont[row_sets[k], j])      # variance of the variances, by object 
        end
    end

    if errors_process == "one per object, per dataset"
        # For each dataset (j), remove NaNs from vector 
        # for j in eachindex(MV)
        #     mask = .!isnan.(prior_μ[j])
        #     prior_μ[j] = prior_μ[j][mask]
        #     prior_σ[j] = prior_σ[j][mask]
        # end

        # For each dataset (j), replace NaNs with 1000
        # for j in eachindex(MV)
        #     replace!(prior_μ[j], NaN=>1000) # As a flag, but will not be in the estimation 
        #     replace!(prior_σ[j], NaN=>1000)
        # end

        return [prior_μ, prior_σ]

    elseif errors_process == "average"
        # For now, return the mean of the means and the mean of the variances 
        prior_μ_means = nan_mean(reduce(hcat, prior_μ))
        prior_σ_means = nan_mean(reduce(hcat, prior_σ))

        return [prior_μ_means, prior_σ_means]
    end
end

function nan_mean(x)
    mask = .!isnan.(x)
    return mean(x[mask])
end

function nan_variance(x)
    mask = .!isnan.(x)
    return var(x[mask])
end

function VAR_fit(y::AbstractArray, p::Int64)
    (T, K) = size(y)
    X = y
    y = transpose(y)
    Y = y[:, p+1:T]
    X = lagmatrix(X, p)'
    β = (Y * X') / (X * X')
    ϵ = Y - β * X
    Σ = ϵ * ϵ' / (T - p * K - 1)
    return diag(Σ)
end

function lagmatrix(x::AbstractArray, p::Int64)
    sk = 1
    T, K = size(x)
    k = K * p + 1
    idx = repeat(1:K, p)
    X = Array{eltype(x)}(undef, (T - p, k))
    # building X (t-1:t-p) allocating data from D matrix - avoid checking bounds
    for j = 1+sk:(sk+K*p)
        for i = 1:(T-p)
            lg = round(Int, ceil((j - sk) / K)) - 1 # create index [0 0 1 1 2 2 ...etc]
            @inbounds X[i, j] = x[i+p-1-lg, idx[j-sk]]
        end
    end
    return X[:, 2:end]
end


#         # Independent Minnesota-Wishart Prior 
# VARIANCE = kron(inv(V.Σ), Matrix(I, T, T))  
# V_post   = inv(inv(V_prior) + Z' * VARIANCE * Z)
# V_post   = nearest_spd(V_post)
#     # show(IOContext(stdout, :limit=>false), MIME"text/plain"(), V_post)        
#     # show(IOContext(stdout, :limit=>false), MIME"text/plain"(), inv(V_prior))
#     # println(all(isapprox(V_post, V_post'; rtol=1e-4)))
# a_post = V_post * (inv(V_prior) * A_prior[:] + Z' * VARIANCE * endogenous[:, 1+lags:end][:])

# Minnesota posterior 
# post_Σ = inv((inv(V_prior) + (kron(inv(V.Σ), X' * X))))  # TODO: This is the Minnesota posterior 
# post_Γ =  post_Σ * (inv(V_prior) * A_prior[:] + kron(inv(V.Σ), X)' * data["endogenous"][:])


#         # # Minnesota-Wishart posterior 
#         # VARIANCE = kron(inv(V.Σ), Matrix(I, T, T))  
#         # post_Σ = inv(inv(V_prior) + Z' * VARIANCE * Z))  
#         # post_Γ =  post_Σ * (inv(V_prior) * A_prior[:] + Z' * VARIANCE * data["endogenous"][:]))




# function var_estimation(endogenous, lags, constant, T)
#         var_Ω = Vector{Vector{Float64}}(undef, length(endogenous))
#         for (b, block) in enumerate(endogenous)
#             replace!(block, NaN=>0.2)  # any small number
#             # Data 
#             X = lagmatrix(block', lags, constant) # copy(V.X) # N by T-1  #TODO: T by N now 
#             Y = view(block, :, 1+lags:size(block,2))'
#             # Z = kron(Matrix(I, state_count, state_count), X)   # stacking the control variables, one per equation 

#             # Variances from the residuals of a least-squares regression
#             A_OLS = inv(X' * X) * (X' * Y)
#             var_Ω[b] = diag(((Y - X * A_OLS)' * (Y - X * A_OLS)) ./ (T - size(block, 1) + 1))  
#         end
#         return vcat(var_Ω...)
# end

