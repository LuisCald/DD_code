# EM Step 
# https://bpb-us-e1.wpmucdn.com/sites.uw.edu/dist/4/10826/files/2023/01/NowcastingBanburaGiannoneReichlin_CorrectedAppendix.pdf
# https://arxiv.org/pdf/1910.03821

# Running the Kalman filter and smoother with current parameters
# %  Description:
# %    EMstep reestimates parameters based on the Estimation Maximization (EM)
# %    algorithm. This is a two-step procedure:
# %    (1) E-step: the expectation of the log-likelihood is calculated using
# %        previous parameter estimates.
# %    (2) M-step: Parameters are re-estimated through the maximisation of
# %        the log-likelihood (maximize result from (1)).
# %

# %
# %  Input:
# %    y:      Series data
# %    A:      Transition matrix
# %    proj:   Observation matrix
# %    Ω:      Covariance for transition equation residuals
# %    Δ:      Covariance for observation matrix residuals
# %    Z_0:    Initial values of factors
# %    V_0:    Initial value of factor covariance matrix
# %    r:      Number of common factors for each block (e.g. vector [1 1 1 1])
# %    p:      Number of lags in transition equation
# %    R_mat:  Estimation structure for quarterly variables (i.e. "tent")
# %    q:      Constraints on loadings
# %    nQ:     Number of quarterly series
# %    i_idio: Indices for monthly variables
# %    blocks: Block structure for each series (i.e. for a series, the structure
# %            [1 0 0 1] indicates loadings on the first and fourth factors)
# %
# %  Output:
# %    proj_new: Updated observation matrix
# %    Δ_new: Updated covariance matrix for residuals of observation matrix
# %    A_new: Updated transition matrix
# %    Ω_new: Updated covariance matrix for residuals for transition matrix
# %    Z_0:   Initial value of state
# %    V_0:   Initial value of covariance matrix
# %    loglik: Log likelihood

# %    x_smoothed: k-by-(nobs+1) matrix, smoothed factor estimates
# %             (i.e. x_smoothed(:,t+1) = Z_t|T)
# %    sigma_smoothed: k-by-k-by-(nobs+1) array, smoothed factor covariance matrices
# %             (i.e. sigma_smoothed(:,:,t+1) = Cov(Z_t|T))
# %    Vsigma_smoothed: k-by-k-by-nobs array, lag 1 factor covariance matrices
# %              (i.e. Cov(Z_t,Z_t-1|T))

#TODO: add controls!
# Initialize

function check_convergence(likelihoods::Vector{Float64}, tol::Float64=1e-4, patience::Int=10)
    # Return false if we don't have enough iterations to evaluate
    if length(likelihoods) < patience + 1
        return false
    end

    # Check if the change in likelihood is below the tolerance over 'patience' iterations
    ll = length(likelihoods)
    recent_changes = [likelihoods[i] - likelihoods[i-1] for i in (ll - patience + 1):ll]
    println(recent_changes)
    if all(abs(change) < tol for change in recent_changes)
        println("Algorithm is stuck: minimal improvement in likelihood.")
        return true
    end

    # Check for oscillating pattern
    diffs = diff(likelihoods[end-(patience):end])
    if all(diffs[i] * diffs[i-1] < 0 for i in 2:length(diffs))
        println("Algorithm is stuck: oscillating likelihood values detected.")
        return true
    end

    return false  # Not stuck if above conditions are not met
end


function run_EM_algorithm(par_final, param_sizes, meas_ind, Σ_ids, model_elements, model_options)
    @unpack y, G, u = model_elements
    @unpack case, number_of_dfs = model_options
    A, B, Ω, Σ                  = matrisize(par_final, param_sizes, case)
    Σ                           = apply_measurement_criteria(A, Σ, meas_ind, Σ_ids, model_options)
    Σ                          .= log.(exp.(Σ) .+ 1)  # softplus transformation
    Ω[diagind(Ω)] = log.(exp.(Ω[diagind(Ω)]) .+ 1)  # softplus transformation 

    A_new = deepcopy(A) # state equation LOM
    B_new = deepcopy(B) # Controls
    Ω_new = deepcopy(Ω) # state equation VCV
    Δ_new = deepcopy(Σ) # measurement equation VCV
    G_new = deepcopy(G) # measurement equation LOM
    y_new = deepcopy(y)

    A_old = deepcopy(A) # state equation LOM
    B_old = deepcopy(B) # Controls
    Ω_old = deepcopy(Ω) # state equation VCV
    Δ_old = deepcopy(Σ) # measurement equation VCV
    G_old = deepcopy(G) # measurement equation LOM
    y_old = deepcopy(y)

    n_iter_ML = 200
    log_D_vec = [] #zeros(n_iter_ML)
    Δ_log     = []
    decrease_count = 0

    for i in 1:n_iter_ML
        # Run Kalman-Filter Smoother
        # log_P, alarm               = prioreval([[A_new..., B_new...], Matrix(Ω_new), diagm(Float64.(Δ_new))], priors, case) # TODO: see why it doesnt work now.
        smoother_output, log_D, alarm             = recurse_kalman_filter(A_new,B_new,Ω_new,Δ_new,G_new,y,u,true)
        println("iteration $i: ", log_D)

        if alarm 
            return A_old, B_old, Ω_old, Δ_old, G_old, log_D_vec, Δ_log
        else
            A_old = deepcopy(A_new)
            B_old = deepcopy(B_new)
            Ω_old = deepcopy(Ω_new)
            Δ_old = deepcopy(Δ_new)
            G_old = deepcopy(G_new)
            y_old = deepcopy(y_new)
        end
        
        push!(log_D_vec, log_D)
        if i > 1 

            # Check if likelihood increased
            if log_D_vec[i] < log_D_vec[i-1]
                decrease_count += 1
            else
                decrease_count = 0
            end

            # Stop if likelihood increased 5 times in a row
            if decrease_count >= 5
                println("Stopping early: Likelihood increased 5 times in a row at iteration $i")
                return A_new, B_new, Ω_new, Δ_new, G_new, log_D_vec, Δ_log
            end
            
            delta_loglik = abs(log_D_vec[i] - log_D_vec[i-1])                                
            avg_loglik = (abs(log_D_vec[i] + log_D_vec[i-1]) + 10^(-3)) / 2 
            push!(Δ_log, delta_loglik / avg_loglik)
            if (delta_loglik / avg_loglik < 1e-6) 
                @info "Convergence reached after $i iterations"
                return A_new, B_new, Ω_new, Δ_new, G_new, log_D_vec, Δ_log
            elseif i == n_iter_ML
                @info "Stopped after $i iterations"
                return A_new, B_new, Ω_new, Δ_new, G_new, log_D_vec, Δ_log
            elseif i % 10 == 0
                # conv_flag = check_convergence(log_D_vec)
                # if conv_flag == true 
                #     @info "Stuck"

                    # Perturb the parameters a bit
                    Δ_new = Δ_new .- 1e-2 .* randn(size(Δ_new))

                # end                
            end
        end

        @unpack x_smoothed, sigma_smoothed, cross_covariances = smoother_output 
        T = size(y, 2)

        # Apply changes 
        y_new = deepcopy(y)
        replace!(y_new, NaN=>0)

        # Convert sigma_smoothed (vector of matrices) to a 3 dimensional object 
        sigma_smoothed    = cat(sigma_smoothed..., dims=3)
        cross_covariances = cat(cross_covariances..., dims=3)

        fₜ   = x_smoothed[:, 2:end]
        fₜ₋₁ = x_smoothed[:, 1:end-1]
        Vₜ   = sigma_smoothed[:, :, 2:end]
        Vₜ₋₁ = sigma_smoothed[:, :, 1:end-1]
        uₜ   = u[:, 2:end]
        
        # Compute components of the E-step
        # E[f_t*f_t' | y_T]
        E_FF = zeros(size(fₜ, 1), size(fₜ, 1))
        for t in 1:T-1
            E_FF += fₜ[:, t] * fₜ[:, t]' + Vₜ[:, :, t]
        end
        
        # E[f_{t-1}*f_{t-1}' | y_T]
        E_FF_lag = zeros(size(E_FF))
        for t in 1:T-1
            E_FF_lag += fₜ₋₁[:, t] * fₜ₋₁[:, t]' + Vₜ₋₁[:, :, t]
        end
        # E_FF = fₜ * fₜ' + sum(Vₜ, dims=3)[:,:] # 2:end because we have t-1 later, so, 1-1 is not possible as an index     
        # E_FF_lag2 = fₜ₋₁ * fₜ₋₁' + sum(Vₜ₋₁, dims=3)[:,:]
        # E[f_t*f_{t-1}' | y_T]
        # E_FF_mix = fₜ * fₜ₋₁' + sum(cross_covariances[:,:,2:end], dims=3)[:,:] 

        E_FF_mix = zeros(size(fₜ, 1), size(fₜ₋₁, 1))
        for t in 1:T-1
            E_FF_mix += fₜ[:, t] * fₜ₋₁[:, t]' + cross_covariances[:, :, t+1]
        end
        
        # Maximization step 
        try
            A_new = E_FF_mix * inv(E_FF_lag) #TODO: breaks here sometimes
        catch ee
            println("here A")
            A_new = E_FF_mix * inv(nearest_spd(E_FF_lag)) # TODO: still doesnt work
        end
        # uu    = uₜ * uₜ' # this is just a diagonal matrix, where each element is equal to size(uₜ, 2) 
        uu    = diagm(ones(size(uₜ, 1)) .* size(uₜ, 2)) # this is just a diagonal matrix, where each element is equal to size(uₜ, 2)
        B_new = ((fₜ - A_new * fₜ₋₁) * uₜ') * inv(uu) # Depends on new 'A_new'
        # n_f = size(fₜ₋₁, 1)
        # n_u = size(uₜ, 1)

        # X = vcat(fₜ₋₁, uₜ)
        # U = zeros(n_f + n_u, n_f + n_u)
        # U[1:n_f, 1:n_f] = sum(Vₜ₋₁, dims=3)[:,:,1]

        # β = inv(X' * X + U) * X' * fₜ
        # B_new = β[:, size(fₜ₋₁, 1)+1:end]
        # println(size(B_new))
        # B_new = fₜ * uₜ' ./ (T-1) 

        Ω_new = diagm(diag((E_FF - A_new * E_FF_mix' - B_new * uu * B_new') ./ (T-1))) # Depends on new 'A_new'
        # Ω_new = diagm(diag((E_FF - A_new * E_FF_mix') ./ (T-1))) # Depends on new 'A_new'
        # Ω_new = (E_FF - A_new * E_FF_mix' - B_new * uu * B_new') ./ (T-1) # Depends on new 'A_new'
        
        # Some parameters  
        n       = size(y, 1)  # Number of observations (rows)
        T       = size(y, 2)  # Number of time periods (columns)
        xₜ_size = size(x_smoothed, 1)

        denom  = zeros(n * xₜ_size, n * xₜ_size)
        num    = zeros(n, xₜ_size) 

        pb = Progress(T, desc = "M Step 4")
        for i in 1:T
            # select_mat   = diagm((!isnan).(y[:, i]))
            select_mat   = (!isnan).(y[:, i])

            if sum(select_mat) != 0
                r_x          = reshape(x_smoothed[:, i], xₜ_size, 1)
                block        = r_x * r_x' + sigma_smoothed[:, :, i] #TODO: 
            
                # Instead of kron(block, select_mat), apply block-wise computation -- since this is a banded matrix
                for j in 1:n
                    if select_mat[j] != 0  # Only process blocks where select_mat is non-zero
                        denom[(j-1) * xₜ_size+1:j * xₜ_size, (j-1) * xₜ_size+1:j * xₜ_size] += block 
                    end
                end

                # denom       += kron(block, select_mat)
                num         += y_new[:, i] * r_x' 
            end

            # nom = nom + ...  E[y_t*f_t' | \Omega_T]
            # y(idx_iM, t) * Zsmooth(bl_idxM(i, :), t+1)' - ...
            # Wt(:, i_idio_i) * (Zsmooth(rp1 + i_idio_ii, t+1) * Zsmooth(bl_idxM(i, :), t+1)' + Vsmooth(rp1 + i_idio_ii, bl_idxM(i, :), t+1));

            if i == T
                println("inverting")
                # invert 'denom' by block
                for j in 1:n
                    if sum(select_mat[j]) != 0
                        # Invert sum of all blocks
                        denom[(j-1) * xₜ_size+1:j * xₜ_size, (j-1) * xₜ_size+1:j * xₜ_size] = inv(denom[(j-1) * xₜ_size+1:j * xₜ_size, (j-1) * xₜ_size+1:j * xₜ_size] .+ 1e-6)
                    end
                end
            end
            # next!(pb)
        end
        G_mat = reshape(denom * num[:], n, xₜ_size) 

        # # Adopting the approach of Matteo Barigozzi and Luciani (2024) and doing this row by row 
        # G_mat       = zeros(n, xₜ_size)
        # denom       = zeros(xₜ_size, xₜ_size)
        # for t in 1:T
        #     denom += x_smoothed[:, t] * x_smoothed[:, t]' + sigma_smoothed[:, :, t]
        # end
        # inv_denom  = inv(denom)
        # # denom       = x_smoothed * x_smoothed' + sum(sigma_smoothed, dims=3)[:,:] #* (!isnan)(y[i, t]) # select_mat 
        # for i in 1:n # Slow, since not column wise 
        #     num         = zeros(xₜ_size, 1)
        #     for t in 1:T
        #         r_x          = reshape(x_smoothed[:, t], xₜ_size, 1)
        #         num         += r_x * y_new[i, t] #* (!isnan)(y[i, t]) # select_mat
        #     end
            
        #     try
        #         G_mat[i, :] = vec(inv_denom * num)
        #     catch ee
        #         println("here2")
        #         G_mat[i, :] = vec(inv(nearest_spd(denom)) * num)
        #     end
        #     next!(pb)
        # end

        Threads.@threads for i in eachindex(G_old)
            G_new[i] = G_mat # Can just have 1 G_new
        end

        # Pre-allocate Δ_temp and other necessary matrices
        Δ_temp      = zeros(n)
        
        # Process observation-by-observation (row-by-column)
        # Make a progress bar 
        pb = Progress(T, desc = "M Step 5")
        for t in 1:T
            E_FF_lag_t = x_smoothed[:, t] * x_smoothed[:, t]' + sigma_smoothed[:,:,t]

            for i in 1:n
                miss  = isnan(y[i, t])
                ## Implementation from Nowcasting code, also with missings 
                r_G = reshape(G_mat[i, :], 1, xₜ_size)

                # # If the observation y[i, t] is not NaN, compute the update
                # residual = y_new[i, t] - (r_G * reshape(x_smoothed[:, t], 7, 1))[1]
                # Δ_temp[i] += (residual * residual)[1]
                
                # # Add the G_new * Vₜ[:,:,t] * G_new' contribution for this observation
                # Δ_temp[i] += (r_G * sigma_smoothed[:, :, t] * r_G')[1]
                
                # # Add the condition matrix part using diagm_Δ_new
                # Δ_temp[i] += miss * Δ_new[i]

                ## Implementation from Banbura and Modugno (2014), but observation by observation: SAME AS ABOVE IN LIKELIHOOD, despite slightly diff. formulation, BUT SLOWER
                # Δ_temp[i] += (y_new[i, t]^2 .- r_G * x_smoothed[:, t] * y_new[i, t] .- y_new[i, t] * x_smoothed[:, t]' * r_G' .+ r_G * E_FF_lag_t * r_G')[1]
                Δ_temp[i] += (y_new[i, t]^2 .-  2 * y_new[i, t] * x_smoothed[:, t]' * r_G' .+ r_G * E_FF_lag_t * r_G')[1] # from Bargozzi 2024
                Δ_temp[i] += miss * Δ_old[i]
            end
            # next!(pb)
        end

        Δ_new = Δ_temp ./ T
    end

    return A_new, B_new, Ω_new, Δ_new, G_new, log_D_vec, Δ_log
end


# # Concatenate lagged states and control variables
# X = vcat(fₜ₋₁, uₜ)  # (k + m) x (T-1)

# # OLS estimation: [A, B] = F * X' * inv(X * X')
# A_B = fₜ * X' * inv(X * X')

# # Extract A and B
# k = size(fₜ, 1)  # Number of state variables
# m = size(uₜ, 1)  # Number of control variables

# A_new = A_B[:, 1:k]  # A is k x k
# B_new = A_B[:, k+1:end]  # B is k x m
# Ω_new = diagm(diag((E_FF - A_new * E_FF_mix') ./ (T-1))) # Depends on new 'A_new'


## Implementation from Banbura and Modugno (2014), matrix form
# Δ_temp      = zeros(n, n)
# Δ_new_mat   = diagm(Float64.(Δ_new))

# for t in 1:T
#     miss    = (isnan).(y[:, t])
#     miss_m  = diagm(miss)
#     Δ_temp += y[:, t] * y[:, t]' - G_new[t] * x_smoothed[:, t] * y[:, t]' - y[:, t] * x_smoothed[:, t]' * G_new[t]' + G_new[t] * E_FF_lag[:, :, t] * G_new[t]'
#     Δ_temp += miss_m * Δ_new_mat * miss_m
#     next!(pb)
# end

# Normalize by the number of time periods
# Δ_new = diag(Δ_temp) ./ T

# Generate new param_vector
# par_final = vectorize(A_new, B_new, Ω_new, Δ_new)

# opttag = "from_mcmc"
# @unpack tmin, tmax = time_params 
# @unpack gdp_series = obs_data 
# @unpack data_sources = func_data
# user_t     = (deepcopy(tmin), deepcopy(tmax))

# dv, _   = reconstruct_data_short(A_new,B_new,Ω_new,Δ_new,G_new, model_elements, obs_data, model_options, time_params, data_sources; reconstruction_to_show=false, dε_smoothed=false)

# within_stat_dict = Dict()
# for (c, k) in enumerate(keys(dv))
#     within_stat_dict[k], dv[k] = export_functional_data(dv[k], k, kind_of_plots, obs_data, func_data, time_params, user_t, model_options, false, true)        
#     if c == length(keys(dv))
#         compare_to_data(dv, func_data, obs_data, user_t, time_params, model_options, kind_of_plots, label)
#         init_path = BASE_PATH
#         dict_path = init_path * "/7_Results/$m_label" * "$tag" * "/$opttag/plots/"
#         export_combined_stat_dict_to_latex(within_stat_dict, [measures..., "copula"], dict_path, label) # [measures..., "copula"]
#         compare_to_external_sources(dv, func_data, obs_data, user_t, time_params, model_options, kind_of_plots, label)
#     end
# end

# A = par_final["matrices"][1]
# B = par_final["matrices"][2]
# Ω = par_final["matrices"][3]
# Σ = par_final["matrices"][4]
# G = par_final["matrices"][5]


function reconstruct_data_short(A_new,B_new,Ω_new,Δ_new,G_new, likeli_vec, Δ_log, model_elements, obs_data, model_options, time_params, data_sources; reconstruction_to_show=false, dε_smoothed=false)
    """Reverse all steps taken in the data preparation."""
    
        @unpack lags, blind_to, case, pre_multiply, estimator, tag    = model_options 
        @unpack y, factor_count, n_less_than_one, u, pcs, βs, trend = model_elements
        
        @unpack gdp_series, df_vec = obs_data
        @unpack tot_periods, tmin, tmax = time_params
        
        number_of_dfs = length(df_vec.data)
        data_names    = df_vec[2]

        # Save the matrices first
        @unpack measures, data_cutoffs = model_options
        println(tag)
        label = "3D_A non-diag"
        m_label   = measures_folder(measures)
        init_path = BASE_PATH
        file_name = init_path * "/7_Results/$m_label" * "$tag" * "/from_mcmc/parameter_vectors/solution" * "$label" * ".jld2" 
        
        JLD2.save(file_name, "matrices", [A_new,B_new, Ω_new, Δ_new, G_new, likeli_vec, Δ_log])

        Plots.scatter(likeli_vec[2:end], msc=:grey, mc=:white, xformatter=:latex, yformatter=:latex, label="", ylabel=L"\textrm{Log\,\,Likelihood}", xlabel=L"\textrm{Iterations}", legend=:best)
        Plots.plot!(likeli_vec[2:end], lc=:blue, label=L"\textrm{Convergence\,\,Path}")
        mkpath(init_path * "/7_Results/$m_label" * "$tag" * "/other_results")
        Plots.savefig(init_path * "/7_Results/$m_label" * "$tag" * "/other_results/likelihood_convergence" * ".pdf")

        Plots.scatter(Δ_log[2:end], msc=:grey, mc=:white, xformatter=:latex, yformatter=:latex, label="", ylabel=L"\textrm{\Delta\,\,\mathcal{l}_{\bar{k}}}", xlabel=L"\textrm{Iterations}", legend=:best)
        Plots.plot!(Δ_log[2:end], lc=:orange, label=L"\textrm{Convergence\,\,Path}")
        Plots.savefig(init_path * "/7_Results/$m_label" * "$tag" * "/other_results/delta_log" * ".pdf")
    
        smoother_output, _, _             = recurse_kalman_filter(A_new,B_new,Ω_new,Δ_new,G_new,y,u,true)
        @unpack x_smoothed, dε_smoothed  = smoother_output
        
        @unpack proj, means, stds, agg_count = model_elements
        @unpack measures                     = model_options
        @unpack grid_pcf, grid_cop           = estimator
        dimension = length(measures)
    
        # Plot the smoothed data
        Plots.plot()
        dts       = QuarterlyDate(tmin["year"], tmin["quarter"]) : Quarter(1) : QuarterlyDate(tmax["year"], tmax["quarter"])
        xaxis     = 1:tot_periods
        init_path = BASE_PATH 
    
        for i in axes(x_smoothed, 1)
            Plots.plot!(
                xaxis, 
                x_smoothed[i, :], 
                xformatter=:latex, 
                yformatter=:latex, 
                xticks=(xaxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_,date) in enumerate(dts[1:20:end])]),
                ylabel = L"\textrm{Factor\,\, Value}",
                label=L"\textrm{Factor \,\, %$(i)}", 
                legend=:best) #, xticks=(collect(1:20:tot_periods), [L"%$(date)"[1:5] * "\$" for (_,date) in enumerate(gdp_series[1:20:end, :time])]))
        end
        mkpath(init_path * "/7_Results/factor_analysis")
        Plots.savefig(init_path * "/7_Results/factor_analysis/distributional_factors $tag.pdf")
        # TODO: undo this later
        # condition = isnan.(means[1])
        # means[1][condition] .= 0 #TODO: maybe take the mean from PSID?
    
        # Reconstruct dataset by dataset and then average  
        ΓF   = (proj * x_smoothed) 
        X    = Vector{Matrix{Float64}}(undef, number_of_dfs)
        T    = size(ΓF, 2)
    
        # Split by object, multiply by stds, reform object 
        ΓF_σ = add_variance(estimator, ΓF, stds, measures)
        X    = add_mean!(X, ΓF_σ, means, data_names, blind_to)
    
        # What kind of trend to add back?
        new_trend = select_trend(trend, "average")
    
        X_dict = Dict()
        X_dict["normal"] = deepcopy(X)
        X_dict["average"] = deepcopy(X)

        # Add trend back for each 
        t = collect(1:tot_periods)
        for j in eachindex(X)
            for t in 1:tot_periods
                # X[j][:, t] = X[j][:, t] .+ βs[:,1,j] .+ βs[:,2,j] .* t 
                # X[j][:, t] = X[j][:, t] .+ βs[:,1,j] .+ βs[:,2,j] .* t .+ βs[:,3,j] .* (2 .* t.^2 ./ tot_periods) # trend[j][:, t] #
            # For copula, weights do not seem to be affected
            X_dict["normal"][j][:, t]  .+= trend[j][:, t] #
            X_dict["average"][j][:, t] .+= new_trend[j] # By adding the nanmean, we can generate estimates for all periods ... to check is how similar are the copula estimates across both
            end
        end
           # Create consensus average. 
    X_array_norm = cat(X_dict["normal"]..., dims=3)
    X_array_avg  = cat(X_dict["average"]..., dims=3)

    X̄_norm            = similar(X_dict["normal"][1])
    X̄_avg             = similar(X_dict["average"][1])
    
    # Take the average of the 3 dimensions, skipping over NaN, element by element
    for i in axes(X_dict["normal"][1], 1)
        for j in axes(X_dict["normal"][1], 2)
            X̄_norm[i, j] = mean(filter(!isnan, X_array_norm[i, j, :]))
            X̄_avg[i, j]  = mean(filter(!isnan, X_array_avg[i, j, :]))
        end
    end

    # Create Dictionary to keep track 
    X̄̄ = Dict() 

    # Create keys via a loop, data sources is already sorted as well as the reconstructed data
    for i in 1:number_of_dfs
        X̄̄[data_sources[i]] = [X_dict["normal"][i], X_dict["average"][i]]
    end
    
    X̄̄["consensus"] = [X̄_norm, X̄_avg]

    # Add the immutable part
    for (k, v) in X̄̄
        X̄̄[k][1]        = add_multidimensional_immutable(estimator, v[1], grid_cop, measures)
        X̄̄[k][2]        = add_multidimensional_immutable(estimator, v[2], grid_cop, measures)
    end

    X̄̄_new = Dict()

    for (k, v) in X̄̄
        X̄̄_new[k]       = Dict()
        X̄̄_new[k]["normal"]  = undo_functional_treatment(estimator, v[1], measures) 
        X̄̄_new[k]["average"]   = undo_functional_treatment(estimator, v[2], measures) 
    end
    
    # n_out_of_bounds  = data_diagnostics(copulas, measures)
    if reconstruction_to_show == false 
        return X̄̄_new, dε_smoothed
    else
        sub_X̄̄= Dict()

        for k in unique(["consensus", reconstruction_to_show])
            sub_X̄̄[k] = X̄̄_new[k]
        end

        return sub_X̄̄, dε_smoothed
    end
end
    

function convert_Δ_new_to_vec(Δ_new, model_options, model_elements)
    @unpack MV = model_elements
    @unpack estimator, measures = model_options
    @unpack grid_pcf = estimator

    dim        = length(measures)
    cop_part, imm_part = retrieve_cop_and_imm_part(estimator, dim)
    cop_id     = cop_part - imm_part


    Δ_to_export = []
    nd          = length(MV)
    partish     = Int(length(Δ_new) / nd)

    # Split Δ_new in seven parts 
    Δ_new_split = [Δ_new[I, :] for I in Iterators.partition(axes(Δ_new, 1), partish)]
    
    j          = 1
    for i in eachindex(MV)
        # Find all rows associated with the copula
        cop_rows = MV[i][1:cop_id, :]

        # Find rows associated non-NaN 
        # cond = all(!isnan, cop_rows, dims=2)[:]
        rows_with_non_nan = findall(row -> any(!isnan, cop_rows[row, :]), axes(cop_rows, 1))
        println(rows_with_non_nan)
        me = [mean(Δ_new_split[i][rows_with_non_nan])]
        println(me)
        append!(Δ_to_export, me)

        start = cop_id
        for k in eachindex(measures)
            # Are the percentile functions observed?
            j += 1
            cond1 = (!isnan).(MV[i][start+1:start+grid_pcf, :])
            if any(cond1)
                me = [mean(Δ_new_split[i][start+1:start+grid_pcf])]
                println(me)
                append!(Δ_to_export, me)
            end
            start += grid_pcf
        end
        j += 1
    end

    filter!(x -> !isnan(x), Δ_to_export)

    return Δ_to_export
end
