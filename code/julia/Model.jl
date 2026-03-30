function set_params(model_elements::StateSpaceModel, time_p::TimeParams, model_options::ModelOptions)
    """Generates the initial draw of parameters and defines their prior distributions."""

    prior_objects, meas_ind = estimate_prior(model_elements, model_options)
    param_vector, param_sizes, Σ_ids, priors = define_parameter_space(model_elements, model_options, prior_objects)

    return param_vector, param_sizes, priors, meas_ind, Σ_ids
end


function zero_out_A(A, n_less_than_one)
    factor_count = size(A, 1)
    n_imp = factor_count - n_less_than_one

    top_A = vec(A[1:n_imp, 1:n_imp])
    bot_A = diag(A)[n_imp+1:factor_count]

    new_A = vcat(top_A, bot_A)

    return new_A
end


function get_rows_of_copulas(G, measures)
    """G = grid_cop, D = dimension """

    # Find number of points that correspond to each copula 
    combos = generate_unique_combinations(measures)
    rows_v = Vector{Vector{Int}}(undef, length(combos))

    start = 0
    for (i, c) in enumerate(combos)
        d = length(c)
        rows_v[i] = (start+1):(start+G^d)  # number of sub copulas for d: binomial(D, d)
        start += G^d
    end

    return rows_v
end


function points_in_copulas(estimator, G, measures)

    dimension = length(measures)

    if typeof(estimator) <: SeriesEstimator
        # Find number of points that correspond to each copula 
        combos = generate_unique_combinations(measures)
        points_in_cop = Vector{Int}(undef, length(combos))

        for (i, c) in enumerate(combos)
            d = length(c)
            points_in_cop[i] = G^d  # number of sub copulas for d: binomial(D, d)
        end

        return points_in_cop
    else
        cop_part, imm_part = retrieve_cop_and_imm_part(estimator, dimension)
        return cop_part - imm_part
    end
end



# function find_observed_objects!(estimator::SeriesEstimator, Σ_ids, MV, cop_id, measures)
#     @unpack grid_cop, grid_pcf = estimator
#     j          = 1

#     # Get the rows for the copula 
#     rows_v = get_rows_of_copulas(grid_cop, measures)

#     for df in eachindex(MV)
#         # Which copulas and subcopulas are observed?
#         for rows in rows_v
#             if all(isnan, MV[df][rows, :])
#                 Σ_ids[j] = NaN
#             end
#             j += 1
#         end

#         start = cop_id

#         for k in eachindex(measures)
#             # Are the percentile functions observed?
#             if all(isnan, MV[df][start+1:start+grid_pcf, :])
#                 Σ_ids[j] = NaN
#             end
#             start += grid_pcf
#             j += 1
#         end
#     end
# end

function find_observed_objects!(estimator::SeriesEstimator, Σ_ids, MV, cop_id, measures)
    @unpack grid_pcf = estimator
    j = 1
    for i in eachindex(MV)
        # Is the copula observed?
        if all(isnan, MV[i][1:cop_id, :])
            Σ_ids[j] = NaN
        end

        start = cop_id
        for k in eachindex(measures)
            # Are the percentile functions observed?
            j += 1
            if all(isnan, MV[i][start+1:start+grid_pcf, :])
                Σ_ids[j] = NaN
            end
            start += grid_pcf
        end
        j += 1
    end

    # For the aggregates, they are observed, so, do nothing
end


function find_observed_objects!(estimator::HistogramEstimator, Σ_ids, MV, cop_id, measures)
    @unpack grid_pcf = estimator
    j = 1
    for i in eachindex(MV)
        # Is the copula observed?
        if all(isnan, MV[i][1:cop_id, :])
            Σ_ids[j] = NaN
        end

        start = cop_id
        for k in eachindex(measures)
            # Are the percentile functions observed?
            j += 1
            if all(isnan, MV[i][start+1:start+grid_pcf, :])
                Σ_ids[j] = NaN
            end
            start += grid_pcf
        end
        j += 1
    end
end


function find_observed_objects!(estimator::KernelEstimator, Σ_ids, MV, cop_id, measures)
    @unpack grid_pcf = estimator
    j = 1
    for i in eachindex(MV)
        # Is the copula observed?
        if all(isnan, MV[i][1:cop_id, :])
            Σ_ids[j] = NaN
        end

        start = cop_id
        for k in eachindex(measures)
            # Are the percentile functions observed?
            j += 1
            if all(isnan, MV[i][start+1:start+grid_pcf, :])
                Σ_ids[j] = NaN
            end
            start += grid_pcf
        end
        j += 1
    end
end

# -----------------------------------------------------------------------------
# NEW FUNCTION: Inverse of u_to_Lcorr 
# -----------------------------------------------------------------------------
function Lcorr_to_u(L::AbstractMatrix)
    n = size(L, 1)
    u_vec = Float64[] # The unconstrained parameter vector u

    k = 1
    for i in 2:n
        c = 1.0 # Initialize cumulative residual variance (c)
        for j in 1:(i-1)
            # 1. Get the partial correlation z = L[i, j] / c
            # This is the z = ζ[k] from the u_to_Lcorr function
            z = L[i, j] / c

            # Check for numerical stability/bounds
            if abs(z) > 1.0 + 1e-10
                @warn "Partial correlation |z| > 1 found. Matrix L is likely invalid."
                z = clamp(z, -1.0, 1.0)
            end

            # 2. Get the unconstrained parameter u = atanh(z)
            u = atanh(z)
            push!(u_vec, u)

            # 3. Update the residual variance 'c'
            c *= sqrt(1.0 - z^2)
            k += 1
        end
        # The diagonal L[i, i] is the final 'c', which is already used implicitly
        # in the construction, so we do not need to check L[i, i] = c, but it should hold.
    end
    return u_vec
end

function u_to_Lcorr(u::AbstractVector, n::Int)
    expected = n * (n - 1) ÷ 2
    @assert length(u) == expected "length(u) must be n*(n-1)/2 = $expected"

    ζ = tanh.(u)  # map to (-1,1)
    L = Matrix{Float64}(I, n, n)
    k = 1
    for i in 2:n
        c = 1.0
        for j in 1:(i-1)
            z = ζ[k]
            @inbounds L[i, j] = z * c
            c *= sqrt(1 - z * z)
            k += 1
        end
        @inbounds L[i, i] = c                # positive diag; ensures diag(R)=1
    end
    return Cholesky(L, 'L', 0)               # R = Matrix(F) gives the correlation matrix
end


function define_parameter_space(model_elements, model_options, prior_objects)
    """Reshapes the parameter vector into 3 matrices and generates the prior, P. Then retains parameters we wish to estimate."""

    @unpack param_vector, priors = prior_objects
    @unpack estimator, case, estimation_object, measures, pre_multiply = model_options
    @unpack n_less_than_one, MV, agg_count = model_elements
    @unpack grid_cop, grid_pcf = estimator

    aggregate_rep = hasproperty(model_options, :aggregate_rep) ? model_options.aggregate_rep : :as_states

    # Set elements in Σ to NaN if the measure is never observed
    dim = length(measures)
    cop_part, imm_part = retrieve_cop_and_imm_part(estimator, dim)
    cop_id = cop_part - imm_part

    if aggregate_rep == :as_inputs
        # ─── :as_inputs ─── 6-block layout: [A, B, Ω_var, Ω_corr, Σ, H] ────
        A = param_vector[1]
        B = param_vector[2]
        Ω_var = param_vector[3]          # factor-only
        Ω_corr_vec = param_vector[4]     # factor-only correlations
        Σ = diag(param_vector[5])        # no aggregate measurement errors
        H = param_vector[6]             # 4 hyperparameters

        Σ_ids = Float64.(collect(1:length(Σ)))
        find_observed_objects!(estimator, Σ_ids, MV, cop_id, measures)

        cond = findall(!isnan, Σ_ids)
        short_Σ = Σ[cond]

        new_param_sizes = Vector(undef, 6)
        for (m, mat) in enumerate([A, B, Ω_var, Ω_corr_vec, short_Σ, H])
            new_param_sizes[m] = size(mat)
        end

        # priors layout: [1]=MvNormal(AB), [2]=Normal(ϕ_F), [3]=LKJ, [4...]=measurement, [end-3:end]=hyperpriors
        n_hyper = length(H)
        priors = [priors[1], priors[2], priors[3], priors[3 .+ cond]..., priors[end-n_hyper+1:end]...]

        return [A[:]; B[:]; Ω_var; Ω_corr_vec; short_Σ[:]; H], new_param_sizes, Σ_ids, priors

    else
        # ─── :as_states (default) ─── 8-block layout: [A, B, C, D, Ω_var, Ω_corr, Σ, H] ─
        A = param_vector[1]
        B = param_vector[2]
        C = param_vector[3]
        D = param_vector[4] # already a vector
        Ω_var = param_vector[5] # a vector
        Ω_corr_vec = param_vector[6] # vector of non-diagonal elements of actual Ω_corr
        Σ = diag(param_vector[7])
        H = param_vector[8] # hyperparameters

        Σ_ids = Float64.(collect(1:length(Σ)))
        find_observed_objects!(estimator, Σ_ids, MV, cop_id, measures)

        # Subset to observed objects
        cond = findall(!isnan, Σ_ids)
        short_Σ = Σ[cond]

        new_param_sizes = Vector(undef, length(param_vector))
        for (m, mat) in enumerate([A, B, C, D, Ω_var, Ω_corr_vec, short_Σ, H])
            new_param_sizes[m] = size(mat)
        end

        # Remove elements from the vector of priors that we don't want to estimate
        # priors = [priors[1], priors[2], priors[3], priors[3 .+ cond]...]
        priors = [priors[1], priors[2], priors[3], priors[4], priors[4 .+ cond]..., priors[end-5:end]...] # keep hyperpriors

        # return [A[:]; B[:]; C[:]; D[:]; Ω[:]; short_Σ[:]], new_param_sizes, Σ_ids, priors # [lb, ub]
        return [A[:]; B[:]; C[:]; D[:]; Ω_var; Ω_corr_vec; short_Σ[:]; H], new_param_sizes, Σ_ids, priors # [lb, ub]
    end
end


# function reconstruct_A(A_vec, factor_count, n_less_than_one)
#     # Construct A from vector to Matrix 
#     new_A = zeros(eltype(A_vec), factor_count, factor_count)
#     n_imp = factor_count - n_less_than_one

#     # Fill  
#     new_A[1:n_imp, 1:n_imp] = reshape(A_vec[1:n_imp*n_imp], (n_imp, n_imp))
#     new_A[n_imp+1:factor_count, n_imp+1:factor_count] = diagm(A_vec[n_imp+1:factor_count])

#     return new_A
# end

#TODO: so here, A comes in with zeros bottomed out (except for the diagonal)
# function matrisize(param_vector, param_sizes, case, estimation_object, factor_count, n_less_than_one) 
#     # For gregors implementation, it requires an nx1 matrix vs. a vec 
#     param_vector = reshape(param_vector, (length(param_vector),)) # TODO: current issue is that my function has to accept a matrix. which is fine but the reshape converts it to a float64 i bleieve 

#     A  = zeros(eltype(param_vector), factor_count, factor_count) 
#     B  = zeros(eltype(param_vector), param_sizes[2][1], param_sizes[2][2]) 
#     Ω  = zeros(eltype(param_vector), param_sizes[3][1], param_sizes[3][1]) 
#     Σ  = zeros(eltype(param_vector), param_sizes[4][1], param_sizes[4][1]) 

#     # Change parameter space based on :case 
#     if case == "diag" 
#         # Convert to matrices  
#         A  .= diagm(param_vector[1:param_sizes[1][1]])  # TODO: works since A has 1 lag
#         B  .= reshape(param_vector[1+param_sizes[1][1]:param_sizes[1][1] + param_sizes[2][1] * param_sizes[2][2]], (param_sizes[2][1], param_sizes[2][2]))
#         Ω  .= diagm(param_vector[1+param_sizes[1][1] + length(B):param_sizes[1][1] + length(B) + param_sizes[3][1]])
#         Σ  .= diagm(param_vector[1+param_sizes[1][1] + length(B) + param_sizes[3][1]:length(param_vector)])

#     elseif case == "A non-diag" && estimation_object == "copulas and percentile functions"
#         # Convert to matrices  
#         A .= reconstruct_A(param_vector[1:param_sizes[1][1]], factor_count, n_less_than_one)
#         # A  .= reshape(param_vector[1:param_sizes[1][1] * param_sizes[1][2]], (param_sizes[1][1], param_sizes[1][2]))  
#         B  .= reshape(param_vector[1+param_sizes[1][1]:param_sizes[1][1] + param_sizes[2][1] * param_sizes[2][2]], (param_sizes[2][1], param_sizes[2][2]))
#         Ω  .= diagm(param_vector[1+param_sizes[1][1]+length(B):param_sizes[1][1]+length(B)+param_sizes[3][1]])
#         Σ  .= diagm(param_vector[1+param_sizes[1][1]+length(B)+param_sizes[3][1]:length(param_vector)])

#     elseif case == "A non-diag" && estimation_object == "levels and percentile functions"
#         # Convert to matrices  
#         A  .= reshape(param_vector[1:param_sizes[1][1] * param_sizes[1][2]], (param_sizes[1][1], param_sizes[1][2]))  
#         B  .= reshape(param_vector[1+length(A):length(A) + param_sizes[2][1] * param_sizes[2][2]], (param_sizes[2][1], param_sizes[2][2]))
#         Ω  .= diagm(param_vector[1+length(A)+length(B):length(A)+length(B)+param_sizes[3][1]])
#         Σ  .= diagm(param_vector[1+length(A)+length(B)+param_sizes[3][1]:length(param_vector)])

#     elseif case == "A, Σ non-diag"
#         A .= reconstruct_A(param_vector[1:param_sizes[1][1]], factor_count, n_less_than_one)
#         # A  .= reshape(param_vector[1:param_sizes[1][1] * param_sizes[1][2]], (param_sizes[1][1], param_sizes[1][2]))
#         B  .= reshape(param_vector[1+param_sizes[1][1]:param_sizes[1][1] + param_sizes[2][1] * param_sizes[2][2]], (param_sizes[2][1], param_sizes[2][2]))
#         Ω  .= diagm(param_vector[1+param_sizes[1][1]+length(B):param_sizes[1][1]+length(B)+param_sizes[3][1]])
#         Σ  .= 1.0 * reshape(param_vector[1+param_sizes[1][1]+length(B)+param_sizes[3][1]:length(param_vector)], (param_sizes[4][1], param_sizes[4][2]))
#         #TODO: implement the case where Σ has cross correlants 
#     end 

#     return A, B, Hermitian(Ω), Hermitian(Σ)  # they need to be positive semi-def. 
# end

function vec_to_lower(v::AbstractVector, n::Int)
    L = zeros(n, n)
    k = 1
    for j in 1:n
        for i in (j+1):n   # rows below the diagonal
            L[i, j] = v[k]
            k += 1
        end
    end
    for i in 1:n
        L[i, i] = 1.0
    end
    return L
end

# function lower_offdiag(L::AbstractMatrix)
#     n = size(L, 1)
#     L_vec = zeros(n * (n - 1) ÷ 2)

#     k = 1
#     for j in 1:n
#         for i in (j+1):n   # rows below the diagonal
#             L_vec[k] = L[i, j]
#             k += 1
#         end
#     end
#     return L_vec
# end

function lower_offdiag(L::AbstractMatrix)
    n = size(L, 1)
    v = Vector{eltype(L)}(undef, n * (n - 1) ÷ 2)
    k = 1
    for j in 1:n-1, i in (j+1):n
        @inbounds v[k] = L[i, j]
        k += 1
    end
    return v
end

function lower_to_vec(L::AbstractMatrix)
    n = size(L, 1)
    L_vec = zeros(n * (n + 1) ÷ 2)

    k = 1
    for j in 1:n
        for i in j:n   # rows below the diagonal
            L_vec[k] = L[i, j]
            k += 1
        end
    end
    return L_vec
end



function matrisize(param_vector, param_sizes; aggregate_rep=:as_states)
    # For gregors implementation, it requires an nx1 matrix vs. a vec
    param_vector = reshape(param_vector, length(param_vector))

    if aggregate_rep == :as_inputs
        # ─── :as_inputs ─── 6-block layout: [A, B, Ω_var, Ω_corr, Σ] (no C, D, no H) ────
        A = zeros(param_sizes[1][1], param_sizes[1][1])
        B = zeros(param_sizes[2][1], param_sizes[2][2])
        Ω_var = zeros(param_sizes[3][1], param_sizes[3][1])
        Σ = zeros(param_sizes[5][1], param_sizes[5][1])

        l_A = length(A)
        l_B = prod(param_sizes[2])
        l_AB = l_A + l_B
        l_Ω_var = param_sizes[3][1]
        l_Ω_corr = param_sizes[4][1]  # number of free elements in lower triangular
        l_Ω = l_Ω_var + l_Ω_corr

        A .= reshape(view(param_vector, 1:l_A), param_sizes[1])
        B .= reshape(view(param_vector, 1+l_A:l_AB), param_sizes[2])
        Ω_var .= diagm(reshape(view(param_vector, 1+l_AB:l_AB+l_Ω_var), l_Ω_var))
        Ω_corr = u_to_Lcorr(view(param_vector, 1+l_AB+l_Ω_var:l_AB+l_Ω_var+l_Ω_corr), param_sizes[3][1])
        Σ .= diagm(view(param_vector, l_AB+l_Ω+1:length(param_vector)))

        return A, B, nothing, nothing, Ω_var, Ω_corr, Hermitian(Σ)

    else
        # ─── :as_states (default) ─── 8-block layout ────────────────────────
        # n_aggs = length(param_sizes[3])  # number of aggregates

        A = zeros(param_sizes[1][1], param_sizes[1][1])
        B = zeros(param_sizes[2][1], param_sizes[2][2])
        C = zeros(param_sizes[3][1], param_sizes[3][2])
        D = zeros(param_sizes[4][1], param_sizes[4][1])
        # Ω = zeros(eltype(param_vector), param_sizes[5][1], param_sizes[5][1])
        Ω_var = zeros(param_sizes[5][1], param_sizes[5][1])
        # Ω_corr = zeros(param_sizes[5][1], param_sizes[5][1])
        Σ = zeros(param_sizes[7][1], param_sizes[7][1])

        l_A = length(A)
        l_B = length(B)
        l_AB = l_A + l_B
        l_C = length(C)
        l_D = param_sizes[4][1]
        l_ABC = l_AB + l_C
        l_ABCD = l_ABC + l_D
        l_Ω_var = param_sizes[5][1]
        l_Ω_corr = param_sizes[6][1]  # number of free elements in lower triangular matrix
        l_Ω = l_Ω_var + l_Ω_corr

        # Change parameter space based on :case
        A .= reshape(view(param_vector, 1:l_A), param_sizes[1])
        B .= reshape(view(param_vector, 1+l_A:l_AB), param_sizes[2])
        C .= reshape(view(param_vector, 1+l_AB:l_ABC), param_sizes[3])
        D .= diagm(view(param_vector, 1+l_ABC:l_ABCD))
        # Ω .= diagm(view(param_vector, 1+l_ABC+param_sizes[4][1]+1:1+l_ABC+param_sizes[4][1]+param_sizes[5][1]))
        Ω_var .= diagm(reshape(view(param_vector, 1+l_ABCD:l_ABCD+l_Ω_var), l_Ω_var))
        Ω_corr = u_to_Lcorr(view(param_vector, 1+l_ABCD+l_Ω_var:l_ABCD+l_Ω_var+l_Ω_corr), param_sizes[5][1]) # maps also "u" to (-1,1)

        # Ω .= Ω_var * Ω_corr * Ω_var

        # Creating Σ, adding zeros for aggregates (perfectly measured)
        # Σ_vec = vcat(view(param_vector, length(A)+length(B)+param_sizes[3][1]+param_sizes[4][1]+1:length(param_vector)), zeros(n_aggs))
        Σ .= diagm(view(param_vector, l_ABCD+l_Ω+1:length(param_vector)))

        # return A, B, C, D, Hermitian(Ω), Hermitian(Σ)  # they need to be positive semi-def.
        return A, B, C, D, Ω_var, Ω_corr, Hermitian(Σ)  # they need to be positive semi-def.
    end
end



function variance_check!(Ω, Σ)
    """For now, this works in the diagonal case. Essentially checking that the params are indeed positive."""
    diag_Ω = diag(Ω)
    diag_Σ = diag(Σ)
    if any(diag_Ω .< 0)
        # println("True, Ω neg.")
        diag_Ω[diag_Ω.<0] = (diag_Ω[diag_Ω.<0]) .^ 2
        Ω[diagind(Ω)] .= diag_Ω
    end
    if any(diag_Σ .< 0)
        # println("True, Σ neg.")
        diag_Σ[diag_Σ.<0] = (diag_Σ[diag_Σ.<0]) .^ 2
        Σ[diagind(Σ)] .= diag_Σ
    end
end

function get_priors(model_elements, model_options, hyper_params)
    """Generates the prior distributions for the parameters."""

    # Define hyperpriors
    @unpack agg_count, factor_count, pcs, u, MV, proj = model_elements
    @unpack estimator, measures, constant, lags, prior, measurement_error, estimation_object, case, errors_process, pre_multiply = model_options
    @unpack grid_pcf, grid_cop = estimator
    @unpack hyperparameters = prior

    aggregate_rep = hasproperty(model_options, :aggregate_rep) ? model_options.aggregate_rep : :as_states

    number_of_dfs = length(MV)
    dimension = length(measures)

    if aggregate_rep == :as_inputs
        # ─── :as_inputs ─── 4 hyperparameters: κ_1, κ_5, σ_ϕ, ν ────────────
        hyper_params[1] = inv_logit(hyper_params[1], 0.2, 0.99)   # κ_1
        hyper_params[2] = inv_logit(hyper_params[2], -0.99, 0.99)  # κ_5

        hyperparameters[1] = 0.05
        hyperparameters[2] = hyper_params[1]
        hyperparameters[6] = hyper_params[2]

        prior_mean, _, _, _, _, V_prior = minnesota_prior(hyperparameters, pcs, u, lags, estimator; aggregate_rep=:as_inputs)

        n_objs = (dimension + 1)

        ρ_F = hyperparameters[6]
        ϕ_F = log(exp(1 - ρ_F^2) - 1)

        σ_ϕ = log(exp(hyper_params[3]) + 1)
        corr_Ω_shape = 1 + log(1 + exp(hyper_params[4]))

        priors = [
            MvNormal(prior_mean, V_prior),           # [1] coefficients of A,B
            Normal(ϕ_F, σ_ϕ),                        # [2] prior for factor persistence
            LKJCholesky(factor_count, corr_Ω_shape),  # [3] factor-only correlations
        ]

        ϕₘ = log(exp(1) - 1)
        for _ in 1:number_of_dfs
            for _ in 1:n_objs
                push!(priors, Normal(ϕₘ, 2))
            end
        end

        # No aggregate measurement error priors

        return priors

    else
        # ─── :as_states (default) ─── 6 hyperparameters ─────────────────────

        # Map hyperparameters to within their Bounds
        # 1. First parameter (e.g., ρ_var) constrained to [0.2, 0.99]
        hyper_params[1] = inv_logit(hyper_params[1], 0.2, 0.99)

        # 2. Second parameter (e.g., ρ_pers_state) constrained to [-0.99, 0.99]
        hyper_params[2] = inv_logit(hyper_params[2], -0.99, 0.99)

        # 3. Third parameter (e.g., ρ_pers_agg) constrained to [-0.99, 0.99]
        hyper_params[3] = inv_logit(hyper_params[3], -0.99, 0.99)


        # Replace hyperparameters
        hyperparameters[1] = 0.05
        hyperparameters[2] = hyper_params[1]
        hyperparameters[6] = hyper_params[2]
        hyperparameters[7] = hyper_params[3]

        # Construct prior hyperparameters
        prior_mean, _, _, _, _, V_prior = minnesota_prior(hyperparameters, pcs, u, lags, estimator)  # TODO: assumes aggs only have 1 lag

        # For measurement error
        n_objs = (dimension + 1)

        # ───  Extract the two persistence hyper-parameters  ─────────────────────────
        ρ_F = hyperparameters[6]        # κ₅  = persistence target for the  r-factor block
        ρ_Y = hyperparameters[7]        # κ₆  = persistence target for the  q-aggregate block

        #  Same log-transform you already used for ρ_F, now applied to both blocks
        ϕ_F = log(exp(1 - ρ_F^2) - 1)    # ⇒ Normal prior on ϕ_F  implies Beta prior on ρ_F
        ϕ_Y = log(exp(1 - ρ_Y^2) - 1)    # ⇒ identical mapping for the aggregate persistence

        # ───  Prior list  ──────────────────────────────────────────────────────────
        σ_ϕ = log(exp(hyper_params[4]) + 1) # These are already corrected
        σ_ϕ_Y = log(exp(hyper_params[5]) + 1)

        corr_Ω_shape = 1 + log(1 + exp(hyper_params[6])) # shape parameter of 2, meaning mild prior towards identity matrix

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

        return priors
    end
end

# Apply the inverse logit transformation to get the constrained values (ρ*)
function inv_logit(z, a, b)
    # Scaled Inverse Logit: a + (b - a) * (1 / (1 + exp(-z)))
    return a + (b - a) * (1.0 / (1.0 + exp(-z)))
end



# TODO: given the new methodology, prioreval should come after the time invariant collapse 
function likeli(model_elements, param_vector, param_sizes, hyperpriors, Σ_ids, model_options; smooth=false)
    """Computes the likelihood of the model given the parameters. 

    Returns smoother estimates if 'smooth'.
    """

    @unpack number_of_dfs = model_options
    @unpack y, G, u, agg_count = model_elements

    aggregate_rep = hasproperty(model_options, :aggregate_rep) ? model_options.aggregate_rep : :as_states

    # Set Priors — dynamic hyperparameter count based on param_sizes
    n_hyper = param_sizes[end][1]
    hyper_params = param_vector[end-n_hyper+1:end]
    other_params = param_vector[1:end-n_hyper]

    # Priors change based on the hyper_params --- hyperpriors do not
    local priors
    try
        priors = get_priors(model_elements, model_options, hyper_params)
        # When there is a cholesky error, return a very low likelihood and an alarm, so that the sampler can move on to the next draw.
    catch e
        # @warn "Error in get_priors: $e. Returning low likelihood and alarm."
        return -1.e12, true
    end
    push!(priors, hyperpriors...)

    # Reshape param_vector into matrices
    A, B, C, D, Ω_var, Ω_corr, Σ = matrisize(other_params, param_sizes; aggregate_rep=aggregate_rep)

    # Reparametization
    if aggregate_rep == :as_inputs
        log_P, alarm = prioreval([[A..., B...], Ω_var, Ω_corr, Matrix(Σ), hyper_params], priors, param_sizes; aggregate_rep=:as_inputs)
    else
        log_P, alarm = prioreval([[A..., B..., C..., diag(D)...], Ω_var, Ω_corr, Matrix(Σ), hyper_params], priors, param_sizes)
    end
    Ω_var[diagind(Ω_var)] = log.(exp.(Ω_var[diagind(Ω_var)]) .+ 1)  # softplus transformation

    # Convert Cholesky factor to correlation matrix
    mat_Ω_corr = Matrix(Ω_corr)

    if !all(isapprox(diag(mat_Ω_corr), ones(size(mat_Ω_corr, 1))))
        return -1.e12, true
    end
    Ω = Ω_var * mat_Ω_corr * Ω_var'  # reconstructing Ω

    if alarm
        # println("not in support")
        log_D = -1.e12
        tot_log_like = log_P + log_D
        return tot_log_like, alarm
    end

    if aggregate_rep == :as_inputs
        # ─── :as_inputs ─── Ω is factor-only; aggregates are exogenous inputs ──
        Σ = apply_measurement_criteria(Σ, Σ_ids, model_options, 0)  # agg_count=0 for measurement

        # Transform parameters related to Σ
        cond = Σ .> 10
        Σ[.!cond] .= log.(exp.(Σ[.!cond]) .+ 1)

        Ω_f = Ω  # entire Ω is factor-only

        if smooth
            smoother_output, log_D, alarm = recurse_kalman_filter_smoother_inputs(A, B, Ω_f, Σ, G, y, u)
            tot_log_like = log_P + log_D
            return smoother_output, tot_log_like, alarm
        else
            local log_D, alarm
            try
                log_D, alarm = recurse_kalman_filter_forward_inputs(A, B, Ω_f, Σ, G, y, u)
            catch e
                log_D = -1.e12
                alarm = true
            end
            tot_log_like = log_P + log_D
            return tot_log_like, alarm
        end
    end

    # Expand to comformable Σ if necessary

    # TODO: IF WE STOP USING SEQUENTIAL, THEN WE HAVE TO UNDO THE DIAG AND CORRECT THE BUTTERFLY EFFECT
    Σ = apply_measurement_criteria(Σ, Σ_ids, model_options, agg_count)

    # Transform parameters related to Σ
    cond = Σ .> 10
    Σ[.!cond] .= log.(exp.(Σ[.!cond]) .+ 1) # Basically, for large numbers, it does not matter if we use log(exp(x) + 1) or log(x + 1)

    # Create the separate Ω_f, Ω_y
    Ω_f = Ω[1:size(A, 1), 1:size(A, 1)]
    Ω_y = Ω[size(A, 1)+1:end, size(A, 1)+1:end]

    # Filter-smoother estimates
    if smooth
        # Generate likelihood of data and smoother output 
        smoother_output, log_D, alarm = recurse_kalman_filter_smoother(A, B, C, D, Ω_f, Ω_y, Σ, G, y)
        tot_log_like = log_P + log_D         # Total likelihood 
        return smoother_output, tot_log_like, alarm

    else
        local log_D, alarm
        try
            log_D, alarm = recurse_kalman_filter_forward(A, B, C, D, Ω_f, Ω_y, Σ, G, y)
        catch e
            log_D = -1.e12
            alarm = true
            # log_D, alarm         = recurse_kalman_filter(A,B,Ω,VΣ̂V,Ĝ,ŷ,u,smooth)
        end

        tot_log_like = log_P + log_D         # Total likelihood 
        return tot_log_like, alarm
    end
end


function recurse_kalman_filter_smoother(A, B, C, D, Ω_f, Ω_y, Σ, G, y)
    """Performs kalman filter recursively, performing a filtering step and updating step thereafter.
    Assumption is: measurements are uncorrelated.

    In terms of timing, everything has been checked so that measurements and aggregates 
    are aligned and that the first estimate is tmin[]
    """
    T = size(y, 2)
    r, q = size(B)
    nₛ = 4r + q
    Tval = eltype(A)

    # --------- constant matrices ------------------------------------
    AI = Matrix{Tval}(I, r, r)

    L = zeros(Tval, nₛ, nₛ)
    @views begin
        L[1:r, 1:r] .= A
        L[1:r, 4r+1:4r+q] .= B
        L[r+1:2r, 1:r] .= AI
        L[2r+1:3r, r+1:2r] .= AI
        L[3r+1:4r, 2r+1:3r] .= AI
        L[4r+1:end, 1:r] .= C
        L[4r+1:end, 4r+1:end] .= D
    end
    Lᵀ = transpose(L)                       # (B)

    bigΩ = zeros(Tval, nₛ, nₛ)
    @views begin
        bigΩ[1:r, 1:r] .= Ω_f
        bigΩ[4r+1:end, 4r+1:end] .= Ω_y
    end

    # --------- initial prior ----------------------------------------
    ixF, ixU = 1:r, 4r+1:4r+q
    Lcore = [A B; C D]
    Qcore = zeros(Tval, r + q, r + q)
    Qcore[1:r, 1:r] .= Ω_f
    Qcore[r+1:end, r+1:end] .= Ω_y
    Pcore = lyapd(Lcore, Qcore)

    P_F = zeros(Tval, nₛ, nₛ)
    P_F[ixF, ixF] .= Pcore[1:r, 1:r]
    P_F[ixU, ixU] .= Pcore[r+1:end, r+1:end]
    P_F[ixF, ixU] .= Pcore[1:r, r+1:end]
    P_F[ixU, ixF] .= Pcore[r+1:end, 1:r]
    P_F[diagind(P_F)[r+1:4r]] .= 1e6

    # Filtered containers  
    x_filtered = zeros(Tval, nₛ, T)
    sigma_filtered = [zeros(Tval, nₛ, nₛ) for _ in 1:T]

    # Updated containers 
    x_updated = zeros(Tval, nₛ, T)
    sigma_updated = [zeros(Tval, nₛ, nₛ) for _ in 1:T]

    # Containers for factors and controls
    x̂ = zeros(Tval, nₛ)
    x_filtered_factors = zeros(Tval, nₛ)
    # x_filtered_controls = zeros(eltype(A), v)

    # Containers for sigma
    sigma_filtered_r = Matrix{Tval}(undef, nₛ, nₛ)
    sigma_filtered_l = Matrix{Tval}(undef, nₛ, nₛ)

    # Model Likelihood 
    logL = 0.0

    # quarterly at the moment 
    for t = 1:T
        # Filtering equations (but better said: prediction equations): uses prior only and C = controls, if any 
        mul!(x_filtered_factors, L, x̂)
        x_filtered[:, t] .= x_filtered_factors # .+ x_filtered_controls


        mul!(sigma_filtered_r, P_F, Lᵀ)
        mul!(sigma_filtered_l, L, sigma_filtered_r)
        sigma_filtered[t] .= sigma_filtered_l .+ bigΩ

        # Update step (also called the forward pass)   
        logL += @views kalman_update!(x_updated, sigma_updated, y[:, t], t, G[t], x_filtered[:, t], sigma_filtered[t], Σ)

        # Update priors 
        copyto!(P_F, sigma_updated[t])
        copy!(x̂, x_updated[:, t])
    end

    alarm = false

    if logL <= -1.e12
        alarm = true
    end

    # Gibbs requires a likelihood over a multivariate dist.
    smoother_results = run_smoother(L, x_filtered, sigma_filtered, x_updated, sigma_updated, T)
    return smoother_results, logL, alarm
end


struct KalmanWorkspace{T}
    g::Vector{T}        # measurement row (n)
    sig_r::Vector{T}        # P_F * g          (n)
    # sig_l::Vector{T}        # g' * P_F        (n)
    K::Vector{T}        # Kalman gain      (n)
    # Kg::Matrix{T}        # K * g'           (n×n)
    # id::Diagonal{T,Vector{T}}               # Identity (n×n)
    # id_Kg::Matrix{T}        # I − K g'         (n×n)
    # id_sig::Matrix{T}        # (I−Kg)P_F        (n×n)
    # eta_buf::Vector{T}        # 1‑element buffer
    # K_eta::Vector{T}        # K * eta         (n)
end

function KalmanWorkspace(n::Int, ::Type{T}=Float64) where {T}
    KalmanWorkspace{T}(
        zeros(T, n),
        zeros(T, n),
        # zeros(T, 1),
        zeros(T, n),
        # zeros(T, n, n),
        # Diagonal{T}(ones(T, n)),
        # zeros(T, n, n),
        # zeros(T, n, n),
        # zeros(T, 1),
        # zeros(T, n)
    )
end

# function recurse_kalman_filter_forward(A, B, C, D, Ω_f, Ω_y, Σ, G, y)
#     # Dimensions
#     r, q = size(B)  # number of factors and controls
#     AI = Matrix{eltype(A)}(I, size(A, 1), size(A, 1))  # Identity matrix

#     # Build transition matrix L
#     L = zeros(eltype(A), 4r + q, 4r + q)
#     @views begin
#         L[1:r, 1:r] .= A
#         L[1:r, 4r+1:4r+q] .= B
#         L[r+1:2r, 1:r] .= AI
#         L[2r+1:3r, r+1:2r] .= AI
#         L[3r+1:4r, 2r+1:3r] .= AI
#         L[4r+1:end, 1:r] .= C
#         L[4r+1:end, 4r+1:end] .= D
#     end

#     # Build covariance matrix big_Ω
#     big_Ω = zeros(eltype(A), 4r + q, 4r + q)
#     @views begin
#         big_Ω[1:r, 1:r] .= Ω_f
#         big_Ω[4r+1:end, 4r+1:end] .= Ω_y
#     end

#     nₛ = size(big_Ω, 1)  # full state size
#     T = size(y, 2)       # number of time steps
#     x̂ = zeros(eltype(A), nₛ)

#     # Initial covariance P_F from Lyapunov solution
#     ixF = 1:r
#     ixU = 4r+1:4r+q

#     Lcore = [A B; C D]
#     Qcore = zeros(eltype(A), r + q, r + q)
#     Qcore[1:r, 1:r] .= Ω_f
#     Qcore[r+1:end, r+1:end] .= Ω_y
#     Pcore = lyapd(Lcore, Qcore)

#     P_F = zeros(eltype(A), 4r + q, 4r + q)
#     P_F[ixF, ixF] .= Pcore[1:r, 1:r]
#     P_F[ixU, ixU] .= Pcore[r+1:end, r+1:end]
#     P_F[ixF, ixU] .= Pcore[1:r, r+1:end]
#     P_F[ixU, ixF] .= Pcore[r+1:end, 1:r]
#     P_F[diagind(P_F)[r+1:4r]] .= 1e6  # placeholder variance

#     # Preallocate
#     x_filtered_factors = zeros(eltype(A), nₛ)
#     sigma_filtered_r = Matrix{eltype(A)}(undef, nₛ, nₛ)
#     sigma_filtered_l = Matrix{eltype(A)}(undef, nₛ, nₛ)
#     sigma_pred = zeros(eltype(A), nₛ, nₛ)
#     # x_updated = zeros(eltype(A), nₛ, T)
#     # sigma_updated = [zeros(eltype(A), nₛ, nₛ) for _ in 1:T]
#     x_updated = zeros(eltype(A), nₛ)
#     sigma_updated = zeros(eltype(A), nₛ, nₛ)
#     logL = 0.0
#     ws = KalmanWorkspace(nₛ, eltype(A))


#     for t = 1:T
#         mul!(x_filtered_factors, L, x̂)
#         mul!(sigma_filtered_r, P_F, L')
#         mul!(sigma_filtered_l, L, sigma_filtered_r)
#         sigma_pred .= sigma_filtered_l .+ big_Ω

#         logL += @views kalman_update_fast!(x_updated, sigma_updated, @views(y[:, t]), t, G[t], x_filtered_factors, sigma_pred, Σ, ws)

#         # copyto!(P_F, sigma_updated[t])
#         # copy!(x̂, x_updated[:, t])
#         copy!(P_F, sigma_updated)
#         copy!(x̂, x_updated)
#     end

#     alarm = logL <= -1e12
#     return logL, alarm
# end

function recurse_kalman_filter_forward(A, B, C, D, Ω_f, Ω_y, Σ, G, y)
    T = size(y, 2)
    r, q = size(B)
    nₛ = 4r + q
    Tval = eltype(A)

    # --------- constant matrices ------------------------------------
    AI = Matrix{Tval}(I, r, r)

    L = zeros(Tval, nₛ, nₛ)
    @views begin
        L[1:r, 1:r] .= A
        L[1:r, 4r+1:4r+q] .= B
        L[r+1:2r, 1:r] .= AI
        L[2r+1:3r, r+1:2r] .= AI
        L[3r+1:4r, 2r+1:3r] .= AI
        L[4r+1:end, 1:r] .= C
        L[4r+1:end, 4r+1:end] .= D
    end
    Lᵀ = transpose(L)                       # (B)

    bigΩ = zeros(Tval, nₛ, nₛ)
    @views begin
        bigΩ[1:r, 1:r] .= Ω_f
        bigΩ[4r+1:end, 4r+1:end] .= Ω_y
    end

    # --------- initial prior ----------------------------------------
    ixF, ixU = 1:r, 4r+1:4r+q
    Lcore = [A B; C D]
    Qcore = zeros(Tval, r + q, r + q)
    Qcore[1:r, 1:r] .= Ω_f
    Qcore[r+1:end, r+1:end] .= Ω_y
    Pcore = lyapd(Lcore, Qcore)

    P_F = zeros(Tval, nₛ, nₛ)
    P_F[ixF, ixF] .= Pcore[1:r, 1:r]
    P_F[ixU, ixU] .= Pcore[r+1:end, r+1:end]
    P_F[ixF, ixU] .= Pcore[1:r, r+1:end]
    P_F[ixU, ixF] .= Pcore[r+1:end, 1:r]
    P_F[diagind(P_F)[r+1:4r]] .= 1e6

    # --------- work buffers -----------------------------------------
    x̂ = zeros(Tval, nₛ)
    x_pred = similar(x̂)
    σ_r = similar(P_F)
    σ_pred = similar(P_F)           # will hold P_{t|t-1} in-place
    x_tmp = similar(x̂)
    σ_tmp = similar(P_F)
    ws = KalmanWorkspace(nₛ, Tval)

    # --------- main loop --------------------------------------------
    logL = zero(Tval)

    @inbounds for t = 1:T
        mul!(x_pred, L, x̂)
        mul!(σ_r, P_F, Lᵀ)
        mul!(σ_pred, L, σ_r)           # σ_pred is P_{t|t-1}

        BLAS.axpy!(1.0, bigΩ, σ_pred)    # σ_pred += bigΩ ; (E)

        # ---- update (C): pass raw pointer to column t without view alloc
        y_col = @view y[:, t]            # a single SubArray per step is OK
        G_t = G[t]                     # assuming G is a Vector{Matrix}


        logL += kalman_update_fast!(
            x_tmp, σ_tmp, y_col, G_t,
            x_pred, σ_pred, Σ, ws)

        # ---- posterior becomes prior (D)
        copy!(x̂, x_tmp)
        copy!(P_F, σ_tmp)
    end

    alarm = logL <= -1e12
    return logL, alarm
end


# ═══════════════════════════════════════════════════════════════════════════════
#  :as_inputs Kalman filters — aggregates as exogenous inputs, state = factors only
# ═══════════════════════════════════════════════════════════════════════════════

function recurse_kalman_filter_smoother_inputs(A, B_input, Ω_f, Σ, G, y, u)
    """Kalman filter-smoother for :as_inputs mode.
    State = [F_t; F_{t-1}; F_{t-2}; F_{t-3}] (4r × 1), no aggregate states.
    Prediction: x_pred = L * x̂ + E * u_t  where u_t are observed aggregates.
    """
    T = size(y, 2)
    r = size(A, 1)
    q = size(B_input, 2)
    nₛ = 4r
    Tval = eltype(A)

    # --------- constant matrices ------------------------------------
    AI = Matrix{Tval}(I, r, r)

    L = zeros(Tval, nₛ, nₛ)
    @views begin
        L[1:r, 1:r] .= A
        L[r+1:2r, 1:r] .= AI
        L[2r+1:3r, r+1:2r] .= AI
        L[3r+1:4r, 2r+1:3r] .= AI
    end
    Lᵀ = transpose(L)

    # Exogenous input matrix E (4r × q): only factors load on aggregates
    E = zeros(Tval, nₛ, q)
    E[1:r, :] .= B_input

    bigΩ = zeros(Tval, nₛ, nₛ)
    bigΩ[1:r, 1:r] .= Ω_f

    # --------- initial prior ----------------------------------------
    Pcore = lyapd(A, Ω_f)
    P_F = zeros(Tval, nₛ, nₛ)
    P_F[1:r, 1:r] .= Pcore
    P_F[diagind(P_F)[r+1:4r]] .= 1e6

    # Filtered containers
    x_filtered = zeros(Tval, nₛ, T)
    sigma_filtered = [zeros(Tval, nₛ, nₛ) for _ in 1:T]

    # Updated containers
    x_updated = zeros(Tval, nₛ, T)
    sigma_updated = [zeros(Tval, nₛ, nₛ) for _ in 1:T]

    # Containers for factors
    x̂ = zeros(Tval, nₛ)
    x_filtered_factors = zeros(Tval, nₛ)

    # Containers for sigma
    sigma_filtered_r = Matrix{Tval}(undef, nₛ, nₛ)
    sigma_filtered_l = Matrix{Tval}(undef, nₛ, nₛ)

    # Exogenous input buffer
    Eu = zeros(Tval, nₛ)

    # Model Likelihood
    logL = 0.0

    for t = 1:T
        # Prediction step: x_pred = L * x̂ + E * u_t
        mul!(x_filtered_factors, L, x̂)
        mul!(Eu, E, view(u, :, t))
        x_filtered_factors .+= Eu
        x_filtered[:, t] .= x_filtered_factors

        mul!(sigma_filtered_r, P_F, Lᵀ)
        mul!(sigma_filtered_l, L, sigma_filtered_r)
        sigma_filtered[t] .= sigma_filtered_l .+ bigΩ

        # Update step
        logL += @views kalman_update!(x_updated, sigma_updated, y[:, t], t, G[t], x_filtered[:, t], sigma_filtered[t], Σ)

        # Update priors
        copyto!(P_F, sigma_updated[t])
        copy!(x̂, x_updated[:, t])
    end

    alarm = logL <= -1.e12

    smoother_results = run_smoother(L, x_filtered, sigma_filtered, x_updated, sigma_updated, T; E=E, u=u)
    return smoother_results, logL, alarm
end


function recurse_kalman_filter_forward_inputs(A, B_input, Ω_f, Σ, G, y, u)
    """Forward-only Kalman filter for :as_inputs mode (no smoother).
    State = [F_t; F_{t-1}; F_{t-2}; F_{t-3}] (4r × 1).
    """
    T = size(y, 2)
    r = size(A, 1)
    q = size(B_input, 2)
    nₛ = 4r
    Tval = eltype(A)

    # --------- constant matrices ------------------------------------
    AI = Matrix{Tval}(I, r, r)

    L = zeros(Tval, nₛ, nₛ)
    @views begin
        L[1:r, 1:r] .= A
        L[r+1:2r, 1:r] .= AI
        L[2r+1:3r, r+1:2r] .= AI
        L[3r+1:4r, 2r+1:3r] .= AI
    end
    Lᵀ = transpose(L)

    E = zeros(Tval, nₛ, q)
    E[1:r, :] .= B_input

    bigΩ = zeros(Tval, nₛ, nₛ)
    bigΩ[1:r, 1:r] .= Ω_f

    # --------- initial prior ----------------------------------------
    Pcore = lyapd(A, Ω_f)
    P_F = zeros(Tval, nₛ, nₛ)
    P_F[1:r, 1:r] .= Pcore
    P_F[diagind(P_F)[r+1:4r]] .= 1e6

    # --------- work buffers -----------------------------------------
    x̂ = zeros(Tval, nₛ)
    x_pred = similar(x̂)
    σ_r = similar(P_F)
    σ_pred = similar(P_F)
    x_tmp = similar(x̂)
    σ_tmp = similar(P_F)
    Eu = zeros(Tval, nₛ)
    ws = KalmanWorkspace(nₛ, Tval)

    # --------- main loop --------------------------------------------
    logL = zero(Tval)

    @inbounds for t = 1:T
        mul!(x_pred, L, x̂)
        mul!(Eu, E, view(u, :, t))
        x_pred .+= Eu

        mul!(σ_r, P_F, Lᵀ)
        mul!(σ_pred, L, σ_r)
        BLAS.axpy!(1.0, bigΩ, σ_pred)

        y_col = @view y[:, t]
        G_t = G[t]

        logL += kalman_update_fast!(
            x_tmp, σ_tmp, y_col, G_t,
            x_pred, σ_pred, Σ, ws)

        copy!(x̂, x_tmp)
        copy!(P_F, σ_tmp)
    end

    alarm = logL <= -1e12
    return logL, alarm
end


function save_estimates!(P_F, x̂, sigma_updated, x_updated)
    P_F .= sigma_updated[:, :, t]
    x̂ .= x_updated[:, t]
    return P_F, x̂
end

function save_estimates(sigma_updated, x_updated)
    P_F = Zygote.Buffer(sigma_updated) # Buffer supports syntax like similar
    x̂ = Zygote.Buffer(x_updated)
    save_estimates!(P_F, x̂, sigma_updated, x_updated)
    return copy(P_F), copy(x̂) # this step makes the Buffer immutable (w/o actually copying)
end



function get_n_of_copulas(estimator, dimension)
    if typeof(estimator) <: SeriesEstimator
        n_of_cops = 0

        for k in 2:dimension
            n_of_cops += binomial(dimension, k)
        end

        return n_of_cops
    else
        return 1
    end
end



# function apply_measurement_criteria(Σ, Σ_ids, model_options, agg_count)
#     """Takes the (grid * number_of_dfs) parameters and extends to fit the actual Σ, accounting for the immutable portions."""

#     @unpack pre_multiply, measurement_error, estimation_object, number_of_dfs, measures, estimator = model_options
#     @unpack grid_pcf, grid_cop = estimator

#     dimension = length(measures)

#     cop_part, imm_part = retrieve_cop_and_imm_part(estimator, dimension)
#     variances = []

#     cop_part, imm_part = retrieve_cop_and_imm_part(estimator, dimension)
#     gp_cop = cop_part - imm_part

#     x = 1 # Σ_ids are of length n_objs * number_of_dfs
#     ii = 1 # diag(Σ) is of length number of objects observed across all dfs

#     for _ in 1:number_of_dfs
#         if isnan(Σ_ids[x])
#             # Append NaNs to the variances, as many NaNs as there are coefficients for this respective copula
#             variances = vcat(variances, repeat([NaN], gp_cop))
#         else
#             # Append the variance parameter assigned to this observed object 
#             variances = vcat(variances, repeat([Σ[ii, ii]], gp_cop))
#             ii += 1
#         end
#         x += 1

#         # for each pcf ....
#         for _ in 1:dimension
#             if isnan(Σ_ids[x])
#                 variances = vcat(variances, repeat([NaN], grid_pcf))
#             else
#                 variances = vcat(variances, repeat([Σ[ii, ii]], grid_pcf))
#                 ii += 1
#             end
#             x += 1
#         end
#     end

#     # Add variances for the aggregates
#     agg_error = [Σ[j, j] for j in (ii):(ii+agg_count-1)]
#     variances = vcat(variances, agg_error)

#     return variances
# end

function apply_measurement_criteria(
    Σ::AbstractMatrix, Σ_ids,
    model_options, agg_count
)
    @unpack number_of_dfs, measures,
    estimator = model_options
    @unpack grid_pcf, grid_cop = estimator

    dim = length(measures)
    cop_part, imm = retrieve_cop_and_imm_part(estimator, dim)
    gp_cop = cop_part - imm           # free copula coeffs / obj
    block_len = gp_cop + dim * grid_pcf  # coeffs per data frame

    total_len = number_of_dfs * block_len + agg_count
    v = Vector{eltype(Σ)}(undef, total_len)

    pos = 1               # write pointer in v
    id_ptr = 1               # pointer in Σ_ids
    Σptr = 1               # pointer on Σ diagonal

    diagΣ = diag(Σ)            # view of diagonal (no alloc)

    @inbounds for _df in 1:number_of_dfs
        # ----- copula part -------------------------------------------
        if isnan(Σ_ids[id_ptr])
            v[pos:pos+gp_cop-1] .= NaN
        else
            val = diagΣ[Σptr]
            v[pos:pos+gp_cop-1] .= val
            Σptr += 1
        end
        pos += gp_cop
        id_ptr += 1

        # ----- pcf part ----------------------------------------------
        for _m in 1:dim
            if isnan(Σ_ids[id_ptr])
                v[pos:pos+grid_pcf-1] .= NaN
            else
                val = diagΣ[Σptr]
                v[pos:pos+grid_pcf-1] .= val
                Σptr += 1
            end
            pos += grid_pcf
            id_ptr += 1
        end
    end

    # ----- aggregates -------------------------------------------------
    @views v[pos:end] .= diagΣ[Σptr:Σptr+agg_count-1]

    return v
end



function run_smoother(L, x_filtered, sigma_filtered, x_updated, sigma_updated, T; E=nothing, u=nothing)
    x_smoothed = copy(x_updated)
    sigma_smoothed = copy(sigma_updated)
    dε_smoothed = zeros(size(x_filtered))

    cross_covariances = [zeros(size(sigma_filtered[1])) for _ in 1:(T-1)]

    for t in (T-1):-1:1
        x_smoothed[:, t], sigma_smoothed[t], cross_covariances[t] = perform_backward_pass(L,
            sigma_filtered[t+1], # was: sigma_updated[:, :, t],
            sigma_updated[t],
            x_updated[:, t],
            x_smoothed[:, t+1],
            x_filtered[:, t+1],
            sigma_smoothed[t+1],
            sigma_smoothed[t]
        )
    end

    # Initial difference
    dε_smoothed[:, 1] = x_smoothed[:, 1] - x_filtered[:, 1]
    for t in 2:T
        if E !== nothing && u !== nothing
            # :as_inputs — subtract exogenous input E * u_t to recover structural shocks
            dε_smoothed[:, t] = x_smoothed[:, t] - L * x_smoothed[:, t-1] - E * view(u, :, t)
        else
            dε_smoothed[:, t] = x_smoothed[:, t] - L * x_smoothed[:, t-1] # - B * u[:, t]
        end
    end

    return SmootherResults(x_filtered, sigma_filtered, x_updated, sigma_updated, x_smoothed, sigma_smoothed, cross_covariances, dε_smoothed)
end

function smoothing_distribution(A, Ω, smoother_results, n_jsd_draws)
    """Returns the smoothing distribution of states.
    From "Backward Simulation Methods for Monte Carlo Statistical Inference", section 1.7
    """

    @unpack x_filtered, sigma_filtered, x_updated, sigma_updated = smoother_results
    n = size(A, 1)
    T = size(x_filtered, 2)

    x_trajectories = zeros(n, T + 1, n_jsd_draws)

    # For each simulation ... 
    for i = 1:n_jsd_draws
        # Draw F_T+1 from the marginal density (from filter)
        x_trajectories[:, T+1, i] = x_filtered[:, end] .+ sigma_filtered[:, :, end] * rand(Normal(), n)
        for t in T:-1:1  # iteration backward from T to 1
            x_trajectories[:, t, i] = perform_backward_simulator(A, Ω, x_updated[:, t], sigma_updated[:, :, t], x_trajectories[:, t+1, i])
        end
    end

    return mean(x_trajectories, dims=3)
end


function perform_backward_simulator(A, Ω, x_updated, sigma_updated, prev_draw)
    """Runs a backward simulator: a standard formulation of the Kalman smoother cannot be used here, since it only provides the marginal densities p(F_t | Y, θ) 
    while for the blocked Gibbs sampler, the joint density p(F | Y, θ) = p(F1,···,FT+1 |Y,θ) is necessary 
    
    For linear-Gaussian models. 
    
    Algorithm:
    1. Draw FN+1 ∼ N(FˆN+1|N, PN+1|N) and iterate the following backward simulation recursion from t = N,...,1

    Ft ∼ MVN(μt, Mt),
    Vt = Pt|t − Pt|tA' (APt|tA' + Ω)−1 * APt|t,
    μt = Ft|t +  L * (Ft+1 - A * Ft|t)
    Lt=Pt|t * A(APt|tA +Ω)^-1
    pg. 15 of Backward Simulation Methods for Monte Carlo Statistical Inference

"""
    # eta 
    comp1 = sigma_updated * A'
    L = comp1 * inv(Ω + A * comp1)
    μ = x_updated .+ L * (prev_draw - A * x_updated)
    V = sigma_updated .- L * A * sigma_updated

    # Draw F_T
    draw = μ + V * rand(Normal(), size(μ, 1))

    return draw
end


function kalman_update_zeros!(x_updated, sigma_updated, measures, t, G, x̂_F, P_F, diag_Σ)
    """Used when data is not transformed ie no G."""

    # Select non-missing  
    # G_u = G[(!isnan).(measures), :]  # Only necessary when not collapsed, but no need for if condition 
    # Σ_u = Σ[(!isnan).(measures), (!isnan).(measures)]
    # Y_u = measures[(!isnan).(measures)]

    # local P, eta, logC 
    # if isempty(Y_u)  #, 1) == 0 
    #     G_u  = 0.0
    #     P    = 0.0

    # measures_valid = (!isnan).(measures)
    # G_u = view(G, measures_valid, :)
    # Σ_u = diag_Σ[measures_valid, measures_valid]
    # Σ_u = diag_Σ[measures_valid]   
    # Y_u = measures[measures_valid]

    local P, eta, logC
    if all(iszero, measures)
        # G_u  = zeros(eltype(G), 0, size(G, 2))
        # P    = zeros(eltype(x̂_F))
        # eta  = zeros(eltype(x̂_F), size(x̂_F)) # forecast error
        logC = 0.0
        x_updated[:, t] .= x̂_F
        # sigma_updated[:, :, t] .= P_F
        sigma_updated[t] .= P_F

    else
        # if isdiagonal(Σ_u) # TODO: with sequential updating, the ordering of the measurements doesnt matter 
        logC = sequential_kalman_update!(x_updated, sigma_updated, measures, t, G, x̂_F, P_F, diag_Σ)
        # else
        # logC = matrix_kalman_update!(x_updated, sigma_updated, Y_u, G_u, Σ_u, x̂_F, P_F, t)
        # end
    end
    return logC
end


function kalman_update_fast!(x_updated, sigma_updated, measures, G, x̂_F, P_F, diag_Σ, ws)
    """Used when data is not transformed ie no G."""

    measures_valid = (!isnan).(measures)

    if !any(measures_valid)
        logC = 0.0
        x_updated .= x̂_F
        sigma_updated .= P_F
    end

    G_u = view(G, measures_valid, :)
    Σ_u = diag_Σ[measures_valid]
    Y_u = measures[measures_valid]

    logC = sequential_kalman_update_fast!(x_updated, sigma_updated, Y_u, G_u, x̂_F, P_F, Σ_u, ws)


    return logC
end

function kalman_update!(x_updated, sigma_updated, measures, t, G, x̂_F, P_F, diag_Σ)
    """Used when data is not transformed ie no G."""

    measures_valid = (!isnan).(measures)
    G_u = view(G, measures_valid, :)
    Σ_u = diag_Σ[measures_valid]
    Y_u = measures[measures_valid]

    local P, eta, logC
    if isempty(Y_u)
        logC = 0.0
        x_updated[:, t] .= x̂_F
        sigma_updated[t] .= P_F
    else
        logC = sequential_kalman_update!(x_updated, sigma_updated, Y_u, t, G_u, x̂_F, P_F, Σ_u)
    end

    return logC
end

function matrix_kalman_update!(x_updated, sigma_updated, Y_u, G_u, Σ_u, x̂_F, P_F, t)
    P = G_u * P_F * G_u' .+ Σ_u
    P = @. 0.5 * (P + P')          # symmetry needed for consideration of the lu/cholesky
    # if isposdef(P) == false 
    #     P .= nearest_spd(P)
    # end 
    eta = Y_u .- G_u * x̂_F # forecast error

    logC = data_likelihood_contribution(eta, P)  # P is returned as inverted P

    # Kalman gain   
    K = P_F * G_u' * inv(P) # Because likelihood contribution inverts P already inplace, this is correct


    # Update Equations 
    comp1 = (I - K * G_u)
    # println(size(sigma_updated[t]))
    @views x_updated[:, t] .= x̂_F .+ K * eta

    if size(Y_u, 1) == 0
        copyto!(comp1, P_F)
    else
        copyto!(comp1, comp1 * P_F * comp1' .+ K * Σ_u * K')
    end #TODO: this could be an issue with copyto!()
    # println(size(comp1))
    @views sigma_updated[t] .= @. 0.5 * (comp1 + comp1')

    if isposdef(sigma_updated[t]) == false
        @views sigma_updated[t] .= nearest_spd(sigma_updated[t])
    end
    return logC
end


function time_varying_collapsed_kalman_update!(x_updated, sigma_updated, measures, t, x̂_F, P_F, Σ)
    """Grouped updating for correlated measurements and time varying Σ."""

    # Select non-missing  
    Y_u = measures[(!isnan).(measures)]

    # Terms 
    P = size(Y_u, 1) == 0 ? 0.0 : P_F .+ Σ
    P = 0.5 .* (P .+ P')          # symmetry needed for consideration of the lu/cholesky
    if isposdef(P) == false && P != 0.0
        P .= nearest_spd(P)
    end

    # Loglikelihood contribution from Model  
    eta = size(Y_u, 1) == 0 ? zeros(eltype(x̂_F), size(x̂_F)) : Y_u .- x̂_F # forecast error
    logC = data_likelihood_contribution(eta, P)  # TODO: P going forward is inverted already 

    # Kalman gain   
    K = P_F * P # Because likelihood contribution inverts P already inplace, this is correct

    # Update Equations 
    comp1 = (I - K)
    x_updated[:, t] .= x̂_F + K * eta

    if size(Y_u, 1) == 0
        copyto!(comp1, P_F)
    else
        copyto!(comp1, comp1 * P_F * comp1' .+ K * Σ * K')
    end
    @views sigma_updated[:, :, t] .= @. 0.5 * (comp1 + comp1')

    if isposdef(view(sigma_updated, :, :, t)) == false
        @views sigma_updated[:, :, t] .= nearest_spd(sigma_updated[:, :, t])
    end

    return logC
end



function perform_backward_pass(L, P_F, sigma_updated, x_updated, x_smooth, x_filtered, sigma_smooth, sigma_smooth_prev)
    """Rauch-Tung-Striebel Two-Pass Smoother: Just the backward pass.
    https://random-walks.org/content/misc/kalman/kalman.html
    The second pass runs backward in time in a sequence from the time 
    t of the last measurement, computing the smoothed state estimate 
    from the intermediate results stored on the forward pass.
    """
    # try
    KS = sigma_updated * L' * inv(lu(P_F))  # Kalman smoothing gain 
    x_s = x_updated .+ KS * (x_smooth .- x_filtered)
    sigma_s = sigma_updated + KS * (sigma_smooth - P_F) * KS'

    # Cross-covariance (Cov(x_t, x_{t-1} | T))
    sigma_cross = KS * sigma_smooth_prev
    return x_s, sigma_s, sigma_cross
    # catch e
    #     println(diag(P_F))
    # end
end


function perform_backward_pass_DK(A, B, D, Ω, Σ, G, y, u, x_filtered, sigma_filtered, x_updated, sigma_updated, T)
    # Durbin and Koopman (2012) -- used to get covariances between factors at time t and t-1, which would otherwise require another state space model with an augmented lag
    # Chapter 4.4, pp.87-91

    x_smoothed = copy(x_updated)
    sigma_smoothed = copy(sigma_updated)
    dε_smoothed = zeros(size(x_filtered))

    cross_covariances = [zeros(size(sigma_filtered[1])) for _ in 1:T]

    C_zero = zeros(eltype(A), size(B, 2), size(A, 1)) # C_{t,t-1 | T} matrix
    LOM = [A B;
        C_zero D]
    r = size(LOM, 1)
    rr = zeros(r, T)
    N = [zeros(r, r) for _ in 1:T]
    L = [zeros(r, r) for _ in 1:T]

    Σ_mat = diagm(Float64.(Σ))

    # one ahead prediction
    sigma_filtered_t_plus1 = LOM * sigma_filtered[T] * LOM' .+ Ω

    for tt = T:-1:2
        # Useful matrices
        mask = (!isnan).(y[:, tt])
        C = G[tt][mask, :]
        CPCR = C * sigma_filtered[tt] * C' + Σ_mat[mask, mask]
        iCPCR = inv(CPCR)
        # CiCPCRC = C' * iCPCR * C
        CiCPCRC = C' * (CPCR \ C)         # Cholesky solve under the hood


        # New
        L[tt] = LOM - LOM * sigma_filtered[tt] * CiCPCRC
        rr[:, tt-1] = C' * iCPCR * (y[mask, tt] - C * x_filtered[:, tt]) + L[tt]' * rr[:, tt]

        # rr[:, tt-1] = C' * iCPCR * (y[mask, tt] - C * x_filtered[:, tt]) + L[tt]' * reshape(rr[:, tt], r, 1)

        # Auxiliary matrices
        N[tt-1] = CiCPCRC + L[tt]' * N[tt] * L[tt]
        # L[tt] = LOM - LOM * sigma_filtered[tt] * CiCPCRC

        # smooth updates
        x_smoothed[:, tt] = x_filtered[:, tt] + sigma_filtered[tt] * rr[:, tt-1]  # xittm = x_filtered
        sigma_smoothed[tt] = sigma_filtered[tt] - sigma_filtered[tt] * N[tt-1] * sigma_filtered[tt]
        if tt == T
            # We need C_{t,t-1 | T} = C'_{t,t+1 | T}
            cross_covariances[tt] = (sigma_filtered[tt] * L[tt]' * (I - N[tt] * sigma_filtered_t_plus1))' # this is C_{t,t+1 | T} transposed
        else
            cross_covariances[tt] = (sigma_filtered[tt] * L[tt]' * (I - N[tt] * sigma_filtered[tt+1]))'
        end
    end
    # t = 1
    mask = .!isnan.(y[:, 1])
    C = G[1][mask, :]
    CPCR = C * sigma_filtered[1] * C' .+ Σ_mat[mask, mask]
    CiCPCRC = C' * (CPCR \ C)                   # no explicit inv
    L[1] = LOM - LOM * sigma_filtered[1] * CiCPCRC

    cross_covariances[1] =
        (sigma_filtered[1] * L[1]' *
         (I - N[1] * sigma_filtered[2]))'

    # cross_covariances[1] = (sigma_filtered[1] * L[1]' * (I - N[1] * sigma_filtered[2]))'


    for t in 2:T
        dε_smoothed[:, t] = x_smoothed[:, t] - LOM * x_smoothed[:, t-1] # - B * u[:, t]
    end
    dε_smoothed[:, 1] = x_smoothed[:, 1] - x_filtered[:, 1]

    return SmootherResults(x_filtered, sigma_filtered, x_updated, sigma_updated, x_smoothed, sigma_smoothed, cross_covariances, dε_smoothed)
end


function data_likelihood_contribution(eta, P)
    """Generates the likelihood contribution p(Y | Θ) for one iteration. """
    # local alarm
    # if sequential == true

    #     logC = ifelse(any(P .< 0), -1.e12, (log(2π) .+ log.(P) + (eta .* eta) .* inv.(P))[1] * -0.5)
    #     # logC = (log(2π) .+ log.(P) + (eta .* eta) .* inv.(P))[1] * -0.5  # white noise case. 
    #     logC = ifelse(isnan(logC), -1.e12, logC)

    # else
    logdet_P, sign_logdet = logabsdet(P)
    # inv_P =  qr(P, Val{true}()) \ Matrix{Float64}(I, size(P)) 
    if sign_logdet < 0
        # println("reversing at ", t)
        logC = -1.e12
        println(logdet_P)
    else
        LinearAlgebra.inv!(lu!(P))  # in-place inverse. alternative to lu! is cholesky!, but I would need a spd matrix, but it doubles the time it takes 
        logC = -(length(measures)*log(2π)+logdet_P.+eta'*(P*eta))[1] .* 0.5  # white noise case. P = inv_P header
        if isnan(logC)
            logC = -1.e12
        end
    end
    # end
    return logC
end

# function joint_likelihood_contribution(A, B, measures, state_trajectories, u, Ω, Σ)
#     """Generates the likelihood contribution for all iterations for the data AND state trajectories, p(X,Y | θ).
#     https://blogs.sas.com/content/iml/2020/07/15/multivariate-normal-log-likelihood.html
#     Conditional Multivariate Normal.
#     """
#     T    = size(measures, 2)
#     n    = size(A, 1)
#     logC = 0.0
#     idtt = Matrix(I, n, n)
#     Γ    = vcat(hcat(A, B), hcat(idtt, zeros(n, size(B, 2))))  # Time-invariant 

#     for t = 1:T
#         if all(measures[:, t] .=== NaN)  
#             logC += 0
#         else
#             # Components for the likelihood
#             a    = vcat(state_trajectories[:, t+1], measures[:, t])
#             b    = vcat(state_trajectories[:, t], u[:, t])
#             Γb   = Γ * b
#             Π_u  = Matrix(blockdiag(sparse(Ω), sparse(Σ[:,:,t])))

#             # mvn = MvNormal(Γb, Π_u)
#             # Check 
#             logdet_Π_u, sign_logdet = logabsdet(Π_u)
#             if sign_logdet < 0
#                 println("reversing at ", t)
#                 logC += -1.e9
#             else
#                 # logC += loglikelihood(mvn, a) 
#                 logC += -0.5 * (n * log(2π) + logdet_Π_u + tr(inv(Π_u) * (a * a') - (a * Γb') - (Γb * a') + (Γb * Γb')))
#                 # logC += -0.5 * (length(measures[t]) * log(2π) + logdet_Π_u + invquad(inv(Symmetric(Π_u)), (a .- Γb)))            
#             end
#         end
#     end
#     return logC
# end


function time_varying_collapse(measures, G, Σ)
    """Collapsing obsevational vector following Jungbacker and Koopman (2008), also found in Koopman (2012) pg.161. 

    Note: An alternative to this approach is to update the state vector one step at a time.
    Note: this only works for time-varying measurement error variances. """

    T = size(measures, 2)
    g = size(G[1], 1)
    n = size(G[1], 2)
    # C   = Vector{Matrix{eltype(Σ)}}(undef, T)
    Σ̂ = zeros(eltype(Σ), n, n, T)

    Σ⁻¹ = inv(Σ)
    ŷ = zeros(eltype(Σ), n, T)

    for t in range(1, stop=T)
        if all(iszero, G[t])
            ŷ[:, t] .= NaN
        else
            m = measures[:, t]
            slice = (!isnan).(m)
            @inbounds G_u = G[t][slice, :]  # Selects respective rows when the measurement is available 
            @inbounds Σ_u = Σ[slice, slice]
            @inbounds Σ⁻¹_u = Σ⁻¹[slice, slice]
            Y_u = m[slice]

            comp1 = G_u' * Σ⁻¹_u  # 31 x 240
            copyto!(comp1, inv(lu(comp1 * G_u)) * comp1)  # 31 x 240 A projection mat. Take Cholesky here?
            ŷ[:, t] = comp1 * Y_u
            Σ̂[:, :, t] = comp1 * Σ_u * comp1'
        end
    end
    return ŷ, Σ̂
end


function pre_multiply_part_correct!(y, G, Σ̂⁻¹²) #TODO: put outside likeli function 
    T = size(y, 2)
    N = size(G[1], 1)

    ŷ = deepcopy(y) # zeros(eltype(Σ̂⁻¹²), N, T)
    Ĝ = deepcopy(G) # Vector{Matrix{eltype(Σ̂⁻¹²)}}(undef, T)

    # Replace NaNs with zero 
    y[isnan.(y)] .= 0.0

    for t in range(1, stop=T)
        Ĝ[t] = Σ̂⁻¹² * G[t]
        ŷ[:, t] = Σ̂⁻¹² * y[:, t]
    end
    return ŷ, Ĝ
end

"""
    pre_multiply_part!(y, G, Σ̂⁻¹², factor_count; n_dist)

Premultiplies the **distributional** rows of the measurement vector `y` and
of each selection–projection matrix `G[t]` by the
pre-whitening matrix `Σ̂⁻¹²  ≡  (Σ̂_dist)^{-1/2}`.

* `y`            : (N_dist+q) × T     stacked distribution + controls  
* `G`            : Vector of length `T`, each element is (N_dist+q) × (r+q)  
* `Σ̂⁻¹²`        : N_dist × N_dist     whitening for the distributional block  
* `factor_count` : r                  (columns 1:r of `G[t]` carry the factors)  
* `n_dist`       : number of distributional rows; default = size(Σ̂⁻¹²,1)

The control rows (indices `n_dist+1 : end`) are left untouched, so they
remain effectively noise-free or carry their own tiny measurement variance.
The function modifies `y` and `G` in place and also returns them for
convenience.
"""
function pre_multiply_part!(y, G, Σ̂⁻¹², factor_count_with_lags)
    T = size(y, 2)
    N = size(y, 1)                       #  N_dist + q
    n_dist = size(Σ̂⁻¹², 1)

    # deep copies so we do not overwrite the original arrays
    ŷ = deepcopy(y)
    Ĝ = deepcopy(G)

    dist_idx = 1:n_dist                  # rows to be pre-whitened

    for t in 1:T
        # which distributional entries are actually observed at time t ?
        obs_dist = .!isnan.(y[dist_idx, t])

        if any(obs_dist)
            # --------   scale the measurement vector
            ŷ[dist_idx[obs_dist], t] =
                Σ̂⁻¹²[obs_dist, obs_dist] * y[dist_idx[obs_dist], t]

            # --------   scale the corresponding rows of G[t]
            Ĝ[t][dist_idx[obs_dist], 1:factor_count_with_lags] =
                Σ̂⁻¹²[obs_dist, obs_dist] * G[t][dist_idx[obs_dist], 1:factor_count_with_lags]
            # rows for the control columns (factor_count+1 : end) are zeros anyway
        else
            ŷ[dist_idx, t] .= NaN
            Ĝ[t][dist_idx, 1:factor_count_with_lags] .= 0
        end
        # control rows (n_dist+1 : N) remain unchanged
    end
    return ŷ, Ĝ
end

# function pre_multiply_part!(y, G, Σ̂⁻¹², factor_count) #TODO: put outside likeli function 
#     T = size(y, 2)
#     N = size(G[1], 1)

#     ŷ = deepcopy(y) # zeros(eltype(Σ̂⁻¹²), N, T)
#     Ĝ = deepcopy(G) # Vector{Matrix{eltype(Σ̂⁻¹²)}}(undef, T)

#     # Replace NaNs with zero 
#     # y[isnan.(y)] .= 0.0

#     for t in range(1, stop=T)
#         condition = (!isnan).(y[:, t])
#         binary_cond = condition[:]

#         if all(iszero, binary_cond) # G[t])
#             ŷ[:, t] .= NaN
#             Ĝ[t] = zeros(eltype(Σ̂⁻¹²), N, factor_count)  # as long as it's a matrix 
#         else
#             m = y[:, t]
#             slice = (!isnan).(m)
#             Ĝ[t][slice, :] = Σ̂⁻¹²[slice, slice] * G[t][slice, :]
#             ŷ[slice, t] = Σ̂⁻¹²[slice, slice] * m[slice]
#         end
#     end
#     return ŷ, Ĝ
# end


function time_invariant_collapse!(measures, G, Σ, proj, number_of_dfs)
    """Collapsing obsevational vector following 
    @misc{jungbacker2015likelihood,
    title={Likelihood-based dynamic factor analysis for measurement and forecasting},
    author={Jungbacker, Borus and Koopman, Siem Jan},
    year={2015},
    publisher={Oxford University Press Oxford, UK}}
  """
    # Params 
    T = size(measures, 2)
    n = size(G[1], 2)

    # Containers 
    # Σ̂     = zeros(eltype(Σ), n, n, T)
    # ŷ     = zeros(eltype(Σ), n, T)
    # Ĝ     = Vector{Matrix{eltype(Σ)}}(undef, T)
    Σ̂ = zeros(n, n, T)
    ŷ = zeros(n, T)
    Ĝ = Vector{Matrix{Float64}}(undef, T)
    # e     = zeros(eltype(Σ), T)

    # Time-invariant Σ̂
    Σ_mat = diagm(Float64.(Σ))
    Σ⁻¹ = inv(Σ_mat)
    Z = repeat(proj, number_of_dfs)
    comp1 = Z' * Σ⁻¹
    Σ̂ = comp1 * Σ_mat * comp1'

    # Σ̂ is not diagonal => make diagonal 
    Vt = svd(Σ̂).Vt
    VΣ̂V = Vt * Σ̂ * Vt'

    # Off-diagonal should be zero 
    f(c) = c.I[1] == c.I[2]
    VΣ̂V[filter(!f, CartesianIndices(size(VΣ̂V)))] .= 0

    # For likelihood: Store generalized least squares (GLS) residual vector for data vector y_t, det(Σ)
    # comp2 = (I - Z * inv(comp1 * Z) * comp1)
    # det_Σ = logabsdet(Σ)[1]
    # det_Σ̂ = logabsdet(Σ̂)[1]
    # LR    =  det_Σ̂ - det_Σ # log ratio

    for t in range(1, stop=T)
        if all(iszero, vec((!isnan).(measures[:, t]))) # G[t])
            ŷ[:, t] .= NaN
            Ĝ[t] = zeros(n, n)  # as long as it's a matrix 
        else
            m = measures[:, t]
            slice = (!isnan).(m)
            Y_u = m[slice]
            Ĝ[t] = Vt * comp1[:, slice] * G[t][slice, :]
            ŷ[:, t] = Vt * comp1[:, slice] * Y_u
        end
    end
    # a = ŷ[:, vec(mapslices(col -> all((!isnan).(col)), ŷ, dims = 1))]
    # println(diag(cov(a')))
    return ŷ, VΣ̂V, Ĝ
end


function sequential_kalman_update!(x_updated, sigma_updated, Y_u, t, G_u, x̂_F, P_F, diag_Σ)
    """Updating the kalman filter one measurement at a time. Only possible when the measurement covariance matrix is diagonal."""
    # Get each object
    # diag_Σ = Σ[diagind(Σ)]  
    x = copy(x̂_F)
    sig = copy(P_F)
    valid_measures = size(G_u, 1)

    # Pre-allocate 
    log_c = 0.0
    P = eltype(x̂_F)[0.0]
    eta = eltype(x̂_F)[0.0]

    K = copy(x̂_F)
    k_len = length(K)
    id = Diagonal{eltype(P_F)}(ones(k_len))
    sig_r = Matrix{eltype(P_F)}(undef, k_len, 1)
    sig_l = Matrix{eltype(P_F)}(undef, 1, 1)
    Kg = Matrix{eltype(P_F)}(undef, k_len, k_len)
    id_sig = Matrix{eltype(P_F)}(undef, k_len, k_len)
    id_Kg = Matrix{eltype(P_F)}(undef, k_len, k_len)
    K_eta = Matrix{eltype(P_F)}(undef, k_len, 1)
    g = Matrix{eltype(P_F)}(undef, k_len, 1)


    # For each measurement ... 
    for i in 1:valid_measures
        # Get the first row 
        @views g .= G_u[i, :]

        mul!(sig_r, sig, g)
        mul!(sig_l, g', sig_r)

        P .= sig_l .+ diag_Σ[i]
        eta .= Y_u[i] .- g' * x # forecast error

        for (e, p) in zip(eta, P)
            try
                log_c += (log(2π).+log(p)+(e*e).*inv(p))[1] * -0.5  # white noise case. # this is the likelihood contribution of 1 time period, given factors and θ
            catch ee
                # println("Error in log_c, $i, $t: ", ee)
                log_c += 0 #-1.e12
            end
        end
        # end

        if isnan(log_c)
            log_c += -1.e12
        end
        mul!(K, sig_r, inv.(P))

        # Store updates and continue 
        mul!(Kg, K, g')
        id_Kg .= (id .- Kg)
        mul!(id_sig, id_Kg, sig)
        copy!(sig, id_sig)

        mul!(K_eta, K, eta)
        x .+= K_eta
    end

    # Store last update 
    x_updated[:, t] .= x
    sigma_updated[t] .= @. 0.5 * (sig + sig')

    # cholesky!(sigma_updated[t], NoPivot(); check=false)

    return log_c
end

function sequential_kalman_update_fast!(x_updated, sigma_updated, Y_u, G_u, x̂_F, P_F, diag_Σ, ws)
    """Updating the kalman filter one measurement at a time. Only possible when the measurement covariance matrix is diagonal."""
    # Get each object
    # diag_Σ = Σ[diagind(Σ)]  
    # x = copy(x̂_F)
    # sig = copy(P_F)
    # valid_measures = size(G_u, 1)

    log_c = 0.0
    # n = length(x̂_F)
    # P = eltype(x̂_F)[0.0]
    # η = eltype(x̂_F)[0.0]

    @inbounds for m in eachindex(Y_u)
        # g ← m-th measurement row (viewed, not copied)
        @views ws.g .= G_u[m, :]

        # innovation covariance vector and scalar variance
        mul!(ws.sig_r, P_F, ws.g)            # sig_r = P_F * g      (n×1)
        σ_g = dot(ws.g, ws.sig_r)            # scalar g' P_F g
        P = σ_g + diag_Σ[m]                # scalar variance
        η = Y_u[m] - dot(ws.g, x̂_F)       # scalar innovation

        try
            log_c += -0.5 * (log(2π) + log(P) + (η * η) / P)
        catch ee
            # println("Error in log_c, $m, $t: ", ee)
            log_c += 0
        end

        invP = 1 / P
        @. ws.K = ws.sig_r * invP            # K = sig_r / P

        # P_F ← P_F − K * sig_r'
        BLAS.ger!(-1.0, ws.K, ws.sig_r, P_F)

        # x̂_F ← x̂_F + K * η   (axpy!: x̂_F += η*K)
        BLAS.axpy!(η, ws.K, x̂_F)
    end

    # # Pre-allocate 
    # log_c = 0.0
    # P = eltype(x̂_F)[0.0]
    # eta = eltype(x̂_F)[0.0]

    # # For each measurement ... 
    # for i in 1:valid_measures
    #     # Get the first row 
    #     @views ws.g .= G_u[i, :]

    #     mul!(ws.sig_r, P_F, ws.g)
    #     mul!(ws.sig_l, ws.g', ws.sig_r)

    #     P .= ws.sig_l .+ diag_Σ[i]
    #     eta .= Y_u[i] .- ws.g' * x̂_F # forecast error

    #     for (e, p) in zip(eta, P)
    #         try
    #             log_c += (log(2π).+log(p)+(e*e).*inv(p))[1] * -0.5  # white noise case. # this is the likelihood contribution of 1 time period, given factors and θ
    #         catch ee
    #             # println("Error in log_c, $i, $t: ", ee)
    #             log_c += 0 #-1.e12
    #         end
    #     end
    #     # end

    #     if isnan(log_c)
    #         log_c += -1.e12
    #     end
    #     mul!(ws.K, ws.sig_r, inv.(P))

    #     # Store updates and continue 
    #     # mul!(ws.Kg, ws.K, ws.g')
    #     # ws.id_Kg .= (ws.id .- ws.Kg)
    #     # mul!(ws.id_sig, ws.id_Kg, P_F)
    #     # copy!(P_F, ws.id_sig)
    #     BLAS.ger!(-1.0, ws.K, ws.sig_r, P_F)            # -1 * K * sig_r' + P_F

    # mul!(ws.K_eta, ws.K, eta)
    # x̂_F .+= ws.K_eta
    # end

    # Store last update 
    x_updated .= x̂_F
    # sigma_updated .= @. 0.5 * (P_F + P_F')
    sigma_updated .= Symmetric(P_F)

    # cholesky!(sigma_updated, NoPivot(); check=false)

    return log_c
end


function isdiagonal(A)
    isdiagonal = all(x -> x == Dual(0, zero(1)), offdiag(A))
    return isdiagonal
end

function offdiag(A)
    [A[ι] for ι in CartesianIndices(A) if ι[1] ≠ ι[2]]
end

# function offdiag(A)
#     println("hi")
#     D = size(A)[1]
#     v = zeros(eltype(A), D*(D-1), 1)
#     println(typeof(v))
#     for i in 1:D
#         for j in 1:(i-1)
#             v[(i-1)*(D-1)+j, :] = A[j,i]
#         end
#         for j in (i+1):D
#             v[(i-1)*(D-1)+j-1, :] = A[j,i]
#         end
#     end
#     v
# end



















# function gen_spd_sigma(A, P_F, Ω)
#     tm_P  = A * P_F * A' .+ Ω 
#     tm_P .= isposdef(tm_P) ? tm_P : nearest_spd(tm_P)   
#     return tm_P
# end

# function gen_spd_P(G_u, P_F, Σ_u)
#     tm_P = G_u * P_F * G_u' .+ Σ_u
#     tm_P .= isposdef(tm_P) ? tm_P : nearest_spd(tm_P)   
# return tm_P
# end

# P = nearest_spd(rand(1000,1000))
# pinv(P)

# # Benchmarking the different inversions. slowest to fastest
# @btime qr(P, Val{true}()) \ Matrix{Float64}(I, size(P))  # 600 ms and 3000 allocations (10x slower, 2900 more allocations)
# @btime I / P;
# @btime inv(P);
# @btime inv(lu(P));
# @btime inv(cholesky(P));











