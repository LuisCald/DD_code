# TODO: For separation strategy 
# https://arxiv.org/pdf/1408.4050.pdf
# https://www.tandfonline.com/doi/pdf/10.1198/016214508000000724?casa_token=Q3nW55JKV_gAAAAA:XRsJ9h7RgMGAcdZTFQBS-5KQrnYV4-VA_7-S7tnhh_iA47TmJ6FfSiLJ0qwFnEswsX0ItcD-_IoJBRkw
# For the variance covariance matrices, we adopt a separation strategy. 
# This means we split Σ into its variances and correlations, Σ = Λ * R * Λ
# When Σ is just a diagonal, R is the identity (perfect correlation), so the prior is just the product of LogNormals x 2
# When Σ is a full matrix, the prior on R = Δ * Q * Δ, where Q  ̃ IW(ν, I) and Δ = diagonal matrix, with Δᵢᵢ = (Qᵢᵢ)^(-1/2)
# To get R, we can just do cov2cor(Σ, sqrt.(diag(Σ))), where the prior is:  det(R)^(−1/2 * (ν+k+1)) * prod(diag(R⁻¹))^ν/2 
# The prior on Σ is then the product of a logNormal * DR * logNormal
# https://stats.stackexchange.com/questions/282380/what-is-the-correct-form-of-metropolis-hasting-step-in-scaled-inverse-wishart-pr 

function set_shock_priors(priors, factor_count, agg_count, n_param)

    # Distributions and draws 
    local Ωf_mode, Ωy_mode, Σ_mode
    try
        Ωf_mode = repeat([StatsBase.mode(priors[2])], factor_count)
        Ωy_mode = repeat([StatsBase.mode(priors[3])], agg_count)
        Σ_mode = [StatsBase.mode(priors[3+i]) for i in 1:n_param]
    catch e
        Ωf_mode = repeat([StatsBase.mode(priors[2].untruncated)], factor_count)
        Ωy_mode = repeat([StatsBase.mode(priors[3].untruncated)], agg_count)
        Σ_mode = [StatsBase.mode(priors[2+i].untruncated) for i in 1:n_param] # For the truncated cauchy
    end

    # Replace the last n_aggs entries of Σ_mode to zero
    Ω_prior = Diagonal(vcat(Ωf_mode, Ωy_mode))
    Σ_prior = Diagonal(Σ_mode)

    return Ω_prior, Σ_prior
end

function correlation_pdf(R, ν)
    k = copy(ν)
    df = (-1 / 2) * (ν + k + 1)
    R⁻¹ = inv(R)
    R = det(R)^df * prod(diag(R⁻¹))^(-ν / 2)
end

function correlation_logpdf(R, ν)
    return log(correlation_pdf(R, ν))
end



# So first we draw from the product distribution 
function LogNormal²(β, ξ)
    """The prior on the variances is the product of logNormals."""
    dgf(n) = (n + 1) / 2
    pd = ProductDistribution(LogNormal(β, ξ), LogNormal(β, ξ))
    rv = Random.rand!(pd, zeros(2, 1000000))
    rv[2, :] = sqrt.(rv[2, :])
    rv = vec(prod(rv, dims=1))

    rv_cdf = ecdf(rv)
    obs = unique(rv)
    # Percentiles 
    cdf_obs = map(x -> rv_cdf(x), sort(obs))

    # Generate probability mass function for next part 
    pmf = Float64[]
    for i in 1:length(obs)
        if i == 1
            append!(pmf, cdf_obs[1])
        else
            append!(pmf, cdf_obs[i] - cdf_obs[i-1])
        end
    end
    # Generate, using MLE, the distribution object 
    d = DiscreteNonParametric(sort(obs), pmf)
    sample = rand(d, 1000000)

    # Make continuous 
    U = kde(sample)
    return U
end

function Distributions.insupport(dist::UnivariateKDE, x::AbstractVector{T}) where T<:Real
    insup = false
    if all(x .> 0) && all(x .< 1e9)
        insup = true
    end
    insup
end

Distributions.logpdf(dist::UnivariateKDE, x::AbstractVector{T}) where T<:Real = sum(log.(pdf(dist, x)))

Distributions.mode(dist::UnivariateKDE, ν::Int, tau²::Real) = tau² / (((ν + 1) * 0.5) + 1)

# We use the notation Σ ∼ SIW (ν, Λ, b, ξ) to refer to this prior. 
# By construction, the σᵢ = δᵢ * sqrt(Qᵢᵢ)  and     Σij = δᵢ * δⱼ * Qᵢⱼ 
# Thus each standard deviation is the product of a log-normal and the square root of a scaled inv-χ2

function estimate_tau²s(MV, pcs, number_of_dfs)
    Σ_tau² = 0.0
    Ω_tau² = 0.0

    for i in 1:number_of_dfs
        cond = vec(mapslices(col -> all((!isnan).(col)), MV[i], dims=1))
        data = MV[i][:, cond]
        var = diag(cov(data'))
        println(var)
        Σ_tau² += mean(var)
    end
    Σ_tau² = Σ_tau² ./ number_of_dfs

    # Estimating scale parameter from the Data for Ω
    cond = vec(mapslices(col -> all((!isnan).(col)), pcs, dims=1))
    data = pcs[:, cond]
    var = diag(cov(data'))
    println(var)
    Ω_tau² = mean(var)
    return Ω_tau², Σ_tau²
end

# TODO: in case we prefer this over the log normals only 
function LogNormal_SqrtScaledInverseChiSquared(ν, tau²)
    dgf(n) = (n + 1) / 2
    pd = ProductDistribution(LogNormal(0, 1), InverseGamma(dgf(ν), dgf(ν) * tau²))
    rv = Random.rand!(pd, zeros(2, 1000000))
    rv[2, :] = sqrt.(rv[2, :])
    rv = vec(prod(rv, dims=1))

    rv_cdf = ecdf(rv)
    obs = unique(rv)
    # Percentiles 
    cdf_obs = map(x -> rv_cdf(x), sort(obs))

    # Generate probability mass function for next part 
    pmf = Float64[]
    for i in 1:length(obs)
        if i == 1
            append!(pmf, cdf_obs[1])
        else
            append!(pmf, cdf_obs[i] - cdf_obs[i-1])
        end
    end
    # Generate, using MLE, the distribution object 
    d = DiscreteNonParametric(sort(obs), pmf)
    sample = rand(d, 1000000)

    # Make continuous 
    U = kde(sample)
    return U
end

function Distributions.insupport(dist::UnivariateKDE, x::AbstractVector{T}) where T<:Real
    insup = false
    if all(x .> 0) && all(x .< 1e9)
        insup = true
    end
    insup
end

Distributions.logpdf(dist::UnivariateKDE, x::AbstractVector{T}) where T<:Real = sum(log.(pdf(dist, x)))

Distributions.mode(dist::UnivariateKDE, ν::Int, tau²::Real) = tau² / (((ν + 1) * 0.5) + 1)