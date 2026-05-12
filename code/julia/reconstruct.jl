"""
Standalone reconstruction of marginal quantile functions and the trivariate
copula density from `*_coefficients_{normal,average}.csv`.

Mirrors `code/python/reconstruct.py`. Self-contained: no dependency on the
model code in this repo. Only CSV.jl + DataFrames.jl + the stdlib.

Usage:

    include("reconstruct.jl")
    using .DistributionalReconstruction

    r = Reconstruction("PSID_coefficients_normal.csv")
    quantile_at(r, "2008-Q3", :consum, [0.1, 0.5, 0.9])
    copula_density_at(r, "2008-Q3", 0.5, 0.5, 0.5)
    copula_density_grid(r, "2008-Q3"; n = 30)

    # Match the published `ciw_*` cell probability masses on a 10x10x10 grid:
    copula_pmf_grid(r, "2008-Q3")   # 10x10x10 array summing to 1

Legendre basis on [0, 1]: Q_o(u) = sqrt(2o + 1) * P_o(2u - 1).
"""

module DistributionalReconstruction

using CSV
using DataFrames
using LinearAlgebra
using Statistics

export Reconstruction,
    FactorMap,
    ProjectionMap,
    load_projection_map,
    quantile_at,
    copula_density_at,
    copula_density_grid,
    copula_pmf_grid,
    quantile_from_row,
    copula_density_from_row,
    Q_m,
    I_m,
    available_dates,
    GRID_COP,
    GRID_PCF,
    MEASURES


# Model constants — match Reconstruction.jl / ModelPrep.jl
const GRID_COP = 12            # 12^3 - 34 = 1694 mutable copula coefs
const GRID_PCF = 12            # 12 marginal coefs per measure
const D = 3
const N_MUTABLE_COP = 1694
const N_QUANTILE = D * GRID_PCF      # 36
const MEASURES = (:consum, :income, :wealth)


# -----------------------------------------------------------------------------
# Orthonormal Legendre basis on [0, 1]
# -----------------------------------------------------------------------------
"""
    legendre_P(n, x)

Standard Legendre polynomial P_n on [-1, 1]. Vectorized.
"""
function legendre_P(n::Int, x)
    n == 0 && return one.(x)
    n == 1 && return float.(x)
    pm2 = one.(x)
    pm1 = float.(x)
    p = similar(pm1)
    for k in 2:n
        @. p = ((2k - 1) * x * pm1 - (k - 1) * pm2) / k
        pm2 = pm1
        pm1 = p
        p = similar(pm1)
    end
    return pm1
end

"""
    Q_m(m, u)

Orthonormal Legendre basis on [0, 1]: Q_m(u) = sqrt(2m+1) * P_m(2u - 1).
"""
Q_m(m::Int, u) = sqrt(2m + 1) .* legendre_P(m, 2 .* u .- 1)

"""
    I_m(m, u)

I_m(u) = integral of Q_m(s) on [0, u]. Closed form from Bonnet's recursion:

    I_0(u) = u
    I_m(u) = (P_{m+1}(2u-1) - P_{m-1}(2u-1)) / (2 * sqrt(2m+1)),  m ≥ 1

The boundary terms at s = 0 cancel because P_n(-1) = (-1)^n.
"""
function I_m(m::Int, u)
    m == 0 && return float.(u)
    x = 2 .* u .- 1
    return (legendre_P(m + 1, x) .- legendre_P(m - 1, x)) ./ (2 * sqrt(2m + 1))
end


# -----------------------------------------------------------------------------
# Index mapping: x1..x1694 -> mutable CartesianIndices in (12, 12, 12)
# -----------------------------------------------------------------------------
"""Iteration order matches Julia's CartesianIndices((12,12,12)): i fastest, then j, then k.
An entry is *immutable* if at least D - 1 = 2 of its 1-indexed components equal 1."""
function build_mutable_indices()
    out = NTuple{3, Int}[]
    for k in 1:GRID_COP, j in 1:GRID_COP, i in 1:GRID_COP
        n_ones = (i == 1) + (j == 1) + (k == 1)
        if n_ones < D - 1
            push!(out, (i, j, k))
        end
    end
    @assert length(out) == N_MUTABLE_COP
    return out
end

const MUTABLE_INDICES = build_mutable_indices()


"""
    unflatten_to_kappa(row)

Take the first 1694 entries of a coefficient row (mutable copula κ) and
build the full (12, 12, 12) tensor. The immutable leading entry
κ[1,1,1] = 1 (uniform marginals); all other immutable entries are 0.

(Indices in the returned tensor are 1-based; polynomial order = index - 1.)
"""
function unflatten_to_kappa(row::AbstractVector{<:Real})
    κ = zeros(GRID_COP, GRID_COP, GRID_COP)
    @inbounds for (n, (i, j, k)) in enumerate(MUTABLE_INDICES)
        κ[i, j, k] = row[n]
    end
    κ[1, 1, 1] = 1.0
    return κ
end


# -----------------------------------------------------------------------------
# Public type
# -----------------------------------------------------------------------------
"""
    Reconstruction(csv_path)

Lazy reader for a `*_coefficients_*.csv` file.
"""
struct Reconstruction
    df::DataFrame
    dates::Vector{String}
    csv_path::String
end

function Reconstruction(csv_path::AbstractString)
    df = CSV.read(csv_path, DataFrame)
    if !("time" in names(df))
        error("Expected a `time` column in $csv_path.")
    end
    n_data = ncol(df) - 1
    expected = N_MUTABLE_COP + N_QUANTILE
    if n_data != expected
        error("$csv_path has $n_data data columns; expected $expected (is this a *_coefficients_*.csv file?).")
    end
    dates = string.(df.time)
    return Reconstruction(df, dates, String(csv_path))
end


available_dates(r::Reconstruction) = copy(r.dates)


function _row_for_date(r::Reconstruction, date::AbstractString)
    idx = findfirst(==(date), r.dates)
    isnothing(idx) && throw(KeyError("Date $date not in file. " *
        "First/last available: $(r.dates[1]) .. $(r.dates[end])"))
    # The first column is `time`; data starts at column 2.
    return Vector{Float64}(r.df[idx, 2:end])
end


function _xi(r::Reconstruction, date::AbstractString, measure::Symbol)
    measure in MEASURES || error("measure must be one of $MEASURES, got $measure")
    m_idx = findfirst(==(measure), MEASURES) - 1   # 0-based
    row = _row_for_date(r, date)
    start = N_MUTABLE_COP + m_idx * GRID_PCF + 1
    return row[start : start + GRID_PCF - 1]
end


_kappa(r::Reconstruction, date::AbstractString) = unflatten_to_kappa(_row_for_date(r, date))


# -----------------------------------------------------------------------------
# Public evaluation methods
# -----------------------------------------------------------------------------
"""
    quantile_at(r, date, measure, u)

Marginal quantile function Ξ⁻¹_{m,t}(u) at quantile points `u ∈ [0, 1]`.
Returns a scalar if `u` is a scalar, otherwise a vector.
"""
function quantile_at(r::Reconstruction, date::AbstractString, measure::Symbol, u)
    ξ = _xi(r, date, measure)
    u_vec = u isa Real ? [Float64(u)] : Vector{Float64}(u)
    basis = hcat([Q_m(o, u_vec) for o in 0:GRID_PCF-1]...)   # (n_u, GRID_PCF)
    out = basis * ξ
    return u isa Real ? out[1] : out
end


"""
    copula_density_at(r, date, u_c, u_y, u_w)

Trivariate copula density dC_t(u_c, u_y, u_w) at one or many points.
Scalar inputs return a scalar; equal-length vectors return a vector
(element-wise).
"""
function copula_density_at(r::Reconstruction, date::AbstractString, u_c, u_y, u_w)
    κ = _kappa(r, date)
    uc = u_c isa Real ? [Float64(u_c)] : Vector{Float64}(u_c)
    uy = u_y isa Real ? [Float64(u_y)] : Vector{Float64}(u_y)
    uw = u_w isa Real ? [Float64(u_w)] : Vector{Float64}(u_w)
    @assert length(uc) == length(uy) == length(uw)
    B_c = reduce(vcat, [reshape(Q_m(o, uc), 1, :) for o in 0:GRID_COP-1])  # (GRID_COP, n)
    B_y = reduce(vcat, [reshape(Q_m(o, uy), 1, :) for o in 0:GRID_COP-1])
    B_w = reduce(vcat, [reshape(Q_m(o, uw), 1, :) for o in 0:GRID_COP-1])
    n = length(uc)
    out = zeros(n)
    @inbounds for a in 1:GRID_COP, b in 1:GRID_COP, c in 1:GRID_COP
        κabc = κ[a, b, c]
        κabc == 0 && continue
        for k in 1:n
            out[k] += κabc * B_c[a, k] * B_y[b, k] * B_w[c, k]
        end
    end
    return u_c isa Real ? out[1] : out
end


"""
    copula_density_grid(r, date; n = 30)

Evaluate dC_t on a regular n×n×n grid of (u_c, u_y, u_w) ∈ (0, 1).
"""
function copula_density_grid(r::Reconstruction, date::AbstractString; n::Int = 30)
    κ = _kappa(r, date)
    u = collect(range(1e-6, 1 - 1e-6; length = n))
    B = reduce(vcat, [reshape(Q_m(o, u), 1, :) for o in 0:GRID_COP-1])   # (GRID_COP, n)
    grid = zeros(n, n, n)
    @inbounds for a in 1:GRID_COP, b in 1:GRID_COP, c in 1:GRID_COP
        κabc = κ[a, b, c]
        κabc == 0 && continue
        @views grid .+= κabc .* (B[a, :] .* reshape(B[b, :], 1, :, 1) .* reshape(B[c, :], 1, 1, :))
    end
    return grid
end


"""
    copula_pmf_grid(r, date; grid = DEFAULT_GRID)

Discrete probability mass on a 10x10x10 grid of (u_c, u_y, u_w) cells with
edges `[1e-6, 0.1, 0.2, ..., 0.9, 0.999999]`. Matches the published
`ciw_<i><j><k>` columns exactly.

Steps (same as `Reconstruction.jl: generate_copula_densities` + `cdf_to_pdf`):
  1) Evaluate the copula CDF at the 10 grid edge points via
        C(x_i, x_j, x_k) = sum_m κ[m] * I_{m1}(x_i) * I_{m2}(x_j) * I_{m3}(x_k)
  2) 3D inclusion-exclusion finite difference → cell probability mass.
  3) Clip negative cells to 0 and renormalize to sum to 1.
"""
const DEFAULT_GRID = collect(0.1:0.1:1.0)
DEFAULT_GRID[end] -= 1e-6

function copula_pmf_grid(r::Reconstruction, date::AbstractString;
        grid::AbstractVector{<:Real} = DEFAULT_GRID)
    κ = _kappa(r, date)
    n = length(grid)
    # I[m+1, ui] = I_m(grid[ui])
    Imat = reduce(vcat, [reshape(I_m(m, grid), 1, :) for m in 0:GRID_COP-1])   # (GRID_COP, n)
    cdf = zeros(n, n, n)
    @inbounds for a in 1:GRID_COP, b in 1:GRID_COP, c in 1:GRID_COP
        κabc = κ[a, b, c]
        κabc == 0 && continue
        @views cdf .+= κabc .* (Imat[a, :] .* reshape(Imat[b, :], 1, :, 1) .* reshape(Imat[c, :], 1, 1, :))
    end
    pmf = zeros(n, n, n)
    @inbounds for i in 1:n, j in 1:n, k in 1:n
        c   = cdf[i, j, k]
        a   = i > 1 ? cdf[i-1, j, k] : 0.0
        b2  = j > 1 ? cdf[i, j-1, k] : 0.0
        d2  = k > 1 ? cdf[i, j, k-1] : 0.0
        ab  = (i > 1 && j > 1) ? cdf[i-1, j-1, k] : 0.0
        ad  = (i > 1 && k > 1) ? cdf[i-1, j, k-1] : 0.0
        bd  = (j > 1 && k > 1) ? cdf[i, j-1, k-1] : 0.0
        abd = (i > 1 && j > 1 && k > 1) ? cdf[i-1, j-1, k-1] : 0.0
        pmf[i, j, k] = c - a - b2 - d2 + ab + ad + bd - abd
    end
    pmf[pmf .< 0] .= 0
    s = sum(pmf)
    s > 0 && (pmf ./= s)
    return pmf
end


# -----------------------------------------------------------------------------
# Row-level evaluators (raw 1730-vector → moments)
# -----------------------------------------------------------------------------
function _extract_xi(row::AbstractVector{<:Real}, measure::Symbol)
    measure in MEASURES || error("measure must be one of $MEASURES, got $measure")
    m_idx = findfirst(==(measure), MEASURES) - 1   # 0-based
    start = N_MUTABLE_COP + m_idx * GRID_PCF + 1
    return row[start : start + GRID_PCF - 1]
end


"""
    quantile_from_row(row, measure, u)

Marginal quantile Ξ⁻¹_{measure}(u) computed from a raw 1730-coefficient row
(handy when `row` is synthesized rather than read from a CSV).
"""
function quantile_from_row(row::AbstractVector{<:Real}, measure::Symbol, u)
    ξ = _extract_xi(row, measure)
    u_vec = u isa Real ? [Float64(u)] : Vector{Float64}(u)
    basis = hcat([Q_m(o, u_vec) for o in 0:GRID_PCF-1]...)
    out = basis * ξ
    return u isa Real ? out[1] : out
end


"""
    copula_density_from_row(row, u_c, u_y, u_w)

Copula density dC(u_c, u_y, u_w) from a raw 1730-coefficient row.
"""
function copula_density_from_row(row::AbstractVector{<:Real}, u_c, u_y, u_w)
    κ = unflatten_to_kappa(row)
    uc = u_c isa Real ? [Float64(u_c)] : Vector{Float64}(u_c)
    uy = u_y isa Real ? [Float64(u_y)] : Vector{Float64}(u_y)
    uw = u_w isa Real ? [Float64(u_w)] : Vector{Float64}(u_w)
    @assert length(uc) == length(uy) == length(uw)
    B_c = reduce(vcat, [reshape(Q_m(o, uc), 1, :) for o in 0:GRID_COP-1])
    B_y = reduce(vcat, [reshape(Q_m(o, uy), 1, :) for o in 0:GRID_COP-1])
    B_w = reduce(vcat, [reshape(Q_m(o, uw), 1, :) for o in 0:GRID_COP-1])
    n = length(uc)
    out = zeros(n)
    @inbounds for a in 1:GRID_COP, b in 1:GRID_COP, c in 1:GRID_COP
        κabc = κ[a, b, c]
        κabc == 0 && continue
        for k in 1:n
            out[k] += κabc * B_c[a, k] * B_y[b, k] * B_w[c, k]
        end
    end
    return u_c isa Real ? out[1] : out
end


# -----------------------------------------------------------------------------
# FactorMap — OLS projection from smoothed factors F_t to coefficient row.
#
#     coef_t  =  α  +  β · F_t  +  ε_t        (one regression per coef)
#
# Useful for counterfactuals: shift a factor, predict the implied coefficient
# row, then pass that row to `quantile_from_row` / `copula_density_from_row`.
# -----------------------------------------------------------------------------
"""
    FactorMap(factors_csv, coefs_csv; n_factors=8)

Learn the OLS map from the published smoothed factors `F_t` to a coefficient
time series. After fitting, `predict(fm, F)` maps any factor vector back to a
1730-element coefficient row.

Rows where any coefficient is `NaN` are dropped from the fit (e.g., PSID has
no consumption coefficients before 1999); the fit is run on the remaining,
fully-observed subsample. Use `n_factors = 43` to use all latent factors in
`smoothed_factors.csv`; the default `8` matches the count of distributional
factors used in the paper's diagnostics.
"""
struct FactorMap
    beta::Matrix{Float64}                # (1 + n_factors, n_coefs)
    y_means::Vector{Float64}             # (n_coefs,)
    y_stds::Vector{Float64}              # (n_coefs,)
    r_squared::Vector{Float64}           # (n_coefs,)
    n_factors::Int
    n_obs_total::Int
    n_obs_used::Int
    dates_used::Vector{String}
end


function FactorMap(factors_csv::AbstractString, coefs_csv::AbstractString;
        n_factors::Int = 8)
    F_df = CSV.read(factors_csv, DataFrame)
    Y_df = CSV.read(coefs_csv, DataFrame)
    for nm in ("time",)
        nm in names(F_df) || error("$factors_csv missing column 'time'")
        nm in names(Y_df) || error("$coefs_csv missing column 'time'")
    end
    factor_cols = ["x$i" for i in 1:n_factors]
    missing_cols = [c for c in factor_cols if !(c in names(F_df))]
    isempty(missing_cols) || error("smoothed_factors.csv missing $(missing_cols)")
    # Align on common dates
    F_df.time = string.(F_df.time)
    Y_df.time = string.(Y_df.time)
    common = sort(collect(intersect(Set(F_df.time), Set(Y_df.time))))
    isempty(common) && error("no overlapping dates between factors and coefs")
    F_df_c = filter(r -> r.time in Set(common), F_df)
    Y_df_c = filter(r -> r.time in Set(common), Y_df)
    sort!(F_df_c, :time); sort!(Y_df_c, :time)
    F = Matrix{Float64}(F_df_c[!, factor_cols])
    Y = Matrix{Float64}(Y_df_c[!, Not(:time)])
    # Drop NaN rows
    valid = .!vec(any(isnan, Y; dims = 2))
    sum(valid) >= n_factors + 2 ||
        error("only $(sum(valid)) clean rows; need at least n_factors+2 = $(n_factors + 2)")
    F_use, Y_use = F[valid, :], Y[valid, :]
    y_means = vec(mean(Y_use; dims = 1))
    y_stds  = vec(std(Y_use; dims = 1, corrected = false))
    safe_stds = ifelse.(y_stds .> 0, y_stds, 1.0)
    Y_std = (Y_use .- y_means') ./ safe_stds'
    X = hcat(ones(size(F_use, 1)), F_use)
    β = X \ Y_std
    Y_std_hat = X * β
    ss_res = vec(sum((Y_std .- Y_std_hat) .^ 2; dims = 1))
    ss_tot = vec(sum(Y_std .^ 2; dims = 1))
    r2 = [t > 0 ? 1 - r / t : NaN for (r, t) in zip(ss_res, ss_tot)]
    dates_used = [d for (d, ok) in zip(common, valid) if ok]
    return FactorMap(β, y_means, y_stds, r2, n_factors, length(common), sum(valid), dates_used)
end


"""
    predict(fm::FactorMap, factors)

Map factor values (length `fm.n_factors`) back to a coefficient row.
Returns a Vector{Float64} of length 1730 in the original CSV's units.
"""
function predict(fm::FactorMap, factors::AbstractVector{<:Real})
    length(factors) == fm.n_factors ||
        error("factors must have length $(fm.n_factors), got $(length(factors))")
    x = vcat(1.0, Vector{Float64}(factors))     # (1 + K,)
    y_std_hat = vec(x' * fm.beta)               # (n_coefs,)
    return y_std_hat .* fm.y_stds .+ fm.y_means
end



# -----------------------------------------------------------------------------
# ProjectionMap — exact reconstruction via the model's PCA loading matrix.
#
#     coef_t = stds ⊙ (Γ · F_dist_t) + means + trend_t
#
# Loaded from the `projection/` folder exported by `Reconstruction.jl` at the
# end of a model run. This is the canonical Distributional_Oil-style map; the
# OLS `FactorMap` above is an approximation that doesn't require Γ.
# -----------------------------------------------------------------------------
"""
    ProjectionMap

Holds the exported projection artifacts for one dataset (`G`, `means`, `stds`,
`trend`). Use `load_projection_map(proj_dir, dataset; trend_kind="normal")` to
construct.
"""
struct ProjectionMap
    G::Matrix{Float64}                   # (n_coefs, n_dist_state)
    means::Vector{Float64}               # (n_coefs,)
    stds::Vector{Float64}                # (n_coefs,) — already expanded per coefficient
    trend::Array{Float64}                # (T, n_coefs) for :normal, or (n_coefs,) for :average
    trend_dates::Union{Nothing, Vector{String}}
    trend_mode::Symbol
    dataset::String
    n_factors::Int
    n_coefs::Int
end


"""
    load_projection_map(proj_dir, dataset; trend_kind="normal")

Read the projection artifacts written by `Reconstruction.jl` for the given
dataset (e.g., `"PSID"`, `"SCF"`, `"CEX"`).

- `proj_dir`: path to the projection folder (typically
  `data/synthetic/projection/` or `<run>/from_mcmc/data/projection/`).
- `trend_kind`: `"normal"` (HP trend, date-anchored — in-sample) or
  `"average"` (time-mean trend — extrapolation-friendly).
"""
function load_projection_map(proj_dir::AbstractString, dataset::AbstractString;
        trend_kind::AbstractString = "normal")
    trend_kind in ("normal", "average") ||
        error("trend_kind must be \"normal\" or \"average\", got $trend_kind")
    G = Matrix{Float64}(CSV.read(joinpath(proj_dir, "$(dataset)_projection.csv"), DataFrame))
    means = Vector{Float64}(CSV.read(joinpath(proj_dir, "$(dataset)_means.csv"), DataFrame).mean)
    stds  = Vector{Float64}(CSV.read(joinpath(proj_dir, "$(dataset)_stds.csv"),  DataFrame).std)
    if trend_kind == "normal"
        tdf = CSV.read(joinpath(proj_dir, "$(dataset)_trend_normal.csv"), DataFrame)
        trend_dates = string.(tdf.time)
        trend = Matrix{Float64}(select(tdf, Not(:time)))      # (T, n_coefs)
        trend_mode = :normal
    else
        trend = Vector{Float64}(CSV.read(joinpath(proj_dir, "$(dataset)_trend_average.csv"),
                                          DataFrame).trend_average)
        trend_dates = nothing
        trend_mode = :average
    end
    return ProjectionMap(G, means, stds, trend, trend_dates, trend_mode,
                         String(dataset), size(G, 2), size(G, 1))
end


"""
    predict(pm::ProjectionMap, factors; date=nothing)

Map factor values (length `pm.n_factors`) to a coefficient row. For
`trend_kind="normal"`, pass `date` (e.g. `"2008-Q3"`) to pick the right trend
row.
"""
function predict(pm::ProjectionMap, factors::AbstractVector{<:Real};
        date::Union{Nothing, AbstractString} = nothing)
    length(factors) == pm.n_factors ||
        error("factors must have length $(pm.n_factors), got $(length(factors))")
    base = pm.stds .* (pm.G * Vector{Float64}(factors)) .+ pm.means
    if pm.trend_mode == :average
        return base .+ pm.trend
    end
    # trend_mode == :normal
    if date === nothing
        return base .+ vec(pm.trend[1, :])
    end
    idx = findfirst(==(String(date)), pm.trend_dates)
    isnothing(idx) && throw(KeyError("date $date not in trend_normal; first/last: " *
        "$(pm.trend_dates[1])..$(pm.trend_dates[end])"))
    return base .+ vec(pm.trend[idx, :])
end


# Shared moment-extractors dispatch on both factor-map types.
quantile_at(pm::ProjectionMap, factors::AbstractVector{<:Real}, measure::Symbol, u;
        date::Union{Nothing, AbstractString} = nothing) =
    quantile_from_row(predict(pm, factors; date), measure, u)

copula_density_at(pm::ProjectionMap, factors::AbstractVector{<:Real}, u_c, u_y, u_w;
        date::Union{Nothing, AbstractString} = nothing) =
    copula_density_from_row(predict(pm, factors; date), u_c, u_y, u_w)

summary(pm::ProjectionMap) =
    "ProjectionMap: dataset=$(pm.dataset), K=$(pm.n_factors), " *
    "n_coefs=$(pm.n_coefs), trend=$(pm.trend_mode)"


# End-to-end convenience: factors → moments (OLS FactorMap)
quantile_at(fm::FactorMap, factors::AbstractVector{<:Real}, measure::Symbol, u) =
    quantile_from_row(predict(fm, factors), measure, u)

copula_density_at(fm::FactorMap, factors::AbstractVector{<:Real}, u_c, u_y, u_w) =
    copula_density_from_row(predict(fm, factors), u_c, u_y, u_w)


function summary(fm::FactorMap)
    med  = median(filter(!isnan, fm.r_squared))
    qs   = quantile(filter(!isnan, fm.r_squared), [0.25, 0.75])
    return "FactorMap: K=$(fm.n_factors), T_used=$(fm.n_obs_used) " *
           "(dropped $(fm.n_obs_total - fm.n_obs_used) NaN rows of $(fm.n_obs_total)), " *
           "R² median=$(round(med, digits = 3)), " *
           "R² P25/P75=($(round(qs[1], digits = 3)), $(round(qs[2], digits = 3)))"
end


end # module
