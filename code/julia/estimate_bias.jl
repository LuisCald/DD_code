"""
Estimate the Legendre polynomial approximation bias by comparing:
  - Truth: direct weighted bin means from the linearized HANK distribution
  - Legendre: fit shifted Legendre polynomials to PSID survey microdata,
              integrate over quintile intervals, compare bin means

Uses the exact functions from DataConstructor.jl and QuadGK for integration.
"""

using CSV, DataFrames, Statistics, Printf, LinearAlgebra, StatsBase, QuadGK

# ═══════════════════════════════════════════════════════════════════════════════
#  Functions copied verbatim from DataConstructor.jl
# ═══════════════════════════════════════════════════════════════════════════════

# DataConstructor.jl line 929
function inverse_hyperbolic_sine(x)
    return log.(x .+ sqrt.(x .^ 2 .+ 1))
end

# DataConstructor.jl line 933
function reverse_inverse_hyperbolic_sine(x)
    return (exp.(2 .* x) .- 1) ./ (2 .* exp.(x))
end

# DataConstructor.jl line 1247
function legendre_polynomial(m, x)
    if m == 0
        return 1.0
    elseif m == 1
        return x
    else
        P_prev_prev = 1.0
        P_prev = x
        P_current = 0.0

        for n in 2:m
            P_current = ((2n - 1) * x * P_prev - (n - 1) * P_prev_prev) / n
            P_prev_prev, P_prev = P_prev, P_current
        end

        return P_current
    end
end

# DataConstructor.jl line 1269
function Q_m(m, x)
    L_m = legendre_polynomial(m, 2x - 1)
    return sqrt(2m + 1) * L_m
end

# DataConstructor.jl line 1215
function series_estimator(data, weights, order)
    Phi = zeros(length(data), order + 1)
    s_weights = cumsum(weights) / sum(weights)

    for i in eachindex(data)
        for j in 0:order
            Phi[i, j+1] = Q_m(j, s_weights[i])
        end
    end

    coefficients = zeros(order + 1)
    for j in 1:order+1
        for i in 1:length(data)
            coefficients[j] += weights[i] * Phi[i, j] * data[i]
        end
        coefficients[j] /= sum(weights)
    end

    return coefficients
end

# DataConstructor.jl line 1121
function scale_to_aggregates(vals, wts, correction)
    tot_data = mean(vals, weights(wts))
    multiplier = abs(correction / tot_data)

    if multiplier > 20
        return fill(NaN, length(vals))
    end

    scaled = vals .* multiplier
    tot_scale = correction - mean(scaled, weights(wts))
    scaled .= scaled .+ tot_scale
    return scaled
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Pipeline matching DataConstructor.jl treat_quantile_functions (line 1195-1202)
# ═══════════════════════════════════════════════════════════════════════════════

function fit_legendre_coefficients(vals::Vector{Float64}, wts::Vector{Float64},
                                   correction::Float64, order::Int)
    # 1. scale_to_aggregates (DataConstructor.jl line 1177)
    scaled = scale_to_aggregates(vals, wts, correction)
    if any(isnan, scaled)
        return nothing
    end

    # 2. Sort by value (DataConstructor.jl line 1196: sort!(non_missing_scaled, rv))
    idx = sortperm(scaled)
    sorted_vals = scaled[idx]
    sorted_wts = wts[idx]

    # 3. IHS transform divided by correction (DataConstructor.jl line 1198)
    t_rv = inverse_hyperbolic_sine(sorted_vals ./ correction)

    # 4. series_estimator (DataConstructor.jl line 1199)
    coefs = series_estimator(t_rv, sorted_wts, order)
    return coefs
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Integration using QuadGK (replaces hand-rolled Gauss-Legendre)
# ═══════════════════════════════════════════════════════════════════════════════

function evaluate_qf(coefs, u)
    val = 0.0
    for j in eachindex(coefs)
        val += coefs[j] * Q_m(j - 1, u)
    end
    return val
end

function legendre_bin_means_quadgk(coefs::Vector{Float64}, avg::Float64;
                                    bin_edges::Vector{Float64}=[0.0, 0.5, 0.9, 1.0])
    n_bins = length(bin_edges) - 1
    means = zeros(n_bins)
    for i in 1:n_bins
        a, b = bin_edges[i], bin_edges[i+1]
        # Integrate ihs_inv(Q(u)) over [a, b] using adaptive quadrature
        integral, _ = quadgk(u -> reverse_inverse_hyperbolic_sine(evaluate_qf(coefs, u)), a, b;
                             rtol=1e-10)
        means[i] = integral * avg / (b - a)
    end
    return means
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════════

data_dir = "/Users/lc/Desktop/BASEtoolbox.jl/bld/baseline_micro_noestim/simulated_data"
econ = "1"

truth_df = CSV.read(joinpath(data_dir, "HANK_truth_$(econ).csv"), DataFrame)
psid_df  = CSV.read(joinpath(data_dir, "HANK_PSID_$(econ).csv"), DataFrame)

bin_edges_q5 = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

order = 11  # Legendre polynomial order (grid_pcf = 12 means order 11)

measures = ["income", "wealth", "consum"]
q_cols = Dict(
    "income" => ["income1q", "income2q", "income3q", "income4q", "income5q"],
    "wealth" => ["wealth1q", "wealth2q", "wealth3q", "wealth4q", "wealth5q"],
    "consum" => ["consum1q", "consum2q", "consum3q", "consum4q", "consum5q"],
)
mean_cols = Dict(
    "income" => "income_per_hh",
    "wealth" => "wealth_per_hh",
    "consum" => "consum_per_hh",
)

psid_times = unique(psid_df.time)

println("=" ^ 80)
println("Legendre Polynomial Bias Estimation (order = $order)")
println("Using DataConstructor.jl functions + QuadGK")
println("=" ^ 80)

for meas in measures
    println("\n─── $meas ───")

    bias_by_q = [Float64[] for _ in 1:5]
    pct_bias_by_q = [Float64[] for _ in 1:5]

    for t in psid_times
        truth_row = truth_df[truth_df.time .== t, :]
        if nrow(truth_row) == 0
            continue
        end

        truth_means = [truth_row[1, c] for c in q_cols[meas]]
        truth_avg = truth_row[1, mean_cols[meas]]

        micro = psid_df[psid_df.time .== t, :]
        vals = Float64.(micro[!, Symbol(meas)])
        wts = Float64.(micro[!, :weight])

        # Fit using exact DataConstructor.jl pipeline
        coefs = fit_legendre_coefficients(vals, wts, truth_avg, order)
        if coefs === nothing
            continue
        end

        # Integrate using QuadGK
        leg_means = legendre_bin_means_quadgk(coefs, truth_avg; bin_edges=bin_edges_q5)

        for q in 1:5
            b = leg_means[q] - truth_means[q]
            pct = 100 * b / abs(truth_means[q])
            push!(bias_by_q[q], b)
            push!(pct_bias_by_q[q], pct)
        end
    end

    @printf("  %-6s %10s %10s %10s %10s\n", "Quint", "Mean Bias", "Median %", "Mean %", "Std %")
    for q in 1:5
        mb = mean(bias_by_q[q])
        med_pct = median(pct_bias_by_q[q])
        mn_pct = mean(pct_bias_by_q[q])
        sd_pct = std(pct_bias_by_q[q])
        @printf("  Q%-5d %10.4f %9.2f%% %9.2f%% %9.2f%%\n", q, mb, med_pct, mn_pct, sd_pct)
    end
end

# ─── Order sweep ─────────────────────────────────────────────────────────

for meas in ["income", "wealth", "consum"]
    println("\n" * "=" ^ 80)
    println("Order sweep for $meas Q5 (mean % bias across survey periods)")
    println("=" ^ 80)

    for ord in [5, 7, 9, 11, 15, 19, 23, 27, 31]
        biases = Float64[]
        for t in psid_times
            truth_row = truth_df[truth_df.time .== t, :]
            if nrow(truth_row) == 0; continue; end

            truth_q5 = truth_row[1, "$(meas)5q"]
            truth_avg = truth_row[1, "$(meas)_per_hh"]

            micro = psid_df[psid_df.time .== t, :]
            vals = Float64.(micro[!, Symbol(meas)])
            wts = Float64.(micro[!, :weight])

            coefs = fit_legendre_coefficients(vals, wts, truth_avg, ord)
            if coefs === nothing; continue; end

            leg_means = legendre_bin_means_quadgk(coefs, truth_avg; bin_edges=bin_edges_q5)
            push!(biases, 100 * (leg_means[5] - truth_q5) / abs(truth_q5))
        end
        @printf("  Order %2d:  mean bias = %+7.2f%%,  median = %+7.2f%%,  std = %5.2f%%\n",
                ord, mean(biases), median(biases), std(biases))
    end
end

println("\nDone.")
