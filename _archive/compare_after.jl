"""
Standalone script: compare HANK_PSID_1 survey data to truth_data_1,
plotting 5 quintiles for each measure. Saves plots with "_after" tag.
"""

using DataFrames, CSV, Dates, PeriodicalDates, Plots

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
data_dir = "/Users/lc/Dropbox/Distributional_Dynamics/5_Code/bld/filtered"
plot_dir = "/Users/lc/Dropbox/Distributional_Dynamics/5_Code/bld/quintile_plots_after"
mkpath(plot_dir)

# -----------------------------------------------------------------------------
# Weighted equal-mass bin means (from testing.jl)
# -----------------------------------------------------------------------------
function weighted_bin_means_by_mass(x::AbstractVector, w::AbstractVector, n::Int)
    @assert length(x) == length(w)
    @assert n >= 1

    keep = trues(length(x))
    if eltype(x) <: Union{Missing,Any}
        keep .&= .!ismissing.(x)
    end
    if eltype(w) <: Union{Missing,Any}
        keep .&= .!ismissing.(w)
    end
    xv = Float64.(x[keep])
    wv = Float64.(w[keep])

    totw = sum(wv)
    totw > 0 || error("weights must sum to > 0")

    if n == 1
        return [sum(xv .* wv) / totw]
    end

    ix = sortperm(xv)
    xs = xv[ix]
    ws = wv[ix]

    target = totw / n
    num = zeros(Float64, n)
    den = zeros(Float64, n)

    bin = 1
    filled = 0.0
    tol = 1e-14 * target

    @inbounds for i in eachindex(xs)
        xi = xs[i]
        wi = ws[i]
        wi <= 0 && continue

        while wi > 0 && bin <= n
            rem = target - filled
            if rem <= tol
                bin += 1
                filled = 0.0
                continue
            end
            take = min(wi, rem)
            num[bin] += xi * take
            den[bin] += take
            wi -= take
            filled += take
            if filled >= target - tol
                bin += 1
                filled = 0.0
            end
        end
    end

    return [den[b] > 0 ? num[b] / den[b] : NaN for b in 1:n]
end

# -----------------------------------------------------------------------------
# Compute survey quintiles by time
# -----------------------------------------------------------------------------
function compute_quintiles_df(df::DataFrame, meas::Symbol)
    out = DataFrame(; time=[], m0=Float64[], m1=Float64[], m2=Float64[],
        m3=Float64[], m4=Float64[], m5=Float64[])
    for t in sort(unique(df.time))
        mask = df.time .== t
        x = df[mask, meas]
        w = df[mask, :weight]

        bins = weighted_bin_means_by_mass(x, w, 5)

        keep = .!ismissing.(x) .& .!ismissing.(w)
        xv = Float64.(x[keep])
        wv = Float64.(w[keep])
        avg = isempty(xv) ? NaN : sum(xv .* wv) / sum(wv)

        push!(out, (time=t, m1=bins[1], m2=bins[2], m3=bins[3],
            m4=bins[4], m5=bins[5], m0=avg))
    end
    sort!(out, :time)
    return out
end

# -----------------------------------------------------------------------------
# Loop over all id_d and generate plots
# -----------------------------------------------------------------------------
# n_dfs = length(shocks_vec)  # from generate_HANK_micro.jl output
for id_d in eachindex(export_idx)
    # Load data
    psid_df = CSV.read(joinpath(data_dir, "HANK_PSID_$(id_d).csv"), DataFrame)
    truth_wide = CSV.read(joinpath(data_dir, "truth_data_$(id_d).csv"), DataFrame)

    # Parse time columns to QuarterlyDate
    psid_df.time = QuarterlyDate.(psid_df.time)
    truth_wide.time = QuarterlyDate.(truth_wide.time)

    # Plot: truth vs PSID survey quintiles (relative to mean)
    # -----------------------------------------------------------------------
    measures = [:consum, :income, :wealth]

    for meas in measures
        # Check columns exist
        truth_q_cols = [Symbol("$(meas)$(q)q") for q in 1:5]
        truth_mean_col = Symbol("$(meas)_per_hh")

        if !all(c -> hasproperty(truth_wide, c), truth_q_cols)
            @warn "Missing truth quintile columns for $meas, skipping"
            continue
        end
        if !hasproperty(truth_wide, truth_mean_col)
            @warn "Missing truth mean column $truth_mean_col, skipping"
            continue
        end
        if !hasproperty(psid_df, meas)
            @warn "Missing $meas in PSID data, skipping"
            continue
        end

        # Truth: quintile / mean (relative to mean)
        truth_mean_vals = Float64.(truth_wide[!, truth_mean_col])

        # Survey quintiles from PSID
        quintiles_df = compute_quintiles_df(psid_df, meas)

        # Map truth means to the survey time grid (handles differing time spans)
        truth_mean_by_time = Dict(truth_wide.time .=> truth_mean_vals)
        truth_mean_on_quintile_grid = [get(truth_mean_by_time, t, NaN) for t in quintiles_df.time]

        for qq in 1:5
            # Truth series
            y_truth = Float64.(truth_wide[!, truth_q_cols[qq]]) ./ truth_mean_vals

            # Survey series
            y_data = Float64.(quintiles_df[!, Symbol("m$(qq)")]) ./ Float64.(quintiles_df.m0) #truth_mean_on_quintile_grid

            p = Plots.plot(
                truth_wide.time, y_truth;
                title="$(meas) quintile $(qq) (relative to mean)",
                label="Truth",
                linewidth=1.5,
                legend=:best,
            )
            Plots.scatter!(
                p, quintiles_df.time, y_data;
                label="PSID Survey",
                markersize=3,
            )

            outpath = joinpath(plot_dir, "HANK_PSID_$(id_d)_$(meas)_q$(qq)_after.png")
            savefig(p, outpath)
            @info "Saved" outpath
        end
    end
end
