# Import
# using Plots
# file_psid = "/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/HANK_PSID_1.csv"
# psid_df = CSV.read(file_psid, DataFrame);
using DataFrames
using CSV
using Dates
using PeriodicalDates
shocks_vec = out[1]
surveys = out[2]
truth_vec = truth_mode == :none ? nothing : out[3]

# -----------------------------------------------------------------------------
# 0) Select which chunk/index to analyze
# -----------------------------------------------------------------------------

# Settings
save_filtered = true
dt_id = 10

@assert truth_vec !== nothing "Need truth_mode != :none to use truth moments in this script."

truth_df = truth_vec[dt_id]
shocks_df = shocks_vec[dt_id]

# Keep a convenient handle to one survey dataset for "data" plotting.
# NOTE: surveys are ordered by the survey_names in your generator script.
psid_df = (surveys === nothing) ? DataFrame() : surveys[3][dt_id]

"""Compute weighted mean in each of `n` equal-mass bins.

This is robust when the series is grid-valued (lots of ties), which is exactly the case
for `:liquid` / `:illiqd` in simulated survey draws.
"""
function weighted_bin_means_by_mass(x::AbstractVector, w::AbstractVector, n::Int)
    @assert length(x) == length(w)
    @assert n >= 1

    # Drop missings defensively (matters for real data).
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

weighted_bin_means_5(x::AbstractVector, w::AbstractVector) = weighted_bin_means_by_mass(x, w, 5)

function compute_quintiles_df(psid_df::DataFrame, meas::Symbol)
    quintiles_df = DataFrame(;
        time=[],
        m0=Float64[],
        m1=Float64[],
        m2=Float64[],
        m3=Float64[],
        m4=Float64[],
        m5=Float64[],
    )

    ts = sort(unique(psid_df.time))
    for t in ts
        mask_t = psid_df.time .== t
        x = psid_df[mask_t, meas]
        w = psid_df[mask_t, :weight]

        mean_by_time = weighted_bin_means_5(x, w)

        # Weighted average (matches the binning)
        keep = .!ismissing.(x) .& .!ismissing.(w)
        xv = Float64.(x[keep])
        wv = Float64.(w[keep])
        mean_meas = isempty(xv) ? NaN : sum(xv .* wv) / sum(wv)

        push!(
            quintiles_df,
            (
                time=t,
                m1=mean_by_time[1],
                m2=mean_by_time[2],
                m3=mean_by_time[3],
                m4=mean_by_time[4],
                m5=mean_by_time[5],
                m0=mean_meas,
            ),
        )
    end

    sort!(quintiles_df, :time)
    return quintiles_df
end

## correction series is now embedded in truth as *_per_hh columns
using Plots

# Build a wide truth table with quintile columns (n=5 only): income1..income5, etc.
function truth_quintiles_wide(truth_df::DataFrame)
    qdf = truth_df[(truth_df.kind.==:binmean).&(truth_df.n.==5), [:time, :var, :bin, :value]]
    # Quintile columns get a trailing `q` to make it explicit they are bin means.
    # Example: income1q, consum5q, wealth3q, ...
    qdf.col = Symbol.(string.(qdf.var) .* string.(qdf.bin) .* "q")
    wide = unstack(qdf, :time, :col, :value)
    sort!(wide, :time)
    return wide
end

"""Build a wide truth table with mean columns named `X_per_hh`.

Uses `kind=:mean` rows from `truth_df` and produces one row per `time` with columns like
`consum_per_hh`, `income_per_hh`, etc.
"""
function truth_means_wide(truth_df::DataFrame)
    mdf = truth_df[truth_df.kind.==:mean, [:time, :var, :value]]
    if nrow(mdf) == 0
        wide = DataFrame(time=QuarterlyDate[])
        return wide
    end
    mdf.col = Symbol.(string.(mdf.var) .* "_per_hh")
    wide = unstack(mdf, :time, :col, :value)
    sort!(wide, :time)
    return wide
end

# Combine quintiles (n=5) with means (X_per_hh) for convenience.
truth_df1 = truth_means_wide(truth_df)
truth_wide_5 = leftjoin(truth_quintiles_wide(truth_df), truth_df1; on=:time)
# if save_filtered
#     out_dir = joinpath("/Users/lc/Dropbox/Distributional_Dynamics/5_Code/bld/", "testing")
#     mkpath(out_dir)

#     # Export wide truth for each index/chunk in `truth_vec`.
#     # (This avoids re-running the simulation just to save all chunks.)
#     if truth_vec === nothing
#         CSV.write(joinpath(out_dir, "truth_data_$(dt_id).csv"), truth_wide_5)
#     else
#         for idx in eachindex(truth_vec)
#             tdf = truth_vec[idx]
#             tdf1 = truth_means_wide(tdf)
#             wide5 = leftjoin(truth_quintiles_wide(tdf), tdf1; on=:time)
#             CSV.write(joinpath(out_dir, "truth_data_$(idx).csv"), wide5)
#         end
#     end
# end

# Optional transform for plotting/exporting moments.
# - :none  -> plot levels
# - :asinh -> robust "signed log"-like transform (safe for negatives/zeros)
# - :log   -> requires strictly positive values (will error otherwise)
transform_kind::Symbol = :none

function apply_transform(x::AbstractVector{<:Real}, kind::Symbol)
    if kind === :none
        return Float64.(x)
    elseif kind === :asinh
        return asinh.(Float64.(x))
    elseif kind === :log
        xf = Float64.(x)
        any(xf .<= 0.0) &&
            error("log transform requires strictly positive values; use transform_kind=:asinh for variables that can be ≤ 0")
        return log.(xf)
    else
        error("Unknown transform_kind=$kind (use :none, :asinh, or :log)")
    end
end

# Export directory for plots
plot_dir = joinpath(@__DIR__, "bld", "quintile_plots")
mkpath(plot_dir)

transform_tag = transform_kind === :none ? "levels" : String(transform_kind)

# Plot series relative to their contemporaneous mean (q/avg).
# - If true: plots (quintile / mean) for both Simulated and Data.
# - If false: plots raw levels.
relative_to_mean::Bool = true
relative_tag = relative_to_mean ? "relmean" : "levels"

# Plot all 5 quintiles for each measure, each in its own figure.
preferred_measures = [:consum, :income, :wealth, :liquid, :illiqd]
available_truth = Set(Symbol.(unique(truth_df.var)))
available_data = Set(Symbol.(propertynames(psid_df)))
measures = [m for m in preferred_measures if (m in available_truth) && (m in available_data)]

for meas in measures
    quintiles_df = compute_quintiles_df(psid_df, meas)

    # Truth mean series for this measure (for relative-to-mean plots)
    truth_mean = truth_df[(truth_df.var.==meas).&(truth_df.kind.==:mean), :]
    sort!(truth_mean, :time)
    mean_map = Dict(truth_mean.time .=> Float64.(truth_mean.value))

    truth_meas = truth_df[(truth_df.var.==meas).&(truth_df.kind.==:binmean).&(truth_df.n.==5), :]
    if nrow(truth_meas) == 0
        @warn "No truth binmeans found for measure" meas
        continue
    end

    # Export m*5 plots: each plot is one quintile moment over time
    for qq in 1:5
        truth_q = truth_meas[truth_meas.bin.==qq, :]
        sort!(truth_q, :time)

        y_sim_levels = Float64.(truth_q[!, "value"])
        y_dat_levels = Float64.(quintiles_df[!, Symbol("m$(qq)")])

        if relative_to_mean
            # Simulated: divide by truth mean at same date
            sim_den = [get(mean_map, t, NaN) for t in truth_q.time]
            y_sim_levels = y_sim_levels ./ sim_den

            # Data: divide by its own weighted mean m0
            dat_den = Float64.(quintiles_df.m0)
            y_dat_levels = y_dat_levels ./ dat_den
        end

        y_sim = apply_transform(y_sim_levels, transform_kind)
        y_dat = apply_transform(y_dat_levels, transform_kind)

        p = Plots.plot(
            truth_q.time,
            y_sim;
            title="$(String(meas)) quintile $(qq) ($(transform_tag), $(relative_tag))",
            label="Simulated",
            legend=:best,
        )
        Plots.scatter!(
            p,
            quintiles_df.time,
            y_dat;
            label="Data",
            markersize=3,
        )

        outpath = joinpath(plot_dir, "$(String(meas))_q$(qq)_$(transform_tag)_$(relative_tag).png")
        savefig(p, outpath)
    end
end

# Plots.plot(
#     truth_df1.time,
#     truth_df1[!, String(meas)*"_per_hh"];
#     label="Simulated Mean $(String(meas))",
# )
# Plots.scatter!(QuarterlyDate.(quintiles_df.time), quintiles_df.m0; label="Data")

# cond = [x for x in QuarterlyDate.(quintiles_df.time) if x in truth_df1.time]
# cond = findall(x -> x in QuarterlyDate.(quintiles_df.time), truth_df1.time)
# ratio = quintiles_df[!, Symbol("m0")] ./ truth_df1[cond, String(meas)*"_per_hh"]

# # Generate 10 directories/folders that are copies of /Users/lc/Dropbox/Distributional_Dynamics/7_Results/consum_and_income_and_wealth HANK
# base = "/Users/lc/Dropbox/Distributional_Dynamics/7_Results/consum_and_income_and_wealth HANK"

# for i = 1:10
#     dst = base * " $(i)"
#     if ispath(dst)
#         @info "Skipping (already exists)" dst
#         continue
#     end
#     cp(base, dst)
#     @info "Created copy" dst
# end
