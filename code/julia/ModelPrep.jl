function X13_seasonality_adjustment!(df_to_des, periods, source)
    # Find rows with at least one observation 
    condition_axes = findall(row -> !all(isnan, row), eachrow(df_to_des))

    qd_vec = [QuarterlyDate(y, q)
              for y in sort(collect(keys(periods)))     # sort so the vector is chronological
              for q in periods[y]]


    yr = minimum(year.(qd_vec))
    qr = periods[yr][1]  # first quarter of the first year

    R"""
    library(x12)                       # x13binary must be installed
    """

    # Loop over these rows
    for i in condition_axes
        indexing = tuple(i, :)

        df_row = df_to_des[indexing...]

        # Find indices of NaNs across columns
        nan_ids = findall(isnan, df_row)

        # Since x12 cannot handle NaNs, we fill them with linear interpolation. We do not use these values anyway, so, its completely fine.
        df_row = fill_between_mean!(df_row)

        # TODO: does not work on my machine due to the R installed on my machine being rosetta based i think
        R"""
        an.error.occured <- FALSE
        tryCatch({
            # 1  Build the ts object (this prints nothing, so no need to silence)
            ts_obj <- ts(
                $(df_row),
                frequency = 4,
                start = c($yr, $qr)
            )

            # 2  Run X‑13ARIMA/SEATS quietly
            invisible(
            capture.output(
                adjusted <- x12(ts_obj)
            )
            )

            # 3  Extract the X‑11 adjusted component
            d <- adjusted@d11
            cat("success row", $i, "\n")

        }, error = function(e) {
            an.error.occured <<- TRUE
            cat("error row", $i, ":", e$message, "\n")
            d <- $(df_row)                 # fallback: keep original series
        })

        """

        @rget d
        # Plots.plot(d, title="X-11 adjusted series for $i", xlabel="Time", ylabel="Value")
        # Plots.plot!(df_row, label="Original series", linestyle=:dash)
        # Plots.savefig("x11_$(source)_$i.pdf")
        df_to_des[i, :] = d
        df_to_des[i, nan_ids] .= NaN
    end

    return df_to_des
end

"""
    fill_between_mean!(y)

In-place linear interpolation of **NaN** gaps in a numeric vector `y`.

* Each contiguous run of `NaN`s is replaced by values on the straight line
  joining the last observed point **before** the gap and the first observed
  point **after** it.
* If the gap touches the start (or end) of the series, that block is
  filled by carrying forward (or backward) the nearest observed value.
"""
function fill_between_mean!(y)
    n = length(y)
    i = 1
    while i ≤ n
        if isnan(y[i])
            # ── start of a NaN block ───────────────────────────────────────────
            start_idx = i
            while i ≤ n && isnan(y[i])
                i += 1
            end
            end_idx = i - 1                     # last NaN position

            # neighbours
            prev_idx = start_idx - 1
            next_idx = i ≤ n ? i : nothing

            prev_val = (prev_idx ≥ 1) ? y[prev_idx] : NaN
            next_val = (next_idx !== nothing) ? y[next_idx] : NaN

            if !isnan(prev_val) && !isnan(next_val)
                # linear interpolation across the gap
                gap_len = end_idx - start_idx + 2      # segments = gap+1
                for (k, idx) in enumerate(start_idx:end_idx)
                    y[idx] = prev_val +
                             (next_val - prev_val) * k / gap_len
                end
            elseif !isnan(prev_val)
                # gap at the end → carry last value forward
                fill!(view(y, start_idx:end_idx), prev_val)
            elseif !isnan(next_val)
                # gap at the beginning → carry first value backward
                fill!(view(y, start_idx:end_idx), next_val)
            end
        else
            i += 1
        end
    end
    return y
end


# x12_object <- new("x12Single", ts = ts($(df[indexing...]), frequency=4, start=c($(time_dict["year"][i]), $(time_dict["quarter"][i]))))

# an.error.occured <- FALSE
# tryCatch( { adjusted_data <- x12(x12_object); print("success") }
#           , error = function(e) {an.error.occured <<- TRUE})
# print(an.error.occured)
# print($i)

# # adjusted_data <- x12(x12_object)
# d <- adjusted_data@x12Output@d11


function find_meas_indices(measures, meas, grid)
    D = length(measures)
    diff = grid + (D - 1) * (grid - 1)

    cop_n = grid^D - diff
    meas_id = findfirst(x -> x == meas, measures)
    pcf_id = cop_n+(meas_id-1)*grid+1:cop_n+meas_id*grid

    return [pcf_id...]
end


function remove_seasonality_from_quarterly_data!(dfs, names, time_dict)
    sim_data = filter(x -> occursin("SimData", x), names)
    datasets_that_need_adjustment = ["SIPP1", "SIPP2", "SIPP3", sim_data..., "CEX"]

    for j in datasets_that_need_adjustment
        try
            df_id = findall(x -> x == j, names)[1]
            X13_seasonality_adjustment!(dfs[df_id], time_dict[df_id], j)
        catch ee
            println(ee)
            println("No dataset $j found.")
        end
        # sipp2_id = findall(x -> x == "SIPP2", names)[1]
    end
end


function estimation_prep(obs_data::ObservedData, model_options::ModelOptions)

    # Unpack data + necessary options 
    @unpack files, agg_data, df_vec, gdp_series = obs_data
    @unpack estimator, number_of_dfs, measures, lags, freq, agg_freq, case, plot_proof, pca_perspective, rm_seasonality, equivalized, data_cutoffs, data_to_mute, tag, agg_lags = model_options

    # Preliminairies
    init_path = DATA_PROCESSING
    dimension = length(measures)

    @info("Extracting observations Tⱼ of the joint distributions, in alphabetical order of datasets' name.")
    @info("Estimation selected: $tag")

    # Collects data for all available years 
    dfs, dfs_unaltered, year_vec, time_dict, freq_type = data_constructor(obs_data, model_options) # TODO: correct freq_type s.t. it makes sense 

    # Remove seasonality of specific datasets, only if OS is not mac
    @info("Removing seasonality from quarterly data.")
    remove_seasonality_from_quarterly_data!(dfs, df_vec.df_names, time_dict)

    # Using time parameters from data and from aggregates, we find an agreeable time frame for the estimation
    time_p = define_timeframe(agg_data, year_vec, freq_type, time_dict, data_cutoffs) # TODO: time_dict needs to be filtered as well  

    cutoff_bounds = align_data_with_timeframe!(dfs, year_vec, time_p, data_cutoffs)

    # Define data intervals # TODO: for the future, run this after data_constructor so that we construct intervals for all data and then subset to the timeframe we want vs. estimating separate intervals for different time frames 
    confidence_intervals, Σ̂⁻¹² = define_data_intervals(df_vec, model_options, init_path, time_p, obs_data)
    println("Intervals defined.")

    # Perform proof of concept reconstruction
    # First, select a df where all measures are observed 

    # Pick the key from confidence intervals corresponding to this id 
    if plot_proof && tag != " Γ all"
        id, OD = df_selector(dfs, df_vec, measures)
        source_of_id = df_vec.df_names[id]

        perform_proof_of_concept_reconstruction(OD, source_of_id, year_vec[id], gdp_series, time_p, freq_type[id], time_dict[id], model_options)
    end


    # Keep functional data for plotting against estimates 
    func_dict = store_functional_data(files, dfs, gdp_series, model_options, time_p)
    confidence_intervals = store_confidence_intervals(files, confidence_intervals, freq, time_p)
    data_sources = sort(collect(keys(files)))

    func_struct = FunctionalData(func_dict, confidence_intervals, year_vec, data_sources)

    @info("Performing data transformations.")

    # Detrending 
    βs = zeros(size(dfs[1], 1), 3, length(dfs))
    fill!(βs, NaN)
    trend = Vector{Matrix{Float64}}(undef, length(dfs))

    for j in eachindex(dfs)
        βs[:, :, j], dfs[j], trend[j] = perform_detrending(dfs[j], time_p, year_vec[j], freq, freq_type[j], time_dict[j])
    end

    # Standardization
    pooled_data, means, stds = perform_standardization(dfs, estimator, dimension, measures)

    # For the perturbation of the projection
    additional_data_blocks = false
    if occursin("PP", tag)
        # Dataset is after the PP
        dataset_names = split(tag[5:end], " ") # there could be multiple datasets
        common_names = intersect(dataset_names, df_vec.df_names)
        println(common_names)
        ids_of_df = Vector{Int8}(undef, length(common_names))

        for i in eachindex(ids_of_df)
            ids_of_df[i] = findall(x -> x == common_names[i], df_vec.df_names)[1]
        end

        # Find the dataset in the pool
        additional_data_blocks = find_data_blocks(dfs, pooled_data, ids_of_df)
    end

    # Dimensionality reduction 
    proj, pcs, _, n_less_than_one = perform_pca(pooled_data, measures, :functional_data, tag; additional_data_blocks)

    # Equation Objects
    MV = Vector{Matrix{Float64}}(undef, number_of_dfs)         # MV = measurement vector.

    # Fill Equation Objects 
    return func_struct, time_p, set_measurements(MV, pcs, means, stds, agg_data, pooled_data, proj, files, time_p, freq_type, time_dict, model_options, n_less_than_one, βs, trend, Σ̂⁻¹², data_to_mute, agg_lags, df_vec, gdp_series)
end


"""Compute number of measurement rows per dataset (block size) for `Σ̂⁻¹²`.

This matches the coefficient-vector length produced by `data_constructor`: copula block
plus `dimension` percentile-function blocks.
"""
function sigma_block_size(dimension::Integer, estimator)
    @unpack grid_pcf, grid_cop = estimator
    immutable = grid_cop + (dimension - 1) * (grid_cop - 1)
    cop_rows = grid_cop^dimension - immutable
    return cop_rows + dimension * grid_pcf
end


"""Extract the relevant per-dataset blocks from a (large) saved `Σ̂⁻¹²`.

`Σ̂⁻¹²_full` is assumed block-diagonal across datasets, in the order given by
`sigma_sources` (parsed from the filename).

This returns a new block-diagonal matrix containing only the blocks for the datasets
in `df_names` (typically 4 here: HANK a–d), in that order.
"""
function extract_sigma_blocks(Σ̂⁻¹²_full::AbstractMatrix,
    sigma_sources,
    df_names,
    number_of_dfs,
    dimension,
    estimator)

    isempty(sigma_sources) && error("sigma_sources is empty; cannot extract sigma blocks")

    K_from_est = sigma_block_size(dimension, estimator)
    K_from_file = size(Σ̂⁻¹²_full, 1) ÷ length(sigma_sources)
    if K_from_file * length(sigma_sources) != size(Σ̂⁻¹²_full, 1)
        error("Sigma size is not divisible by number of sources")
    end
    if K_from_est != K_from_file
        @warn "Sigma block size mismatch (estimator vs file); using file-implied block size" K_from_est K_from_file
    end
    K = K_from_file

    # Map df names to sigma source tokens
    function target_source(df_name::AbstractString)
        if occursin("HANK a", df_name)
            return "PSID"
        elseif occursin("HANK b", df_name)
            return "CPS"
        elseif occursin("HANK c", df_name)
            return "CEX"
        elseif occursin("HANK d", df_name)
            return "SCF"
        elseif occursin("HANK e", df_name)
            return "SIPP1"
        end
        for s in sigma_sources
            if occursin(s, df_name)
                return s
            end
        end
        return ""
    end

    n_keep = min(Int(number_of_dfs), length(df_names))
    out = zeros(eltype(Σ̂⁻¹²_full), K * n_keep, K * n_keep)

    for j in 1:n_keep
        tok = target_source(df_names[j])
        tok == "" && error("Could not map df name to sigma source token: $(df_names[j])")
        idx = findfirst(==(tok), sigma_sources)
        idx === nothing && error("Token $(tok) not found in sigma_sources")

        src_r = (idx-1)*K+1:idx*K
        dst_r = (j-1)*K+1:j*K
        @views out[dst_r, dst_r] .= Σ̂⁻¹²_full[src_r, src_r]
    end

    return out
end


function find_data_blocks(dfs, pool, ids_of_df)
    additional_data_blocks = Vector{Matrix{Float64}}(undef, length(ids_of_df))

    for (x, id) in enumerate(ids_of_df)
        i = 0
        for j in 1:(id-1)
            i += size(dfs[j], 2)
        end

        # The additional block 
        additional_data_block = pool[:, i+1:i+size(dfs[id], 2)]

        # Clean the data by replacing the NaNs with 0 
        additional_data_block[isnan.(additional_data_block)] .= 0

        additional_data_blocks[x] = additional_data_block
    end

    return additional_data_blocks
end


function df_selector(dfs, df_vec, measures)

    # Find the df with all the measures by checking whether one of the columns has zero NaNs
    # for (j, df) in enumerate(dfs)
    #     if any(col -> !any(isnan, col), eachcol(df))
    #         return j, deepcopy(df)
    #     end
    # end
    if measures == ["income", "wealth"]
        id_of_df = findall(x -> x == "SCF", df_vec.df_names)[1]
        return id_of_df, deepcopy(dfs[id_of_df])

    elseif measures == ["consum", "income", "wealth"]
        id_of_df = findall(x -> x == "PSID", df_vec.df_names)[1]
        return id_of_df, deepcopy(dfs[id_of_df])

    elseif measures == ["consum", "income"] || measures == ["consum", "income", "liquid"]
        id_of_df = findall(x -> occursin("CEX", x), df_vec.df_names)[1]
        return id_of_df, deepcopy(dfs[id_of_df])
    end
end

function moving_average(A::AbstractArray, m::Int)
    out = similar(A)
    R = CartesianIndices(A)
    Ifirst, Ilast = first(R), last(R)
    I1 = m ÷ 2 * oneunit(Ifirst)
    for I in R
        n, s = 0, zero(eltype(out))
        for J in max(Ifirst, I - I1):min(Ilast, I + I1)
            s += A[J]
            n += 1
        end
        out[I] = s / n
    end
    return out
end


function generate_correct_indices(year_vec, freq_type, freq, tmin, time_dict)
    ind_of_obs = year_vec .- tmin["year"]
    correct_indices = []

    # To do this properly, the index must go to the proper .. I need the actual dates like year and quarter
    if freq_type == "quarter"
        for (k, yr) in enumerate(unique(year_vec))
            for quarter in sort(time_dict[yr])
                push!(correct_indices, (freq * (yr - tmin["year"])) + (quarter - 1) + 1)
            end
        end
        correct_indices .= correct_indices .- (tmin["quarter"] .- 1) # Reason: think of an 1967 observation with 1967 being the start of the estimation 

    elseif freq_type == "year"
        correct_indices = [freq * ind_of_obs[i] for i in eachindex(ind_of_obs)]  # Relocates index based on frequency. For quarterly freq, I assume annual data corresponds to the 4th quarter.  
        correct_indices .= correct_indices .- freq .+ 1
    end
    return correct_indices
end


function perform_detrending(df, time_p, year_vec, freq, freq_type, time_dict; return_only_data=false)
    @unpack tot_years, tmin, tmax, tot_periods = time_p

    # Time index by dataset
    t = collect(1:1:tot_periods)
    TT = hcat(ones(tot_periods), t, 2 .* t .^ 2 ./ tot_periods)
    # TT               = hcat(ones(tot_periods), t)
    interp_series = zeros(tot_periods)
    trend = zeros(size(df, 1), tot_periods)

    # Create a matrix of NaN of size trend to store the coefficients
    fill!(trend, NaN)
    βs = zeros(size(df, 1), size(TT, 2))

    # First, interpolate series 
    correct_indices = generate_correct_indices(year_vec, freq_type, freq, tmin, time_dict)
    # correct_indices  = generate_correct_indices(year_vec[4], freq_type[4], freq, tmin, time_dict[4])

    # Subset to rows that are not completely NaN 
    condition = findall(row -> !all(isnan, row), eachrow(df))
    # Interpolate, store trend and keep residuals
    for i in condition
        # Linear Interpolation 
        mask = (!isnan).(df[i, :])
        interp_linear = linear_interpolation(correct_indices[mask], df[i, mask], extrapolation_bc=Line())
        interp_series .= convert.(Float64, interp_linear.(t))

        # Find the index of the first non-NaN column
        first_non_nan_index = findfirst(x -> all((!isnan).(x)), df[i, :])

        # Extract Trend 
        trend[i, correct_indices[first_non_nan_index]:end] .= HP(interp_series[correct_indices[first_non_nan_index]:end], 1600) # 
        df[i, mask] .= df[i, mask] .- trend[i, correct_indices[mask]]
        βs[i, :] .= 1 # HP(interp_series, 6000)


        # # Slice the matrix to skip the first NaN columns
        # trimmed_int = interp_series[correct_indices[first_non_nan_index]:end]
        # trimmed_TT  = TT[correct_indices[first_non_nan_index]:end, :]


        # βs[i, :]         = (trimmed_TT' * trimmed_TT) \ (trimmed_TT' * trimmed_int)
        # df[i, mask]      = df[i, mask] .- TT[correct_indices[mask], :] * βs[i, :]

        # ols              = lm(TT, interp_series)
        # df[i, mask]      = residuals(ols)[correct_indices[mask]]
        # βs[i, :]        .= coef(ols) 
        # if mod(i, 10) == 0
        #     Plots.plot(1:tot_periods, interp_series)
        #     Plots.plot!(correct_indices, df[i, :])
        #     Plots.plot!(correct_indices[mask], TT[correct_indices[mask], :] * βs[i, :])
        #     # Plots.plot!(1:tot_periods, trend[i, :])
        #     Plots.savefig("test$i" * ".png") 
        # end
    end
    if return_only_data
        return df
    else
        return βs, df, trend
    end
end


function identify_grid(estimator, grid_choice)
    est_tag = typeof(estimator) <: SeriesEstimator ? "_series" : typeof(estimator) <: KernelEstimator ? "_kernel" : typeof(estimator) <: HistogramEstimator ? "_hist" : "_notagyet"

    local grid_granularity
    if grid_choice == 10
        grid_granularity = "_deciles"
    elseif grid_choice == 5
        grid_granularity = "_quintiles"
    elseif grid_choice == 100
        grid_granularity = "_percentiles"
    elseif grid_choice == 20
        grid_granularity = "_ventiles"
    end

    return grid_granularity * est_tag
end


function define_data_intervals(df_vec, model_options, init_path, time_p, obs_data)
    """Estimate the confidence intervals of the data via a bootstrap procedure. It is a bit independent of the estimator, but for the 
    series estimator, it will be of a higher granularity. Because of this higher granularity and simplicity, I also use the beta kernel 
    to generate the copulas when a series estimator is specified."""

    @unpack measures, estimator, errors_process, freq, equivalized, bottom_coded, blind_to, data_cutoffs, tag = model_options
    @unpack time_dict, freq_type, year_vec, tmin, tmax, tot_periods = time_p
    @unpack gdp_series, agg_data = obs_data
    @unpack grid_pcf, grid_cop = estimator

    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid, integral_cop_grid = estimator
    end

    grid_choice = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf

    # For labels and tagging 
    end_year = data_cutoffs["end"] != "" ? data_cutoffs["end"][1:4] : "all"

    sources = df_vec.df_names  # e.g., SCF, PSID, etc.
    data_label = data_tag(sources) # e.g., CEX_CPS_CPS2_PSID_SCF

    grid_tag = identify_grid(estimator, grid_choice) # e.g., _deciles_series
    m_label = measures_folder(measures) # e.g., consum_income_wealth
    dimension = length(measures)

    gdp_series[!, "date"] = QuarterlyDate.(gdp_series[!, "time"])
    filter!(row -> row.date >= QuarterlyDate(tmin["year"], tmin["quarter"]), gdp_series)
    filter!(row -> row.date <= QuarterlyDate(tmax["year"], tmax["quarter"]), gdp_series)

    # Define series for the confidence intervals --- are based on the estimation results 
    objects = sort(["quantiles", "levels", "shares"])

    local series
    if grid_choice == 10 || grid_choice == 100 || grid_choice == 20
        series = ["bottom50", "next40", "top10"]
    elseif grid_choice == 5
        series = ["bottom40", "next40", "top20"]
    end

    # Loop over data sources and generate intervals
    confidence_intervals = Dict()
    Σ̂⁻¹²ⱼ = Vector{Matrix{Float64}}(undef, length(sources))
    max_draws = 999

    # File to save the sigma matrix. Only higher orders require a different computation
    ci_tag = tag == " higher order15" ? tag : ""
    sigma_file_name = init_path * "/noise_distributions/sigma_" * m_label * grid_tag * "_" * data_label * "_$end_year" * ci_tag * ".jld2"
    sigma_exists = isfile(sigma_file_name)

    # Threads.@threads 
    for j in axes(df_vec[1], 1) # TODO: threading possible when I don't compute the correlations. Why? it doesn't use RCall
        # Data name + container
        source = sources[j]
        ci_source = occursin("SCF", source) ? "SCF" : source
        confidence_intervals[source] = Dict()

        # File name for the confidence intervals and noise distributions 
        ci_file_name = init_path * "/confidence_intervals/ci_draws_" * m_label * grid_tag * "_" * ci_source * "_$end_year" * ci_tag * ".jld2"
        noise_file_name = init_path * "/noise_distributions/noise_draws_" * m_label * grid_tag * "_" * ci_source * "_$end_year" * ci_tag * ".jld2"

        println(ci_file_name)
        println(noise_file_name)

        # Check if the files exist
        ci_file_exists = isfile(ci_file_name)
        noise_file_exists = isfile(noise_file_name)

        if ci_file_exists && noise_file_exists
            # Read in confidence intervals
            sub_boot_dict = jldopen(ci_file_name, "r")["ci"]

            # only keep the lower and upper bounds
            confidence_intervals[source]["ci_u"], confidence_intervals[source]["ci_l"] = construct_confidence_intervals(sub_boot_dict, 0.025, 0.975, measures, year_vec[j], estimator)

            if !sigma_exists
                DCT_boot = jldopen(noise_file_name, "r")["noise"]
                Σ̂⁻¹²ⱼ[j] = transform_DCT_boot(DCT_boot, time_p, year_vec[j], freq, freq_type[j], time_dict[j], estimator, dimension, measures, source)
            end

        else
            df = df_vec[1][j]
            draws = occursin("CPS", source) ? Int(round(max_draws * 0.10, digits=0)) : max_draws # CPS is a large dataset --- tight intervals anyway 
            draws = occursin("SIPP", source) ? Int(round(draws * 0.50, digits=0)) : draws # SIPP is a quarterly dataset
            # draws = occursin("HANK", source) ? Int(round(draws * 0.50, digits=0)) : draws # HANK is a large simulated dataset
            data, _ = select_data(df, measures, equivalized, bottom_coded, blind_to, source)

            # HANK special case: empirical (data-based) noise is disabled (we don't add noise from data).
            # Keep `hank_noise_draws = nothing` so we don't require/loading noise files.
            local hank_noise_draws
            hank_noise_draws = nothing

            # Generate confidence intervals 
            # For HANK, also estimate Monte Carlo (finite-sample) noise by bootstrapping the simulated microdata
            # within each period, then combine it with empirical measurement noise.
            sub_boot_dict, DCT_boot = estimate_confidence_intervals!(data, objects, series, year_vec[j], time_dict[j], freq_type[j], estimator, measures, draws, gdp_series, source; noise_draws=hank_noise_draws, sd_scale=0.0)

            # Save the sub_boot_dict and DCT_boot using jld2 
            JLD2.save(ci_file_name, "ci", sub_boot_dict)
            JLD2.save(noise_file_name, "noise", DCT_boot)

            # only keep the lower and upper bounds
            confidence_intervals[source]["ci_u"], confidence_intervals[source]["ci_l"] = construct_confidence_intervals(sub_boot_dict, 0.025, 0.975, measures, year_vec[j], estimator)

            Σ̂⁻¹²ⱼ[j] = transform_DCT_boot(DCT_boot, time_p, year_vec[j], freq, freq_type[j], time_dict[j], estimator, dimension, measures, source)
        end
    end

    local Σ̂⁻¹²
    if sigma_exists
        Σ̂⁻¹² = jldopen(sigma_file_name, "r")["sigma"]
    else
        println("Sigma did not exist!")
        Σ̂⁻¹² = cat(Σ̂⁻¹²ⱼ...; dims=(1, 2))
        JLD2.save(sigma_file_name, "sigma", Σ̂⁻¹²)
    end

    return confidence_intervals, Σ̂⁻¹²
end


function compute_Kullback_Leibler_divergence(P, Q; compare_to_average=false)
    # Ensure P and Q are valid probability distributions
    # if nansum(P) ≈ 1.0 && nansum(Q) ≈ 1.0 || compare_to_average 
    # Avoid division by zero or log of zero by only considering elements where P(i) > 0 and Q(i) > 0
    mask = .!isnan.(P) .& .!isnan.(Q) .& (P .> 0) .& (Q .> 0)

    # Compute the KL divergence
    D_KL = dot(P[mask], log.(P[mask] ./ Q[mask]))
    return D_KL
    # else
    #     println(nansum(P))
    #     println(nansum(Q))
    #     error("P and Q must be valid probability distributions.")
    # end
end


"""Map a simulated HANK dataset name to its target empirical dataset name."""
function hank_target_dataset(df_name::AbstractString)
    if occursin("HANK a", df_name)
        return "PSID"
    elseif occursin("HANK b", df_name)
        return "CPS"
    elseif occursin("HANK c", df_name)
        return "CEX"
    elseif occursin("HANK d", df_name)
        return "SCF"
    elseif occursin("HANK e", df_name)
        return "SIPP1"
    else
        return ""
    end
end



function transform_DCT_boot(DCT_boot, time_p, year_vec, freq, freq_type, time_dict, estimator, dimension, measures, source)
    # STEP 1: detrend 
    n_bootstraps = axes(DCT_boot, 3)

    if occursin("CEX", source) || occursin("SIPP", source)
        # Threads.@threads 
        for b in n_bootstraps
            X13_seasonality_adjustment!(DCT_boot[:, :, b], time_dict, source)
            DCT_boot[:, :, b] = perform_detrending(DCT_boot[:, :, b], time_p, year_vec, freq, freq_type, time_dict; return_only_data=true)
            DCT_boot[:, :, b] = perform_standardization([DCT_boot[:, :, b]], estimator, dimension, measures)[1]
        end
    else
        Threads.@threads for b in n_bootstraps
            DCT_boot[:, :, b] = perform_detrending(DCT_boot[:, :, b], time_p, year_vec, freq, freq_type, time_dict; return_only_data=true)
            DCT_boot[:, :, b] = perform_standardization([DCT_boot[:, :, b]], estimator, dimension, measures)[1]
        end
    end

    # Demean across bootstraps 
    DCT_boot .= DCT_boot .- mean(DCT_boot, dims=3)

    # VCV estimation 
    sizes = size(DCT_boot)
    DCT_boot_reshaped = reshape(DCT_boot, (sizes[1], sizes[2] * sizes[3]))
    Σ = nancov(DCT_boot_reshaped, dims=2) # time-invariant VCV of the measurements 
    replace!(Σ, -0.0 => 0.0)

    # Find indices along diagonal that are 0.0 
    zero_inds = find_zero_diagonal_indices(Σ)

    for (v, w) in zero_inds
        Σ[v, w] = 1 #[v == 0 ? v + .1 : v for v in diag(Σ[k])]
    end

    Σ = sqrt(nearest_spd(inv(Σ)))

    for (v, w) in zero_inds
        Σ[v, w] = 0 #[v == 0 ? v + .1 : v for v in diag(Σ[k])]
    end

    return Σ
end

# correct_draws      = size(DCT_boot, 3)
# σ_means[source]    = zeros(n_objs, correct_draws) 
# σ_vars[source]     = zeros(n_objs, correct_draws)

# for b in n_bootstraps

#     # Create MV object
#     start = 1  # necessary to move from dataset to dataset  

#     MV_b   = order_measures(DCT_boot_copy[k][:,:,b], year_vec[j], tot_periods, freq, freq_type[j], tmin, time_dict[j], start) 
#     Σ      = estimate_measurement_VCV([MV_b], dimension, grid, errors_process, time_p)        

#     σ_means[k][:, b]        = Σ[1][1] # there is a frequency mismatch
#     σ_vars[k][:, b]         = Σ[2][1]

#     # For state error, not dataset or object specific 
#     # ω[b]              = Ω
#     # ω_means[b]        = mean(Ω)
#     # ω_vars[b]         = var(Ω)
# end
# Create diagonal block matrix, consisting of the different Σ̂⁻¹²[j] 


# For each SCF set ... 
# for set in [set1, set2]
# for k in sources
#     DCT_boot_copy = deepcopy(DCT_boot)
#     j             = findall(x -> x == k, sources)[1]
#     n_bootstraps  = axes(DCT_boot_copy[k], 3)



# # Plot the distribution from the bootstraps
# s_init_path = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
# path        = s_init_path * "/7_Results/$m_label/other_results/noise_distributions/measurement_noise/" 
# obj_tag     = ["copula", measures...]

# for k in keys(DCT_boot)
#     for i in axes(σ_means[k], 1)
#         try
#             p         = barhist(σ_means[k][i, :]) 
#             Plots.savefig(p, path * k * "_" * obj_tag[i] * "_mean_of_variances_" * ".pdf")
#         catch ee
#             println(ee)
#             println("Nothing to plot")
#         end
#     end

#     for i in axes(σ_vars[k], 1)
#         try 
#             p         = barhist(σ_vars[k][i, :]) 
#             Plots.savefig(p, path * k * "_" * obj_tag[i] * "_variance_of_variances_" * ".pdf")
#         catch ee
#             println(ee)
#             println("Nothing to plot")
#         end
#     end
# end

#     # # Process noise 
#     # path      = s_init_path * "/7_Results/$m_label/other_results/noise_distributions/process_noise/"
#     # obj_tag   = ["copula", measures...]

#     # p         = barhist(ω_means) 
#     # Plots.savefig(p, path * "joint_distribution" * "_mean_of_variances_" * ".pdf")

#     # p         = barhist(ω_vars) 
#     # Plots.savefig(p, path * "joint_distribution" * "_variance_of_variances_" * ".pdf")
# JLD2.save(noise_file_name, "noise", [σ_means, σ_vars, Σ̂⁻¹²])
# JLD2.save(ci_file_name, "ci", confidence_intervals)


function three_dim_standardization!(mat, grid, dimension)
    mat .= mat .- mean(mat, dims=3)
    n_cop_rows = grid^dimension - (grid + (dimension - 1) * (grid - 1))
    mat[1:n_cop_rows, :, :] .= mat[1:n_cop_rows, :, :] ./ nanstd(mat[1:n_cop_rows, :, :])

    # Divide the percentile function rows by the standard deviation of the percentile function rows, in sets of size grid 
    for i in 1:dimension
        mat[n_cop_rows+(i-1)*grid+1:n_cop_rows+i*grid, :, :] .= mat[n_cop_rows+(i-1)*grid+1:n_cop_rows+i*grid, :, :] ./ nanstd(mat[n_cop_rows+(i-1)*grid+1:n_cop_rows+i*grid, :, :])
    end
end


function find_zero_diagonal_indices(matrix)
    zero_indices = []
    for i in 1:min(size(matrix)...)
        if matrix[i, i] == 0.0
            push!(zero_indices, (i, i))
        end
    end
    return zero_indices
end


function find_completely_informed_periods(dfs, freq, time_p)
    @unpack year_vec, tmin, tmax = time_p
    # first, combine year vector and dfs 
    n_of_dfs = length(dfs)
    comb_dfs = Vector{Matrix{Float64}}(undef, n_of_dfs)
    comp_years = Vector{Vector{Float64}}(undef, n_of_dfs)

    hack = 0.01
    hack_size = 0.01
    for (j, df) in enumerate(dfs)
        years_of_df = year_vec[j]' .+ hack
        dated_df = vcat(years_of_df, df)  # hack is for dealing with duplicate years and tagging the dataframe 
        println("test")
        valid_columns = [all(!isnan.(col)) for col in eachcol(comb_dfs[j])]
        println(dated_df[:, valid_columns][1, :])
        println("test")
        println(dated_df[:, mapslices(col -> all((!isnan).(col)), comb_dfs[j], dims=1)[:]][1, :])
        comp_years[j] = dated_df[:, mapslices(col -> all((!isnan).(col)), comb_dfs[j], dims=1)[:]][1, :]  # Only complete data, actual years 
        # dated_df[1, :]    = round.((dated_df[1, :] .- tmin["year"]) * freq, digits=2)
        hack += hack_size
    end

    # Get years for which complete data is observed (for OLS estimation in prior)
    complete_years = Int.(floor.(vcat(comp_years...)))

    return complete_years
end

_nanfunc(f, A, ::Colon) = f(filter(!isnan, A))
_nanfunc(f, A, dims) = mapslices(a -> _nanfunc(f, a, :), A, dims=dims)
nanfunc(f, A; dims=:) = _nanfunc(f, A, dims)


function retrieve_cop_and_imm_part(estimator, dimension)
    @unpack grid_pcf, grid_cop = estimator

    cop_part = grid_cop^dimension
    imm_part = grid_cop + (dimension - 1) * (grid_cop - 1)

    return cop_part, imm_part
end

function perform_standardization(dfs, estimator, dimension, measures)
    """We demean the data at the grid point level, but divide by σₒ,
    which only varies at the object level (copulae, percentile functions).
    """
    @unpack grid_pcf, grid_cop = estimator

    # We store the means, one vector per dataset and stds, one per object 
    n_dfs = length(dfs)
    means = Vector{Vector{Float64}}(undef, n_dfs)

    # compute means for each observation and demean 
    for j in eachindex(dfs)
        means[j] = nanfunc(mean, dfs[j], dims=2)[:]
        dfs[j] .= dfs[j] .- means[j]
    end

    # Subset to columns where everything is observed
    pool = hcat(dfs...)

    # Compute standard deviations for each object 
    cop_part, imm_part = retrieve_cop_and_imm_part(estimator, dimension)
    # n_cops             = get_n_of_copulas(estimator, dimension)
    # rows_v             = get_rows_of_copulas(grid_cop, measures)

    # Copula stuff 
    n_objs = 1 + dimension
    cop_rows = cop_part - imm_part
    copula_data = pool[1:cop_rows, :]

    # Pcfs 
    pcfs_data = pool[cop_rows+1:end, :]
    pcf_partition = size(pcfs_data, 1) ÷ dimension
    pcfs = [pcfs_data[I, :] for I in Iterators.partition(axes(pcfs_data, 1), pcf_partition)]

    # standard deviations 
    stds = Vector{Float64}(undef, n_objs)

    stds[1] = std(filter(!isnan, copula_data))
    copula_data .= copula_data ./ stds[1]

    start = 2
    for o in start:n_objs
        stds[o] = std(filter(!isnan, pcfs[o-1]))
        pcfs[o-1] .= pcfs[o-1] ./ stds[o]
    end

    transformed_data = vcat(copula_data, vcat(pcfs...))

    return transformed_data, means, stds
end

function define_timeframe(agg_data, year_vec, freq_type, time_dict, data_cutoffs)
    """This defines the time bounds for the estimates."""

    # Dictionairies for the time bounds 
    min = Dict()
    max = Dict()
    years = unique(vcat(year_vec...))

    # Add time columns to aggregate data 
    try
        agg_data[!, "time"] = QuarterlyDate.(agg_data[!, "time"])
    catch
        agg_data[!, "time"] = QuarterlyDate.(agg_data[!, "year"], agg_data[!, "quarter"])
    end
    agg_data[!, "year"] = Dates.year.(agg_data[!, "time"])
    agg_data[!, "quarter"] = Dates.quarterofyear.(agg_data[!, "time"])

    # Find what the minimums, maximums are for the aggregates 
    agg_min_yr = minimum(agg_data[!, "year"])
    agg_max_yr = maximum(agg_data[!, "year"])

    # Find what the minimums, maximums are for the data 
    meas_min_yr = (minimum ∘ vcat)(year_vec...)
    meas_max_yr = (maximum ∘ vcat)(year_vec...)

    # Now, focusing on the minimum quarter  
    meas_min_q = []

    for j in eachindex(time_dict)
        if meas_min_yr in collect(keys(time_dict[j]))
            if freq_type[j] == "year"
                push!(meas_min_q, 4)   # in this setup, the quarter can be anything 

            elseif freq_type[j] == "quarter"
                push!(meas_min_q, minimum(time_dict[j][meas_min_yr]))
            end
        else
            continue
        end
    end

    # Find the minimum of meas_min_q and the index of the minimum
    meas_min_q = minimum(meas_min_q)
    source_min = findfirst(x -> x == meas_min_q, meas_min_q)

    # Maximum quarter 
    meas_max_q = []

    for j in eachindex(time_dict)
        if meas_max_yr in collect(keys(time_dict[j]))
            if freq_type[j] == "year"
                push!(meas_max_q, 4)

            elseif freq_type[j] == "quarter"
                push!(meas_max_q, maximum(time_dict[j][meas_max_yr]))

            end
        else
            continue
        end
    end

    # Find the maximum of meas_max_q and the index of the maximum
    meas_max_q = maximum(meas_max_q)
    source_max = findfirst(x -> x == meas_max_q, meas_max_q)

    # What is the first and last quarter that appears in the aggregate data ?
    year_data = filter(row -> row.year == agg_min_yr, agg_data)
    agg_min_qtr = minimum(year_data[:, "quarter"])

    year_data = filter(row -> row.year == agg_max_yr, agg_data)
    agg_max_qtr = maximum(year_data[:, "quarter"])

    # Print the current time bounds 
    @info("The aggregates run from $agg_min_yr, $agg_min_qtr to $agg_max_yr, $agg_max_qtr")
    @info("The data runs from $meas_min_yr, $meas_min_q to $meas_max_yr, $meas_max_q")

    # First case: if agg min year is equal to the measurement year, we check whether they also start in the same quarter. If so, we shift to the next data year. And from that year we take the quarter from the data as the minimum (we assume the agg data is balanced) 
    # If the measurement quarter is already larger than the agg quarter, then we can ensure that there are aggregate lags to use for the measurements. 
    if meas_min_yr == agg_min_yr
        if agg_min_qtr < meas_min_q  # agg earlier than data 
            min["year"] = meas_min_yr
            min["quarter"] = meas_min_q

        elseif agg_min_qtr >= meas_min_q # agg later than data
            years_over = sort(years[years.>=agg_min_yr+1])  # We take first year over, to give aggregates cushion. Otherwise, we can't regress on the lags 
            min["year"] = first(years_over)
            o_source_min = find_source(time_dict, min["year"])

            # Finding minimum 
            if freq_type[o_source_min] == "year"
                min["quarter"] = 4

            elseif freq_type[o_source_min] == "quarter"
                # Find the source that contains this new min["year"], using time_dict
                min["quarter"] = minimum(time_dict[o_source_min][min["year"]]) #TODO: a bit lucky that it works, since we designate as min. an entirely new year, which may not even exist 
            end
        end

        # Second case: if the agg min year is larger than the measurement year, then we take first the year that is equal or larger to the agg min. 
        # If larger, we take all info from the measurements, since the aggregates can back it up 
        # If equal, we need to check that the quarter of the aggs is less than the measurement quarter. Else, we shift it by the next data year and take the minimum quarter from the data.
    elseif meas_min_yr < agg_min_yr
        years_over = sort(years[years.>=agg_min_yr])  # there is a lag here, so no need to do +1 ie, go to next year 
        meas_min_y = first(years_over)
        o_source_min = find_source(time_dict, meas_min_y)

        local meas_min_q2
        if freq_type[o_source_min] == "year"
            meas_min_q2 = 4
        elseif freq_type[o_source_min] == "quarter"
            meas_min_q2 = minimum(time_dict[o_source_min][meas_min_y])
        end

        if meas_min_y > agg_min_yr
            min["year"] = meas_min_y
            min["quarter"] = meas_min_q2
        elseif meas_min_y == agg_min_yr
            # If they are the same year, we need to make sure there's cushion 
            if agg_min_qtr < meas_min_q2
                min["year"] = meas_min_y
                min["quarter"] = meas_min_q2
            elseif agg_min_qtr >= meas_min_q2
                years_over = sort(years[years.>=agg_min_yr+1])
                min["year"] = first(years_over)
                o_source_min = find_source(time_dict, min["year"])

                # Finding minimum
                if freq_type[o_source_min] == "year"
                    min["quarter"] = 4
                elseif freq_type[o_source_min] == "quarter"
                    min["quarter"] = minimum(time_dict[o_source_min][min["year"]])
                end
            end
        end
    elseif meas_min_yr > agg_min_yr
        min["year"] = meas_min_yr
        min["quarter"] = meas_min_q
    end

    # Defining the max bounds. At the moment, we just take the max of the data.
    max["year"] = agg_max_yr
    max["quarter"] = agg_max_qtr

    # if meas_max_yr < agg_max_yr  
    #     max["year"]    = meas_max_yr  # if 'agg_max_yr', this implies that there will be some periods where we observe agg data but no measures i.e., out of sample performance 
    #     max["quarter"] = meas_max_q 

    # elseif meas_max_yr > agg_max_yr  
    #     max["year"]    = agg_max_yr  
    #     max["quarter"] = agg_max_qtr 

    # elseif meas_max_yr == agg_max_yr      
    #     if agg_max_qtr == meas_max_q
    #         max["year"]    = meas_max_yr  
    #         max["quarter"] = meas_max_q 
    #     elseif agg_max_qtr < meas_max_q
    #         max["year"]    = agg_max_yr  
    #         max["quarter"] = agg_max_qtr
    #     elseif agg_max_qtr > meas_max_q
    #         max["year"]    = agg_max_yr  
    #         max["quarter"] = meas_max_q 
    #     end
    # end
    println(min, max)
    @info("Agreed upon estimation period: $(min["year"]), $(min["quarter"]) to $(max["year"]), $(max["quarter"])")

    # Filter, in-place, 'year_vec'
    year_vec = [filter!(x -> x .>= min["year"], year_vec[j]) for j in eachindex(year_vec)]  # removing all years from year_vec that are outside the estimation 

    # Filter in place 'time_dict'
    # for j in eachindex(year_vec)
    #     keys_to_remove = [key for key in keys(time_dict[j]) if !(key in year_vec[j])]
    #     for key in keys_to_remove
    #         pop!(time_dict[j], key)
    #     end
    # end

    # # Extract from the min year the quarter list and ensure to remove all the quarters below the minimum 
    # for j in eachindex(year_vec)
    #     if min["year"] in keys(time_dict[j]) && freq_type[j] == "quarter"
    #         filter!(x->x .>= min["quarter"], time_dict[j][min["year"]])
    #     else
    #         @info("No quarters to subset to.")
    #     end
    # end

    # Correct the min and max for user imposed dates 
    # if data_cutoffs["end"] != ""
    #     end_year = parse(Int, data_cutoffs["end"][1:4])
    #     end_qtr  = parse(Int, data_cutoffs["end"][end])

    #     if max["year"] >= end_year
    #         max["year"] = end_year
    #         if max["quarter"] >= end_qtr
    #             max["quarter"] = end_qtr
    #         end
    #     end
    # end

    # Periods in the estimation  
    n_dates = QuarterlyDate(min["year"], min["quarter"]):Quarter(1):QuarterlyDate(max["year"], max["quarter"]) #TODO: for the moment, no extrapolation on the lower end 
    tot_years = max["year"] - min["year"] + 1
    tot_periods = length(n_dates)

    # n_dates     = QuarterlyDate(min["year"], min["quarter"]) : Quarter(1) : QuarterlyDate(max["year"], max["quarter"]) 
    # tot_years   = max["year"] - min["year"] + 1
    # tot_periods = length(n_dates)

    return TimeParams(year_vec, min, max, tot_years, tot_periods, time_dict, freq_type)
end


function find_source(time_dict, year)
    dfs_with_year = []
    for (j, dict) in time_dict
        if year in collect(keys(dict))
            push!(dfs_with_year, j)
        end
    end
    # Among the dfs, find the one with the lowest quarter
    quarters_of_df = [minimum(time_dict[j][year]) for j in dfs_with_year]
    source = dfs_with_year[findfirst(x -> x == minimum(quarters_of_df), quarters_of_df)]

    return source
end


function align_data_with_timeframe!(dfs, year_vec, time_p, data_cutoffs)
    @unpack time_dict, tmin, tmax = time_p

    # Find year of data_cutoffs 
    beg_yr = tmin["year"]
    beg_qtr = tmin["quarter"]
    end_year = data_cutoffs["end"] != "" ? parse(Int, data_cutoffs["end"][1:4]) : tmax["year"]
    end_qtr = data_cutoffs["end"] != "" ? parse(Int, data_cutoffs["end"][end]) : tmax["quarter"]

    n_inside_lb = Vector(undef, length(dfs))
    n_inside_ub = Vector(undef, length(dfs))
    cutoff_bounds = Vector{Vector{Int}}(undef, length(dfs))

    # This step removes observations from above, according to 'data_cutoffs'
    for j in eachindex(dfs)
        println(j)
        # Identify the upper bound of the data and the estimation upper bound wanted 
        ub_cutoff = [key for key in keys(time_dict[j]) if !(key <= end_year)] # keys that are above 2016 e.g.
        lb_cutoff = [key for key in keys(time_dict[j]) if !(key >= beg_yr)] # keys that are below 1980 e.g.

        # Count how many observations that is over the year
        total_obs_ub = 0
        total_obs_lb = 0

        for k in keys(time_dict[j])
            if k in ub_cutoff
                total_obs_ub += length(time_dict[j][k])
            elseif k in lb_cutoff
                total_obs_lb += length(time_dict[j][k])
            else
                continue
            end
        end

        # Count how many observations that are over within the year e.g., I have 4 quarters in 2016, but i want to subset to just Q1 and Q2
        if end_year in keys(time_dict[j])
            total_obs_ub += sum(time_dict[j][end_year] .> end_qtr)

            # Remove from time_dict the respective quarter that is over 
            time_dict[j][end_year] = filter(x -> x .<= end_qtr, time_dict[j][end_year])
            n_inside_ub[j] = length(time_dict[j][end_year])

            if n_inside_ub[j] == 0
                pop!(time_dict[j], end_year)
            end
        end


        if beg_yr in keys(time_dict[j])
            total_obs_lb += sum(time_dict[j][beg_yr] .< beg_qtr)

            # Remove from time_dict the respective quarter that is over 
            time_dict[j][beg_yr] = filter(x -> x .>= beg_qtr, time_dict[j][beg_yr])
            n_inside_lb[j] = length(time_dict[j][beg_yr])

            if n_inside_lb[j] == 0
                pop!(time_dict[j], beg_yr)
            end
        end

        for key in vcat(lb_cutoff, ub_cutoff)
            pop!(time_dict[j], key)
        end

        # Subset the dfs from above, removing all observations above upper bound 
        dfs[j] = dfs[j][:, 1+total_obs_lb:end-total_obs_ub]
        cutoff_bounds[j] = [total_obs_lb, total_obs_ub]
    end

    # Correct year_vec as well 
    year_vec = [filter!(x -> x .<= end_year, year_vec[j]) for j in eachindex(year_vec)]

    for j in eachindex(year_vec)
        if end_year in year_vec[j]
            # Remove all years associated with end_year
            filter!(x -> x != end_year, year_vec[j])

            # Append the correct amount of "year" = quarters to the year_vec
            append!(year_vec[j], repeat([end_year], n_inside_ub[j]))
        end
    end

    return cutoff_bounds
end

# for j in eachindex(year_vec)
#     keys_to_remove = [key for key in keys(time_dict[j]) if !(key in year_vec[j])]
#     for key in keys_to_remove
#         pop!(time_dict[j], key)
#     end
# end

# # Extract from the min year the quarter list and ensure to remove all the quarters below the minimum 
# for j in eachindex(year_vec)
#     if min["year"] in keys(time_dict[j]) && freq_type[j] == "quarter"
#         filter!(x->x .>= min["quarter"], time_dict[j][min["year"]])
#     else
#         @info("No quarters to subset to.")
#     end
# end


# TODO: ideally, the user should be able to set the dimension and call only 1 function, estimate_topologies 
function collect_topologies(orig_data::XLSX.XLSXFile, measures, grid, type)
    years = create_list_of_years(orig_data)  # for the 1 dimensional case, this will be an issue in the PSID 
    mat = length(measures) == 1 ? get_shares_levels(orig_data, measures, years, grid, type) : estimate_topologies(orig_data, measures, years, grid, type)
    return mat, years
end


function set_measurements(MV, pcs, means, stds, agg_data, pool, proj, files, time_p, freq_type, time_dict, model_options, n_less_than_one, βs, trend, Σ̂⁻¹², data_to_mute, agg_lags, df_vec, gdp_series)
    @unpack freq, agg_freq, number_of_dfs, lags, case, measures, pre_multiply, pca_perspective, tag, plot_proof, best_aggs, estimator = model_options
    @unpack tot_years, year_vec, tmin, tmax, tot_periods = time_p
    @unpack integral_pcf_grid, integral_cop_grid, grid_pcf, grid_cop = estimator


    # Order Measures
    start = 1  # necessary to move from dataset to dataset  
    for (j, _) in enumerate(sort!(collect(keys(files))))
        MV[j] = order_measures(pool, year_vec[j], tot_periods, freq, freq_type[j], tmin, time_dict[j], start)
        start += length(year_vec[j])  #TODO: very important parameter 
    end

    # Collect aggregate data
    u_proj, controls, agg_count = select_aggregates(agg_data, measures, tot_periods, tmin, tmax, agg_lags, tag, best_aggs)

    # Create law of motion from states to measurements 
    n_dfs = length(MV)

    # mute observations that correspond to 'data_to_mute'

    # Find the time index that corresponds to 'data_to_mute'
    dts = QuarterlyDate(tmin["year"], tmin["quarter"]):Quarter(1):QuarterlyDate(tmax["year"], tmax["quarter"])
    if typeof(data_to_mute) == Dict{String,QuarterlyDate}
        tag_check = [occursin(m, tag) for m in measures]
        if sum(tag_check) != 0
            # Case of removing wealth/income/other observations during the housing cycle

            # First, find wealth indices
            meas_to_remove = tag[end-5:end]
            meas_to_remove_loc = findfirst(x -> x == meas_to_remove, measures)
            cop_end = grid_cop^length(measures) - (grid_cop + (length(measures) - 1) * (grid_cop - 1))
            meas_to_remove_ids = (cop_end+grid_pcf*(meas_to_remove_loc-1)+1):(cop_end+grid_pcf*meas_to_remove_loc)

            if data_to_mute["begin"] != ""
                begin_id = findall(x -> x == data_to_mute["begin"], dts)
                end_id = findall(x -> x == data_to_mute["end"], dts)

                for j in eachindex(MV)
                    if isempty(end_id)
                        MV[j][meas_to_remove_ids, begin_id[1]:end] .= NaN
                    else
                        MV[j][meas_to_remove_ids, begin_id[1]:end_id[1]] .= NaN
                    end
                end
            end
        else
            # General case
            if data_to_mute["begin"] != ""
                begin_id = findall(x -> x == data_to_mute["begin"], dts)
                end_id = findall(x -> x == data_to_mute["end"], dts)

                for j in eachindex(MV)
                    if isempty(end_id)
                        MV[j][:, begin_id[1]:end] .= NaN
                    else
                        MV[j][:, begin_id[1]:end_id[1]] .= NaN
                    end
                end
            end
        end
    elseif typeof(data_to_mute) == Dict{String,Vector{QuarterlyDate}}
        # For specific datasets e.g., muting CEX
        tag_check = [occursin(m, tag) for m in measures]
        if sum(tag_check) != 0
            # Case of removing wealth/income/other observations during the housing cycle

            # First, find wealth indices
            meas_to_remove = tag[end-5:end]
            meas_to_remove_loc = findfirst(x -> x == meas_to_remove, measures)
            cop_end = grid_cop^length(measures) - (grid_cop + (length(measures) - 1) * (grid_cop - 1))
            meas_to_remove_ids = (cop_end+grid_pcf*(meas_to_remove_loc-1)+1):(cop_end+grid_pcf*meas_to_remove_loc)

            if data_to_mute["begin"] != ""
                begin_id = findall(x -> x == data_to_mute["begin"], dts)
                end_id = findall(x -> x == data_to_mute["end"], dts)

                for j in eachindex(MV)
                    if isempty(end_id)
                        MV[j][meas_to_remove_ids, begin_id[1]:end] .= NaN
                    else
                        MV[j][meas_to_remove_ids, begin_id[1]:end_id[1]] .= NaN
                    end
                end
            end
        else
            for df_name in collect(keys(data_to_mute))
                j = findfirst(x -> x == df_name, df_vec.df_names)
                for dt in data_to_mute[df_name]
                    date_id = findall(x -> x == dt, dts)[1]
                    MV[j][:, date_id] .= NaN
                end
            end
        end
    end

    if occursin(" Γ all", tag)
        if tag == " Γ all"
            pcs, proj = run_TW_algorithm(MV, 0.8, "second tallest")
        elseif tag == " Γ all 85"
            pcs, proj = run_TW_algorithm(MV, 0.85, "second tallest")
        end
        println("The size of the projection matrix for TW is $(size(proj))")

        if plot_proof
            perform_proof_of_concept_Γ_comparison(MV, df_vec, gdp_series, time_p, freq_type, time_dict, model_options, means, stds, trend)
        end
    end
    # ───────────────────────────────────────────────
    # 0.  Dimensions
    # ----------------------------------------------
    factor_count = size(proj, 2)                                 # # factors
    q = size(controls, 1)                             # # controls

    # First, for each MV, I want to define a G
    Gⱼ = [zeros(size(proj, 1), factor_count * 4) for _ in 1:n_dfs]  # each Gⱼ is (meas × state)

    # Now, depending on the df name and measure, we multiply it by 1/4, 1/2, or 1
    dimension = length(measures)  # number of objects, e.g., 1 for copula, 2 for copula and pcf, etc.
    cop_part, imm_part = retrieve_cop_and_imm_part(estimator, dimension)
    cop_rows = cop_part - imm_part
    id_dict = Dict(measures[i] => (cop_rows+1)+grid_pcf*(i-1):(cop_rows)+i*grid_pcf for i in eachindex(measures))  # mapping measures to their ids
    id_dict["copula"] = 1:cop_rows  # copula rows are always the first rows
    cop_proj = proj[id_dict["copula"], :]

    for (j, df_name) in enumerate(df_vec.df_names)
        if df_name == "SCF" || df_name == "PSID" || occursin("CPS", df_name)
            # For these, income, consumption, wealth  is annual --- we have a timestamp for wealth
            for m in measures
                id_tup = tuple(id_dict[m], :)
                sub_proj = proj[id_tup...]
                if m == "consum" || m == "income"
                    Gⱼ[j][id_tup...] .= 1 / 4 .* hcat(sub_proj, [sub_proj for _ in 1:3]...)  # Each row has the projection, zeros for the lags, and zeros for the controls
                elseif m == "wealth"
                    Gⱼ[j][id_tup...] .= hcat(sub_proj, [zeros(size(sub_proj)) for _ in 1:3]...)  # Each row has the projection, zeros for the lags, and zeros for the controls
                end
            end
            # Dealing with copula
            Gⱼ[j][1:cop_rows, :] .= 1 / 4 .* hcat(cop_proj, [cop_proj for _ in 1:3]...)  # Copula rows are scaled by 1/4

        elseif occursin("SIPP", df_name) || occursin("CEX", df_name) || occursin("HANK", df_name)
            # These datasets are quarterly/high-frequency
            for m in measures
                id_tup = tuple(id_dict[m], :)
                sub_proj = proj[id_tup...]
                Gⱼ[j][id_dict[m], :] .= hcat(sub_proj, [zeros(size(sub_proj)) for _ in 1:3]...)  # Each row has the projection, zeros for the lags, and zeros for the controls
            end
            Gⱼ[j][1:cop_rows, :] .= hcat(cop_proj, [zeros(size(cop_proj)) for _ in 1:3]...)  # Copula rows are scaled by 1/4
        end
    end

    # Defining components of the measurement equation
    y = vcat(MV...)
    N_dist = size(y, 1)                     # distributional measurements
    y = vcat(y, controls) # TODO: timing?
    meas_count = N_dist + q

    # ───────────────────────────────────────────────
    # 1.  Stacked projection Γ̃
    # ----------------------------------------------
    # proj_dist = repeat(proj, n_dfs)                   # N_dist × r   (distributional loadings)
    proj_dist = vcat(Gⱼ...)

    # ───────────────────────────────────────────────
    # 2.  Pre–allocate measurement matrices
    # ----------------------------------------------
    G = [zeros(meas_count, factor_count * 4 + q) for _ in 1:tot_periods]  # each G[t] is (meas × state)

    # ───────────────────────────────────────────────
    # 3.  Fill each G[t]
    # ----------------------------------------------
    for t in 1:tot_periods
        # 3a. Selector for distributional block
        mask_dist = (!isnan).(y[1:N_dist, t])
        H_dist = Diagonal(mask_dist)                 # N_dist × N_dist

        # 3b. Upper block: H_t * proj_dist     (loads on F_t, zeros on Y_t)
        G[t][1:N_dist, 1:factor_count*4] .= H_dist * proj_dist

        # 3c. Lower block: I_q on the Y-columns (columns factor_count+1 : factor_count+q)
        for i in 1:q
            G[t][N_dist+i, factor_count*4+i] = 1.0            # pick out Y_t element i exactly
        end
    end


    # stacked_proj = repeat(proj, n_dfs)
    # factor_count = size(proj, 2)
    # meas_count = size(y, 1)
    # G = [zeros(meas_count, factor_count) for _ in 1:tot_periods]   # The selection matrix * stacked projection 

    # # Measurement equation: Part of Term #1 - Selection matrix times the projection matrix 
    # H = Diagonal{Float64}(ones(meas_count))
    # for t in eachindex(G)
    #     condition = (!isnan).(y[:, t])[:]
    #     copyto!(H, Diagonal(condition))
    #     @views G[t] = H * stacked_proj
    # end

    if pre_multiply
        ŷ, Ĝ = pre_multiply_part!(y, G, Σ̂⁻¹², factor_count * 4)
        return StateSpaceModel(MV, ŷ, Ĝ, Gⱼ, controls, pcs, means, stds, proj, u_proj, factor_count, agg_count, n_less_than_one, βs, trend, Σ̂⁻¹²)
    else
        return StateSpaceModel(MV, y, G, Gⱼ, controls, pcs, means, stds, proj, u_proj, factor_count, agg_count, n_less_than_one, βs, trend, Σ̂⁻¹²)
    end
end


function perform_pca(pool, measures, type, tag; additional_data_blocks=false, best_aggs=false)
    local M, proj, pcs
    # if pca_perspective == "bayesian"
    #     M = MultivariateStats.fit(PPCA, data_matrix; maxiter=3000, method=:bayes)  # probabilistic PCA 
    #     # TODO: for some reason, this only works for :ml, will have to think about it 
    # else
    if type == :aggs
        # Create new data matrix with information on the lags 
        r_max = 30 # just some high number 
        MOdim = n_factors(pool', r_max) # pool' is a T x N matrix
        println("The maximum number of aggregate factors with other estimator is $MOdim")
        println(tag)

        M = MultivariateStats.fit(PCA, pool; pratio=0.95, method=:svd, mean=0)  # TODO: unfix this!
        pcs = MultivariateStats.transform(M, pool)
        Mdim = tag == " less AF" ? 3 : tag == " more AF" ? 15 : tag == " all AF" ? 30 : tag == " less DF and AF" ? 3 : occursin("HANK", tag) ? 9 : 11

        λ = sqrt.(principalvars(M))
        pcs_s = pcs ./ λ
        proj = projection(M) * diagm(λ)
        println("The size of the projection matrix for the aggs. is $(size(proj))")

        # Import jld2
        if best_aggs
            m_label = measures_folder(measures)
            init_path = BASE_PATH
            file_name = init_path * "/7_Results/$(m_label)$(tag)/other_results/R2_dict.jld2"
            best_aggs = jldopen(file_name, "r")["aggs"][1:Mdim]
            pcs_s = pcs_s[best_aggs, :]
            proj = proj[:, best_aggs]
        else
            pcs_s = pcs_s[1:Mdim, :]
            proj = proj[:, 1:Mdim]
        end

        return proj, pcs_s, M

    elseif type == :functional_data
        data_matrix = pool[:, mapslices(col -> all((!isnan).(col)), pool, dims=1)[:]]  # Only complete data

        # Only complete data
        pr = variation_selector(measures, tag)
        local M
        if tag == " 6 factors"
            M = MultivariateStats.fit(PCA, data_matrix; maxoutdim=6, method=:svd) # mean=0
        elseif tag == " 7 factors"
            M = MultivariateStats.fit(PCA, data_matrix; maxoutdim=7, method=:svd) # mean=0
        elseif occursin("HANK full", tag)
            M = MultivariateStats.fit(PCA, data_matrix; maxoutdim=8, method=:svd) # mean=0
        elseif occursin("HANK", tag)
            # M = MultivariateStats.fit(PCA, data_matrix; maxoutdim=8, method=:svd) # mean=0
            M = MultivariateStats.fit(PCA, data_matrix; pratio=pr, method=:svd) # mean=0
        else
            M = MultivariateStats.fit(PCA, data_matrix; pratio=pr, method=:svd) # mean=0
        end
        # Why do I do this? X = P * F, but if we want F = pcs to have unit variance, we must divide it by the square root of the eigenvalues of the covariance matrix of X.
        # But to maintain equality of the two sides, we must multiply P by the square root of the eigenvalues of the covariance matrix of X.
        pcs = MultivariateStats.transform(M, data_matrix) # predict() produces the same
        λ = sqrt.(principalvars(M))
        proj = projection(M) * diagm(λ)
        pcs_s = pcs ./ λ

        # Compute difference between observed and predicted values
        # rr = data_matrix - proj * pcs  # This is the residuals, which we will use to augment the projection matrix
        # M = MultivariateStats.fit(PCA, rr; method=:svd)

        if additional_data_blocks != false
            if tag == " PP CEX every 4 years"
                file_name = BASE_PATH * "/7_Results/consum_and_income_and_wealth PP CEX_annual/other_results/CEX_block.jld2"
                additional_data_blocks = jldopen(file_name, "r")["block"]
            end
            println(size(additional_data_blocks[1]))

            # Partial out variation from original data block 
            for g in eachindex(additional_data_blocks)
                A = (I - proj * inv(proj' * proj) * proj') * additional_data_blocks[g]

                AM = MultivariateStats.fit(PCA, A; pratio=pr, method=:svd) # mean=0
                Aλ = sqrt.(principalvars(AM))

                # Extract loadings corresponding to the first factor
                Aproj = projection(AM) * diagm(Aλ)

                # Augment the projection matrix with the loadings of the first factor 
                proj = hcat(proj, Aproj[:, 1])
                pcs_A = MultivariateStats.transform(AM, A)[1, 1:size(pcs_s, 2)]
                pcs_s = vcat(pcs_s, pcs_A')
            end
        end


        # proj  = projection(M) * sqrt.(Diagonal(pcs * pcs'))
        # pcs_s = pcs ./ sqrt.(diag(pcs * pcs')) 

        # Seeing how many factors are less than 1%
        variance_explained = principalvars(M) ./ tprincipalvar(M)
        n_less_than_one = length(findall(variance_explained .< 0.01))

        return proj, pcs_s, M, n_less_than_one
    end
    # end
end

function variation_selector(measures, tag)
    if tag == " additional factors"
        return 0.99
    elseif tag == " less factors" || tag == " less DF and AF"
        return 0.5
    else
        return 0.99 # TODO: new baseline
    end
end


function order_measures(pool, year_vec, tot_periods, freq, freq_type, tmin, time_dict, start)
    # Time index by dataset
    t = collect(1:1:tot_periods)
    # ind_of_obs       = year_vec .- tmin["year"] .+ 1                                        
    # correct_indices  = [freq * ind_of_obs[i] for i in eachindex(ind_of_obs)]  # Relocates index based on frequency. For quarterly freq, I assume annual data corresponds to the 4th quarter.  
    # correct_indices .= correct_indices .- freq .+ 1    

    correct_indices = []
    # Getting correct indices for interpolation
    ind_of_obs = year_vec .- tmin["year"]

    # To do this properly, the index must go to the proper .. I need the actual dates like year and quarter
    if freq_type == "quarter"
        # for (k, ind) in enumerate(unique(ind_of_obs))
        for (k, yr) in enumerate(unique(year_vec))
            for quarter in sort(time_dict[yr])
                push!(correct_indices, freq * (yr - tmin["year"]) + (quarter - 1) + 1)
            end
        end

        # I have to shift the indices by the amount of quarters in the first year 
        correct_indices .= correct_indices .- (tmin["quarter"] - 1)

    elseif freq_type == "year"
        correct_indices = [freq * ind_of_obs[i] for i in eachindex(ind_of_obs)]  # Relocates index based on frequency. For quarterly freq, I assume annual data corresponds to the 4th quarter.  
        correct_indices .= correct_indices .- freq .+ 1
    end

    # Fill matrix with density, percentile function measurements 
    meas_vec = zeros(size(pool, 1), tot_periods)  # Factors by Time 
    meas_vec[:, correct_indices] = copy(pool[:, start:(start+length(year_vec)-1)])   # the -1 is there because we'd otherwise exceed the length of the data 

    # Fill rest of indices with NaN
    ind_of_missing = setdiff(t, correct_indices)
    meas_vec[:, ind_of_missing] .= NaN

    return meas_vec
end


function select_aggregates(aggregates, measures, tot_periods, tmin, tmax, agg_lags, tag, best_aggs; until=false)

    aggregates[:, "date"] = QuarterlyDate.(aggregates[:, "time"])

    # # Since aggregates require a lag
    # q_correction = tmin["quarter"] > 1 ? tmin["quarter"] - 1 : 4
    # y_correction = tmin["quarter"] > 1 ? tmin["year"] : tmin["year"] - 1

    # Calculate the quarter correction
    q_correction = tmin["quarter"] > agg_lags ? tmin["quarter"] - agg_lags : 4 - (agg_lags - tmin["quarter"]) % 4 # eg,

    # Calculate the year correction
    y_correction = tmin["quarter"] > agg_lags ? tmin["year"] : trunc(tmin["year"] - ((agg_lags + 3 - 1) / 4) + 1)

    filter!(row -> row.date >= QuarterlyDate(y_correction, q_correction), aggregates)
    # a = filter(row -> row.date >= QuarterlyDate(y_correction, q_correction), aggregates)
    # tot_rows_lb = nrow(a)

    # dropped_rows_lb = tot_rows - tot_rows_lb

    # filter!(row -> row.date <= QuarterlyDate(tmax["year"], tmax["quarter"]), aggregates)
    local b
    if until != false
        b = filter(row -> row.date <= QuarterlyDate(until["year"], until["quarter"]), aggregates)
    else
        b = filter(row -> row.date <= QuarterlyDate(tmax["year"], tmax["quarter"]), aggregates)
    end

    dropped_rows_ub = nrow(aggregates) - nrow(b)

    data_only_aggs = select(aggregates, Not(["date", "time", "year", "quarter"]))

    # @assert size(data_only_aggs, 1) == tot_periods + 1

    # Fill matrix with density, percentile function measurements 
    data_only_aggs = copy(transpose(Matrix{Float64}(data_only_aggs)))  # removes date column, K x T

    if occursin("HANK", tag)
        # HANK: shocks are independent AR(1)s, lags are redundant — use contemporaneous only
        super_aggs = data_only_aggs[:, 1:end-agg_lags]
        standardize_aggs!(super_aggs)

        u_proj, agg_pcs, _ = perform_pca(super_aggs, measures, :aggs, tag; best_aggs=best_aggs)

        agg_pcs = agg_pcs[:, 1:end-dropped_rows_ub]
        agg_count = size(agg_pcs, 1)
        return u_proj, agg_pcs, agg_count
    else
        super_aggs = data_only_aggs[:, 1:end-agg_lags] # N by T

        for i in 1:agg_lags
            super_aggs = vcat(super_aggs, data_only_aggs[:, i+1:end-agg_lags+i])
        end

        standardize_aggs!(super_aggs)

        u_proj, agg_pcs, _ = perform_pca(super_aggs, measures, :aggs, tag; best_aggs=best_aggs)

        agg_pcs = agg_pcs[:, 1:end-dropped_rows_ub]
        agg_count = size(agg_pcs, 1)
        return u_proj, agg_pcs, agg_count
    end
end

# function impose_stationarity(aggregates)
# """For each series, performs hypothesis tests until L differences induces a stationary series."""
#     lag_vector  = zeros(Int64, size(aggregates, 1))
#     Δ₁_row      = size(aggregates, 1)
#     Δ₁_col      = size(aggregates, 2) - 1
#     new_agg_mat = zeros(Δ₁_row, Δ₁_col)

#     # See how many differences make it stationary  
#     for k in axes(aggregates, 1)
#         remove_time_seasonal_variation!(aggregates[k,:])
#         new_agg_mat[k, :] .= diff(aggregates[k, :])
#     end
#     #     p_val = 1
#     #     while p_val > .05
#     #         # Hypothesis test
#     #         p_val = aug_df_test(aggregates[k, :], lag_vector[k]+1) 
#     #         if p_val <= .05
#     #             nothing
#     #         else 
#     #             aggregates[k, :] .= diff(aggregates[k, :])
#     #         end 
#     #     end
#     # end
#     return new_agg_mat 
# end

# function remove_time_seasonal_variation!(series)
#     T       = length(series)
#     t_s     = collect(1:T)
#     season  = vcat(repeat([1:4], Int(floor(T/4)))...)
#     diff    = T - length(season)
#     diff_m  = mod(diff, 4)
#     cseason = diff < 0 ? season[1:end+diff] : vcat(season, collect(1:diff_m))
#     data    = DataFrame(s=series, t=t_s, t2=t_s.^2, season=cseason)
#     ols     = lm(@formula(s ~ t + t2 + season), data)  # constant is there automatically 
#     series .= series - GLM.predict(ols)  # keep residual variation 
# end

# function aug_df_test(series, l)
#     ADFResults = ADFTest(series, :none, l)
#     @unpack stat = ADFResults

#     # Check p-values. H₁ = stationarity  
#     p_val = adf_pv_aux(stat, :none)
#     return p_val 
# end

# function standardize_dct(dfs, grid, measures) 
#     """Standardize DCT coefficients."""
#     # Containers, params 
#     std_df_vec  = Vector{Matrix{Float64}}(undef, length(dfs))  # standardized df vector 
#     cop_grid    = grid^length(measures)
#     vec_pcf     = Vector(undef, length(measures))

#     # For each df ...
#     for (j, df) in enumerate(dfs) 
#         T        = size(df, 2)
#         data_cop = view(df, 1:cop_grid, :)
#         start    = cop_grid + 1
#         stretch  = copy(grid)
#         s_bar    = zeros(T)

#         # For each measure ...
#         for i in 1:length(measures)
#             vec_pcf[i] = view(df, start:cop_grid + stretch, :) 
#             start     += grid 
#             stretch   += grid 
#         end 

#         # For the copula 
#         for t in 1:T
#             s_bar[t] = norm(data_cop[:, t])
#             std_df_vec[j] = data_cop ./ mean(s_bar) 
#         end
#         std_df_vec[j] = data_cop ./ mean(s_bar) 

#         # For the percentile functions 
#         for pcf in vec_pcf
#                 for t in 1:T
#                     s_bar[t] = (sqrt ∘ sum)(pcf[:, t].^2)
#                 end
#                 std_df_vec[j] = vcat(std_df_vec[j], pcf ./ mean(s_bar))
#         end
#     end
#     return std_df_vec
# end

function varimax(A; gamma=1.0, minit=20, maxit=1000, reltol=1e-12)
    # Get the sizes of input matrix
    d, m = size(A)

    # If there is only one vector, then do nothing.
    if m == 1
        return A
    end

    if d == m && rank(A) == d
        return Matrix{Float64}(I, d, m)
    end

    # Warm up step: start with a good initial orthogonal matrix T by SVD and QR
    T = Matrix{Float64}(I, m, m)
    B = A * T
    L, _, M = svd(A' * (d * B .^ 3 - gamma * B * Diagonal(sum(B .^ 2, dims=1)[:])))
    T = L * M'
    if norm(T - Matrix{Float64}(I, m, m)) < reltol
        T, _ = qr(randn(m, m)).Q
        B = A * T
    end

    # Iteration step: get better T to maximize the objective (as described in Factor Analysis book)
    D = 0
    for k in 1:maxit
        Dold = D
        L, s, M = svd(A' * (d * B .^ 3 - gamma * B * Diagonal(sum(B .^ 2, dims=1)[:])))
        T = L * M'
        D = sum(s)
        B = A * T
        if (abs(D - Dold) / D < reltol) && k >= minit
            break
        end
    end

    # Adjust the sign of each rotated vector such that the maximum absolute value is positive.
    for i in 1:m
        if abs(maximum(B[:, i])) < abs(minimum(B[:, i]))
            B[:, i] .= -B[:, i]
        end
    end

    return B
end


function standardize_aggs!(agg_data)
    # Standardization -- here, each row is a time series so N by T
    agg_data .= (agg_data .- mean(agg_data, dims=2)) ./ std(agg_data, dims=2)
end


function generate_ecdf(rv, weight=false)
    # Estimate cdf and generate percentiles
    rv = convert.(Float64, rv)
    local rv_cdf
    if weight != false
        rv_cdf = ecdf(rv, weights=weight)
    end
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
    pmf[pmf.<0] .= 0  # there's some round off error garbage 
    d = DiscreteNonParametric(sort(obs), pmf)
    return d
end


function q_ignore_nan(v, p)
    w = v[isfinite.(v)]              # keeps finite; drops NaN/Inf
    return isempty(w) ? NaN : quantile(w, p)
end


function construct_confidence_intervals(sub_boot_dict, lb, ub, measures, years, estimator)

    # q_grid, cop_grid = pick_grid_for_confidence_intervals(estimator)

    @unpack grid_pcf, grid_cop = estimator
    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid, integral_cop_grid = estimator
    end
    grid_choice_pcf = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf
    grid_choice_cop = typeof(estimator) <: SeriesEstimator ? integral_cop_grid : grid_cop

    dim = length(measures)
    l_yrs = length(years)

    ci_u = Dict()
    ci_b = Dict()

    for meas in measures
        ci_u[meas] = Dict()
        ci_b[meas] = Dict()

        ci_u[meas]["quantiles"] = zeros(grid_choice_pcf, l_yrs)
        ci_b[meas]["quantiles"] = zeros(grid_choice_pcf, l_yrs)
        fill!(ci_u[meas]["quantiles"], NaN)
        fill!(ci_b[meas]["quantiles"], NaN)
    end


    for (m, meas) in enumerate(measures)
        for yr in 1:l_yrs
            for o in ["quantiles"]
                U = size(ci_b[meas][o], 1)
                for i in 1:U
                    vv = sub_boot_dict[meas][o][i, :, yr]
                    try
                        ci_b[meas][o][i, yr] = q_ignore_nan(vv, lb)
                        ci_u[meas][o][i, yr] = q_ignore_nan(vv, ub)
                    catch e
                        nothing
                    end
                end
            end
        end
    end

    ci_u["copula"] = zeros(grid_choice_cop^dim, l_yrs)
    ci_b["copula"] = zeros(grid_choice_cop^dim, l_yrs)

    for yr in 1:l_yrs

        # Copula stuff
        for i in axes(sub_boot_dict["copula"], 1)
            try
                vv = sub_boot_dict["copula"][i, :, yr]
                ci_b["copula"][i, yr] = q_ignore_nan(vv, lb)
                ci_u["copula"][i, yr] = q_ignore_nan(vv, ub)
            catch e
                ci_b["copula"][i, yr] = NaN
                ci_u["copula"][i, yr] = NaN
            end
        end
    end

    return ci_u, ci_b
end


function reduce_array(mat)

    # Initialize the new array
    sizes = size(mat)
    new_array = zeros(sizes[1], Int(sizes[2] / 5), sizes[3])

    # Iterate and average every 5 columns
    for i in axes(new_array, 2)
        start_col = (i - 1) * 5 + 1
        end_col = i * 5
        new_array[:, i, :] = mean(mat[:, start_col:end_col, :], dims=2)
    end
    return new_array
end

function reduce_vector(vec)

    # Initialize the new array
    new_vec = zeros(Int(length(vec) / 5))

    # Iterate and average every 5 columns
    for i in eachindex(new_vec)
        start_col = (i - 1) * 5 + 1
        end_col = i * 5
        new_vec[i] = mean(vec[start_col:end_col])
    end
    return new_vec
end


function get_param_vector(measures, kind_of_plots, label, data_cutoffs, tag)
    init_path = BASE_PATH
    m_label = measures_folder(measures)
    end_year = data_cutoffs["end"] != "" ? data_cutoffs["end"][1:4] : "all"

    local θ_cop
    if kind_of_plots == :mcmc
        cop_nd_sol_path = init_path * "/7_Results/$m_label" * "$tag" * "/from_mcmc/parameter_vectors/solution$label" * "_$end_year" * ".csv"
        θ_cop = vec(Matrix(CSV.read(cop_nd_sol_path, DataFrame, header=0)))
    else
        cop_nd_sol_path = init_path * "/7_Results/$m_label" * "$tag" * "/from_optimization/parameter_vectors/solution$label" * "_$end_year" * ".csv"
        θ_cop = vec(Matrix(CSV.read(cop_nd_sol_path, DataFrame, header=0)))
    end
    return θ_cop
end


# function adf_pv_aux(adf_stat, deterministic)
#     # helper function for p-value computation
#     #
#     # based on James G. MacKinnon, "Approximate Asymptotic Distribution Functions for
#     # Unit-Root and Cointegration Tests", Journal of Business & Economic Statistics,
#     # Vol. 12, No. 2 (Apr., 1994), pp. 167-176, http://www.jstor.org/stable/1391481
#     #
#     # A slightly more accurate approach could be based on an algorithm described in
#     # James G. MacKinnon, "Numerical Distribution Functions for Unit Root and Cointegration
#     # Tests", Journal of Applied Econometrics, Vol. 11, No. 6 (Nov.-Dec., 1996),
#     # pp. 601-618, http://www.jstor.org/stable/2285154
#     # and made available under GPL as Fortran routines (used e.g. by R's urca package) at
#     # http://qed.econ.queensu.ca/pub/faculty/mackinnon/numdist/

#     pv_coeff_smallp = [
#         0.6344 1.2378 0.032496 -1.04 -19.04
#         2.1659 1.4412 0.038269 -1.61 -18.83
#         3.2512 1.6047 0.049588 -2.89 -16.18
#         4.0003 1.6580 0.048288 -3.21 -17.17
#     ]

#     pv_coeff_largep = [
#         0.4797 0.93557 -0.06999  0.033066 -1.04 Inf
#         1.7339 0.93202 -0.12745 -0.010368 -1.61 2.74
#         2.5261 0.61654 -0.37956 -0.060285 -2.89 0.70
#         3.0778 0.49529 -0.41477 -0.059359 -3.21 0.54
#     ]

#     if deterministic == :none
#         tab_row = 1
#     elseif deterministic == :constant
#         tab_row = 2
#     elseif deterministic == :trend
#         tab_row = 3
#     elseif deterministic == :squared_trend
#         tab_row = 4
#     else
#         throw(ArgumentError("deterministic = $(deterministic) is invalid"))
#     end

#     if adf_stat < pv_coeff_smallp[tab_row, 5]
#         aux_var = -Inf
#     elseif adf_stat > pv_coeff_largep[tab_row, 6]
#         aux_var = Inf
#     else
#         if adf_stat < pv_coeff_smallp[tab_row, 4]
#             aux_var = pv_coeff_smallp[tab_row, 1] + pv_coeff_smallp[tab_row, 2] * adf_stat +
#                         pv_coeff_smallp[tab_row, 3] * (adf_stat^2)
#         else
#             aux_var = pv_coeff_largep[tab_row, 1] + pv_coeff_largep[tab_row, 2] * adf_stat +
#                         pv_coeff_largep[tab_row, 3] * (adf_stat^2) +
#                         pv_coeff_largep[tab_row, 4] * (adf_stat^3)
#         end

#     end
# end


function n_factors(X, r_max; include_plot::Int=0, τ::Float64=0.5)
    # Estimates the number of relevant factors for dataset X following Freyaldenhoven (2021)
    #
    # INPUT
    # X: (Txn) matrix of data
    # r_max: upper bound on number of factors
    # include_plot(optional): If set to 1 includes illustrative figures)
    # τ(optional): default is 0.5
    #
    # OUTPUT
    # FR: Estimate for the number of factors

    T, n = size(X)
    z = round(Int, min(0.7 * n^τ * sqrt(log(log(n))), n))

    # SVD decomposition
    _, d, V = svd(X ./ sqrt(T))

    # Factor Loadings i.e., Projection matrix
    Lambda = V[:, 1:r_max] * Diagonal(d[1:r_max])

    # Sorting
    sorted = sort(abs.(Lambda), dims=1, rev=true)
    largest_z = sorted[1:z, :]

    error_part = d[r_max:end]
    estimate_variance = sum(error_part .^ 2) / n

    Shat = zeros(r_max)
    T2 = zeros(r_max)

    for k in 1:r_max
        Shat[k] = ((largest_z[:, k]' * largest_z[:, k] / z) / sqrt(Lambda[:, k]' * Lambda[:, k] / n))
        T2[k] = (Lambda[:, k]' * Lambda[:, k]) * Shat[k]^2
    end

    incl_mock_T2 = [estimate_variance * n; T2]
    T2_ratio = incl_mock_T2[1:r_max] ./ T2[1:r_max]
    FR = findall(x -> x == maximum(T2_ratio), T2_ratio)[1]
    FR -= 1

    # Determine 'maxoutdim'
    eig_X = d[1:min(n, T)] .^ 2
    MOdim = sum(incl_mock_T2[1:r_max] .> eig_X[1:r_max])
    println("Max dimension outputed: ", MOdim)
    println("FR: ", FR)

    init_path = BASE_PATH

    if include_plot == 1
        eig_X = d[1:min(n, T)] .^ 2

        plot1 = Plots.scatter(0:r_max, [1; Shat .^ 2], ylab=L"\textrm{Value}", ms=5, markerstrokewidth=3, mc=:orange, guidefontsize=12, tickfontsize=12, legendefontsize=12, xformatter=:latex, yformatter=:latex, xlabel=L"k", label=L"\hat{S}^2", legend=:topright)
        Plots.savefig(init_path * "/7_Results/factor_analysis/factor_analysis_plot1.pdf")

        plot2 = Plots.plot(0:r_max, [estimate_variance * n; eig_X[1:r_max]], lw=2, guidefontsize=12, tickfontsize=12, legendefontsize=12, xformatter=:latex, yformatter=:latex, xlabel=L"k", ylabel=L"\textrm{Value}", label=[L"\hat{\Upsilon}^0_k" L"\hat{\Upsilon}^2_k"], legend=:topright)
        Plots.scatter!([0], [0], label=" ", ms=0, mc=:white, msc=:white)
        Plots.plot!(0:r_max, incl_mock_T2, seriestype=:line, lw=2, label=L"\hat{\Upsilon}^2_k")

        Plots.savefig(init_path * "7_Results/factor_analysis/factor_analysis_plot2.pdf")

        # Determine 'maxoutdim'
        MOdim = sum(incl_mock_T2[1:r_max] .> eig_X[1:r_max])
        println("Max dimension outputed:", MOdim)

        M = MultivariateStats.fit(PCA, Matrix(X' ./ sqrt(T)); maxoutdim=MOdim, method=:svd, mean=0)

        # Plot the factors 
        pcs = MultivariateStats.transform(M, Matrix(X' ./ sqrt(T)))
        saxis = 1:size(pcs, 2)

        ls_vec = [:solid, :solid, :solid, :dash, :dot, :dashdot, :solid, :solid, :solid, :solid, :solid, :solid, :solid, :solid, :solid, :solid]
        lc_vec = [:black, :red, :blue, :green, :orange, :purple, :brown, :pink, :cyan, :magenta, :yellow, :gray]
        lo_vec = [1.0, 0.9, 0.8, 0.6, 0.4, 0.2, [0.2 for _ in 1:size(pcs, 2)-6]...]
        Plots.plot()
        for i in axes(pcs, 1)
            Plots.plot!(saxis, pcs[i, :], xformatter=:latex, yformatter=:latex, legend_columns=3, label=L"\textrm{Factor \,\, %$(i)}", linealpha=lo_vec[i], ls=ls_vec[i], lc=lc_vec[i], ylab=L"\textrm{Factor\,\, Value}", xticks=(saxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(fred_qd_stationary[1:20:end, :daten])]))
        end
        Plots.savefig(init_path * "7_Results/factor_analysis/agg_factors.pdf")
    end

    return MOdim
end

