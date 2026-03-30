function total_grid_points(N::Int, G::Int)
    # Initialize total points count
    total_points = 0

    # Add up points for all dimensions from 1 to N
    for k in 2:N
        # The number of k-dimensional sub-copulas is binomial(N, k)
        # Each such sub-copula has G^k points
        total_points += binomial(N, k) * G^k
    end

    return total_points
end


"""Infer the source ordering embedded in a saved sigma filename.

Expected patterns include e.g.
`..._series_CEX_and_CPS_and_PSID_all.jld2` or
`..._quintiles_series_CEX_and_CPS2_and_SIPP1_all.jld2`.

Returns an empty vector if the pattern is not recognized.
"""
function sigma_sources_from_path(path::AbstractString)
    base = basename(path)
    namestr = replace(base, r"^.*_series_" => "", r"_all\.jld2$" => "")
    df_map = split(namestr, "_and_")
    return df_map
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


"""Compute the standard deviation of `v` ignoring NaNs."""
function _nanstd(v::AbstractVector{<:Real})
    s = 0.0
    s2 = 0.0
    n = 0
    @inbounds for x in v
        if !isnan(x)
            n += 1
            s += x
            s2 += x * x
        end
    end
    n <= 1 && return NaN
    μ = s / n
    # numerically safe-ish variance
    var = max(0.0, (s2 - n * μ * μ) / (n - 1))
    return sqrt(var)
end




"""Perturb coefficients with independent N(0, sd^2) noise (per coefficient and period).

Accepted `sd` inputs:

- K×T×D array of *draws* (e.g. empirical bootstrap coefficient draws). In this case we
  estimate a time-invariant per-coefficient SD by demeaning within each period across draws
  and pooling the within-period variance across time.
- K×T or K×1 array of SDs (used directly / pooled across time).
"""
@inline function _nanstd_view(x)
    n = 0
    mean = 0.0
    m2 = 0.0
    @inbounds for v in x
        if !isnan(v)
            n += 1
            δ = v - mean
            mean += δ / n
            m2 += δ * (v - mean)
        end
    end
    return n > 1 ? sqrt(m2 / (n - 1)) : NaN
end

# Carry/back-fill NaNs in-place for a 1D row; returns whether the row had any non-NaN.
function _carry_fill_row!(row)
    T = length(row)

    first_idx = 0
    @inbounds for t in 1:T
        if !isnan(row[t])
            first_idx = t
            break
        end
    end
    first_idx == 0 && return (false, 0)

    n_filled = 0
    first_val = row[first_idx]

    @inbounds for t in 1:(first_idx-1)
        if isnan(row[t])
            row[t] = first_val
            n_filled += 1
        end
    end

    @inbounds for t in (first_idx+1):T
        if isnan(row[t])
            row[t] = row[t-1]
            n_filled += 1
        end
    end

    return (true, n_filled)
end

# Build K×T SD matrix from either draws (K×Tn×Dn) or precomputed SDs (K×T or K×1).
function _build_sds(sd, K::Int, T::Int)
    Tn = size(sd, 2)
    @assert T <= Tn "by_period=true requires sd time dimension ≥ size(X,2)"

    sds = fill(NaN, K, Tn)

    @inbounds for k in 1:K, t in 1:Tn
        sds[k, t] = _nanstd_view(@view sd[k, t, :])
    end

    n_filled = 0
    n_all_missing = 0
    @inbounds for k in 1:K
        row = @view sds[k, :]
        ok, filled = _carry_fill_row!(row)
        if ok
            n_filled += filled
        else
            n_all_missing += 1
        end
    end

    n_filled > 0 && @info "Filled $n_filled missing SD entries by carry/back-fill"
    n_all_missing > 0 && @info "$n_all_missing coefficients had no SD information at any time"

    return sds
end

# --- main ------------------------------------------------------------------

function perturb_coefficients_sd!(
    X::AbstractMatrix,
    sd;
    rng=Random.default_rng(),
    scale::Real=1.0,
    by_period::Bool=true,
)
    by_period || error("by_period=false is no longer supported; use period-by-period SDs")

    K, T = size(X)
    @assert size(sd, 1) == K

    sds = _build_sds(sd, K, T)

    n_no_sd = 0
    @inbounds @views for t in 1:T
        xcol = X[:, t]
        for k in 1:K
            x = xcol[k]
            if !isnan(x)
                s = sds[k, t]
                if !isnan(s) && s > 0
                    xcol[k] = x + scale * s * randn(rng)
                else
                    n_no_sd += 1
                end
            end
        end
    end

    n_no_sd > 0 && @info "$n_no_sd coefficients had non-positive/NaN SD; left unperturbed"
    return X
end


function data_constructor(obs_data::ObservedData, model_options)
    """Prepares data, of any frequency, for estimation."""

    # sigma_invhalf_path=nothing, sigma_sources=nothing, rng=Random.default_rng(), sigma_noise_scale::Real=1.0

    @unpack files, df_vec, gdp_series = obs_data
    @unpack estimator, number_of_dfs, measures, information, equivalized, bottom_coded, blind_to, tag = model_options

    @unpack grid_pcf, grid_cop = estimator

    # Pre-allocate
    dfs = Vector{Matrix{Float64}}(undef, number_of_dfs)
    time_vec = Vector{Vector}(undef, number_of_dfs)
    freq_type = Vector{String}(undef, number_of_dfs)

    # Extract year and quarter from time variable
    gdp_series[!, "time"] = QuarterlyDate.(gdp_series[!, "time"])
    gdp_series[!, "year"] = Dates.year.(gdp_series[!, "time"])
    gdp_series[!, "quarter"] = Dates.quarterofyear.(gdp_series[!, "time"])

    # Create dictionary that stores the years and quarters of the sample data
    time_dict = Dict()

    un_perturbed_dfs = Vector{Matrix{Float64}}(undef, number_of_dfs)

    # For each dataset ... 
    for (j, df) in enumerate(df_vec.data)
        df_name = df_vec.df_names[j]
        time_dict[j] = Dict()

        # Quickly: rename "consumption" to "consum" if it exists
        if "consumption" ∈ names(df)
            rename!(df, "consumption" => "consum")
        end

        data, non_missing_cols = select_data(df, measures, equivalized, bottom_coded, blind_to, df_name)

        # Stores years of measurements based on 'information' criteria  
        time_vec[j], freq_type[j] = define_time_of_measures(data)

        # Pre-allocate 
        T = length(time_vec[j]) # number of periods - observed
        D = length(measures)

        local copulas, pcfs, immutable # sub_copulas
        if typeof(estimator) <: HistogramEstimator
            immutable = grid_cop + (D - 1) * (grid_cop - 1)
            copulas = zeros(grid_cop^D - immutable, T)  # filled with NaN later
            pcfs = [zeros(grid_pcf, T) for _ in 1:D]

        elseif typeof(estimator) <: KernelEstimator
            immutable = grid_cop + (D - 1) * (grid_cop - 1) # TODO: kernel has immutable?
            copulas = zeros(grid_cop^D - immutable, T)  # filled with NaN later
            pcfs = [zeros(grid_pcf, T) for _ in 1:D]

        elseif typeof(estimator) <: SeriesEstimator
            # In this case, we need to have separate containers for the sub-copulas and the larger copula. Why? That way the coefficients of a slice, where only the slice could be observed, is not comparable to the coefficients of a slice of the larger dimensional object 
            # Also, with this estimator, there is no need for immutables 

            # Based on the number of dimensions, generate a container for the sub-copulas
            # sub_copula_l  = D < 3 ? 0 : grid_cop^(D-1) * D # the order of the subcopulas should (1,2) (1,3) (2,3), (1,2,3) -- in a long vector form
            # nc            = total_grid_points(D, grid_cop) # total number of points in the copula + sub-copulas
            immutable = grid_cop + (D - 1) * (grid_cop - 1) # TODO: kernel has immutable?
            copulas = zeros((grid_cop)^D - immutable, T)  # filled with NaN later
            pcfs = [zeros(grid_pcf, T) for _ in 1:D]
        end

        # Initialize counter of all periods 
        q = 0
        years = unique(time_vec[j])
        ft = freq_type[j]

        # Per time period ...    
        for yr in years
            year_data = filter(row -> row.year == yr, data)
            correction = filter(row -> row.year == yr, gdp_series) # TODO: this only works for annual and quarterly data       

            # returns observed periods in the year
            time_dict[j][yr] = get_periods(year_data, ft)

            for (p, actual_period) in enumerate(time_dict[j][yr])
                obs_meas = String[]
                period_data = filter(row -> row[ft] == actual_period, year_data)

                s = 1
                # Per measure ...  
                for (m, meas) in enumerate(measures)
                    # Check which measures are observed for each time period 
                    if meas ∈ names(period_data)
                        if any(.!isnan.(period_data[:, meas]))
                            push!(obs_meas, meas)

                            # For the respective macro series, extract the correction data for the respective quarter 
                            local correction_q, correction_agg
                            if ft == "year"
                                correction_q = filter(row -> row[ft] == Int.(round(actual_period, digits=0)), correction)[:, meas*"_per_hh"][end]
                            else
                                correction_q = filter(row -> row[ft] == Int.(round(actual_period, digits=0)), correction)[:, meas*"_per_hh"]
                            end

                            # Get percentile functions for each period
                            pcfs[m][:, q+p] = get_percentile_functions!(period_data, meas, non_missing_cols, obs_meas, correction_q, df_name, estimator)
                        else
                            pcfs[m][:, q+p] .= NaN
                        end
                    else
                        pcfs[m][:, q+p] .= NaN
                    end
                end

                if df_name != "SCF" #TODO: automate this s.t. any dataset with "impnum" will be treated (in a loop)
                    copulas[:, q+p] .= get_copulas(period_data, measures, obs_meas, estimator) # TODO: of course, a subcopula is not guaranteed to be observed -> use length of obs_meas
                else
                    temp_cop = zeros(size(copulas[:, q+p]))

                    for i in 1:5
                        period_dataᵢ = filter(row -> row.impnum == i, period_data)
                        # copulas[:, q + p]   .+= get_copulas(period_dataᵢ, measures, obs_meas, estimator)
                        temp_cop .+= get_copulas(period_dataᵢ, measures, obs_meas, estimator)
                    end

                    copulas[:, q+p] .= temp_cop ./ 5
                end
            end
            # Shift the time dimension as we move year by year, period by period
            q += length(time_dict[j][yr])
        end
        # Concatenate topologies 
        dfs[j] = vcat(copulas, vcat(pcfs...))

        un_perturbed_dfs[j] = deepcopy(dfs[j])

        # If requested, inject empirical-style coefficient noise into simulated HANK datasets.
        # Priority: if SDs-from-draws are provided for the target dataset, use those.
        # Otherwise fall back to the dataset-specific block from Σ̂⁻¹² (block diagonal across sources).
        # if occursin("HANK", df_name)
        #     target = hank_target_dataset(df_name)
        #     noise_df_files = infer_noise_paths(files, model_options)

        #     # Import noise file 
        #     noise_df = jldopen(noise_df_files[target], "r")["noise"]

        #     # perturb coefficients
        #     dfs[j] = perturb_coefficients_sd!(dfs[j], noise_df)
        # end
    end

    return dfs, un_perturbed_dfs, time_vec, time_dict, freq_type
end


function infer_noise_paths(files::Dict{String,String}, model_options::ModelOptions)
    isempty(files) && return Dict{String,String}()

    any_path = first(values(files))
    base_dir = dirname(any_path)                      # .../2_Data_processing
    noise_dir = joinpath(base_dir, "noise_distributions")

    m_label = _noise_measures_label(model_options.measures)
    grid_tag = _noise_grid_tag(model_options.estimator)
    end_year = get(model_options.data_cutoffs, "end", "") != "" ? String(model_options.data_cutoffs["end"])[1:4] : "all"
    ci_tag = model_options.tag == " higher order15" ? model_options.tag : ""

    out = Dict{String,String}()
    for src in infer_empirical_sources(files)
        ci_source = occursin("SCF", src) ? "SCF" : src
        fname = "noise_draws_" * m_label * grid_tag * "_" * ci_source * "_" * end_year * ci_tag * ".jld2"
        out[src] = joinpath(noise_dir, fname)
    end

    return out
end

"""Infer which empirical datasets are implied by the provided `files` mapping.

This is primarily intended for HANK runs where keys look like "HANK a 2" but
file basenames look like "HANK_PSID_2.csv".
"""
function infer_empirical_sources(files::Dict{String,String})
    sources = Set{String}()
    known = Set(["PSID", "CPS", "CEX", "SCF", "CPS2", "SIPP1", "SIPP2", "SIPP3"])

    for (k, p) in files
        if k in known
            push!(sources, k)
        end

        b = basename(p)
        if occursin("HANK_", b)
            occursin("HANK_PSID", b) && push!(sources, "PSID")
            occursin("HANK_CPS_", b) && push!(sources, "CPS")
            occursin("HANK_CEX", b) && push!(sources, "CEX")
            occursin("HANK_SCF", b) && push!(sources, "SCF")
            occursin("HANK_CPS2", b) && push!(sources, "CPS2")
            occursin("HANK_SIPP1", b) && push!(sources, "SIPP1")
            occursin("HANK_SIPP2", b) && push!(sources, "SIPP2")
            occursin("HANK_SIPP3", b) && push!(sources, "SIPP3")
        end
    end

    return sort!(collect(sources))
end



function _noise_grid_tag(estimator)
    grid_choice = typeof(estimator) <: SeriesEstimator ? estimator.integral_pcf_grid : estimator.grid_pcf
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
    else
        # Fallback (rare): encode raw grid_choice
        grid_granularity = "_grid$(grid_choice)"
    end

    return grid_granularity * est_tag
end



function _noise_measures_label(measures::AbstractVector{<:AbstractString})
    label = ""
    sorted_measures = sort(collect(measures))
    for (i, meas) in enumerate(sorted_measures)
        if i < length(sorted_measures)
            label *= String(meas) * "_and_"
        else
            label *= String(meas)
        end
    end
    return label
end



function generate_mesh(grid_points, dimension)
    # Generate the Cartesian product of [1/grid, 2/grid, ..., grid/grid] across 'dimension' dimensions
    grid_points = Iterators.product(fill(grid_points, dimension)...)

    # Use a list comprehension to construct the mesh matrix directly
    mesh = [collect(point) for point in grid_points]

    return hcat(mesh...)' # #TODO: maybe i could make it faster
end



# Error: ArgumentError("unable to check bounds for indices of type Missing")
# Error: MethodError(StatsAPI.weights, (Union{Missing, Float64}[18191.949, 23


function get_periods(year_data, freq_type_of_data)
    # Return unique values in the year column of the data
    if freq_type_of_data == "year"
        exact_periods_in_year = sort(unique(year_data[:, :year]))
    elseif freq_type_of_data == "quarter"
        exact_periods_in_year = sort(unique(year_data[:, :quarter]))
    elseif freq_type_of_data == "month"
        exact_periods_in_year = sort(unique(year_data[:, :month]))
    else
        error("Data frequency not recognized.")
    end
    return exact_periods_in_year
end

# TODO: if some data has missings within the year, then drop those before estimation 
function get_percentile_functions!(period_data, meas, non_missing_cols, obs_meas, correction, key, estimator)

    # If the measure is always missing ... 
    if meas ∉ non_missing_cols
        return NaN
    end

    # If the measure is missing this period 
    if meas ∉ obs_meas
        return NaN
    end


    pcf = perform_percentile_corrections!(period_data, meas, key, estimator, true, correction)

    return pcf

end

#TODO: not a good function per se. 
function define_time_of_measures(data)
    # Get unique years
    year_vec = Int.(sort(unique(data[:, :year])))  # all possible years for this df 

    local freq_type
    # A string that tells me the frequency of the data 
    if "quarter" ∈ names(data) && "year" ∈ names(data) && "month" ∉ names(data)
        freq_type = "quarter"
    elseif "month" ∈ names(data) && "year" ∈ names(data)
        freq_type = "month"
    elseif "year" ∈ names(data) && "quarter" ∉ names(data) && "month" ∉ names(data)
        freq_type = "year"
    else
        error("Data frequency not recognized.")
    end

    @info("The frequency of the data is $freq_type" * "ly in some way")

    # For each year, see how many periods there are and fill in a vector with that count 
    periods_in_year = Vector{Int64}(undef, length(year_vec))

    for (i, yr) in enumerate(year_vec)
        # get data for the year 
        year_data = filter(row -> row.year == yr, data)
        periods_in_year[i] = length(unique(year_data[:, freq_type]))
    end

    # Now repeat year_vec by the number of periods in each year
    time_vec = vcat(fill.(year_vec, periods_in_year)...)[:] # repeat(year_vec, periods_in_year)

    return time_vec, freq_type
end


function select_data(df, measures, equivalized, bottom_coded, blind_to, key)

    # From CSV import
    non_missing_cols = clean_data!(df, measures, equivalized, bottom_coded, blind_to, key)

    return df, non_missing_cols
end



# function estimate_levels(year_data, levels_vec, meas, non_missing_cols, grid, gdp)
#     if meas ∉ non_missing_cols
#         levels_vec .= NaN
#         return dct(levels_vec)
#     end

#     if all(isnan.(year_data[:, meas]))
#         levels_vec .= NaN 
#         return dct(levels_vec)
#     end

#     try 
#         assign_quantile_groups!(year_data, meas, grid)  #TODO: I no longer need to worry about this 
#     catch ex
#         error("Failed to assign quantile groups: $ex")
#     end

#     non_missing        = filter(meas => !isnan, year_data)
#     # total              = wsum(non_missing[:, meas], non_missing[:, :weight])  # weighted sum 
#     for g in 1:grid 
#         filt_d         = filter(row -> row[meas * "_quantile"] == float(g), non_missing)  
#         levels_vec[g]  = wsum(filt_d[:, meas], weights(filt_d[:, :weight])) / gdp
#     end

#     return dct(levels_vec)
# end

# # Find the index of the weight column 
# weight_col = findfirst(isequal("weight"), names(period_data))

# # Use RCall to generate copulas 
# R_data   = Matrix(period_data[!, [measures..., "weight"]])

# # Repeat Observations based on the size of the weights 
# R_data   = vcat([repeat(R_data[i,:]', inner=(max.(Int_NaN.(ceil.(R_data[i,3])), 0),1)) for i in axes(R_data,1)]...)
# R_data   = R_data[:, 1:2]  # drop the weight column
# mesh     = create_mesh(grid, dimension)

# R"""
# ranked  <- apply($R_data, 2, rank) / (nrow($R_data) + 1)
# kde.fit <- kdecopula::kdecop(ranked) # perhaps try method = "MR", 
# # cop_den <- kdecopula::dkdecop($mesh, kde.fit)  # density 
# cop_dis <- kdecopula::pkdecop($mesh, kde.fit) # distribution
# """
# @rget cop_dis 

# # convert distribution to density 
# # cop_den = copula_to_density(cop_dis, cop_size, grid)

# cop .= dct(reshape(cop_dis, cop_size...))


function estimate_copula!(::Type{<:HistogramEstimator}, cop, cop_ind, period_data, obs_measures, grid, grid_type, meas_v, scale)
    # Assign quantiles 
    period_data = assign_quantile_groups_for_copula!(period_data, obs_measures, grid, grid_type)

    # Subset to non-missing 
    if "quarter" in names(period_data)
        # select!(period_data, [meas_v..., measures..., "quarter", "year", "weight"]) # for smoothing
        select!(period_data, [meas_v..., "quarter", "year", "weight"])

    else
        # select!(period_data, [meas_v..., measures..., "year", "weight"]) # for smoothing 
        select!(period_data, [meas_v..., "year", "weight"])
    end

    filter!(row -> !any(isnan(x) for x in row), period_data)

    cop[cop_ind...] .= Array(prop(freqtable(period_data, tuple(meas_v...)..., weights=period_data.weight))) # slow because i have to index 

    @assert sum(cop[cop_ind...] .< 1) == length(cop[cop_ind...]) && sum(cop[cop_ind...] .>= 0) == length(cop[cop_ind...]) && sum(cop[cop_ind...]) ≈ 1
    cop[cop_ind...] .= dct(cop[cop_ind...]) ./ scale
end


function estimate_copula!(::Type{<:FakeHistEstimator}, cop, cop_ind, period_data, obs_measures, grid, grid_type, meas_v, scale)
    # Assign quantiles 
    # println(grid, grid_type)
    # println(size(cop))
    # println(cop_ind)
    period_data = assign_quantile_groups_for_copula!(period_data, obs_measures, grid, grid_type)

    # Subset to non-missing 
    if "quarter" in names(period_data)
        # select!(period_data, [meas_v..., measures..., "quarter", "year", "weight"]) # for smoothing
        select!(period_data, [meas_v..., "quarter", "year", "weight"])

    else
        # select!(period_data, [meas_v..., measures..., "year", "weight"]) # for smoothing 
        select!(period_data, [meas_v..., "year", "weight"])
    end

    filter!(row -> !any(isnan(x) for x in row), period_data)

    cop[cop_ind...] .= Array(prop(freqtable(period_data, tuple(meas_v...)..., weights=period_data.weight))) # slow because i have to index 

    @assert sum(cop[cop_ind...] .< 1) == length(cop[cop_ind...]) && sum(cop[cop_ind...] .>= 0) == length(cop[cop_ind...]) && sum(cop[cop_ind...]) ≈ 1
    cop[cop_ind...] .= dct(cop[cop_ind...]) ./ scale
end



function estimate_copula!(::Type{<:KernelEstimator}, cop, cop_ind, period_data, obs_measures, grid, grid_type, meas_v, scale)
    # Assign quantiles 
    period_data = assign_quantile_groups_for_copula!(period_data, obs_measures, grid, grid_type)

    # Subset to non-missing 
    if "quarter" in names(period_data)
        # select!(period_data, [meas_v..., measures..., "quarter", "year", "weight"]) # for smoothing
        select!(period_data, [meas_v..., "quarter", "year", "weight"])

    else
        # select!(period_data, [meas_v..., measures..., "year", "weight"]) # for smoothing 
        select!(period_data, [meas_v..., "year", "weight"])
    end

    filter!(row -> !any(isnan(x) for x in row), period_data)

    cop[cop_ind...] .= beta_copula_estimator(period_data, obs_measures, grid, grid_type)
    cop[cop_ind...] .= dct(cop[cop_ind...]) ./ scale
end


function estimate_copula!(::Type{<:FakeKernelEstimator}, cop, cop_ind, period_data, obs_measures, grid, grid_type, meas_v, scale)
    # Assign quantiles 
    period_data = assign_quantile_groups_for_copula!(period_data, obs_measures, grid, grid_type)

    # Subset to non-missing 
    if "quarter" in names(period_data)
        # select!(period_data, [meas_v..., measures..., "quarter", "year", "weight"]) # for smoothing
        select!(period_data, [meas_v..., "quarter", "year", "weight"])

    else
        # select!(period_data, [meas_v..., measures..., "year", "weight"]) # for smoothing 
        select!(period_data, [meas_v..., "year", "weight"])
    end

    filter!(row -> !any(isnan(x) for x in row), period_data)

    cop[cop_ind...] .= beta_copula_estimator(period_data, obs_measures, grid, grid_type)
    cop[cop_ind...] .= dct(cop[cop_ind...]) ./ scale
end


function select_grid_points(grid)

    # local grid_points
    # if grid_type == "uniform"
    interval = 1 / grid
    grid_points = collect(interval:interval:1)

    # Add an epsilon to the bottom and subtract an epsilon from the top
    grid_points[end] = grid_points[end] - 1e-6

    return grid_points
end


# function beta_copula_estimator(period_data, obs_measures, grid, grid_type)

#     grid_points = select_grid_points(grid, grid_type)
#     X = select(period_data, obs_measures)

#     # Some cleaning 
#     X = coalesce.(X, NaN)
#     X = filter(row -> !any(isnan(x) for x in row), X)

#     # Make grid 
#     lg = length(grid_points)
#     dimension = length(obs_measures)
#     mesh_grid = generate_mesh(grid_points, dimension)

#     R"""
#     EC = copula::C.n($mesh_grid, X = $X, smoothing = "beta")        
#     """
#     @rget EC

#     # Reshape arbitrarily sized matrix
#     ocop_size = tuple([grid for i in 1:dimension]...)
#     mat_dist = reshape(EC, ocop_size...)

#     # Correct the matrix s.t. it is a proper copula
#     mat_dist .= mat_dist ./ maximum(mat_dist)

#     return mat_dist
# end


## Removed from get_copulas
# mesh        = generate_mesh(grid_points, obs_dim) 
# period_data = assign_quantile_groups_for_copula!(period_data, obs_measures, grid, grid_points)

# sel_cols = period_data[:, meas_v] ./ (grid + 1)
# w        = period_data[:, "weight"]

# println(sel_cols[1:5,:])

# R"""
# vcop <- rvinecopulib::vinecop($sel_cols,
#                 var_types = rep("c", $obs_dim), # continuous variables
#                 nonpar_method = "linear", # local likelihood estimatin of order 1 -- best estimator from simulations
#                 mult = 1, # greater than 1 = more smooth
#                 selcrit = "aic", # criterion for family selection -- doesnt matter for us 
#                 weights = $w,
#                 presel = TRUE, # pre-select families that better represent the data -- for us, it's always non-parametric 
#                 trunc_lvl = Inf, # no truncation on the trees -> looks at all pairs of variables vs. making an assumption on certain pairs
#                 tree_crit = "tau",
#                 )
# cop_dens <- rvinecopulib::dvinecop($mesh, vcop, cores = 1)            
# """
# @rget cop_dens

# # Normalize the copula --- 'cop_dens' is a vector
# cop_dens        .= cop_dens ./ sum(cop_dens)
# cop_dens         = reshape(cop_dens, ocop_size...)
# cop[cop_ind...] .= dct(cop_dens) ./ scale


function get_copulas(period_data, measures, obs_measures, estimator; with_immutable=false)
    """Copulas are created when all dimensions are observed. Sub-copulas are created with partial observability.
    For series estimators, we generate both the copula and the sub-copulas.
    For dim=1, copula does not exist, so returns NaNs."""

    @unpack grid_cop = estimator

    # Timeless parameters 
    dimension = length(measures)
    cop_size = tuple([grid_cop for i in 1:dimension]...) # doesnt affect SeriesEstimator
    cop = fill!(Array{Float64}(undef, cop_size...), NaN)

    # Based on observed 
    obs_dim = length(obs_measures)
    meas_v = obs_measures .* "_quantile"
    cop_ind = copula_case(obs_measures, measures)
    scale = sqrt(grid_cop)^(dimension - obs_dim)

    # All is observed 
    if obs_dim <= 1
        return NaN
    else
        if typeof(estimator) <: HistogramEstimator || typeof(estimator) <: KernelEstimator #TODO: kernel not yet estimated 
            # At least two measures are observed
            estimate_copula!(typeof(estimator), cop, cop_ind, period_data, obs_measures, grid_cop, grid_type_cop, meas_v, scale)
            # println(size(cop))
            if with_immutable
                return cop[:]
            else
                return remove_immutable(cop)
            end

        elseif typeof(estimator) <: SeriesEstimator # tracking sub-copulas as well 
            cop_weights = series_approximate_copula(period_data, obs_measures, estimator)

            # Fill cop based on case 
            cop[cop_ind...] = cop_weights

            # Return a vector of the copula
            if with_immutable
                return cop[:]
            else
                return remove_immutable(cop)
            end
        end
    end
end




function generate_unique_combinations(V)

    # Determine the length of the tuple
    n = length(V)


    # Create an empty array to store all combinations
    all_combinations = []

    # Generate combinations for each size from 1 to n (inclusive)
    for k in 2:n # I don't want to generate the 1-combinations
        append!(all_combinations, collect(combinations(V, k)))
    end

    a = unique(sort.(all_combinations))

    return a
end


# unique(sort.(all_combinations))

function series_approximate_copula(period_data, obs_measures, estimator)

    @unpack grid_cop = estimator

    cols_to_keep = vcat(obs_measures, "weight")

    # Some data cleaning
    X1 = select(period_data, cols_to_keep) #TODO: may not work with datasets without weights
    X2 = coalesce.(X1, NaN)
    for m in cols_to_keep
        filter!(m => !isnan, X2)
    end

    # Separate data from weights
    W = Vector{Float64}(X2.weight)
    # select!(X2, obs_measures)

    # Convert to Matrix for fast BLAS path
    R = Matrix{Float64}(X2[:, obs_measures])

    # Rank the data
    n = size(R, 1)
    for j in axes(R, 2)
        R[:, j] .= rankdata(@view(R[:, j])) ./ (n + 1)
    end

    # Series estimate the copula — vectorized BLAS-based version
    # cop_coefs = get_copula_coefficients(X2, W, grid_cop - 1)
    cop_coefs = get_copula_coefficients_fast(R, W, grid_cop - 1)

    return cop_coefs[:]
end


function hist_series_copula_estimator(X, obs_measures, grid_cop, granularity_cop, grid_type)
    period_data = assign_quantile_groups_for_copula!(X, obs_measures, granularity_cop, grid_type)
    meas_v = obs_measures .* "_quantile"

    # Subset to non-missing 
    if "quarter" in names(period_data)
        # select!(period_data, [meas_v..., measures..., "quarter", "year", "weight"]) # for smoothing
        select!(period_data, [meas_v..., "quarter", "year", "weight"])

    else
        # select!(period_data, [meas_v..., measures..., "year", "weight"]) # for smoothing 
        select!(period_data, [meas_v..., "year", "weight"])
    end

    filter!(row -> !any(isnan(x) for x in row), period_data)

    mat_dist = Array(prop(freqtable(period_data, tuple(meas_v...)..., weights=period_data.weight))) # slow because i have to index 

    # Approximate the copula
    lc = length(obs_measures)
    dom_1 = [1.0, 0.0]
    A_plan, _ = define_copula_plan(grid_cop, granularity_cop, lc, dom_1)

    return chebyshev_weights(mat_dist, A_plan)[:]
end


# function beta_series_copula_estimator(X3, obs_measures, grid_cop, granularity_cop)

#     X = Matrix(select(X3, obs_measures))
#     lc = length(obs_measures)
#     dom_1 = [1.0, 0.0]
#     A_plan, mesh_g = define_copula_plan(grid_cop, granularity_cop, lc, dom_1)

#     R"""
#     EC = copula::C.n($mesh_g, X = $X, smoothing = "beta") 
#     """
#     @rget EC # A Vector{Float64}

#     # Reshape arbitrarily sized matrix
#     ocop_size = tuple([granularity_cop for i in 1:lc]...)
#     mat_dist = reshape(EC, ocop_size...)

#     # Correct the matrix s.t. it is a proper copula
#     mat_dist .= mat_dist ./ maximum(mat_dist)

#     return chebyshev_weights(mat_dist, A_plan)[:]
# end


do_nothing(x) = x

function define_copula_plan(order_cop, granularity, lc, dom_1)
    tens_p = ntuple(_ -> nodes(granularity, :chebyshev_nodes, dom_1), lc) #TODO: change this to incorporate expand_grid 

    g = Grid(tens_p)
    O = ntuple(_ -> order_cop - 1, lc) # Minus 1 because the number of coefficients given is always the order + 1
    dom = hcat([dom_1 for _ in 1:lc]...)
    A_plan = CApproxPlan(g, O, dom) # last argument is a function that does nothing, can be changed if grid points need to be re-mapped            
    mesh_g = generate_mesh(tens_p[1].points, lc)

    return A_plan, mesh_g
end


function logit_transform_func(x)
    return log.(x ./ (1 .- x))
end



function inverse_hyperbolic_sine(x)
    return log.(x .+ sqrt.(x .^ 2 .+ 1))
end

function reverse_inverse_hyperbolic_sine(x)
    return (exp.(2 .* x) .- 1) ./ (2 .* exp.(x))
end


# a = DataFrame(rand(10000, 3), :auto)
# a[!, "a_quantile"] = rand(1:10, 10000)
# a[!, "b_quantile"] = rand(1:10, 10000)
# a[!, "c_quantile"] = rand(1:10, 10000)

# b = Array(prop(freqtable(a, tuple(["a_quantile", "b_quantile", "c_quantile"]...)...)))  
# b_T = treat_copula(b, true) 

# # case 1
# c = Array(prop(freqtable(a, tuple(["a_quantile", "b_quantile"]...)...)))
# c_T = treat_copula(c, true) ./ sqrt(10)

# # case 2
# c = Array(prop(freqtable(a, tuple(["a_quantile", "c_quantile"]...)...)))
# c_T = treat_copula(c) ./ sqrt(10)

# # case 3
# c = Array(prop(freqtable(a, tuple(["b_quantile", "c_quantile"]...)...)))
# c_T = treat_copula(c) ./ sqrt(10)

# # Test
# d = zeros(10,10,10)

# # case 1: missing last
# d[:, :, 1] .= c_T
# d = zeros(10,10,10)

# # case 2: missing second 
# d[:, 1, :] .= c_T
# d = zeros(10,10,10)

# # case 3: missing first 
# d[1, :, :] .= c_T
# d = zeros(10,10,10)

# f(c)         = sum((==(1)).(c.I)) >= length(c.I) - 1
# to_keep      = d[filter(!f, CartesianIndices(size(d)))]

# measures  = ["consum", "income", "wealth"]
# test_meas = ["consum", "income"]
# test_id   = findall(x -> x ∈ measures, test_meas)
# test_id_not = findall(x -> x ∉ test_meas, measures)
# colons = ntuple(_ -> Colon(), 3)  # Create a tuple of colons

# modified_indices = Base.setindex(modified_indices, 1, 2)

function copula_case(obs_meas, measures)
    dim = length(measures)
    colons = ntuple(_ -> Colon(), dim)

    # Get id's 
    id_wrt_measures = findall(x -> x ∈ measures, obs_meas)
    id_of_missings = findall(x -> x ∉ obs_meas, measures)

    if length(id_of_missings) == 0
        return colons
    else
        modified_indices = Base.setindex(colons, 1, id_of_missings[1])

        for i in 2:length(id_of_missings)
            modified_indices = Base.setindex(modified_indices, 1, id_of_missings[i])
        end

        # # Get the cartesian indices 
        # size_of_cop = tuple([grid for _ in 1:dim]...)
        # get_cartesian_indices

        return modified_indices
    end
end


function copula_case_validation(id_of_observed, measures)
    dim = length(measures)
    colons = ntuple(_ -> Colon(), dim)

    if length(id_of_observed) == dim
        return colons
    else
        # Create indices 
        id_of_missings = filter(x -> x ∉ id_of_observed, 1:dim)
        modified_indices = Base.setindex(colons, 1, id_of_missings[1])

        for i in 2:length(id_of_missings)
            modified_indices = Base.setindex(modified_indices, 1, id_of_missings[i])
        end

        # # Get the cartesian indices 
        # size_of_cop = tuple([grid for _ in 1:dim]...)
        # get_cartesian_indices

        return modified_indices
    end
end


function copula_to_density(cop_dis, cop_size, grid)
    cop_dis = reshape(cop_dis, cop_size...)
    cop_den = zeros(grid, grid)
    for h in 1:grid
        for j in 1:grid
            if h == 1 && j == 1
                cop_den[h, j] = cop_dis[h, j]
            elseif h == 1 && j != 1
                cop_den[h, j] = cop_dis[h, j] - cop_dis[h, j-1]
            elseif h != 1 && j == 1
                cop_den[h, j] = cop_dis[h, j] - cop_dis[h-1, j]
            else
                cop_den[h, j] = cop_dis[h, j] - cop_dis[h-1, j] - cop_dis[h, j-1] + cop_dis[h-1, j-1]
            end
        end
    end
    return cop_den
end

function logit_sum(dimension)
    if dimension == 1
        return NaN
    elseif dimension == 2
        return [-475, -520] # -490 is the mean
    elseif dimension == 3
        return [-7160, -7280] # -7210 is the mean
    end
end

function remove_immutable(A)
    """This removes the components corresponding to the first index of every dimension s.t. 
    the correct components are removed."""

    # Extract all cartesian indices that have D-1 1s in them 
    f(c) = sum((==(1)).(c.I)) >= length(c.I) - 1
    to_keep = A[filter(!f, CartesianIndices(size(A)))]

    return to_keep
end

function store_immutable(A)
    """This removes the components corresponding to the first index of every dimension s.t. 
    the correct components are removed."""

    f(c) = sum((==(1)).(c.I)) >= length(c.I) - 1
    to_store = A[filter(f, CartesianIndices(size(A)))]

    return to_store
end


function assign_quantile_groups_for_copula!(df, obs_meas, grid, grid_type)
    local grid_points
    if grid_type == "uniform"
        interval = 1 / grid
        grid_points = collect(interval:interval:1)

    elseif grid_type == "chebyshev"
        grid_points = chebyshev_nodes(grid)

    end


    # Remove rows where weight is missing
    df = coalesce.(df, NaN)
    filter!("weight" => !isnan, df)

    for rv in obs_meas
        sort!(df, [Symbol.(rv), :id])
        df[!, :running_weight] = cumsum(df[!, :weight])
        df[!, :running_weight] = df[!, :running_weight] ./ maximum(df[!, :running_weight])

        df[!, rv*"_quantile"] .= NaN
        for d in 1:grid
            upper_cond = (df[:, :running_weight] .<= grid_points[d])
            if d == 1
                df[upper_cond, rv*"_quantile"] .= d
            elseif d > 1 && d <= grid
                lower_cond = (df[!, :running_weight] .> grid_points[d-1])
                df[(lower_cond.&upper_cond).||(df[:, :running_weight].>grid_points[end]), rv*"_quantile"] .= d
            end
        end
        select!(df, Not(:running_weight))
    end
    return df
end

function scale_to_aggregates(non_missing, rv, correction)

    # Get weighted mean of the measure
    tot_data = mean(non_missing[!, rv], weights(non_missing[:, :weight]))

    # # See how that sum compares to the aggregate
    # println(typeof(correction))
    # println(size(correction))
    # println(correction)

    # println(typeof(tot_data))
    # println(size(tot_data))
    # println(tot_data)
    multiplier = abs.(correction / tot_data)

    # TODO: this is endogenous to the series given i.e., make sure the series is really aligned to the measure 
    if multiplier[1] > 20
        non_missing[!, rv] .= NaN
        return non_missing
    end

    # println("$rv : $multiplier")
    # Multiply the data by this multiplier
    non_missing[!, rv] .= non_missing[!, rv] .* multiplier

    # To correct the average 
    tot_scale = correction .- mean(non_missing[!, rv], weights(non_missing[:, :weight])) # should be zero if data is all positive 
    non_missing[!, rv] .= non_missing[:, rv] .+ tot_scale # in this way, total is correct, the average is corrected and ranks don't change. 
    return non_missing
end


function assign_quantile_groups!(non_missing, rv, grid, grid_points)
    sort!(non_missing, [Symbol.(rv), :id])

    non_missing[!, :running_weight] = cumsum(non_missing[!, :weight])
    non_missing[!, :running_weight] = non_missing[!, :running_weight] ./ maximum(non_missing[!, :running_weight])
    non_missing[:, rv*"_quantile"] .= NaN


    for d in 1:grid
        upper_cond = (non_missing[:, :running_weight] .<= grid_points[d])
        if d == 1
            non_missing[upper_cond, rv*"_quantile"] .= d
        elseif d > 1 && d <= grid
            lower_cond = (non_missing[!, :running_weight] .> grid_points[d-1])
            non_missing[(lower_cond.&upper_cond).||(non_missing[:, :running_weight].>grid_points[end]), rv*"_quantile"] .= d
        end
    end

    select!(non_missing, Not(:running_weight))
end

function treat_quantile_functions(non_missing, rv, grid, grid_points, correction, estimator)

    q_vec = fill(NaN, grid)
    non_missing_scaled = scale_to_aggregates(non_missing, rv, correction)

    if typeof(estimator) <: HistogramEstimator
        assign_quantile_groups!(non_missing_scaled, rv, grid, grid_points)

        if typeof(non_missing_scaled) != DataFrame
            nothing
        else
            for i in 1:grid
                data_q = filter(x -> x[rv*"_quantile"] == i, non_missing_scaled)
                q_vec[i] = mean(data_q[:, rv], weights(data_q[:, :weight])) # conditional (weighted) mean
            end
        end

        pcf = dct(inverse_hyperbolic_sine(q_vec ./ correction)) # the inverse hyperbolic sine is good for the tails of the DCT vector -> doesnt overwhelm PCA

        return pcf

    elseif typeof(estimator) <: SeriesEstimator
        sort!(non_missing_scaled, rv)

        t_rv = inverse_hyperbolic_sine(non_missing_scaled[:, rv] ./ correction)
        coefs = series_estimator(t_rv, non_missing_scaled[:, :weight], grid - 1) # -1 because the order starts at zero and otherwise, we'd have one extra coefficient

        return coefs
    end
end

function legendre_polynomials(x, order)
    P = [1, x]
    for n in 2:order
        Pn = ((2n - 1) * x * P[end] - (n - 1) * P[end-1]) / n
        push!(P, Pn)
    end
    return P
end


function series_estimator(data, weights, order)
    # Estimate Phi using the legendre polynomials
    Phi = zeros(length(data), order + 1)

    s_weights = cumsum(weights) / sum(weights)

    for i in eachindex(data)
        for j in 0:order
            Phi[i, j+1] = Q_m(j, s_weights[i]) # Q_m places [0,1] data internally to [-1, 1]
        end
    end

    # Incorporate weights into the Phi matrix
    W = Diagonal(sqrt.(weights))  # Use square root of weights for correct weighting
    Phi_weighted = W * Phi

    # coefficients = (transpose(Phi_weighted) * Phi_weighted) \ (transpose(Phi_weighted) * (W * data))

    # More stable and efficient way to estimate the coefficients
    coefficients = zeros(order + 1)
    for j in 1:order+1
        for i in 1:length(data)
            coefficients[j] += weights[i] * Phi[i, j] * data[i]
        end
        coefficients[j] /= sum(weights)
    end

    return coefficients

end

# Define the orthonormal shifted Legendre polynomials as in the paper 
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



# Define the Legendre polynomial of degree m
function Q_m(m, x)
    L_m = legendre_polynomial(m, 2x - 1)

    return sqrt(2m + 1) * L_m
end



# Estimate the copula coefficients ρ_m
function estimate_phi(R, W, m)
    n, d = size(R)
    rho_m = 0.0

    if all(iszero, m)
        return 1.0
    end

    # If d-1 elements of m are zero, then the product is zero
    if sum(m .== 0) >= d - 1
        return 0.0
    end

    # Threads.@threads  # removed: called from within threaded bootstrap loop
    for i in 1:n
        product = 1.0

        for j in 1:d
            product *= Q_m(m[j], R[i, j])
        end

        rho_m += W[i] * product
    end

    return rho_m / sum(W)
end

# Helper function to compute ranks
function rankdata(a)
    order = sortperm(a)
    ranks = similar(order)
    ranks[order] .= 1:length(a)
    return ranks
end

# Construct the N-th order estimator for copula density
function get_copula_coefficients(X, W, N)
    d = size(X, 2)

    ranges = [(0:N) for j in 1:d] # TODO: can be made more flexible
    cl = length(collect(Iterators.product(ranges...)))
    c_N = 0.0
    rho_m = zeros(cl)

    # Threads.@threads  # removed: called from within threaded bootstrap loop
    for (xx, m) in collect(enumerate(Iterators.product(ranges...)))
        tup_m = Tuple(m)
        rho_m[xx] = estimate_phi(X, W, tup_m)
    end

    return rho_m
end

# Vectorized copula coefficient estimation using BLAS matrix operations.
# Precomputes basis matrices once, then uses matrix multiplies instead of
# looping over each multi-index with separate N_obs passes.
# For D=3, N=11: 12 BLAS dgemm calls vs 26M scalar Q_m evaluations.
function get_copula_coefficients_fast(R_input, W, N)
    R = R_input isa Matrix ? R_input : Matrix{Float64}(R_input)
    n, d = size(R)
    K = N + 1  # number of basis functions per dimension

    # Step 1: Precompute basis matrices Φ_d[i, k] = Q_m(k-1, R[i, d])
    Phi = [Matrix{Float64}(undef, n, K) for _ in 1:d]
    for dd in 1:d
        for i in 1:n
            @inbounds for k in 1:K
                Phi[dd][i, k] = Q_m(k - 1, R[i, dd])
            end
        end
    end

    # Step 2: Weighted first-dimension basis: WΦ₁[i, k] = w[i] / sum(W) * Φ₁[i, k]
    w_sum = sum(W)
    WPhi1 = Phi[1] .* (W ./ w_sum)

    if d == 2
        # ρ[k1, k2] = WΦ₁' * Φ₂  — single BLAS dgemm call
        rho_mat = WPhi1' * Phi[2]
        rho_m = vec(rho_mat)

    elseif d == 3
        # For each k3: ρ[:,:,k3] = (WΦ₁ .* Φ₃[:,k3])' * Φ₂
        rho_arr = Array{Float64}(undef, K, K, K)
        for k3 in 1:K
            WPhi1_k3 = WPhi1 .* @view(Phi[3][:, k3])
            rho_arr[:, :, k3] = WPhi1_k3' * Phi[2]
        end
        rho_m = vec(rho_arr)
    else
        error("get_copula_coefficients_fast: only d=2 and d=3 supported, got d=$d")
    end

    # Step 3: Apply special cases
    # (0,0,...,0) → 1.0 (normalization)
    rho_m[1] = 1.0

    # Marginal coefficients (d-1 indices are zero, one is non-zero) → 0.0
    for (xx, m) in enumerate(Iterators.product([(0:N) for _ in 1:d]...))
        if sum(m .== 0) >= d - 1 && !all(iszero, m)
            rho_m[xx] = 0.0
        end
    end

    return rho_m
end


inverse_hyperbolic_sine(x) = log.(x .+ sqrt.(x .^ 2 .+ 1))
scale_to_interval(x_values) = 2 * (x_values .- minimum(x_values)) / (maximum(x_values) - minimum(x_values)) .- 1



function generate_percentiles(rv, weight, method, grid_size)
    # Estimate cdf and generate percentiles 
    obs = unique(rv)
    local rv_cdf, U
    if method == "ecdf"
        rv_cdf = ecdf(rv, weights=weight) # #TODO: technically i should not only do unique above, i should aggregate the weights per unique observation as well

    elseif method == "kde"
        # h      = silverman_bandwidth(rv, weight)
        # x_grid = minimum(rv):2000:maximum(rv)
        # rv_cdf = kernel_distribution(rv, x_grid, h, weight)
        U = kde(rv, weights=weight, npoints=8192)
    end

    # Percentiles 
    cdf_obs = method == "ecdf" ? map(x -> rv_cdf(x), sort(obs)) : cumsum(pdf(U, sort(obs))) ./ sum(pdf(U, sort(obs)))

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
    if any(pmf .< 0)
        # println("Warning: PMF has negative values")
        # println(pmf[pmf .< 0])

        # scale s.t. it is positive
        pmf = pmf .+ abs(minimum(pmf))

        # scale s.t. it sums to 1
        pmf = pmf ./ sum(pmf)
    end

    # Σ_p = sum(pmf)
    # if Σ_p ≉ 1
    #     println(Σ_p)
    #     if round(Σ_p, digits=4) == 1
    #         pmf = pmf ./ sum(pmf)
    #     else
    #         println("Warning: PMF does not sum to 1")
    #     end
    # end
    d = DiscreteNonParametric(sort(obs), pmf)

    # Sampling from the distribution object # out = rand(d,100) # Test it: plot(sort(out))
    # Generate percentile functions 
    # grid      = collect(0.05:1/grid_size:.95)
    # grid[end] = .99
    rv_vals = zeros(length(grid_size), 1)

    # Inverse transform 
    for (i, v) in enumerate(grid_size)
        rv_vals[i] = quantile(d, v)  # Gives you the income value for some percentile from distribution "d" 
    end

    return rv_vals
end


function chebyshev_nodes(n)
    return [0.5 * (1 + cos((2i - 1) * π / (2n))) for i in n:-1:1]
end


function perform_percentile_corrections!(df, rv, df_name, estimator, pcf=false, correction=false)
    non_missing = filter(rv => !isnan, df)
    non_missing = coalesce.(non_missing, NaN)
    filter!("weight" => !isnan, non_missing)

    @unpack grid_pcf = estimator

    grid_points = zeros(grid_pcf) # unused by series estimator

    if typeof(estimator) <: HistogramEstimator
        @unpack grid_type_pcf, = estimator
        if grid_type_pcf == "uniform"
            interval = 1 / grid_pcf
            grid_points = collect(interval:interval:1) # in the series estimator, basically we choose some grid, could be fine or not 

        # # To avoid the upper bound
        # if grid_pcf <= 20
        #     grid_points[end]          = .99
        # else
        #     grid_points[end]          = typeof(estimator) <: SeriesEstimator || typeof(estimator) <: KernelEstimator ? (1 + grid_points[end-1]) / 2 : 1
        # end

        elseif grid_type_pcf == "chebyshev"
            grid_points = nodes(grid_pcf, :chebyshev_nodes, [1.0, 0.0]).points #TODO: for pcf, many chebyshev nodes create a mass near end points -> not good for approximation
        end
    end

    # Sort columns by the measure and id  and creating column of running weights 
    if df_name != "SCF"
        if pcf != false
            return treat_quantile_functions(non_missing, rv, grid_pcf, grid_points, correction, estimator)
        end
    elseif df_name == "SCF"
        pcf_I = zeros(grid_pcf, 5)

        # # Must loop over all implicates 
        if pcf != false
            for i in 1:5
                non_missingᵢ = filter(row -> row.impnum == i, non_missing)
                pcf_I[:, i] .= treat_quantile_functions(non_missingᵢ, rv, grid_pcf, grid_points, correction, estimator)

                # return treat_quantile_functions(non_missingᵢ, rv, grid_pcf, grid_points, correction, estimator)
            end
        end
        # Return the mean of the means
        pcf = mean(pcf_I, dims=2)

        return pcf
    end
end


function assign_quantile_groups_for_bootstrap(df, rv, estimator)

    @unpack grid_pcf, grid_type_pcf = estimator

    grid_points = select_grid_points(grid_pcf)

    non_missing = filter(rv => !isnan, df)
    non_missing = coalesce.(non_missing, NaN)
    filter!("weight" => !isnan, non_missing)

    # Sort columns by the measure and id  and creating column of running weights 
    sort!(non_missing, [Symbol.(rv), :id])
    non_missing[!, :running_weight] = cumsum(non_missing[!, :weight])
    non_missing[!, :running_weight] = non_missing[!, :running_weight] ./ maximum(non_missing[!, :running_weight])

    # Generate new variable 
    non_missing[:, rv*"_quantile"] .= NaN

    for d in 1:grid_pcf
        upper_cond = (non_missing[:, :running_weight] .<= grid_points[d])
        if d == 1
            non_missing[upper_cond, rv*"_quantile"] .= d
        elseif d > 1 && d <= grid_pcf
            lower_cond = (non_missing[!, :running_weight] .> grid_points[d-1])
            non_missing[(lower_cond.&upper_cond).||(non_missing[:, :running_weight].>grid_points[end]), rv*"_quantile"] .= d
        end
    end

    # Drop 'running_weight' to not interfere with subsequent runs (will be different each time if missings is measure dependent)
    select!(non_missing, Not(:running_weight))

    return non_missing, flag
end



function clean_data!(df, measures, equivalized, bottom_coded, blind_to, key)
    # First check if data has the measure. Subset to the columns of choice
    local non_missing_cols
    if all([col in names(df) for col in measures])
        @info("All columns in measures are present in the $key")
        non_missing_cols = copy(measures)
    else
        # Identify the non-missing columns 
        @info("At least one column in measures is not present in the $key")
        non_missing_cols = [col for col in measures if col ∈ names(df)]
    end

    if length(blind_to) != 0 && key ∈ keys(blind_to)
        meas_to_blind = blind_to[key]

        @info("We are blind to $(blind_to[key]) in the $(key)")

        # remove the strings from the vector
        non_missing_cols = filter(x -> !(x in meas_to_blind), non_missing_cols)
    end

    # For the random-perturbantions later 
    Random.seed!(123456789)

    # Minor cleaning  
    # relevant_cols = setdiff(names(df), ["year", "quarter", "month", "time", "id"])

    for c in non_missing_cols
        replace!(df[!, c], missing => NaN)

        # Convert to Float for future calculations   
        df[!, c] = convert.(Float64, df[!, c])

        scaled_perturbations = rand(length(df[:, c])) .* 0.01

        # Equivalize step 
        if c != "weight" && c != "hhequiv"
            if equivalized == true
                df[!, c] = (df[:, c] ./ df[:, "hhequiv"]) .+ scaled_perturbations
            else
                # Perturb data randomly a bit to avoid ties -> uniform distribution 
                df[!, c] = df[:, c] + scaled_perturbations
            end
        end
    end

    # Bottom code step
    if !isempty(bottom_coded)
        for c in bottom_coded
            df[df[!, c].<0, c] .= 0
            @assert minimum(df[:, c][(!isnan).(df[:, c])]) >= 0.0
        end
    end

    return non_missing_cols
end




