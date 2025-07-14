# cd("/Users/lc/Dropbox/Distributional_Dynamics/5_Code")
include("DistributionalDynamics.jl")


function run_interval_estimation(obs_data::ObservedData, method_options::MethodOptions, intervals_for=:all)

    @unpack grid, measures, lags, freq, agg_freq, case, pca_perspective, pca_method, plot_proof, rm_seasonality, data_cutoffs = method_options
    @unpack files, agg_data, df_vec, gdp_series = obs_data

    number_of_dfs = length(df_vec.data)
    init_path = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) * "/2_Data_processing" : pwd()
    dimension = length(measures)

    @info("Taking raw data and placing them in the correct containers, sorted in alphabetical order of datasets' name.")

    # Collects data for all available years 
    dfs, year_vec, time_dict, freq_type = data_constructor(obs_data, method_options)

    # Using time parameters from data and from aggregates, we find an agreeable time frame for the estimation
    time_p = define_timeframe(agg_data, year_vec, freq_type, time_dict, data_cutoffs) # TODO: time_dict needs to be filtered as well  

    cutoff_bounds = align_data_with_timeframe!(dfs, year_vec, time_p, data_cutoffs)

    # Define data intervals for the specific dataset or all 
    confidence_intervals = define_data_intervals_only(df_vec, method_options, init_path, time_p, obs_data, intervals_for)

    for j in collect(keys(confidence_intervals))
        # if occursin("SCF", j)
        #     k = findall(x -> x == "SCF", data_sources)[1]
        # else
        k = findall(x -> x == j, df_vec.df_names)[1]
        # end
        confidence_intervals[j] = interval_time_correction(confidence_intervals[j], year_vec[k], time_p, k, cutoff_bounds[k])

        init_path = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
        save_path = init_path * "/confidence_intervals/ci_deciles_$j" * "npimp.jld2"
        JLD2.save(save_path, "ci", confidence_intervals[j])
    end
end


function define_data_intervals_only(df_vec, method_options, init_path, time_p, obs_data, intervals_for)
    @unpack measures, equivalized, bottom_coded, blind_to, grid, errors_process, scf_ci, freq = method_options
    @unpack time_dict, freq_type, year_vec, tmin, tmax, tot_periods = time_p
    @unpack gdp_series, agg_data = obs_data

    id_for_intervals = Dict()
    new_df_dict = Dict()
    for d in intervals_for
        id_for_intervals[d] = findall(x -> x == d, df_vec.df_names)[1]
        new_df_dict[d] = df_vec.data[id_for_intervals[d]]
    end

    # For labels and tagging 
    sources = copy(intervals_for) #df_vec.df_names  # e.g., SCF, PSID, etc.
    data_label = data_tag(sources)
    grid_granularity = identify_grid(grid)
    m_label = measures_folder(measures)
    dimension = length(measures)
    ci_exists = isfile(init_path * "/confidence_intervals/ci_" * m_label * grid_granularity * "_" * data_label * scf_ci * ".jld2")
    noise_exists = isfile(init_path * "/noise_distributions/noise_" * m_label * grid_granularity * "_" * data_label * scf_ci * ".jld2")

    gdp_series[!, "date"] = QuarterlyDate.(gdp_series[!, "time"])
    filter!(row -> row.date >= QuarterlyDate(tmin["year"], tmin["quarter"]), gdp_series)
    filter!(row -> row.date <= QuarterlyDate(tmax["year"], tmax["quarter"]), gdp_series)

    # For the SCF
    interval_tags = ["npimp", "np"]

    # For the estimation of measurement and process noise 
    DCT_boot = Dict()

    if ci_exists && noise_exists
        # Read in confidence intervals
        f = jldopen(init_path * "/confidence_intervals/ci_" * m_label * grid_granularity * "_" * data_label * scf_ci * ".jld2", "r")
        vw = jldopen(init_path * "/noise_distributions/noise_" * m_label * grid_granularity * "_" * data_label * scf_ci * ".jld2", "r")
        return f["ci"], vw["noise"]
    else
        objects = sort(["quantiles", "levels", "shares"])

        # Define series for the confidence intervals --- are based on the estimation results 
        local series
        if grid == 10 || grid == 100 || grid == 20
            series = ["bottom50", "next40", "top10"]
        elseif grid == 5
            series = ["bottom40", "next40", "top20"]
        end

        # Loop over data sources and generate intervals
        confidence_intervals = Dict()
        draws = 999
        for source in collect(keys(new_df_dict))
            # if source != "SCF"
            confidence_intervals[source] = Dict()
            data, _ = select_data(new_df_dict[source], measures, equivalized, bottom_coded, blind_to, source)

            # Generate non-parametric confidence intervals 
            j = id_for_intervals[source]
            sub_boot_dict, DCT_boot[source] = bootstrap_inverse_transform!(data, objects, series, year_vec[j], time_dict[j], freq_type[j], grid, measures, draws, gdp_series, source, scf_ci)

            # only keep the lower and upper bounds
            confidence_intervals[source]["ci_u"], confidence_intervals[source]["ci_l"] = construct_confidence_intervals(sub_boot_dict, 0.005, 0.995, objects, measures, length(year_vec[j]), source)
        end
        JLD2.save(init_path * "/confidence_intervals/ci_" * m_label * grid_granularity * "_" * data_label * scf_ci * ".jld2", "ci", confidence_intervals)
        return confidence_intervals
    end
end

run_interval_estimation(obs_data, method_options, ["CEX"])
