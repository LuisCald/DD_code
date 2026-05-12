# Forecasting
# What to do?
# We take the estimated series, using all of the data, etc and then we forecast the series using the estimated parameters.
# We then remove a few observations from the estimation sample, but using the same parameters. We then forecast again 
function create_stat_dict(data_sources, measures)
    """This function creates a dictionary that stores the MSE for each observation that is removed. """
    # Create dictionary to store results 
    stat_dict = Dict()

    # Loop through all data sources 
    for source in data_sources
        stat_dict[source] = Dict()
        for meas in measures
            stat_dict[source][meas] = Dict()
        end
        stat_dict[source]["copula"] = Dict()
    end
    return stat_dict
end


function perform_forecast(dataset_choice, par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, filtering_criteria, obs_data, model_options, model_elements, time_params, user_params, func_data, type, forecast_type)
    if forecast_type[1] == "iterative"
        return perform_iterative_forecast(dataset_choice, par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, obs_data, model_options, model_elements, time_params, func_data, type)

    elseif forecast_type[1] == "clumpy"
        perform_clumpy_forecast(dataset_choice, par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, filtering_criteria, obs_data, model_options, model_elements, time_params, user_params, func_data, type)

    elseif forecast_type[1] == "extensive"
        perform_extensive_forecast(dataset_choice, par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, filtering_criteria, obs_data, model_options, model_elements, time_params, user_params, func_data, type, forecast_type)
    end
end

function perform_iterative_table_forecast(par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, obs_data, model_options, model_elements, time_params, func_data, statistic)
    """This loops through all observations in some dataset and removes the corresponding one observation, then estimates the series again. """

    @unpack MV = model_elements
    @unpack data_sources, func_dict = func_data
    @unpack case, measures, std_method, grid = model_options

    dimension = length(measures)
    label = "$dimension" * "D" * "_$case"
    orig_MV = deepcopy(MV) # for later 


    # First, get the original reconstruction 
    smoother_output, _ = likeli(model_elements, par_final, param_sizes, hyperpriors, Σ_ids, model_options; smooth=true)
    @unpack x_smoothed = smoother_output
    @unpack proj = model_elements
    dv, _ = reconstruct_data(par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources)

    # Compute the statistic of choice 
    local stat_dict_r
    if statistic == :mse
        stat_dict_r = compute_forecast_mse(dv, func_dict, model_elements, orig_MV, MV, data_sources, measures, grid, dimension, par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, obs_data, model_options, time_params)
    end

    # Try to present everything in a table format for latex, where we have 4 columns: release date, copula, income pcf, wealth pcf 
    @unpack time_dict = time_params

    # Create a table for each dataset
    full_stat_dict = Dict()

    # Create a column for each object 
    for (j, k) in enumerate(data_sources)
        println(k)
        stat_df = DataFrame()

        # create date column
        rel_vec = []
        for y in sort!(collect(keys(time_dict[j])))
            for i in time_dict[j][y]
                push!(rel_vec, "$y" * "Q" * "$i")
            end
        end
        stat_df[!, :release_date] = rel_vec
        println(length(rel_vec))
        println(length(stat_dict_r[k]["copula"]))

        # Create copula column 
        stat_df[!, "copula"] = stat_dict_r[k]["copula"]

        # Create columns for the pcfs
        for m in measures
            stat_df[!, m] = stat_dict_r[k][m]
        end

        # Add to the full dataframe
        full_stat_dict[k] = stat_df
    end

    return full_stat_dict
end


function compute_forecast_mse(dv, func_dict, model_elements, orig_MV, MV, data_sources, measures, grid, dimension, par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, obs_data, model_options, time_params)
    N(x) = length(x)

    # Create dictionary to store results of the MSE for the "all-data" case
    stat_dict = create_stat_dict(data_sources, measures)
    cop_ids = Dict()
    pcf_ids = Dict()
    id_sets = [I for I in Iterators.partition(1:grid*dimension, grid)]

    data_cop = Dict()
    pcfs = Dict()

    for (j, k) in enumerate(sort!(data_sources)) # TODO: includes consensus 
        println(k)
        # Get indices of the observations 
        pcfs[k] = Dict()
        pcf_ids[k] = Dict()
        cop_ids[k] = findall(x -> all(.!isnan.(x)), eachcol(func_dict[k]["copulas"]["data"]))
        data_cop[k] = func_dict[k]["copulas"]["data"][:, cop_ids[k]]

        for (n, m) in enumerate(measures)
            # Calculate squared errors for pcfs => reconstruction minus data 
            try
                pcf_ids[k][m] = findall(x -> !all(isnan.(x)), eachcol(func_dict[k][m]["quantiles"]["data"]))
            catch ee
                println(ee)
                pcf_ids[k][m] = []
            end
            pcfs[k][m] = func_dict[k][m]["quantiles"]["data"][:, pcf_ids[k][m]]
            stat_dict[k][m] = N(pcf_ids[k][m]) != 0 ? sum((dv[k][2][id_sets[n], pcf_ids[k][m]] .- pcfs[k][m]) .^ 2, dims=1) : NaN
            # stat_dict[k][m] = N(pcf_ids[k][m]) != 0 ? sum((dv[k][2][id_sets[n], pcf_ids[k][m]] .- pcfs[k][m]).^2) / N(pcf_ids[k][m]) : NaN
        end

        # Calculate TOTAL squared errors for copula => reconstruction minus data 
        # stat_dict[k]["copula"]  = N(cop_ids[k]) != 0 ? sum((reshape(dv[k][1][:, :, cop_ids[k]], (grid^dimension, N(cop_ids[k]))) .- data_cop[k]).^2) / N(cop_ids[k]) : NaN 
        stat_dict[k]["copula"] = N(cop_ids[k]) != 0 ? sum((reshape(dv[k][1][:, :, cop_ids[k]], (grid^dimension, N(cop_ids[k]))) .- data_cop[k]) .^ 2, dims=1) : NaN
    end

    # Now remove one observation at a time, for each dataset, calculate the MSE and store it in the dictionary
    stat_dict_r = Dict()

    for (j, k) in enumerate(sort!(data_sources)) # TODO: includes consensus 
        stat_dict_r[k] = Dict()
        indices = findall(x -> !all(isnan.(x)), eachcol(MV[j]))
        # indices         = pcf_ids[k][m]

        # Fill Containers
        stat_dict_r[k]["copula"] = Vector{Float64}(undef, N(indices))

        for (n, m) in enumerate(measures)
            stat_dict_r[k][m] = Vector{Float64}(undef, N(indices))
        end

        # for each observation 
        # indices = findall(x -> !all(isnan.(x)), eachcol(MV[j]))

        for (t, id) in enumerate(indices)
            MV = deepcopy(orig_MV)
            # Remove the observation in MV for this dataset 
            MV[j][:, id] .= NaN
            @pack! model_elements = MV

            # Reconstruction with less data
            dv_r, _ = reconstruct_data(par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources)

            # Calculate squared errors for each object relative to total MSE => reconstruction minus data
            # cop_recon                   = reshape(dv_r[k][1][:, :, cop_ids[k]], (grid^dimension, N(cop_ids[k])))
            local cop_recon
            try
                cop_recon = reshape(dv_r[k][1][:, :, cop_ids[k][t]], (grid^dimension, 1))
            catch ee
                cop_recon = []
            end
            # stat_dict_r[k]["copula"][i] = N(cop_ids[k]) != 0 ? 100 * (1 - ((sum((cop_recon .- data_cop[k]).^2) / N(cop_ids[k])) / stat_dict[k]["copula"])) : NaN
            # stat_dict_r[k]["copula"][t] = N(cop_ids[k]) != 0 ? 100 * (1 - (sum((cop_recon .- data_cop[k]).^2) / stat_dict[k]["copula"])) : NaN
            stat_dict_r[k]["copula"][t] = cop_recon != [] ? 100 * (1 - (sum((cop_recon .- data_cop[k][:, t]) .^ 2) ./ stat_dict[k]["copula"][1, t])) : NaN

            # We only look at quantiles
            for (n, m) in enumerate(measures)
                # stat_dict_r[k][m][i] = N(pcf_ids[k][m]) != 0 ? 100 * (1 - ((sum((dv_r[k][2][id_sets[n], pcf_ids[k][m]] .- pcfs[j][n]).^2) / N(pcf_ids[k][m])) / stat_dict[k][m])) : NaN
                try
                    stat_dict_r[k][m][t] = N(pcf_ids[k][m]) != 0 ? 100 * (1 - (sum((dv_r[k][2][id_sets[n], pcf_ids[k][m][t]] .- pcfs[k][m][:, t]) .^ 2) ./ stat_dict[k][m][1, t])) : NaN # TODO: check first N()
                catch ee
                    stat_dict_r[k][m][t] = NaN #TODO: issue is I use all indices
                end
            end
        end
    end
    return stat_dict_r
end



function perform_clumpy_forecast(dataset_choice, par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, how_much, obs_data, model_options, model_elements, time_params, user_params, func_data, type)
    """This function only removes the last 3 observations, all measures, from the SCF. """
    # Take estimated series that used all data 
    @unpack MV = model_elements
    @unpack data_sources = func_data
    @unpack case, measures, std_method, grid = model_options

    dimension = length(measures)
    label = "$dimension" * "D" * "_$case"
    orig_MV = deepcopy(MV) # for later 


    # First, get the original reconstruction 
    smoother_output, _ = likeli(model_elements, par_final, param_sizes, hyperpriors, Σ_ids, model_options; smooth=true)
    @unpack x_smoothed = smoother_output
    factors_all_obs = deepcopy(x_smoothed)
    dv, _ = reconstruct_data(par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources)

    # Now remove the last few observations from the estimation sample
    # obs_removed = Vector{Matrix{Float64}}(undef, 1)

    # Find the "j" that corresponds to the SCF 
    data_id = findfirst(x -> x == dataset_choice, data_sources)

    # Find the indices in MV where at least something is observed 
    indices = findall(x -> !all(isnan.(x)), eachcol(MV[data_id]))

    # Retain the observations that we will remove
    # obs_removed[1] = deepcopy(MV[data_id][:, indices[end-how_much:end]])

    # Set the last X wealth observations in SCF to missing
    MV[data_id][:, indices[end-how_much:end]] .= NaN

    # Define new measurement vector 
    @pack! model_elements = MV

    # Estimate the new series 
    # smoother_output, _               = likeli(model_elements, par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, model_options, true)
    # @unpack x_smoothed               = smoother_output
    # @unpack proj                     = model_elements
    # @unpack means, stds, agg_count   = model_elements

    # Reconstruction with less data
    dv_r, _ = reconstruct_data(par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources)

    # Repack the original MV
    MV = orig_MV
    @pack! model_elements = MV
    # Compare construction with less data with original estimation 
    @info("Plotting forecasts")
    setup_forecasts(dv, dv_r, obs_data, func_data, time_params, user_params, model_options, type, label, "clumpy")

end

function perform_extensive_forecast(dataset_choice, par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, filtering_criteria, obs_data, model_options, model_elements, time_params, user_params, func_data, type, forecast_type)
    """This forecast drops all wealth observations from all data that correspond to the last X quarters."""

    # Take estimated series that used all data 
    @unpack MV, y = model_elements
    @unpack data_sources = func_data
    @unpack case, measures, estimator = model_options
    @unpack integral_pcf_grid, integral_cop_grid = estimator
    @unpack tmin, tmax, time_dict = time_params
    @unpack df_vec = obs_data

    dimension = length(measures)
    label = "$dimension" * "D" * "_$case"
    orig_y = deepcopy(y) # for later 


    # First, get the original reconstruction 
    smoother_output, _ = likeli(model_elements, par_final, param_sizes, hyperpriors, Σ_ids, model_options; smooth=true)
    @unpack x_smoothed = smoother_output
    # factors_all_obs     = deepcopy(x_smoothed)
    dv, _ = reconstruct_data(par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources)

    # find indices related to wealth 
    local data_indices
    if forecast_type[2] == "income_only"
        data_indices = find_wealth_indices(measures, grid) # needs to be changed s.t. it refers to 'y'
    elseif forecast_type[2] == "data_only"
        # Find index associated with the given 'dataset_choice'
        data_id = findall(x -> x == dataset_choice, df_vec.df_names)[1]

        # Based on grid, find the indices associated with the dataset
        n = size(MV[1], 1)
        data_indices = collect((data_id-1)*n+1:data_id*n)
        println((data_id-1)*n+1:data_id*n)
    else
        n = size(y, 1)
        data_indices = collect(1:n) #weird, but necessary 
    end

    local start_ind, end_ind
    if collect(keys(filtering_criteria)) == ["dates"]
        # Find indices associated with the dates 
        dates_of_est = QuarterlyDate(tmin["year"], tmin["quarter"]):Quarter(1):QuarterlyDate(tmax["year"], tmax["quarter"])
        start_ind = findall(dates_of_est .== QuarterlyDate(filtering_criteria["dates"][1]))[1]
        end_ind = findall(dates_of_est .== QuarterlyDate(filtering_criteria["dates"][2]))[1]
        println(start_ind, end_ind)
        y[data_indices, start_ind:end_ind] .= NaN # removes X quarters of data, independent of which dataset it is or object 
    elseif collect(keys(filtering_criteria)) == ["periods"]
        T = size(y, 2)
        how_much = filtering_criteria["periods"]
        start_ind = T - (how_much - 1)
        end_ind = T
        y[data_indices, start_ind:end_ind] .= NaN # removes X quarters of data, independent of which dataset it is or object 
    elseif collect(keys(filtering_criteria)) == ["periods_to_remove"]
        dates_of_est = QuarterlyDate(tmin["year"], tmin["quarter"]):Quarter(1):QuarterlyDate(tmax["year"], tmax["quarter"])

        # Find indices associated with the dates
        all_indices = []
        for period in filtering_criteria["periods_to_remove"]
            ind = findall(dates_of_est .== period)[1]
            push!(all_indices, ind)
        end
        # Sort 
        all_indices = sort(all_indices)

        # Mute observations associated with dataset 
        y[data_indices, all_indices] .= NaN # removes X quarters of data, independent of which dataset it is or object 
    end
    # Define new measurement vector 

    @pack! model_elements = y

    # Estimate the new series 
    # smoother_output, _               = likeli(model_elements, par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, model_options, true)
    # @unpack x_smoothed               = smoother_output
    # @unpack proj                     = model_elements
    # @unpack means, stds, agg_count   = model_elements

    # Reconstruction with less data
    dv_r, _ = reconstruct_data(par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources)


    # TODO: Forecast the series using the estimated parameters, but check whether we have aggs that run further enough 

    # Quickly calculate the MSE for the factors and observations 
    # T = size(x_smoothed, 2)

    # Calculate squared errors for each data point
    # squared_errors = (factors_all_obs .- x_smoothed).^2

    # Calculate the sum of squared errors
    # sum_squared_errors = sum(squared_errors)

    # Calculate the Mean Squared Error: distance from 
    # mse1 = sum_squared_errors / T

    # Now, compute the difference between obs_removed and the estimated series 
    # recon_data =  (proj * x_smoothed)[:, end-(how_much-1):end]

    # # Split by object, multiply by stds, reform object 
    # ΓF_σ = add_variance(recon_data, stds, grid, dimension)
    # X    = add_mean!(X, ΓF_σ, means, data_names, blind_to)

    # # Calculate squared errors for each data point
    # y              = reduce(vcat, obs_removed)
    # y_r            = reduce(vcat, X)
    # squared_errors = (y .- y_r).^2

    # # Calculate the sum of squared errors
    # sum_squared_errors = sum(squared_errors)

    # # Calculate the Mean Squared Error
    # mse2 = sum_squared_errors / T

    # Compare construction with less data with original estimation 
    @info("Plotting forecasts")
    setup_forecasts(dataset_choice, dv, dv_r, "normal", obs_data, func_data, time_params, user_params, model_options, type, label, forecast_type, filtering_criteria)

    # Pack the orig MV
    y = deepcopy(orig_y)
    @pack! model_elements = y
end

function find_grid_points_for_series(grid)
    if grid == 10
        return [[1:5]..., [6:9]..., 10]
    elseif grid == 20
        return [[1:10]..., [11:18]..., [19:20]...]
    elseif grid == 5
        return [[1:2]..., [3:4]..., [5]...]
    end
end

function find_wealth_indices(measures, grid)
    D = length(measures)
    diff = grid + (D - 1) * (grid - 1)

    cop_n = grid^D - diff
    wealth_id = findfirst(x -> x == "wealth", measures)
    wealth_pcf_id = cop_n+(wealth_id-1)*grid+1:cop_n+wealth_id*grid

    return [1:cop_n..., wealth_pcf_id...]
end


function setup_forecasts(dataset_choice, dv, dv_r, ty, obs_data, func_data, time_params, user_params, model_options, type, label, forecast_type, filtering_criteria)
    """Plot forecasts for a select number of series."""
    # Place all data into dictionaries of the same structure 
    avg_series = Dict()
    avg_series["recon"] = Dict()
    avg_series["missing"] = Dict()

    # Create New dictionaries with dataset_choice only 
    dv2 = Dict(dataset_choice => dv[dataset_choice])
    dv_r2 = Dict(dataset_choice => dv_r[dataset_choice])


    for k in [dataset_choice]
        dv2[k][ty], avg_series["recon"][k] = export_functional_data(dv2[k][ty], ty, k, type, obs_data, func_data, time_params, user_params, model_options, false, false, :forecast)
    end

    for k in [dataset_choice]
        dv_r2[k][ty], avg_series["missing"][k] = export_functional_data(dv_r2[k][ty], ty, k, type, obs_data, func_data, time_params, user_params, model_options, false, false, :forecast)
    end


    # Plot the forecasts
    plot_forecasts(dv2, dv_r2, ty, func_data, obs_data, avg_series, user_params, time_params, model_options, type, label, forecast_type[1], filtering_criteria)
end


function plot_forecasts(data_dt, data_dt_r, ty, func_data, obs_data, avg_series, user_params, time_params, model_options, type, label, forecast_type, filtering_criteria)
    @unpack measures, lags, tag = model_options
    @unpack estimator = model_options
    @unpack integral_pcf_grid, integral_cop_grid = estimator
    @unpack tmin, tmax, time_dict = time_params
    @unpack func_dict, year_vec, confidence_intervals = func_data
    @unpack gdp_series, df_vec = obs_data


    data_sources = setdiff(sort(collect(keys(data_dt))), ["consensus"])
    smin, smax = user_params
    m_label = measures_folder(measures)
    dates = Dict()
    base_jump, end_jump = find_subset_frame(smin, smax, tmin, tmax)

    # Path situation
    folder = type == :optimization ? "from_optimization" : "from_mcmc"
    init_path = BASE_PATH
    path = init_path * "/7_Results/$m_label" * tag * "/$folder/plots/"
    mkpath(path)

    # Common annual series 
    top = Dict("consensus" => Dict())
    bot = Dict("consensus" => Dict())
    mid = Dict("consensus" => Dict())

    # Same, but for the forecast 
    top_r = Dict("consensus" => Dict())
    bot_r = Dict("consensus" => Dict())
    mid_r = Dict("consensus" => Dict())

    for dt in [top, bot, mid, top_r, bot_r, mid_r]
        for source in data_sources
            dt[source] = Dict()
            dt[source*"_obs"] = Dict()
            dt[source]["lb"] = Dict()
            dt[source]["ub"] = Dict()
            for obj in ["shares", "quantiles"]
                dt[source][obj] = Dict()
                dt[source*"_obs"][obj] = Dict()
                dt[source]["lb"][obj] = Dict()
                dt[source]["ub"][obj] = Dict()
                dt["consensus"][obj] = Dict()
            end
        end
        for obj in ["shares", "quantiles"]
            dt["consensus"][obj] = Dict()
        end
    end

    # Defining the series 
    series = Dict()
    if integral_pcf_grid == 10 || integral_pcf_grid == 20 || integral_pcf_grid == 100
        series["bot"] = "bottom50"
        series["mid"] = "next40"
        series["top"] = "top10"
    elseif integral_pcf_grid == 5
        series["bot"] = "bottom40"
        series["mid"] = "next40"
        series["top"] = "top20"
    end

    gps = find_grid_points_for_series(integral_pcf_grid)

    # Trends 
    # gdp_series[!, "date"]    = QuarterlyDate.(gdp_series[!, "time"])
    # correction_names         = [meas * "_per_hh" for meas in measures]
    # select_series            = select(gdp_series, ["date", correction_names...])

    dates = Dict()

    for meas in measures
        M = uppercasefirst(meas)
        M = meas == "consum" ? "Consumption" : M

        # Align dates, create xaxis
        dates[meas] = Dict()
        dates[meas]["tmin"] = Dict()
        dates[meas]["tmax"] = Dict()

        dates[meas]["tmin"]["year"] = smin["year"]
        dates[meas]["tmax"]["year"] = smax["year"]

        dates[meas]["tmin"]["quarter"] = smin["quarter"]
        dates[meas]["tmax"]["quarter"] = smax["quarter"]

        @info(dates[meas])
        for obj in ["quantiles"]
            # Adjust frequencies, Confine to dates
            for source in data_sources

                if obj == "quantiles"
                    # correction first 
                    corrected_ci_l = confidence_intervals[source]["ci_l"][meas][obj] ./ avg_series["recon"][source][base_jump:end-end_jump, meas*"_per_hh"]'
                    corrected_ci_u = confidence_intervals[source]["ci_u"][meas][obj] ./ avg_series["recon"][source][base_jump:end-end_jump, meas*"_per_hh"]'

                    # subset to asked timeframe 
                    top[source]["lb"][obj][meas] = subset_to_cutoff(corrected_ci_l[gps[3], :], dates[meas], tmin, tmax, obj, 1)
                    mid[source]["lb"][obj][meas] = subset_to_cutoff(sum(corrected_ci_l[gps[2], :], dims=1)' ./ length(gps[2]), dates[meas], tmin, tmax, obj, 1)
                    bot[source]["lb"][obj][meas] = subset_to_cutoff(sum(corrected_ci_l[gps[1], :], dims=1)' ./ length(gps[1]), dates[meas], tmin, tmax, obj, 1)

                    top[source]["ub"][obj][meas] = subset_to_cutoff(corrected_ci_u[gps[3], :], dates[meas], tmin, tmax, obj, 1)
                    mid[source]["ub"][obj][meas] = subset_to_cutoff(sum(corrected_ci_u[gps[2], :], dims=1)' ./ length(gps[2]), dates[meas], tmin, tmax, obj, 1)
                    bot[source]["ub"][obj][meas] = subset_to_cutoff(sum(corrected_ci_u[gps[1], :], dims=1)' ./ length(gps[1]), dates[meas], tmin, tmax, obj, 1)

                elseif obj == "shares"
                    top[source]["lb"][obj][meas] = subset_to_cutoff(confidence_intervals[source]["ci_l"][meas][obj][3, :], dates[meas], tmin, tmax, obj)  # TODO: generalize grid
                    mid[source]["lb"][obj][meas] = subset_to_cutoff(confidence_intervals[source]["ci_l"][meas][obj][2, :], dates[meas], tmin, tmax, obj)
                    bot[source]["lb"][obj][meas] = subset_to_cutoff(confidence_intervals[source]["ci_l"][meas][obj][1, :], dates[meas], tmin, tmax, obj)

                    top[source]["ub"][obj][meas] = subset_to_cutoff(confidence_intervals[source]["ci_u"][meas][obj][3, :], dates[meas], tmin, tmax, obj)
                    mid[source]["ub"][obj][meas] = subset_to_cutoff(confidence_intervals[source]["ci_u"][meas][obj][2, :], dates[meas], tmin, tmax, obj)
                    bot[source]["ub"][obj][meas] = subset_to_cutoff(confidence_intervals[source]["ci_u"][meas][obj][1, :], dates[meas], tmin, tmax, obj)
                end

                # For the reconstructions 
                top[source][obj][meas] = subset_to_cutoff(data_dt[source][ty][meas][obj]["common series"][series["top"]], dates[meas], tmin, tmax, obj, avg_series["recon"][source][base_jump:end-end_jump, meas*"_per_hh"])
                bot[source][obj][meas] = subset_to_cutoff(data_dt[source][ty][meas][obj]["common series"][series["bot"]], dates[meas], tmin, tmax, obj, avg_series["recon"][source][base_jump:end-end_jump, meas*"_per_hh"])
                mid[source][obj][meas] = subset_to_cutoff(data_dt[source][ty][meas][obj]["common series"][series["mid"]], dates[meas], tmin, tmax, obj, avg_series["recon"][source][base_jump:end-end_jump, meas*"_per_hh"])

                top_r[source][obj][meas] = subset_to_cutoff(data_dt_r[source][ty][meas][obj]["common series"][series["top"]], dates[meas], tmin, tmax, obj, avg_series["missing"][source][base_jump:end-end_jump, meas*"_per_hh"])
                bot_r[source][obj][meas] = subset_to_cutoff(data_dt_r[source][ty][meas][obj]["common series"][series["bot"]], dates[meas], tmin, tmax, obj, avg_series["missing"][source][base_jump:end-end_jump, meas*"_per_hh"])
                mid_r[source][obj][meas] = subset_to_cutoff(data_dt_r[source][ty][meas][obj]["common series"][series["mid"]], dates[meas], tmin, tmax, obj, avg_series["missing"][source][base_jump:end-end_jump, meas*"_per_hh"])

                # Observations, but only for original reconstruction 
                if obj == "quantiles"
                    top[source*"_obs"][obj][meas] = func_dict[source][meas][obj]["common series"][series["top"]] ./ avg_series["recon"][source][base_jump:end-end_jump, meas*"_per_hh"]
                    bot[source*"_obs"][obj][meas] = func_dict[source][meas][obj]["common series"][series["bot"]] ./ avg_series["recon"][source][base_jump:end-end_jump, meas*"_per_hh"]
                    mid[source*"_obs"][obj][meas] = func_dict[source][meas][obj]["common series"][series["mid"]] ./ avg_series["recon"][source][base_jump:end-end_jump, meas*"_per_hh"]

                    top[source*"_obs"][obj][meas] = subset_to_cutoff(top[source*"_obs"][obj][meas], dates[meas], tmin, tmax, obj, 1)
                    bot[source*"_obs"][obj][meas] = subset_to_cutoff(bot[source*"_obs"][obj][meas], dates[meas], tmin, tmax, obj, 1)
                    mid[source*"_obs"][obj][meas] = subset_to_cutoff(mid[source*"_obs"][obj][meas], dates[meas], tmin, tmax, obj, 1)

                else
                    top[source*"_obs"][obj][meas] = subset_to_cutoff(func_dict[source][meas][obj]["common series"][series["top"]], dates[meas], tmin, tmax, obj)
                    bot[source*"_obs"][obj][meas] = subset_to_cutoff(func_dict[source][meas][obj]["common series"][series["bot"]], dates[meas], tmin, tmax, obj)
                    mid[source*"_obs"][obj][meas] = subset_to_cutoff(func_dict[source][meas][obj]["common series"][series["mid"]], dates[meas], tmin, tmax, obj)
                end

            end

            # Plots of both data, 2 lines per plot (we take the first quarter for the year)
            xaxis = QuarterlyDate(dates[meas]["tmin"]["year"], dates[meas]["tmin"]["quarter"]):Quarter(1):QuarterlyDate(dates[meas]["tmax"]["year"], dates[meas]["tmax"]["quarter"])
            dist_dict = Dict("top" => [top, top_r], "mid" => [mid, mid_r], "bot" => [bot, bot_r])

            @info("Building reconstruction plots, comparing reconstructions to data: $obj")
            # Top
            for (ds, dist) in dist_dict
                for (j, source) in enumerate(data_sources)

                    plot_name = occursin("CEX", source) ? "CEX" : source

                    # Get the R2 of the data on the two Reconstructions
                    cond = .!isnan.(top[source*"_obs"][obj][meas]) # indices between 1-238 where obs. are observed 
                    c_data = dist[1][source*"_obs"][obj][meas][cond]

                    if c_data != []
                        # Actual indices of observations 
                        p_cond = findall(cond .!= 0)
                        start = length(p_cond) == 0 ? 1 : p_cond[1]
                        stop = length(p_cond) == 0 ? 1 : p_cond[end]
                        p_xaxis = collect(xaxis[start:stop])

                        # Estimates from when the first observation begins 
                        p_data = dist[1][source][obj][meas][start:stop]
                        p_data2 = dist[2][source][obj][meas][start:stop]
                        p_cond = cond[start:stop]

                        # Correlate construction with forecast 
                        ρ = round(cor(p_data, p_data2), digits=2)

                        # Find indices corresponding to missing periods
                        local id_for_missing, id_for_data, within_data_missing_ids, within_data_obs_ids, range_missing
                        if collect(keys(filtering_criteria))[1] == "dates"
                            start2 = findall(x -> x == QuarterlyDate(filtering_criteria["dates"][1]), p_xaxis)[1]
                            stop2 = findall(x -> x == QuarterlyDate(filtering_criteria["dates"][2]), p_xaxis)[1]
                            range_missing = collect(start2:stop2)

                            # Find the id of missing for observations only 
                            source_id = findall(x -> x == source, df_vec.df_names)[1]
                            years_df = year.([QuarterlyDate(filtering_criteria["dates"][1]), QuarterlyDate(filtering_criteria["dates"][2])])
                            years_df = collect(years_df[1]:years_df[2])

                            miss_indices = []
                            obs_indices = []
                            within_data_missing_ids = []
                            within_data_obs_ids = []

                            for (j, period) in enumerate(xaxis[cond])
                                ind = findall(x -> x == period, p_xaxis)[1]
                                if year(period) ∈ years_df
                                    push!(miss_indices, ind)
                                    push!(within_data_missing_ids, j)
                                else
                                    push!(obs_indices, ind)
                                    push!(within_data_obs_ids, j)
                                end
                            end

                            id_for_missing = sort(miss_indices)
                            id_for_data = sort(obs_indices)
                            within_data_missing_ids = sort(within_data_missing_ids)
                            within_data_obs_ids = sort(within_data_obs_ids)

                        elseif collect(keys(filtering_criteria))[1] == "periods_to_remove"
                            # Find the id of missing for observations only 
                            source_id = findall(x -> x == source, df_vec.df_names)[1]
                            years_df = year.(filtering_criteria["periods_to_remove"])

                            miss_indices = []
                            obs_indices = []
                            within_data_missing_ids = []
                            within_data_obs_ids = []

                            for (j, period) in enumerate(xaxis[cond]) #enumerate(sort(collect(keys(time_dict[source_id]))))
                                ind = findall(x -> x == period, p_xaxis)[1]

                                if year(period) in years_df
                                    push!(miss_indices, ind)
                                    push!(within_data_missing_ids, j)
                                else
                                    push!(obs_indices, ind)
                                    push!(within_data_obs_ids, j)
                                end
                            end

                            # Sort 
                            id_for_missing = sort(miss_indices)
                            id_for_data = sort(obs_indices)
                            within_data_missing_ids = sort(within_data_missing_ids)
                            within_data_obs_ids = sort(within_data_obs_ids)
                        end


                        Plots.plot()
                        intd = plot_name == "CEX" || occursin("SIPP", plot_name) ? 20 : 40

                        # if data_bounds != false
                        Plots.plot!(p_xaxis,
                            p_data,
                            ylabel=obj == "quantiles" ? L"\textrm{%$(M)\, \, rel.\,  to\,\,   average}" : L"\textrm{%$(M)\,\,%$(series[ds]) \,\,%$(obj)}",
                            xformatter=:latex,
                            yformatter=:latex,
                            xticks=(p_xaxis[1:intd:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(p_xaxis[1:intd:end])]),
                            legend=:best,
                            xtickfontsize=14,
                            ytickfontsize=14,
                            legendfontsize=10,
                            guidefontsize=14,
                            lc=:red,
                            lw=4,
                            # yerror = (dist[1][source][obj][meas] - dist[1][source]["lb"][obj][meas], dist[1][source]["ub"][obj][meas] - dist[1][source][obj][meas]),
                            # label=L"\textrm{Baseline}"
                            label=""
                        )

                        # define difference condition 
                        # diff_cond = findall(dist[1][source][obj][meas] .- dist[2][source][obj][meas] .!= 0)

                        Plots.plot!(p_xaxis,
                            p_data2,
                            lc=:blue,
                            ls=:dot,
                            la=0.5,
                            lw=4,
                            # label=L"\textrm{Less\,\,data}"
                            label=""
                        )

                        if plot_name == "CEX"
                            Plots.plot!(xaxis[cond],
                                dist[1][source]["lb"][obj][meas][cond],
                                fillrange=dist[1][source]["ub"][obj][meas][cond],
                                fillalpha=0.1,
                                fillcolor=:red,
                                la=0.0,
                                lw=4, dpi=500,
                                label=""
                            )
                            Plots.scatter!(xaxis[cond], # findall(!isnan, top[source * "_obs"][obj][meas])
                                c_data,
                                marker=:dot,
                                markercolor=:black,
                                ms=5,
                                lc=:black,
                                la=0.5,
                                lw=4,
                                label=""
                                # label= ds == "bot" ? L"\textrm{Observation}" : ""
                            )
                        else
                            # The observations 
                            Plots.scatter!(xaxis[cond], # findall(!isnan, top[source * "_obs"][obj][meas])
                                c_data,
                                marker=:dot,
                                markercolor=:black,
                                ms=5,
                                lc=:black,
                                la=0.5,
                                lw=4,
                                yerror=(c_data - dist[1][source]["lb"][obj][meas][cond], dist[1][source]["ub"][obj][meas][cond] - c_data),
                                label=""
                                # label= ds == "bot" ? L"\textrm{Observation}" : ""
                            )
                        end

                        # Missing Observations 
                        Plots.scatter!(p_xaxis[id_for_missing], # findall(!isnan, top[source * "_obs"][obj][meas])
                            c_data[within_data_missing_ids],
                            marker=:dot,
                            mc=:white, msc=:black, msw=5,
                            la=0.5,
                            lw=4, dpi=500,
                            label=""
                        )

                        ymin = nanminimum(p_data) - 0.05 * nanminimum(p_data)
                        ymax = nanmaximum(p_data) + 0.05 * nanmaximum(p_data)
                        if collect(keys(filtering_criteria))[1] == "periods_to_remove"
                            nothing
                        else
                            Plots.plot!(
                                p_xaxis[range_missing],
                                repeat([ymin], length(p_xaxis[range_missing])),
                                fillrange=ymax,
                                fillalpha=0.3,
                                color=:gray,
                                label=""
                            )
                        end

                        Plots.savefig(path * "/$meas/" * "/forecasts/$(forecast_type)/all_data/" * plot_name * "_" * meas * "_" * obj * "_" * series[ds] * "_" * label * ".pdf")
                    end
                end
            end


            # # mid 
            # for (j, source) in enumerate(data_sources)
            #     ρ = round(cor(mid[source][obj][meas], mid_r[source][obj][meas]) , digits=2)

            #     # Get the R2 of the data on the two Reconstructions
            #     c_data  = mid[source * "_obs"][obj][meas][.!isnan.(mid[source * "_obs"][obj][meas])]
            #     c_recon = mid[source][obj][meas][.!isnan.(mid[source * "_obs"][obj][meas])]
            #     c_forec = mid_r[source][obj][meas][.!isnan.(mid[source * "_obs"][obj][meas])]
            #     R2        = round(1 - sum((c_data .- c_recon).^2) / sum((c_data .- mean(c_data)).^2), digits=2) 
            #     R2_r      = round(1 - sum((c_data .- c_forec).^2) / sum((c_data .- mean(c_data)).^2), digits=2) 

            #     Plots.plot()
            #         # if data_bounds != false
            #             Plots.plot!(xaxis,
            #                 mid[source][obj][meas],
            #                 xlabel = L"\textrm{Year}",
            #                 ylabel = L"\textrm{%$(M)}\,\,\textrm{%$(obj)}",
            #                 xformatter=:latex, 
            #                 yformatter=:latex, 
            #                 xticks=(xaxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_,date) in enumerate(xaxis[1:20:end])]),
            #                 legend=:outertopright,
            #                 lw=2,
            #                 # yerror = (mid[source][obj][meas] - mid[source]["lb"][obj][meas], mid[source]["ub"][obj][meas] - mid[source][obj][meas]),
            #                 label=L"\textrm{%$(source)}\,\,\textrm{Reconstruction}"
            #             )

            #             Plots.plot!(xaxis,
            #             mid_r[source][obj][meas],
            #             xlabel = L"\textrm{Year}",
            #             ylabel = L"\textrm{%$(M)}\,\,\textrm{%$(obj)}",
            #             xformatter=:latex, 
            #             yformatter=:latex, 
            #             xticks=(xaxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_,date) in enumerate(xaxis[1:20:end])]),
            #             legend=:outertopright,
            #             lw=2,
            #             # yerror = (mid[source][obj][meas] - mid[source]["lb"][obj][meas], mid[source]["ub"][obj][meas] - mid[source][obj][meas]),
            #             label=L"\textrm{%$(source)}\,\,\textrm{Forecast}"
            #         )

            #         Plots.plot!(xaxis,
            #         mid_r[source][obj][meas],
            #         lw=0,
            #         linecolor=:white,
            #         # yerror = (top[source][obj][meas] - top[source]["lb"][obj][meas], top[source]["ub"][obj][meas] - top[source][obj][meas]),
            #         label=L"\textrm{corr.: %$(ρ)}", 
            #         )


            #     Plots.scatter!(xaxis[findall(!isnan, mid[source * "_obs"][obj][meas])],
            #         c_data,
            #         marker=:diamond,
            #         markercolor=:black,
            #         lc=:black,
            #         la=0.5,
            #         lw=2,
            #         label=L"\textrm{%$(source)}\,\,\textrm{obs}"
            #     )

            #     Plots.plot!(xaxis,
            #     mid_r[source][obj][meas],
            #     lw=0,
            #     linecolor=:white,
            #     # yerror = (top[source][obj][meas] - top[source]["lb"][obj][meas], top[source]["ub"][obj][meas] - top[source][obj][meas]),
            #     label=L"\textrm{model-to-forecast\,\,corr.: %$(ρ)}", 
            #     )

            #     Plots.plot!(xaxis,
            #     mid_r[source][obj][meas],
            #     lw=0,
            #     linecolor=:white,
            #     # yerror = (top[source][obj][meas] - top[source]["lb"][obj][meas], top[source]["ub"][obj][meas] - top[source][obj][meas]),
            #     label=L"R^2\,\,\textrm{data-to-model: %$(R2)}", 
            #     )

            #     Plots.plot!(xaxis,
            #     mid_r[source][obj][meas],
            #     lw=0,
            #     linecolor=:white,
            #     # yerror = (top[source][obj][meas] - top[source]["lb"][obj][meas], top[source]["ub"][obj][meas] - top[source][obj][meas]),
            #     label=L"R^2\,\,\textrm{data-to-forecast: %$(R2_r)}", 
            #     )

            #     Plots.savefig(path * "/$meas/" * "/forecasts/$(forecast_type)/all_data/" * source * "_" * meas * "_" * obj * "_" * series["mid"] * "_" * label * ".pdf")
            # end

            # # bot  
            # for (j, source) in enumerate(data_sources)
            #     ρ = round(cor(bot[source][obj][meas], bot_r[source][obj][meas]) , digits=2)

            #     # Get the R2 of the data on the two Reconstructions
            #     c_data  = bot[source * "_obs"][obj][meas][.!isnan.(bot[source * "_obs"][obj][meas])]
            #     c_recon = bot[source][obj][meas][.!isnan.(bot[source * "_obs"][obj][meas])]
            #     c_forec = bot_r[source][obj][meas][.!isnan.(bot[source * "_obs"][obj][meas])]
            #     R2        = round(1 - sum((c_data .- c_recon).^2) / sum((c_data .- mean(c_data)).^2), digits=2) 
            #     R2_r      = round(1 - sum((c_data .- c_forec).^2) / sum((c_data .- mean(c_data)).^2), digits=2) 

            #     Plots.plot()
            #     # if data_bounds != false
            #         Plots.plot!(xaxis,
            #             bot[source][obj][meas],
            #             xlabel = L"\textrm{Year}",
            #             ylabel = L"\textrm{%$(M)}\,\,\textrm{%$(obj)}",
            #             xformatter=:latex, 
            #             yformatter=:latex, 
            #             xticks=(xaxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_,date) in enumerate(xaxis[1:20:end])]),
            #             legend=:outertopright,
            #             lw=2,
            #             # yerror = (bot[source][obj][meas] - bot[source]["lb"][obj][meas], bot[source]["ub"][obj][meas] - bot[source][obj][meas]),
            #             label=L"\textrm{%$(source)}\,\,\textrm{Reconstruction}"
            #         )

            #         Plots.plot!(xaxis,
            #         bot_r[source][obj][meas],
            #         xlabel = L"\textrm{Year}",
            #         ylabel = L"\textrm{%$(M)}\,\,\textrm{%$(obj)}",
            #         xformatter=:latex, 
            #         yformatter=:latex, 
            #         xticks=(xaxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_,date) in enumerate(xaxis[1:20:end])]),
            #         legend=:outertopright,
            #         lw=2,
            #         # yerror = (bot[source][obj][meas] - bot[source]["lb"][obj][meas], bot[source]["ub"][obj][meas] - bot[source][obj][meas]),
            #         label=L"\textrm{%$(source)}\,\,\textrm{Forecast}"
            #     )

            #     Plots.plot!(xaxis,
            #     bot_r[source][obj][meas],
            #     lw=0,
            #     linecolor=:white,
            #     # yerror = (top[source][obj][meas] - top[source]["lb"][obj][meas], top[source]["ub"][obj][meas] - top[source][obj][meas]),
            #     label=L"\textrm{corr.: %$(ρ)}", 
            #     )


            #     Plots.scatter!(xaxis[findall(!isnan, bot[source * "_obs"][obj][meas])],
            #         c_data,
            #         marker=:diamond,
            #         markercolor=:black,
            #         lc=:black,
            #         la=0.5,
            #         lw=2,
            #         label=L"\textrm{%$(source)}\,\,\textrm{obs}"
            #     )

            #     Plots.plot!(xaxis,
            #     bot_r[source][obj][meas],
            #     lw=0,
            #     linecolor=:white,
            #     # yerror = (top[source][obj][meas] - top[source]["lb"][obj][meas], top[source]["ub"][obj][meas] - top[source][obj][meas]),
            #     label=L"\textrm{model-to-forecast\,\,corr.: %$(ρ)}", 
            #     )

            #     Plots.plot!(xaxis,
            #     bot_r[source][obj][meas],
            #     lw=0,
            #     linecolor=:white,
            #     # yerror = (top[source][obj][meas] - top[source]["lb"][obj][meas], top[source]["ub"][obj][meas] - top[source][obj][meas]),
            #     label=L"R^2\,\,\textrm{data-to-model: %$(R2)}", 
            #     )

            #     Plots.plot!(xaxis,
            #     bot_r[source][obj][meas],
            #     lw=0,
            #     linecolor=:white,
            #     # yerror = (top[source][obj][meas] - top[source]["lb"][obj][meas], top[source]["ub"][obj][meas] - top[source][obj][meas]),
            #     label=L"R^2\,\,\textrm{data-to-forecast: %$(R2_r)}", 
            #     )

            #     Plots.savefig(path * "/$meas/" * "/forecasts/$(forecast_type)/all_data/" * source * "_" * meas * "_" * obj * "_" * series["bot"] * "_" * label * ".pdf")
            # end
        end
    end
end



function perform_iterative_forecast(dataset_choice, par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, obs_data, model_options, model_elements, time_params, func_data, opttag)
    @unpack y, MV = model_elements
    @unpack data_sources, func_dict = func_data
    @unpack case, measures, estimator = model_options
    @unpack tmin, tmax = time_params

    @unpack grid_cop, grid_pcf = estimator

    user_t = (deepcopy(tmin), deepcopy(tmax))
    data_id = findfirst(x -> x == dataset_choice, data_sources)

    optfolder = opttag == :optimization ? "from_optimization" : "from_mcmc"

    dimension = length(measures)
    diff = grid_cop + (dimension - 1) * (grid_cop - 1)
    n = grid_cop^dimension + grid_pcf * dimension - diff
    label = "$dimension" * "D" * "_$case"
    orig_y = deepcopy(y) # for later 
    N(x) = length(x)
    cop_tup = ntuple(_ -> (:), dimension)

    # First, get the original reconstruction 
    smoother_output, _ = likeli(model_elements, par_final, param_sizes, hyperpriors, Σ_ids, model_options; smooth=true)
    @unpack x_smoothed = smoother_output
    @unpack proj = model_elements
    dv, _ = reconstruct_data(par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources)
    predictions_dict = Dict()

    dv[dataset_choice]["normal"], avg_series = export_functional_data(dv[dataset_choice]["normal"], "normal", dataset_choice, opttag, obs_data, func_data, time_params, user_t, model_options, false, false, :forecast)


    # I run this only to get the object, so I can fill it with the predictions later 
    # predictions_dict[dataset_choice] = export_functional_data(dv[dataset_choice], dataset_choice, opttag, obs_data, func_data, time_params, user_t, model_options, false, false, :forecast)     
    predictions_dict[dataset_choice] = deepcopy(dv[dataset_choice]["normal"])

    # Second, generate the predictions, removing one observation at a time, keep prediction corresponding to that removed observation 
    # Now remove one observation at a time, for each dataset, calculate the MSE and store it in the dictionary

    for k in [dataset_choice] #enumerate(sort!(data_sources)) # TODO: includes consensus 
        # Indices of observed data
        indices = findall(x -> !all(isnan.(x)), eachcol(MV[data_id]))
        start = (data_id - 1) * n + 1
        data_end = data_id * n
        data_indices = start:data_end

        for (t, id) in enumerate(indices)
            # Copy the original observations
            y = deepcopy(orig_y)

            # Remove one observation in MV for this dataset 
            y[data_indices, id] .= NaN

            # Define this MV as our new measurement vector 
            @pack! model_elements = y

            # Reconstruct data with new MV
            dv_r, _ = reconstruct_data(par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources)

            # for (c, k) in enumerate(keys(dv)) 
            # Relative to the average
            dv_r[k]["normal"], _ = export_functional_data(dv_r[k]["normal"], "normal", k, opttag, obs_data, func_data, time_params, user_t, model_options, false, false, :forecast)

            # Now, I just want to store the prediction corresponding to the removed observation
            for m in measures
                for o in ["quantiles"]
                    predictions_dict[dataset_choice][m][o]["data"][:, id] = dv_r[k]["normal"][m][o]["data"][:, id]
                    for s in ["next40", "top10", "bottom50"]
                        predictions_dict[dataset_choice][m][o]["common series"][s][id] = dv_r[k]["normal"][m][o]["common series"][s][id]
                    end
                end
            end
            predictions_dict[dataset_choice]["copulas"]["data"][cop_tup..., id] = dv_r[k]["normal"]["copulas"]["data"][cop_tup..., id]
        end

        # Now drop everything that isnt in indices 
        ids_to_exclude = setdiff(1:size(MV[data_id], 2), indices)

        for m in measures
            for o in ["quantiles"]
                predictions_dict[dataset_choice][m][o]["data"][:, ids_to_exclude] .= NaN
                for s in ["next40", "top10", "bottom50"]
                    predictions_dict[dataset_choice][m][o]["common series"][s][ids_to_exclude] .= NaN
                end
            end
        end
        predictions_dict[dataset_choice]["copulas"]["data"][cop_tup..., ids_to_exclude] .= NaN
    end

    # Now we toss this into the plotting function and plot the estimated model predictions, the data and the predictions from removing the observation 
    # dv[dataset_choice] = export_functional_data(dv[dataset_choice], dataset_choice, opttag, obs_data, func_data, time_params, user_t, model_options, false, false, :forecast)     

    plot_iterative_forecasts(dataset_choice, dv, predictions_dict, func_data, obs_data, user_t, time_params, model_options, opttag, label, avg_series) # select_series could go in here, but it's already relative to the average  
    y = deepcopy(orig_y)
    @pack! model_elements = y
end


function plot_iterative_forecasts(dataset_choice, data_dt, data_dt_r, func_data, obs_data, user_params, time_params, model_options, type, label, avg_series)
    @unpack measures, lags, estimator, tag = model_options
    @unpack tmin, tmax = time_params
    @unpack func_dict, year_vec, confidence_intervals = func_data
    @unpack gdp_series = obs_data

    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid = estimator
    else
        @unpack grid_pcf, grid_cop = estimator
    end
    grid_choice = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf

    data_sources = [dataset_choice]
    smin, smax = user_params
    m_label = measures_folder(measures)
    dates = Dict()

    # Path situation
    folder = type == :optimization ? "from_optimization" : "from_mcmc"
    init_path = BASE_PATH
    path = init_path * "/7_Results/$m_label" * tag * "/$folder/plots/"
    mkpath(path)

    # # Common annual series
    top = Dict()
    bot = Dict()
    mid = Dict()

    # Same, but for the forecast 
    top_r = Dict()
    bot_r = Dict()
    mid_r = Dict()

    for dt in [top, bot, mid, top_r, bot_r, mid_r]
        for source in data_sources
            dt[source] = Dict()
            dt[source*"_obs"] = Dict()
            dt[source]["lb"] = Dict()
            dt[source]["ub"] = Dict()

            for obj in ["shares", "quantiles"]
                dt[source][obj] = Dict()
                dt[source*"_obs"][obj] = Dict()
                dt[source]["lb"][obj] = Dict()
                dt[source]["ub"][obj] = Dict()
            end
        end
    end

    # Defining the series 
    series = Dict()
    if grid_choice == 10 || grid_choice == 20 || grid_choice == 100
        series["bot"] = "bottom50"
        series["mid"] = "next40"
        series["top"] = "top10"
    elseif grid_choice == 5
        series["bot"] = "bottom40"
        series["mid"] = "next40"
        series["top"] = "top20"
    end

    gps = find_grid_points_for_series(grid_choice)

    # This is to correct the observed data  
    # gdp_series[!, "date"]    = QuarterlyDate.(gdp_series[!, "time"])
    # correction_names         = [meas * "_per_hh" for meas in measures]
    # select_series            = select(gdp_series, ["date", correction_names...])

    # filter!(row -> row.date >= QuarterlyDate(tmin["year"], tmin["quarter"]), select_series)
    # filter!(row -> row.date <= QuarterlyDate(tmax["year"], tmax["quarter"]), select_series)


    dates = Dict()

    for meas in measures
        M = uppercasefirst(meas)

        # Align dates, create xaxis
        dates[meas] = Dict()
        dates[meas]["tmin"] = Dict()
        dates[meas]["tmax"] = Dict()

        dates[meas]["tmin"]["year"] = smin["year"]
        dates[meas]["tmax"]["year"] = smax["year"]

        dates[meas]["tmin"]["quarter"] = smin["quarter"]
        dates[meas]["tmax"]["quarter"] = smax["quarter"]

        @info(dates[meas])
        for obj in ["quantiles"] #TODO: levels?

            # Adjust frequencies, Confine to dates
            for source in data_sources

                # For the confidence intervals
                # if obj == "quantiles"
                top[source]["lb"][obj][meas] = subset_to_cutoff(confidence_intervals[source]["ci_l"][meas][obj][gps[3], :], dates[meas], tmin, tmax, obj, avg_series[!, meas*"_per_hh"])  # TODO: generalize grid

                # Average of the next40 and average of the bottom50 
                mid[source]["lb"][obj][meas] = subset_to_cutoff(sum(confidence_intervals[source]["ci_l"][meas][obj][gps[2], :], dims=1)' ./ length(gps[2]), dates[meas], tmin, tmax, obj, avg_series[!, meas*"_per_hh"])
                bot[source]["lb"][obj][meas] = subset_to_cutoff(sum(confidence_intervals[source]["ci_l"][meas][obj][gps[1], :], dims=1)' ./ length(gps[1]), dates[meas], tmin, tmax, obj, avg_series[!, meas*"_per_hh"])

                top[source]["ub"][obj][meas] = subset_to_cutoff(confidence_intervals[source]["ci_u"][meas][obj][gps[3], :], dates[meas], tmin, tmax, obj, avg_series[!, meas*"_per_hh"])
                mid[source]["ub"][obj][meas] = subset_to_cutoff(sum(confidence_intervals[source]["ci_u"][meas][obj][gps[2], :], dims=1)' ./ length(gps[2]), dates[meas], tmin, tmax, obj, avg_series[!, meas*"_per_hh"])
                bot[source]["ub"][obj][meas] = subset_to_cutoff(sum(confidence_intervals[source]["ci_u"][meas][obj][gps[1], :], dims=1)' ./ length(gps[1]), dates[meas], tmin, tmax, obj, avg_series[!, meas*"_per_hh"])

                # elseif obj == "shares"
                #     top[source]["lb"][obj][meas]         = subset_to_cutoff(confidence_intervals[source]["ci_l"][meas][obj][3, :], dates[meas], tmin, tmax, obj, false)  # TODO: generalize grid
                #     mid[source]["lb"][obj][meas]         = subset_to_cutoff(confidence_intervals[source]["ci_l"][meas][obj][2, :], dates[meas], tmin, tmax, obj, false)
                #     bot[source]["lb"][obj][meas]         = subset_to_cutoff(confidence_intervals[source]["ci_l"][meas][obj][1, :], dates[meas], tmin, tmax, obj, false)

                #     top[source]["ub"][obj][meas]         = subset_to_cutoff(confidence_intervals[source]["ci_u"][meas][obj][3, :], dates[meas], tmin, tmax, obj, false)
                #     mid[source]["ub"][obj][meas]         = subset_to_cutoff(confidence_intervals[source]["ci_u"][meas][obj][2, :], dates[meas], tmin, tmax, obj, false)
                #     bot[source]["ub"][obj][meas]         = subset_to_cutoff(confidence_intervals[source]["ci_u"][meas][obj][1, :], dates[meas], tmin, tmax, obj, false)
                # end

                # For the reconstructions 
                top[source][obj][meas] = subset_to_cutoff(data_dt[source]["normal"][meas][obj]["common series"][series["top"]], dates[meas], tmin, tmax, obj, avg_series[!, meas*"_per_hh"])
                bot[source][obj][meas] = subset_to_cutoff(data_dt[source]["normal"][meas][obj]["common series"][series["bot"]], dates[meas], tmin, tmax, obj, avg_series[!, meas*"_per_hh"])
                mid[source][obj][meas] = subset_to_cutoff(data_dt[source]["normal"][meas][obj]["common series"][series["mid"]], dates[meas], tmin, tmax, obj, avg_series[!, meas*"_per_hh"])

                top_r[source][obj][meas] = subset_to_cutoff(data_dt_r[source][meas][obj]["common series"][series["top"]], dates[meas], tmin, tmax, obj, avg_series[!, meas*"_per_hh"])
                bot_r[source][obj][meas] = subset_to_cutoff(data_dt_r[source][meas][obj]["common series"][series["bot"]], dates[meas], tmin, tmax, obj, avg_series[!, meas*"_per_hh"])
                mid_r[source][obj][meas] = subset_to_cutoff(data_dt_r[source][meas][obj]["common series"][series["mid"]], dates[meas], tmin, tmax, obj, avg_series[!, meas*"_per_hh"])

                # For the observations
                # if obj == "quantiles"
                top[source*"_obs"][obj][meas] = func_dict[source][meas][obj]["common series"][series["top"]] ./ avg_series[!, meas*"_per_hh"]
                bot[source*"_obs"][obj][meas] = func_dict[source][meas][obj]["common series"][series["bot"]] ./ avg_series[!, meas*"_per_hh"]
                mid[source*"_obs"][obj][meas] = func_dict[source][meas][obj]["common series"][series["mid"]] ./ avg_series[!, meas*"_per_hh"]

                top[source*"_obs"][obj][meas] = subset_to_cutoff(top[source*"_obs"][obj][meas], dates[meas], tmin, tmax, obj, 1)
                bot[source*"_obs"][obj][meas] = subset_to_cutoff(bot[source*"_obs"][obj][meas], dates[meas], tmin, tmax, obj, 1)
                mid[source*"_obs"][obj][meas] = subset_to_cutoff(mid[source*"_obs"][obj][meas], dates[meas], tmin, tmax, obj, 1)

                # else
                #     top[source * "_obs"][obj][meas]          = subset_to_cutoff(func_dict[source][meas][obj]["common series"][series["top"]], dates[meas], tmin, tmax, obj, false)
                #     bot[source * "_obs"][obj][meas]          = subset_to_cutoff(func_dict[source][meas][obj]["common series"][series["bot"]], dates[meas], tmin, tmax, obj, false)
                #     mid[source * "_obs"][obj][meas]          = subset_to_cutoff(func_dict[source][meas][obj]["common series"][series["mid"]], dates[meas], tmin, tmax, obj, false)    
                # end

            end

            # Plots of both data, 2 lines per plot (we take the first quarter for the year)
            xaxis = QuarterlyDate(dates[meas]["tmin"]["year"], dates[meas]["tmin"]["quarter"]):Quarter(1):QuarterlyDate(dates[meas]["tmax"]["year"], dates[meas]["tmax"]["quarter"])
            dist_dict = Dict("top" => [top, top_r], "bottom" => [bot, bot_r], "middle" => [mid, mid_r])
            @info("Building reconstruction plots, comparing reconstructions to data: $obj")

            # Top
            for (ds, dist) in dist_dict
                for (j, source) in enumerate(data_sources)
                    Plots.plot()
                    plot_name = occursin("CEX", source) ? "CEX" : source

                    # Reconstruction 
                    Plots.plot!(xaxis,
                        dist[1][source][obj][meas],
                        # xlabel = L"\textrm{Year}",
                        ylabel=M == "Consum" ? L"\textrm{Consumption\, \, rel.\,  to\,\,  average}" : L"\textrm{%$(M)\, \, rel.\,  to\,\,  average}",
                        xformatter=:latex,
                        yformatter=:latex,
                        xticks=(xaxis[1:40:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(xaxis[1:40:end])]),
                        legend=ds == "bottom" ? :best : false,
                        lc=:red, #select_color(plot_name),
                        lw=4,
                        xtickfontsize=14,
                        ytickfontsize=14,
                        legendfontsize=10,
                        guidefontsize=14,
                        label=L"\textrm{Baseline}" #L"\textrm{%$(plot_name)}\,\,\textrm{Approximation}"
                    )

                    # Predictions 
                    Plots.scatter!(xaxis,
                        dist[2][source][obj][meas],
                        marker=:square,
                        markercolor=:black,
                        label=L"\textrm{Prediction}"
                    )

                    # Observations
                    Plots.scatter!(xaxis,
                        dist[1][source*"_obs"][obj][meas],
                        marker=:circle,
                        ms=4,
                        markercolor=:black,
                        lc=:black,
                        la=0.5,
                        lw=2,
                        yerror=(dist[1][source*"_obs"][obj][meas] - dist[1][source]["lb"][obj][meas], dist[1][source]["ub"][obj][meas] - dist[1][source*"_obs"][obj][meas]),
                        label=L"\textrm{Observation}"
                    )

                    Plots.savefig(path * "/$meas/" * "/forecasts/iterative/" * plot_name * "_" * meas * "_" * obj * "_" * ds * "_" * label * ".pdf")
                end
            end
        end
    end
end