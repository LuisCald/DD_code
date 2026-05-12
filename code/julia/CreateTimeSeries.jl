function create_time_series_dictionary(data_vector, estimator, measures)

    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid = estimator
    else
        @unpack grid_pcf = estimator
    end
    grid_choice = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf

    dimension = length(measures)

    # Organize data into dictionaries 
    data_dict = Dict("copulas" => Dict())
    data_dict["copulas"]["data"] = data_vector[1]

    # Levels data 
    quantiles = data_vector[2]
    levels = data_vector[3]
    shares = data_vector[4]

    # By measure ...
    indices = [I for I in Iterators.partition(axes(data_vector[2], 1), grid_choice)]

    for (i, meas) in enumerate(measures)
        data_dict[meas] = Dict("levels" => Dict(), "quantiles" => Dict(), "shares" => Dict())
        data_dict[meas]["levels"]["data"] = levels[indices[i], :]
        data_dict[meas]["shares"]["data"] = shares[indices[i], :]
        data_dict[meas]["quantiles"]["data"] = quantiles[indices[i], :]

        data_dict[meas]["levels"]["common series"] = Dict()
        data_dict[meas]["shares"]["common series"] = Dict()
        data_dict[meas]["quantiles"]["common series"] = Dict()

        if grid_choice == 10 # deciles        
            data_dict[meas]["levels"]["common series"]["top10"] = levels[indices[i], :][10, :]
            data_dict[meas]["levels"]["common series"]["next40"] = sum(levels[indices[i], :][6:9, :], dims=1)'
            data_dict[meas]["levels"]["common series"]["bottom50"] = sum(levels[indices[i], :][1:5, :], dims=1)'

            data_dict[meas]["shares"]["common series"]["top10"] = shares[indices[i], :][10, :]
            data_dict[meas]["shares"]["common series"]["next40"] = sum(shares[indices[i], :][6:9, :], dims=1)'
            data_dict[meas]["shares"]["common series"]["bottom50"] = sum(shares[indices[i], :][1:5, :], dims=1)'

            data_dict[meas]["quantiles"]["common series"]["top10"] = quantiles[indices[i], :][10, :]
            data_dict[meas]["quantiles"]["common series"]["next40"] = sum(quantiles[indices[i], :][6:9, :], dims=1)' ./ 4
            data_dict[meas]["quantiles"]["common series"]["bottom50"] = sum(quantiles[indices[i], :][1:5, :], dims=1)' ./ 5
        elseif grid_choice == 5
            data_dict[meas]["levels"]["common series"]["top20"] = levels[indices[i], :][5, :]
            data_dict[meas]["levels"]["common series"]["next40"] = sum(levels[indices[i], :][3:4, :], dims=1)'
            data_dict[meas]["levels"]["common series"]["bottom40"] = sum(levels[indices[i], :][1:2, :], dims=1)'

            data_dict[meas]["shares"]["common series"]["top20"] = shares[indices[i], :][5, :]
            data_dict[meas]["shares"]["common series"]["next40"] = sum(shares[indices[i], :][3:4, :], dims=1)'
            data_dict[meas]["shares"]["common series"]["bottom40"] = sum(shares[indices[i], :][1:2, :], dims=1)'

            data_dict[meas]["quantiles"]["common series"]["top20"] = quantiles[indices[i], :][5, :]
            data_dict[meas]["quantiles"]["common series"]["next40"] = sum(quantiles[indices[i], :][3:4, :], dims=1)' ./ 2
            data_dict[meas]["quantiles"]["common series"]["bottom40"] = sum(quantiles[indices[i], :][1:2, :], dims=1)' ./ 2

        elseif grid_choice == 100
            data_dict[meas]["levels"]["common series"]["top10"] = sum(levels[indices[i], :][91:100, :], dims=1)'
            data_dict[meas]["levels"]["common series"]["next40"] = sum(levels[indices[i], :][51:90, :], dims=1)'
            data_dict[meas]["levels"]["common series"]["bottom50"] = sum(levels[indices[i], :][1:50, :], dims=1)'

            data_dict[meas]["shares"]["common series"]["top10"] = sum(shares[indices[i], :][91:100, :], dims=1)'
            data_dict[meas]["shares"]["common series"]["next40"] = sum(shares[indices[i], :][51:90, :], dims=1)'
            data_dict[meas]["shares"]["common series"]["bottom50"] = sum(shares[indices[i], :][1:50, :], dims=1)'

            data_dict[meas]["quantiles"]["common series"]["top10"] = sum(quantiles[indices[i], :][91:100, :], dims=1)' ./ 10
            data_dict[meas]["quantiles"]["common series"]["next40"] = sum(quantiles[indices[i], :][51:90, :], dims=1)' ./ 40
            data_dict[meas]["quantiles"]["common series"]["bottom50"] = sum(quantiles[indices[i], :][1:50, :], dims=1)' ./ 50
        elseif grid_choice == 20
            data_dict[meas]["levels"]["common series"]["top10"] = sum(levels[indices[i], :][19:20, :], dims=1)'
            data_dict[meas]["levels"]["common series"]["next40"] = sum(levels[indices[i], :][11:18, :], dims=1)'
            data_dict[meas]["levels"]["common series"]["bottom50"] = sum(levels[indices[i], :][1:10, :], dims=1)'

            data_dict[meas]["shares"]["common series"]["top10"] = sum(shares[indices[i], :][19:20, :], dims=1)'
            data_dict[meas]["shares"]["common series"]["next40"] = sum(shares[indices[i], :][11:18, :], dims=1)'
            data_dict[meas]["shares"]["common series"]["bottom50"] = sum(shares[indices[i], :][1:10, :], dims=1)'

            data_dict[meas]["quantiles"]["common series"]["top10"] = sum(quantiles[indices[i], :][19:20, :], dims=1)' ./ 2
            data_dict[meas]["quantiles"]["common series"]["next40"] = sum(quantiles[indices[i], :][11:18, :], dims=1)' ./ 8
            data_dict[meas]["quantiles"]["common series"]["bottom50"] = sum(quantiles[indices[i], :][1:10, :], dims=1)' ./ 10
        else
            error("Grid size not supported")
        end
    end

    return data_dict
end


function export_raw_data(data_dict, estimator, source, measures, time_p, tag)
    @unpack tmin, tmax = time_p
    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid, grid_pcf = estimator
    else
        @unpack grid_pcf = estimator
    end
    grid_choice = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf
    grid_points = select_grid_points(grid_choice)


    # Loop through each measure in the dictionary
    for meas in measures
        # Get the data and variable names for the current measure
        levels = data_dict[meas]["levels"]["data"]'
        shares = data_dict[meas]["shares"]["data"]'
        quants = data_dict[meas]["quantiles"]["data"]'

        # Put together 
        col_names = vcat(["$meas" * "_quants_$i" for i in grid_points], ["$meas" * "_levels_$i" for i in grid_points], ["$meas" * "_shares_$i" for i in grid_points])
        df = DataFrame(hcat(quants, levels, shares), Symbol.(col_names))
        df[!, "time"] = QuarterlyDate(tmin["year"], tmin["quarter"]):Quarter(1):QuarterlyDate(tmax["year"], tmax["quarter"])

        # Create DataFrame for the current measure
        m_label = measures_folder(measures)
        init_path = BASE_PATH
        path = init_path * "/7_Results/$m_label" * tag * "/other_results/raw_data/"
        mkpath(path)
        f_cols = filter(e -> e != "time", names(df))
        subset!(df, (f_cols .=> ByRow(x -> !(x isa Number && isnan(x))))...)

        CSV.write(path * source * "_" * meas * ".csv", select(df, :time, :))
    end
end


function generate_specific_plots(data_dict, ty, func_data, data_name, time_params, timeframe, model_options, type, select_series, gdp_series, posterior_bounds=false)
    @unpack measures, case, number_of_dfs, estimator, equivalized, bottom_coded, compare_to_other_est, tag, freq = model_options
    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid, integral_cop_grid = estimator
    else
        @unpack grid_pcf, grid_cop = estimator
    end
    # grid_choice_pcf = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf
    grid_choice_cop = typeof(estimator) <: SeriesEstimator ? integral_cop_grid : grid_cop

    equiv = equivalized == true ? "eq" : ""
    botcod = !isempty(bottom_coded) ? "bc" : ""
    label = "$case" * "_$equiv" * "$botcod"

    # Wrestling with time
    @unpack tmin, tmax = time_params  # m = model
    smin, smax = timeframe    # what user wants 

    @assert smin["year"] >= tmin["year"]
    @assert smax["year"] <= tmax["year"]

    # Find keys associated with objects  
    # keys1 = filter(key -> key in ["copulas"], keys(data_dict))

    # Create new dictionaries. one for copulas, other for pcfs  
    # d1    = Dict(key => data_dict[key] for key in keys1)

    # # Generate plots!
    # if number_of_dfs !=1
    #     generate_copula_plots(d1, func_data, data_name, tmin, tmax, grid_choice_cop, type, measures, label, posterior_bounds)
    # end

    # if occursin("HANK full", tag)
    #     within_stat_dict = generate_quantiles_shares_levels_HANK_full(data_dict, ty, func_data, data_name, smin, smax, tmin, tmax, estimator, label, type, measures, time_params, select_series, gdp_series, posterior_bounds, compare_to_other_est, tag, freq)
    if occursin("HANK", tag)
        within_stat_dict = generate_quantiles_shares_levels_HANK(data_dict, ty, func_data, data_name, smin, smax, tmin, tmax, estimator, label, type, measures, time_params, select_series, gdp_series, posterior_bounds, compare_to_other_est, tag, freq)
    else
        within_stat_dict = generate_quantiles_shares_levels_plots(data_dict, ty, func_data, data_name, smin, smax, tmin, tmax, estimator, label, type, measures, time_params, select_series, gdp_series, posterior_bounds, compare_to_other_est, tag)
    end

    return within_stat_dict
end


function generate_copula_plots(data_dict, func_data, data_name, tmin, tmax, grid, type, measures, label, posterior_bounds)
    @unpack func_dict, year_vec, data_sources = func_data

    # Preliminaries
    m_label = measures_folder(measures)
    dimension = length(measures)
    dts = QuarterlyDate(tmin["year"], tmin["quarter"]):Quarter(1):QuarterlyDate(tmax["year"], tmax["quarter"])
    folder = type == :optimization ? "from_optimization" : "from_mcmc"
    init_path = BASE_PATH
    path = init_path * "/7_Results/$m_label" * tag * "/$folder/copulas/"
    mkpath(path)

    # Plotting copulas for 3 separate years
    axis = [1:grid]  # percentile groups 


    # TODO: this is tricky because the data is not always observed and they're observed in different quarters even, so "Q4" does not work. Q3 snags SCF, Q2 snags PSID
    # TODO: so for now I just pick to show the SCF 
    # We want 1970, 1990, 2010
    the_years = []
    for i in [1971, 1992, 2010] # this has to change since the SCF will not always observe every dimension + for these set of years 
        id = findfirst(x -> x == QuarterlyDate("$i-Q3"), dts)
        push!(the_years, id)
    end
    colons = ntuple(_ -> (:), dimension)
    # yMin      = 0
    # yMax      = 35

    # Generate side by side copulas  
    for y in the_years
        ci = Dict()
        if dimension == 2
            measure1 = uppercasefirst(measures[1])
            measure2 = uppercasefirst(measures[2])
            m1 = uppercasefirst(measures[1])[1]
            m2 = uppercasefirst(measures[2])[1]

            Plots.surface(
                axis,
                axis,
                data_dict["copulas"]["data"][colons..., y],  # copula 
                xlabel=L"\textrm{%$(measure1)}",
                ylabel=L"\textrm{%$(measure2)}",
                zlabel=L"dC(%$(m1), %$(m2))",
                xformatter=:latex,
                yformatter=:latex,
                zformatter=:latex,
                legend=false,
                camera=(30, 10),
                size=(400, 400),
                color=:winter,
                dpi=500,
                display_option=Plots.GR.OPTION_SHADED_MESH)

            Plots.savefig(path * "copula_" * "$(dts[y])"[1:4] * "_" * label * ".pdf")

            if data_name == "SCF"
                Plots.surface(
                    axis,
                    axis,
                    reshape(func_dict[data_name]["copulas"]["data"][:, y], (grid, grid)),  # copula 
                    xlabel=L"\textrm{%$(measure1)}",
                    ylabel=L"\textrm{%$(measure2)}",
                    zlabel=L"dC(%$(m1), %$(m2))",
                    xformatter=:latex,
                    yformatter=:latex,
                    zformatter=:latex,
                    legend=false,
                    camera=(30, 10),
                    color=:winter,
                    size=(400, 400),
                    dpi=500,
                    display_option=Plots.GR.OPTION_SHADED_MESH)

                Plots.savefig(path * data_name * "_copula_" * "$(dts[y])"[1:4] * "_reconstructed_" * label * ".pdf")
            end
        end
    end
end

function find_subset_frame(smin, smax, tmin, tmax)
    base_jump = QuarterlyDate(smin["year"], smin["quarter"]) - QuarterlyDate(tmin["year"], tmin["quarter"]) + Quarter(1)
    end_jump = QuarterlyDate(tmax["year"], tmax["quarter"]) - QuarterlyDate(smax["year"], smax["quarter"])

    return base_jump.value, end_jump.value
end

function get_obs_meas(func_dict, data_name, measures; top="top10")
    """ With slight abuse """

    # Find some random series and check if there are only NaNs
    obs_meas = []
    for m in measures
        if all(isnan.(func_dict[data_name][m]["quantiles"]["common series"][top])) == false
            push!(obs_meas, m)
        end
    end
    return obs_meas
end

function fill_data_dicts(dd, measures, pushup, pushdown, estimator)
    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid = estimator
    else
        @unpack grid_pcf = estimator
    end
    grid_choice = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf

    top_k = grid_choice == 10 ? "top10" : grid_choice == 5 ? "top20" : grid_choice == 100 ? "top10" : grid_choice == 20 ? "top10" : error("Grid size not supported")
    mid_k = "next40"
    bot_k = grid_choice == 10 ? "bottom50" : grid_choice == 5 ? "bottom40" : grid_choice == 100 ? "bottom50" : grid_choice == 20 ? "bottom50" : error("Grid size not supported")

    data_dict = Dict()
    # data_dict["copulas"] = Dict()

    pushup = pushup >= Quarter(0) ? pushup + Quarter(1) : pushup

    if pushup.value >= 0 && pushdown.value >= 0
        # data_dict["copulas"]["data"] = dd["copulas"]["data"][:, pushup.value:end-pushdown.value] 

        for meas in measures
            data_dict[meas] = Dict()
            for o in ["levels", "shares", "quantiles"]
                data_dict[meas][o] = Dict()
                data_dict[meas][o]["data"] = dd[meas][o]["data"][:, pushup.value:end-pushdown.value]
                data_dict[meas][o]["common series"] = Dict()
                data_dict[meas][o]["common series"][top_k] = dd[meas][o]["common series"][top_k][pushup.value:end-pushdown.value]
                data_dict[meas][o]["common series"][mid_k] = dd[meas][o]["common series"][mid_k][pushup.value:end-pushdown.value]
                data_dict[meas][o]["common series"][bot_k] = dd[meas][o]["common series"][bot_k][pushup.value:end-pushdown.value]
            end
        end
    elseif pushup.value >= 0 && pushdown.value < 0
        # Add NaNs to the end
        # data_dict["copulas"]["data"] = hcat(dd["copulas"]["data"][:, pushup.value:end], fill!(Array{Float64}(undef, (size(dd["copulas"]["data"], 1), abs(pushdown.value)), NaN)

        for meas in measures
            data_dict[meas] = Dict()
            for o in ["levels", "shares", "quantiles"]
                data_dict[meas][o] = Dict()
                data_dict[meas][o]["data"] = hcat(dd[meas][o]["data"][:, pushup.value:end], fill!(Array{Float64}(undef, size(dd[meas][o]["data"], 1), abs(pushdown.value)), NaN))
                data_dict[meas][o]["common series"] = Dict()
                data_dict[meas][o]["common series"][top_k] = vcat(dd[meas][o]["common series"][top_k][pushup.value:end], fill!(Array{Float64}(undef, abs(pushdown.value)), NaN))
                data_dict[meas][o]["common series"][mid_k] = vcat(dd[meas][o]["common series"][mid_k][pushup.value:end], fill!(Array{Float64}(undef, abs(pushdown.value)), NaN))
                data_dict[meas][o]["common series"][bot_k] = vcat(dd[meas][o]["common series"][bot_k][pushup.value:end], fill!(Array{Float64}(undef, abs(pushdown.value)), NaN))
            end
        end

        # println("wrap")
        # println(pushdown.value)
        # println(pushup.value)
        # println("wrap")

        # pushup       = QuarterlyDate(tmin["year"], tmin["quarter"]) - min_time # Tricky part is if min_time is later than tmin ... I guess i could add NaNs for missing periods. Makes sense!
        # pushdown     = max_time - QuarterlyDate(tmax["year"], tmax["quarter"]) # max_time solved similarly

    elseif pushup.value < 0 && pushdown.value < 0
        # Add NaNs to both beginning and end
        # data_dict["copulas"]["data"] = hcat(fill!(Array{Float64}(undef, (size(dd["copulas"]["data"], 1), abs(pushup.value)), dd["copulas"]["data"], fill!(Array{Float64}(undef, (size(dd["copulas"]["data"], 1), abs(pushdown.value)), NaN)))

        for meas in measures
            data_dict[meas] = Dict()
            for o in ["levels", "shares", "quantiles"]
                data_dict[meas][o] = Dict()
                data_dict[meas][o]["data"] = hcat(fill!(Array{Float64}(undef, size(dd[meas][o]["data"], 1), abs(pushup.value)), NaN), dd[meas][o]["data"], fill!(Array{Float64}(undef, size(dd[meas][o]["data"], 1), abs(pushdown.value)), NaN))
                data_dict[meas][o]["common series"] = Dict()
                data_dict[meas][o]["common series"][top_k] = vcat(fill!(Array{Float64}(undef, abs(pushup.value)), NaN), dd[meas][o]["common series"][top_k], fill!(Array{Float64}(undef, abs(pushdown.value)), NaN))
                data_dict[meas][o]["common series"][mid_k] = vcat(fill!(Array{Float64}(undef, abs(pushup.value)), NaN), dd[meas][o]["common series"][mid_k], fill!(Array{Float64}(undef, abs(pushdown.value)), NaN))
                data_dict[meas][o]["common series"][bot_k] = vcat(fill!(Array{Float64}(undef, abs(pushup.value)), NaN), dd[meas][o]["common series"][bot_k], fill!(Array{Float64}(undef, abs(pushdown.value)), NaN))
            end
        end

    elseif pushup.value < 0 && pushdown.value >= 0
        # Add NaNs to the beginning 
        # data_dict["copulas"]["data"] = hcat(fill!(Array{Float64}(undef, (size(dd["copulas"]["data"], 1), abs(pushup.value)), NaN), dd["copulas"]["data"][:, 1:end-pushdown.value])

        for meas in measures
            data_dict[meas] = Dict()
            for o in ["levels", "shares", "quantiles"]
                data_dict[meas][o] = Dict()
                data_dict[meas][o]["data"] = hcat(fill!(Array{Float64}(undef, size(dd[meas][o]["data"], 1), abs(pushup.value)), NaN), dd[meas][o]["data"])[1:end-pushdown.value]
                data_dict[meas][o]["common series"] = Dict()
                data_dict[meas][o]["common series"][top_k] = vcat(fill!(Array{Float64}(undef, abs(pushup.value)), NaN), dd[meas][o]["common series"][top_k])[1:end-pushdown.value]
                data_dict[meas][o]["common series"][mid_k] = vcat(fill!(Array{Float64}(undef, abs(pushup.value)), NaN), dd[meas][o]["common series"][mid_k])[1:end-pushdown.value]
                data_dict[meas][o]["common series"][bot_k] = vcat(fill!(Array{Float64}(undef, abs(pushup.value)), NaN), dd[meas][o]["common series"][bot_k])[1:end-pushdown.value]
            end
        end
    end

    return data_dict
end

# function get_scf_data(time_p, measures, grid)
#     @unpack year_vec = time_p

#     # Download SCF estimates 
#     init_path               = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
#     scf_estimates           = CSV.read(init_path * "/7_Results/consum_and_income_and_wealth full CEX/from_mcmc/data/SCF_functional_data_A non-diag_.csv", DataFrame)

#     # Decompose this into copulas, pcfs, levels and shares
#     T          = nrow(scf_estimates)
#     D          = length(measures)
#     colons     = ntuple(_ -> (:), D) # TODO: needs to be corrected for 4 dimensional case 
#     cop_size   = tuple([grid for i in 1:D]...)

#     cop_mat_size   = tuple([grid for i in 1:D]..., T)


#     copulas    = fill!(Array{Float64}(undef, cop_mat_size...), NaN) 
#     cop_cols   = filter(x -> occursin("ciw", x), names(scf_estimates)) # TODO: reflect initials of measures 
#     cop_est    = Matrix(select(scf_estimates, cop_cols))
#     data_pcf   = Matrix(transpose(Matrix(select(scf_estimates, filter(x -> occursin("quantiles", x), names(scf_estimates))))))
#     levels     = Matrix(transpose(Matrix(select(scf_estimates, filter(x -> occursin("levels", x), names(scf_estimates))))))
#     shares     = Matrix(transpose(Matrix(select(scf_estimates, filter(x -> occursin("shares", x), names(scf_estimates))))))

#     for t in 1:T
#         copulas[colons..., t] = reshape(cop_est[t, :], cop_size) 
#     end

#     # Putting it all together 
#     scf_estimates_dict                          = create_time_series_dictionary([copulas, data_pcf, levels, shares], grid, sort(measures))
#     scf_estimates_dict["copulas"]["data"]       = Matrix(transpose(cop_est))

#     @unpack tmin, tmax = time_p

#     pushup       = QuarterlyDate(tmin["year"], tmin["quarter"]) - QuarterlyDate(1962, 3) + Quarter(1)
#     pushdown     = QuarterlyDate(2021, 4) - QuarterlyDate(tmax["year"], tmax["quarter"])

#     scf_est_dict  = fill_data_dicts(scf_estimates_dict, measures, pushup, pushdown)

#     return scf_est_dict
# end


function get_estimates_for_comparison(data_name, ty, time_p, measures, estimator)
    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid, integral_cop_grid = estimator
    else
        @unpack grid_pcf, grid_cop = estimator
    end

    grid_choice_pcf = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf
    grid_choice_cop = typeof(estimator) <: SeriesEstimator ? integral_cop_grid : grid_cop

    plot_name = occursin("CEX", data_name) ? "CEX" : data_name

    if data_name != "SCF" && plot_name != "CEX"
        return Dict()
    else
        @unpack year_vec = time_p
        init_path = BASE_PATH
        objects = sort(["quantiles", "levels", "shares"])
        meas_folder = measures_folder(measures)

        # Download CEX_all confidence intervals
        # cex_all                     = jldopen(init_path * "/2_Data_processing/confidence_intervals/ci_deciles_CEX_allnpimp" * ".jld2", "r")
        # cex_confidence_intervals    = cex_all["ci"]

        # Download cex_all data
        # cex_all_data                = jldopen(init_path * "/2_Data_processing/data_consum_and_income_and_wealth_deciles_CEX_allnpimp" * ".jld2", "r")
        # cex_all_data                = cex_all_data["CEX"]
        # measures                    = setdiff(collect(keys(cex_all_data)), ["copulas"])

        estimates_dict = Dict()

        # Download cex_all estimates 
        local est_tags
        if plot_name == "CEX"
            est_tags = ["every 4 years"] # CEX for 'every 4 years', CEX_all for 'additional factors', "PP CEX every 4 years"
        elseif plot_name == "SCF"
            # #TODO: for some reason, each vector of tags cannot overlap ... so e.g., "additional factors", "less factors" cannot be used for "SCF" and "CEX" at the same time
            est_tags = [
                "excluding recent 20 quarters",
                "excluding housing cycle",
                "excluding housing cycle wealth",
                "excluding housing cycle short",
                # "PP CEX excluding housing cycle",
                # "PP CEX excluding housing cycle short",
                # "PP CEX excluding recent 20 quarters",
                "6 factors",
                "7 factors",
                "less DF and AF",
                "less AF",
                "more AF",
                # "Γ estimated",
                # "Γ all",
                # "Γ all 85",
                # "PP CEX",
                # "PP SCF",
                # "PP CEX SCF"
            ]
        elseif plot_name == "HANK"
            est_tags = ["HANK full"]
        else
            est_tags = [""]
        end

        data_tag = ty == "normal" ? "" : "_detrended"

        for tag in est_tags
            local estimates
            # if occursin("every 4 years", tag) && plot_name == "CEX"
            #     estimates = CSV.read(init_path * "/7_Results/$meas_folder " * tag * "/from_mcmc/data/CEX_functional_data" * data_tag * "_A non-diag_.csv", DataFrame)
            # elseif !occursin("every 4 years", tag) && plot_name == "CEX"
            #     estimates = CSV.read(init_path * "/7_Results/$meas_folder " * tag * "/from_mcmc/data/CEX_all_functional_data" * data_tag * "_A non-diag_.csv", DataFrame)
            if plot_name == "CEX"
                estimates = CSV.read(init_path * "/7_Results/$meas_folder " * tag * "/from_mcmc/data/CEX_functional_data" * data_tag * "_A non-diag_.csv", DataFrame)
            elseif tag == "" && plot_name == "SCF"
                estimates = CSV.read(init_path * "/7_Results/$meas_folder" * tag * "/from_mcmc/data/$(plot_name)_functional_data" * data_tag * "_A non-diag_.csv", DataFrame)
            else
                estimates = CSV.read(init_path * "/7_Results/$meas_folder " * tag * "/from_mcmc/data/$(plot_name)_functional_data" * data_tag * "_A non-diag_.csv", DataFrame)
            end


            # Decompose this into copulas, pcfs, levels and shares
            T = nrow(estimates)
            D = length(measures)
            colons = ntuple(_ -> (:), D) # TODO: needs to be corrected for 4 dimensional case 
            cop_size = tuple([grid_choice_cop for i in 1:D]...)

            cop_mat_size = tuple([grid_choice_cop for i in 1:D]..., T)

            copulas = fill!(Array{Float64}(undef, cop_mat_size...), NaN)
            # cop_cols   = filter(x -> occursin("ciw", x), names(estimates))
            # cop_est    = Matrix(select(estimates, cop_cols))

            data_pcf = Matrix(transpose(Matrix(select(estimates, filter(x -> occursin("quantiles", x), names(estimates))))))
            levels = Matrix(transpose(Matrix(select(estimates, filter(x -> occursin("levels", x), names(estimates))))))
            shares = Matrix(transpose(Matrix(select(estimates, filter(x -> occursin("shares", x), names(estimates))))))

            # for t in 1:T
            #     copulas[colons..., t] = reshape(cop_est[t, :], cop_size) 
            # end

            # Putting it all together 
            estimates_dict2 = create_time_series_dictionary([copulas, data_pcf, levels, shares], estimator, sort(measures))
            # estimates_dict2["copulas"]["data"]       = Matrix(transpose(cop_est)) # TODO: check this line out

            @unpack tmin, tmax = time_p

            # Filter cex_all_data to match the overall estimation time frame. cex_all_data is already in quarters and is from 1962Q3 to 2021Q4. cex_all_estimates as well.
            min_time = minimum(estimates[!, "time"]) # the first period of the estimation of the outside estimates 
            max_time = maximum(estimates[!, "time"]) # the last period of the estimation of the outside estimates

            pushup = QuarterlyDate(tmin["year"], tmin["quarter"]) - QuarterlyDate(min_time) # Tricky part is if min_time is later than tmin ... I guess i could add NaNs for missing periods. Makes sense!
            pushdown = QuarterlyDate(max_time) - QuarterlyDate(tmax["year"], tmax["quarter"]) # max_time solved similarly

            estimates_dict[tag] = fill_data_dicts(estimates_dict2, measures, pushup, pushdown, estimator)
            # println(tag)
            # println(min_time)
            # println(max_time)
            # println(pushup)
            # println(pushdown)
        end

        return estimates_dict
    end
end

function select_linestyle(source)
    if source == "WID"
        return :dot
    elseif source == "DFA"
        return :dash
    end
end


function select_color(source; ext=false)
    a = theme_palette(:auto)

    if occursin("SCF", source)
        return a[1]
    elseif source == "PSID"
        return a[2]
    elseif source == "CEX"
        return a[3]
    elseif source == "CPS"
        return a[4]
    elseif source == "CPS2"
        return a[5]
    elseif source == "consensus"
        return a[6]
    elseif source == "WID"
        return :green
    elseif source == "DFA"
        return :blue
    elseif source == "SIPP1"
        return a[7]
    elseif source == "SIPP2"
        return a[8]
    elseif source == "SIPP3"
        return a[9]
    elseif source == "HANK"
        return a[10]
    end
end


function generate_quantiles_shares_levels_plots(data_dict, ty, func_data, data_name, smin, smax, tmin, tmax, estimator, label, type, measures, time_params, select_series, gdp_series, posterior_bounds, compare_to_other_est, tag)

    @unpack func_dict, confidence_intervals, year_vec, data_sources = func_data
    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid, integral_cop_grid = estimator
    else
        @unpack grid_pcf, grid_cop = estimator
    end

    grid_choice_pcf = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf
    grid_choice_cop = typeof(estimator) <: SeriesEstimator ? integral_cop_grid : grid_cop

    local est_dict
    if compare_to_other_est
        est_dict = get_estimates_for_comparison(data_name, ty, time_params, measures, estimator)
    end

    # Import CEX intervals 
    # Dates for CEX all: 1984Q1 to 2021Q4


    # Subset data for the reconstruction
    base_jump, end_jump = find_subset_frame(smin, smax, tmin, tmax)
    m_label = measures_folder(measures)
    folder = type == :optimization ? "from_optimization" : "from_mcmc"
    init_path = BASE_PATH
    path = init_path * "/7_Results/$m_label" * "$tag" * "/$folder/plots/"
    mkpath(path * "correlations/")
    for m in measures
        mkpath(path * "$m/quantiles_levels/")
    end

    local bot, mid, top
    if grid_choice_pcf == 10 || grid_choice_pcf == 100 || grid_choice_pcf == 20
        bot, mid, top = "bottom50", "next40", "top10"
    elseif grid_choice_pcf == 5
        bot, mid, top = "bottom40", "next40", "top20"
    end

    # Get the observed measures
    local obs_meas
    if data_name == "consensus"
        obs_meas = measures
    else
        obs_meas = get_obs_meas(func_dict, data_name, measures, top=top)
    end

    within_stat_dict = Dict()
    plot_name = occursin("CEX", data_name) ? "CEX" : data_name

    # Generate the within statistic for the copulas 
    dimension = length(measures)

    if data_name != "consensus" && ty == "normal" # ty == "average" doesn't work because it has no comparison to the data
        within_stat_dict["copula"] = compute_copula_within_stat(data_dict, confidence_intervals, base_jump, end_jump, data_name, dimension, grid_choice_cop)
    end

    local correlations_dict, correlations_dict_cycle
    if compare_to_other_est && (plot_name == "SCF" || plot_name == "CEX")
        correlations_dict = Dict()
        correlations_dict_cycle = Dict()
        models = collect(keys(est_dict))
        # segments = ["bottom", "middle", "top"]

        # For correlations table 
        for meas in obs_meas
            correlations_dict[meas] = Dict()
            correlations_dict_cycle[meas] = Dict()

            # Example correlation data (replace these with your actual correlations)
            for m in models
                correlations_dict[meas][m] = Dict("bottom" => [NaN, NaN], "middle" => [NaN, NaN], "top" => [NaN, NaN])
                correlations_dict_cycle[meas][m] = Dict("bottom" => [NaN, NaN], "middle" => [NaN, NaN], "top" => [NaN, NaN])
            end
        end
    end

    for meas in obs_meas # TODO: this issue here is that not all measures are observed ... ofc, we can use the reconstructed data but not the confidence intervals

        # All quantiles #TODO: I will hold off on this for now
        qu = data_dict[meas]["quantiles"]["data"][:, base_jump:end-end_jump] ./ select_series[base_jump:end-end_jump, meas*"_per_hh"]'

        local lb_qu, ub_qu
        if posterior_bounds != false
            lb_qu = posterior_bounds["lb"][meas]["quantiles"]["data"][:, base_jump:end-end_jump] ./ select_series[base_jump:end-end_jump, meas*"_per_hh"]'
            ub_qu = posterior_bounds["ub"][meas]["quantiles"]["data"][:, base_jump:end-end_jump] ./ select_series[base_jump:end-end_jump, meas*"_per_hh"]'
        end

        # For the (o)bservations 
        qu_o = Vector{Any}(undef, 3)

        if data_name != "consensus"
            # All quantiles 
            qu_o[1] = func_dict[data_name][meas]["quantiles"]["data"][:, base_jump:end-end_jump] ./ select_series[base_jump:end-end_jump, meas*"_per_hh"]'
            qu_o[2] = confidence_intervals[data_name]["ci_l"][meas]["quantiles"][:, base_jump:end-end_jump] ./ select_series[base_jump:end-end_jump, meas*"_per_hh"]'
            qu_o[3] = confidence_intervals[data_name]["ci_u"][meas]["quantiles"][:, base_jump:end-end_jump] ./ select_series[base_jump:end-end_jump, meas*"_per_hh"]'
        end


        local qu_outside_est
        if compare_to_other_est
            if plot_name == "CEX" || plot_name == "SCF"
                qu_outside_est = Dict()
                for tag in collect(keys(est_dict))
                    # Generate the 'select_series' for these estimations. 'gdp_series' is of the length of the baseline estimation ... other tags are stuffed with NaNs to accomodate this length
                    tag_select_series = generate_average_series(est_dict[tag], gdp_series, measures)

                    qu_outside_est[tag] = est_dict[tag][meas]["quantiles"]["data"][:, base_jump:end-end_jump] ./ tag_select_series[base_jump:end-end_jump, meas*"_per_hh"]'
                end

            end
        end

        # Plots 
        dts = QuarterlyDate(smin["year"], smin["quarter"]):Quarter(1):QuarterlyDate(smax["year"], smax["quarter"])
        markers = [:diamond, :utriangle, :circle, :star5, :cross, :diamond, :utriangle, :circle, :star5, :cross, :diamond, :utriangle, :circle, :star5, :cross, :diamond, :utriangle, :circle, :star5, :cross] #collect(range(0.3, 1, length=length(data_sources)))

        # x-axis
        xaxis = collect(1:length(dts))
        M = uppercasefirst(meas)

        # Important tag 
        detrended_or_not = ty == "average" ? "detrended_" : ""

        # All possible line-styles
        line_styles = [:solid :dash :dot :dashdot :dashdotdot :solid :dash :dot :dashdot :dashdotdot] # has to be a matrix

        within_stat_dict[meas] = Dict()

        # Doing quantiles now 
        label_quantiles = label_qs(grid_choice_pcf) # e.g., L"\textrm{%$(i)th}"
        sequences = define_sequences(grid_choice_pcf) # gets indices of each group
        dist_dict = Dict("bottom" => [sequences[1]], "middle" => [sequences[2]], "top" => [sequences[3]])
        log_qu = log_transformation(deepcopy(qu))


        for (obj, dist) in dist_dict
            within_stat_dict[meas][obj] = Dict()

            cond = data_name != "consensus" ? findall(!isnan, qu_o[1][dist[1][1], :]) : [1] # .!isnan.(all_lv_o[1][:, j])
            s_axis = xaxis[cond[1]:cond[end]] # start at the first observation 
            s_data = qu[:, cond[1]:cond[end]]
            s_log_data = log_qu[:, cond[1]:cond[end]]
            s_dts = dts[cond[1]:cond[end]]

            Plots.plot()
            # Plotting estimates
            for (lsⱼ, j) in enumerate(dist[1])
                Plots.plot!(s_axis,
                    s_data[j, :],
                    ylabel=L"\textrm{%$(M)\, \, rel.\,  to\,\,  average}",
                    lc=obj != "top" ? palette(:glasbey_bw_n256)[j] : :red,
                    xformatter=:latex,
                    yformatter=:latex,
                    xtickfontsize=14,
                    ytickfontsize=14,
                    legendfontsize=10,
                    guidefontsize=14,
                    xticks=(s_axis[1:40:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(s_dts[1:40:end])]),
                    legend=:best,
                    label=obj != "top" ? label_quantiles[:, j][1] : L"\textrm{Model}",
                    lw=obj == "top" ? 4 : 2, dpi=500, ls=obj == "top" ? :solid : line_styles[:, lsⱼ][1],
                )
            end

            # Plotting data to go along with estimates
            if data_name != "consensus" && ty == "normal"
                for j in dist[1]
                    c_data = Vector{Any}(undef, 3)
                    for i in 1:3
                        c_data[i] = qu_o[i][j, :][cond] # c_data = complete, All non-NaN 
                    end

                    # See how many points fall within the confidence intervals
                    r_data = qu[j, :][cond] # estimates that correspond to the indices of the data points 
                    num = count(c_data[2] .<= r_data .<= c_data[3])
                    den = length(r_data)
                    within_stat = floor(Int, (num ./ den) * 100)

                    within_stat_dict[meas][obj][label_quantiles[:, j][1]] = "$num" * "/" * "$den"

                    Plots.scatter!(xaxis[cond],
                        c_data[1],
                        marker=:utriangle, #markers[j],
                        markercolor=:black,
                        ms=5,
                        la=0.5,
                        lw=2, dpi=500,
                        label=L"\textrm{Data}",
                        yerror=(c_data[1] - c_data[2], c_data[3] - c_data[1]),
                    )
                    Plots.plot!([], [], ls=:dash, lc=:black, la=0.0, label=L"\textrm{Within\, \,  bounds: %$(within_stat)\%}",)

                end
            end

            Plots.savefig(path * "$meas/" * "quantiles_levels/" * plot_name * "_$meas" * "_$obj" * "_quantiles_" * detrended_or_not * label * ".pdf")

            # Plotting for the different aggregate groups now 
            # if obj != "top"
            cond = data_name != "consensus" ? findall(!isnan, qu_o[1][dist[1][1], :]) : [1] # .!isnan.(all_lv_o[1][:, j])
            s_axis = xaxis[cond[1]:cond[end]] # start at the first observation 

            # Aggregation
            s_data = vec(sum(qu[dist[1], cond[1]:cond[end]], dims=1)' ./ length(dist[1]))
            s_dts = dts[cond[1]:cond[end]]

            # Ordinary plots
            intd = plot_name == "CEX" || occursin("SIPP", plot_name) ? 20 : 40
            Plots.plot(s_axis,
                s_data,
                ylabel=M == "Consum" ? L"\textrm{Consumption\, \, rel.\,  to\,\, average}" : L"\textrm{%$(M)\, \, rel.\,  to\,\, average}",
                lc=:red, #select_color(plot_name),
                xformatter=:latex,
                yformatter=:latex,
                xticks=(s_axis[1:intd:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(s_dts[1:intd:end])]),
                legend=:best,
                linewidth=4,
                xtickfontsize=14,
                ytickfontsize=14,
                legendfontsize=10,
                guidefontsize=14,
                # label=L"\textrm{Model}",
                label="",
                dpi=500, ls=:solid,
            )

            if data_name != "consensus" && ty == "normal"

                c_data = Vector{Any}(undef, 3) # data + intervals = 3
                for i in 1:3
                    c_data[i] = vec(sum(qu_o[i][dist[1], cond], dims=1)' ./ length(dist[1]))  # c_data = complete, All non-NaN 
                end

                # See how many points fall within the confidence intervals
                r_data = vec(sum(qu[dist[1], cond], dims=1)' ./ length(dist[1])) # estimates that correspond to these data points 
                within_stat = floor(Int, (count(c_data[2] .<= r_data .<= c_data[3]) ./ length(r_data)) * 100)

                if plot_name == "CEX"
                    Plots.plot!(xaxis[cond],
                        c_data[2],
                        fillrange=c_data[3],
                        fillalpha=0.2,
                        fillcolor=:red,
                        la=0.0,
                        lc=:white,
                        lw=4, dpi=500,
                        label="",
                    )
                    Plots.scatter!(xaxis[cond],
                        c_data[1],
                        marker=markers[2],
                        ms=5,
                        markercolor=:black,
                        la=0.5,
                        lw=2, dpi=500,
                        label="",
                    )
                else
                    Plots.scatter!(xaxis[cond],
                        c_data[1],
                        marker=markers[2],
                        ms=5,
                        markercolor=:black,
                        la=0.5,
                        lw=2, dpi=500,
                        # label=L"\textrm{Data}",
                        label="",
                        yerror=(c_data[1] - c_data[2], c_data[3] - c_data[1]),
                    )
                end

                Plots.plot!([], [], ls=:dash, lc=:black, la=0.0, label=L"\textrm{Within\, \,  bounds: %$(within_stat)\%}",)
            end

            Plots.savefig(path * "$meas/" * "quantiles_levels/" * plot_name * "_$meas" * "_$obj" * "_quantilegroup_" * detrended_or_not * label * ".pdf")

            # Log-HP now
            cond = data_name != "consensus" ? findall(!isnan, qu_o[1][dist[1][1], :]) : [1] # .!isnan.(all_lv_o[1][:, j])
            s_axis = xaxis[cond[1]:cond[end]] # start at the first observation 
            s_data = log_transformation(vec(sum(qu[dist[1], cond[1]:cond[end]], dims=1)' ./ length(dist[1])))

            local s_data_cycle
            try
                # s_data_cycle = s_data .- HP(s_data, 1600)
                s_data_cycle = 100 .* HP(s_data .- HP(s_data, 1600), 6)
            catch ee
                s_data_cycle = copy(s_data)
            end

            s_dts = dts[cond[1]:cond[end]]

            Plots.plot(s_axis,
                s_data_cycle, # to have it in percent
                ylabel=M == "Consum" ? L"\% \Delta \textrm{\,\, Consumption}" : L"\% \Delta \textrm{\,\,%$(M)}",
                lc=:red, #select_color(plot_name),
                xformatter=:latex,
                yformatter=:latex,
                xticks=(s_axis[1:40:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(s_dts[1:40:end])]),
                legend=:best,
                linewidth=4,
                xtickfontsize=14,
                ytickfontsize=14,
                legendfontsize=10,
                guidefontsize=14,
                label=L"\textrm{Model}",
                dpi=500, ls=:solid,
            )

            # No data counterpart because you cannot HP data

            Plots.savefig(path * "$meas/" * "quantiles_levels/" * plot_name * "_$meas" * "_$obj" * "_loghp_quantilegroup_" * detrended_or_not * label * ".pdf")


            if compare_to_other_est
                if plot_name == "CEX" || plot_name == "SCF"
                    for esttag in collect(keys(qu_outside_est))
                        cond_cex = findall(!isnan, qu_o[1][dist[1][1], :]) # indices of observations used in current estimation    
                        cond = plot_name != "consensus" ? findall(!isnan, qu_o[1][dist[1][1], :]) : [1]
                        s_axis = xaxis[cond[1]:cond[end]] # start at the first observation 
                        s_data = qu_outside_est[esttag][:, cond[1]:cond[end]]
                        s_dts = dts[cond[1]:cond[end]]

                        # First, plot OUTSIDE estimates
                        Plots.plot()
                        intd = plot_name == "CEX" || occursin("SIPP", plot_name) ? 20 : 40
                        for (lsⱼ, j) in enumerate(dist[1])
                            Plots.plot!(s_axis,
                                s_data[j, :],
                                ylabel=M == "Consum" ? L"\textrm{Consumption\, \, rel.\,  to\,\,  average}" : L"\textrm{%$(M)\, \, rel.\,  to\,\,  average}",
                                lc=obj != "top" ? palette(:glasbey_bw_n256)[j] : select_color(plot_name),
                                xformatter=:latex,
                                yformatter=:latex,
                                xticks=(s_axis[1:intd:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(s_dts[1:intd:end])]),
                                legend=:best,
                                xtickfontsize=14,
                                ytickfontsize=14,
                                legendfontsize=10,
                                guidefontsize=14,
                                label=label_quantiles[:, j][1],
                                lw=4, dpi=500, ls=line_styles[:, lsⱼ][1],
                            )
                        end

                        # Next the idea is to plot the confidence intervals/data, which is at some kind of yearly 
                        local cond_ci, non_overlap_ids, overlap_ids
                        if ty == "normal" # For average, we currently do not have a comparison to the data
                            for j in dist[1]
                                c_data = Vector{Any}(undef, 3)

                                for i in eachindex(c_data)
                                    c_data[i] = qu_o[i][j, cond_cex]
                                end

                                # See how many points fall within the confidence intervals
                                r_data = qu_outside_est[esttag][j, cond_cex]
                                within_stat = floor(Int, (count(c_data[2] .<= r_data .<= c_data[3]) ./ length(r_data)) * 100)

                                # Correlate estimates with only the data that didn't go into the model
                                other_model_tags = [
                                    "6 factors",
                                    "7 factors",
                                    "less DF and AF",
                                    "less AF",
                                    "more AF",
                                    "Γ estimated",
                                    "Γ all",
                                    "Γ all 85",
                                    "PP CEX",
                                    "PP SCF",
                                    "PP CEX SCF"
                                ]

                                local ρ, housing_ids, last20_ids
                                if occursin("every 4 years", esttag)
                                    # Idea here is we want to get the last observation the same across both estimations -> move data point in starved model to the right 1
                                    # non_overlap_ids = [j for (i, j) in enumerate(cond_cex) if i % 4 != 2]
                                    # overlap_ids = [j for (i, j) in enumerate(cond_cex) if i % 4 == 2]
                                    non_overlap_ids = [j for (i, j) in enumerate(cond_cex) if (i - 1) % 16 ≥ 4]
                                    overlap_ids = [j for (i, j) in enumerate(cond_cex) if (i - 1) % 16 < 4]

                                    ρ = round(cor(qu_o[1][j, non_overlap_ids], qu_outside_est[esttag][j, non_overlap_ids]), digits=2) # correlate data with estimates 
                                elseif occursin("excluding housing cycle", esttag)
                                    # s_dts           = dts[cond[1]:cond[end]]
                                    housing_ids = findall(x -> x >= QuarterlyDate(2004, 1) && x <= QuarterlyDate(2009, 4), dts)
                                    non_overlap_ids = [i for i in cond_cex if i ∈ housing_ids]
                                    overlap_ids = [i for i in cond_cex if i ∉ housing_ids]
                                    ρ = round(cor(qu_o[1][j, non_overlap_ids], qu_outside_est[esttag][j, non_overlap_ids]), digits=2) # correlate data with estimates 
                                elseif occursin("excluding housing cycle short", esttag)
                                    housing_ids = findall(x -> x >= QuarterlyDate(2007, 4) && x <= QuarterlyDate(2011, 4), dts)
                                    non_overlap_ids = [i for i in cond_cex if i ∈ housing_ids]
                                    overlap_ids = [i for i in cond_cex if i ∉ housing_ids]
                                    ρ = non_overlap_ids == [] ? NaN : round(cor(qu_o[1][j, non_overlap_ids], qu_outside_est[esttag][j, non_overlap_ids]), digits=2) # correlate data with estimates 
                                elseif occursin("excluding recent 20 quarters", esttag)
                                    last20_ids = findall(x -> x >= QuarterlyDate(2020, 1) && x <= QuarterlyDate(2024, 4), dts)
                                    non_overlap_ids = [i for i in cond_cex if i ∈ last20_ids]
                                    overlap_ids = [i for i in cond_cex if i ∉ last20_ids]
                                    ρ = round(cor(qu_o[1][j, non_overlap_ids], qu_outside_est[esttag][j, non_overlap_ids]), digits=2) # correlate data with estimates 
                                elseif esttag ∈ other_model_tags
                                    non_overlap_ids = [i for i in cond_cex if i ∈ []]
                                    overlap_ids = [i for i in cond_cex] # basically all data points
                                    ρ = round(cor(qu_o[1][j, overlap_ids], qu_outside_est[esttag][j, overlap_ids]), digits=2) # correlate everything since all overlap 
                                end

                                # Plotting the confidence intervals 
                                Plots.plot!(xaxis[cond_cex],
                                    c_data[2],
                                    fillrange=c_data[3],
                                    fillalpha=0.2,
                                    fillcolor=obj != "top" ? palette(:glasbey_bw_n256)[j] : select_color(plot_name),
                                    la=0.0,
                                    lc=:white,
                                    lw=4, dpi=500,
                                    label=L"\textrm{Within\, \,  bounds: %$(within_stat)\%}"
                                )

                                # Plotting the CEX data 
                                Plots.scatter!(xaxis[cond_cex],
                                    c_data[1],
                                    marker=markers[j],
                                    markercolor=:black,
                                    markersize=5,
                                    la=0.5,
                                    lw=4, dpi=500,
                                    label=L"\textrm{Corr. \,\, with\,\, data: %$(ρ)}",
                                )
                            end
                        end

                        Plots.savefig(path * "$meas/" * "quantiles_levels/" * plot_name * "_$meas" * "_$obj" * "_quantiles_" * esttag * "_" * detrended_or_not * label * ".pdf")

                        # Do it for all groups at once
                        if obj == "top"
                            # grid_choice_pcf = 20
                            sequences = define_sequences(grid_choice_pcf)
                            dist_dict = Dict("bottom" => [sequences[1]], "middle" => [sequences[2]], "top" => [sequences[3]])

                            # legend_labels     = [L"\textrm{Bottom\, 50 \,--\,\,}", L"\textrm{Next\,\, 40}", L"\textrm{Top\,\, 10}"]
                            s_data_all = qu[:, cond[1]:cond[end]]
                            plot_tag = split(esttag, " ")
                            plot_tag = join(plot_tag, "\\,\\,")

                            file_name = init_path * "/7_Results/$m_label" * "$tag" * "/$folder/data/" * plot_name * "_bounds.jld2"
                            local lb, ub
                            if (plot_name == "SCF" || plot_name == "CEX") && isfile(file_name) && ty == "normal"
                                # Now plots with the posterior bounds 

                                # Get relevant data and subset
                                bounds = jldopen(file_name, "r")
                                lb = bounds["lb"][meas][:, cond[1]:cond[end]]
                                ub = bounds["ub"][meas][:, cond[1]:cond[end]]
                            end

                            seg_dict = Dict("bottom" => 1, "middle" => 2, "top" => 3)

                            for (k, g) in dist_dict
                                n_lines = length(collect(g...))
                                outside_est = vec(sum(s_data[g..., :], dims=1)' ./ n_lines) # outside estimates 
                                cex_all_est = vec(sum(s_data_all[g..., :], dims=1)' ./ n_lines) # conditional on all data 

                                Plots.plot()

                                Plots.plot!(s_axis,
                                    cex_all_est,
                                    lc=:red,
                                    # label=L"\textrm{Baseline}",
                                    label="",
                                    lw=4, dpi=500, ls=:solid,
                                    legend=false
                                )

                                Plots.plot!(s_axis,
                                    outside_est,
                                    ylabel=M == "Consum" ? L"\textrm{Consumption\, \, rel.\,  to\,\,  average}" : L"\textrm{%$(M)\, \, rel.\,  to\,\,  average}",
                                    lc=:blue,
                                    xformatter=:latex,
                                    yformatter=:latex,
                                    lw=4,
                                    la=0.5,
                                    xtickfontsize=14,
                                    ytickfontsize=14,
                                    legendfontsize=10,
                                    guidefontsize=14,
                                    xticks=(s_axis[1:40:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(s_dts[1:40:end])]),
                                    legend=:best, # BUG: with twinx(), legend=:outertopright widens the plot in a weird way
                                    # label=esttag == "7 factors" ? L"\textrm{7 \,\, Factors}" : esttag == "6 factors" ? L"\textrm{6 \,\, Factors}" : esttag == "less factors" ? L"\textrm{Less \,\, Factors}" : esttag == "higher order15" ? L"\textrm{Higher \,\, Order}" : esttag == "less AF" ? L"\textrm{Less \,\,Agg.\,\, Factors}" : esttag == "less DF and AF" ? L"\textrm{Compact\,\, Model}" : esttag == "Γ estimated" ? L"\textrm{Γ\,\, estimated}" : esttag == "Γ all" ? L"\textrm{Γ-10}" : esttag == "more AF" ? L"\textrm{More\,\,Agg.\,\,Factors}" : esttag == "Γ all 85" ? L"\textrm{Γ-12}" : esttag == "PP CEX" ? L"\textrm{\Gamma_{aug}\,\, CEX}" : esttag == "PP SCF" ? L"\textrm{\Gamma_{aug}\,\, SCF}" : esttag == "PP CEX SCF" ? L"\textrm{\Gamma_{aug}\,\, CEX-SCF}" : L"\textrm{Less \,\, Data}",
                                    label="",
                                    dpi=500, ls=:dot,
                                )


                                # Perform the within stat thing, once for the data and once for the cycle
                                if ty == "normal"
                                    r_data_all = sum(qu[g..., cond_cex], dims=1)' ./ n_lines
                                    c_int_ub(x) = sum(qu_o[3][g..., x], dims=1)' ./ n_lines
                                    c_int_lb(x) = sum(qu_o[2][g..., x], dims=1)' ./ n_lines
                                    c_data(x) = sum(qu_o[1][g..., x], dims=1)' ./ n_lines
                                    r_data(x) = sum(qu_outside_est[esttag][g..., x], dims=1)' ./ n_lines

                                    # Statistics
                                    within_stat = floor(Int, (count(c_int_lb(cond_cex) .<= r_data(cond_cex) .<= c_int_ub(cond_cex)) ./ length(r_data(cond_cex))) * 100)

                                    other_model_tags = [
                                        "6 factors",
                                        "7 factors",
                                        "less DF and AF",
                                        "less AF",
                                        "more AF",
                                        "Γ estimated",
                                        "Γ all",
                                        "Γ all 85"
                                    ]

                                    if occursin("every 4 years", esttag) || occursin("excluding housing cycle", esttag) || occursin("excluding housing cycle short", esttag) || occursin("excluding recent 20 quarters", esttag)
                                        local ids_to_use
                                        if occursin("excluding recent 20 quarters", esttag) # since there is no extrapolation, I must stop it at the last observation
                                            ids_to_use = findall(x -> x >= QuarterlyDate(2020, 1) && x <= QuarterlyDate(2024, 4), dts)
                                            filter!(x -> x >= cond[1], ids_to_use)
                                            filter!(x -> x <= cond[end], ids_to_use)
                                            ids_to_use = (ids_to_use .- cond[1]) .+ 1
                                            # ids_to_use = collect(ids_to_use[1]:ids_to_use[end]) # Since the other esttags are in a single interval. every 4 years is multiple intervals.
                                        elseif occursin("excluding housing cycle", esttag)  # this one I need to just stop at the end of the housing cycle
                                            ids_to_use = findall(x -> x >= QuarterlyDate(2004, 1) && x <= QuarterlyDate(2009, 4), dts)
                                            ids_to_use = (ids_to_use .- cond[1]) .+ 1
                                        elseif occursin("excluding housing cycle short", esttag) # this one I need to just stop at the end of the housing cycle
                                            ids_to_use = findall(x -> x >= QuarterlyDate(2007, 4) && x <= QuarterlyDate(2011, 4), dts)
                                            ids_to_use = (ids_to_use .- cond[1]) .+ 1
                                        elseif occursin("every 4 years", esttag)
                                            ids_to_use = collect(1:length(dts))
                                            filter!(x -> x ∉ overlap_ids, ids_to_use)
                                            filter!(x -> x >= cond[1], ids_to_use)
                                            filter!(x -> x <= cond[end], ids_to_use)
                                            ids_to_use = (ids_to_use .- cond[1]) .+ 1
                                        end
                                        # last20_ids      = findall(x -> x >= QuarterlyDate(2020, 1) && x <= QuarterlyDate(2024, 4), dts)  
                                        # ids_to_use = esttag == "excluding housing cycle" ? housing_ids : esttag == "excluding recent 20 quarters" ? last20_ids : esttag == "every 4 years" ? periods_of_int : []
                                        # println(esttag)
                                        # println(meas)
                                        # println(ids_to_use)
                                        # println(cex_all_est)
                                        # println(outside_est)
                                        ρ1 = round(nancor(cex_all_est, outside_est), digits=2) # correlate estimates with estimates of full CEX model
                                        ρ2 = round(nancor(cex_all_est[ids_to_use], outside_est[ids_to_use]), digits=2) # correlate estimates with estimates of full CEX model

                                        # Correlation for entire series and subsection
                                        println([ρ1, ρ2])
                                        correlations_dict[meas][esttag][k] .= [ρ1, ρ2]
                                    end

                                    # Only for plotting the additional correlation label
                                    #     Plots.plot!(xaxis[cond_cex],
                                    #         c_data(cond_cex),
                                    #         la=0.0,
                                    #         label=L"\textrm{Corr. \,\, of\,\, (A)\,\,with\,\,(C): %$(ρ)}",
                                    #     )

                                    #     Plots.plot!(xaxis[cond_cex],
                                    #     c_data(cond_cex),
                                    #     la=0.0,
                                    #     label=L"\textrm{Corr. \,\, between \,\, (A)\,\,and \,\, (B): %$(ρ3)}",
                                    # )

                                    # First plot everything as is 
                                    if plot_name == "CEX"
                                        Plots.plot!(xaxis[cond_cex],
                                            c_int_lb(cond_cex),
                                            fillrange=c_int_ub(cond_cex),
                                            la=0.0,
                                            fillalpha=0.2,
                                            fillcolor=:red,
                                            dpi=500,
                                            label=""
                                        )

                                        Plots.scatter!(xaxis[cond_cex],
                                            c_data(cond_cex),
                                            marker=:dot,
                                            markercolor=:black,
                                            markersize=5,
                                            dpi=500,
                                            label=""
                                            # label=g == 1 ? L"\textrm{Data}" : "",
                                            # yerror=(c_int_ub(cond_cex) - c_data(cond_cex), c_data(cond_cex) - c_int_lb(cond_cex)),
                                        )

                                        if non_overlap_ids != []
                                            Plots.scatter!(xaxis[non_overlap_ids],
                                                c_data(non_overlap_ids),
                                                mc=:white, msc=:black, msw=3, #(5, :white, stroke(1, :black)),
                                                dpi=500,
                                                label="",
                                                # label=g == 1 ? L"\textrm{Missing\,\, data}" : "",
                                            )
                                        end
                                    else
                                        Plots.scatter!(xaxis[cond_cex],
                                            c_data(cond_cex),
                                            marker=:dot,
                                            markercolor=:black,
                                            markersize=5,
                                            la=0.5,
                                            lw=4, dpi=500,
                                            label="",
                                            # label=g == 1 ? L"\textrm{Data}" : "",
                                            yerror=(c_int_ub(cond_cex) - c_data(cond_cex), c_data(cond_cex) - c_int_lb(cond_cex)),
                                        )

                                        if non_overlap_ids != []
                                            Plots.scatter!(xaxis[non_overlap_ids],
                                                c_data(non_overlap_ids),
                                                mc=:white, msc=:black, msw=3, #(5, :white, stroke(1, :black)),
                                                la=0.5,
                                                lw=4, dpi=500,
                                                label="",
                                                # label=g == 1 ? L"\textrm{Missing\,\, data}" : "",
                                            )
                                        end
                                    end
                                end

                                Plots.savefig(path * "$meas/" * "quantiles_levels/" * plot_name * "_$(meas)_" * "$(k)_" * "$(esttag)_" * detrended_or_not * label * ".pdf")

                                if (plot_name == "SCF" || plot_name == "CEX") && isfile(file_name) && ty == "normal"

                                    Plots.plot!(s_axis,
                                        lb[seg_dict[k], :],
                                        fillrange=ub[seg_dict[k], :],
                                        fillalpha=0.2,
                                        fillcolor=:red,
                                        la=0.0,
                                        lc=:white,
                                        lw=4, dpi=500,
                                        label=L"\textrm{Credible\,\, Int.}"
                                    )
                                    Plots.savefig(path * "$meas/" * "quantiles_levels/" * plot_name * "_$(meas)_" * "$(k)_" * "$(esttag)_" * detrended_or_not * "w_bounds_" * label * ".pdf")
                                end
                            end
                        end
                    end

                    # Compare to other estimates, but now in logs and HP 
                    println("Second comparison")

                    for esttag in collect(keys(qu_outside_est))
                        cond_cex = findall(!isnan, qu_o[1][dist[1][1], :]) # indices of observations used in current estimation    
                        cond = plot_name != "consensus" ? cond_cex : [1]
                        s_axis = xaxis[cond[1]:cond[end]] # start at the first observation 
                        s_dts = dts[cond[1]:cond[end]]
                        s_data = qu_outside_est[esttag][:, cond[1]:cond[end]]

                        log_s_data = log_transformation(s_data)

                        local s_data_cycle
                        try
                            for j in axes(s_data_cycle)
                                s_data_cycle[j, :] = 100 .* HP(log_s_data[j, :] .- HP(log_s_data[j, :], 1600), 6)
                            end
                        catch ee
                            s_data_cycle = copy(log_s_data)
                        end

                        # First, plot estimates 
                        Plots.plot()
                        for (lsⱼ, j) in enumerate(dist[1])
                            Plots.plot!(s_axis,
                                s_data_cycle[j, :],
                                ylabel=M == "Consum" ? L"\% \Delta \textrm{\,\,Consumption}" : L"\% \Delta \textrm{\,\, %$(M)}",
                                lc=obj != "top" ? palette(:glasbey_bw_n256)[j] : select_color(plot_name),
                                xformatter=:latex,
                                yformatter=:latex,
                                xtickfontsize=14,
                                ytickfontsize=14,
                                legendfontsize=10,
                                guidefontsize=14,
                                xticks=(s_axis[1:40:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(s_dts[1:40:end])]),
                                legend=:best,
                                label=label_quantiles[:, j][1],
                                lw=4, dpi=500, ls=line_styles[:, lsⱼ][1],
                            )
                        end

                        # Next the idea is to plot the confidence intervals, which is annual 
                        # The data is also annual 
                        # for j in dist[1]
                        #     # Correlate estimates with only the data that didn't go into the model
                        #     other_model_tags = ["6 factors", "higher order15", "less DF and AF", "less AF", "more AF", "Γ estimated", "Γ all", ""]

                        # local non_overlap_ids, overlap_ids, ids_in_question
                        # local ρ
                        # if esttag == "every 4 years"
                        #     # Idea here is we want to get the last observation the same across both estimations -> move data point in starved model to the right 1
                        #     non_overlap_ids = [j for (i, j) in enumerate(cond_cex) if i % 4 != 2]
                        #     overlap_ids     = [j for (i, j) in enumerate(cond_cex) if i % 4 == 2]
                        #     ρ               = round(cor(qu_o[1][j, non_overlap_ids], qu_outside_est[esttag][j, non_overlap_ids]), digits=2) # correlate data with estimates 
                        # elseif esttag == "excluding housing cycle"
                        #     ids_in_question     = findall(x -> x >= QuarterlyDate(2004, 1) && x <= QuarterlyDate(2009, 4), dts)
                        #     non_overlap_ids = [i for i in cond_cex if i ∈ ids_in_question]
                        #     overlap_ids     = [i for i in cond_cex if i ∉ ids_in_question]
                        #     ρ               = round(cor(qu_o[1][j, non_overlap_ids], qu_outside_est[esttag][j, non_overlap_ids]), digits=2) # correlate data with estimates 
                        # elseif esttag == "excluding recent 20 quarters"
                        #     ids_in_question      = findall(x -> x >= QuarterlyDate(2020, 1) && x <= QuarterlyDate(2024, 4), dts)  
                        #     non_overlap_ids = [i for i in cond_cex if i ∈ ids_in_question]
                        #     overlap_ids     = [i for i in cond_cex if i ∉ ids_in_question]
                        #     ρ               = round(cor(qu_o[1][j, non_overlap_ids], qu_outside_est[esttag][j, non_overlap_ids]), digits=2) # correlate data with estimates 
                        # elseif esttag ∈ other_model_tags
                        #     non_overlap_ids = [i for i in cond_cex if i ∈ []]
                        #     overlap_ids     = [i for i in cond_cex] # basically all data points
                        #     ρ               = round(cor(qu_o[1][j, overlap_ids], qu_outside_est[esttag][j, overlap_ids]), digits=2) # correlate everything since all overlap 
                        # end

                        # c_data = Vector{Any}(undef, 3)

                        # for i in eachindex(c_data)
                        #     c_data[i]  = log_transformation(qu_o[i][j, cond_cex])
                        #     c_data[i]  = 100 * HP(c_data[i] .- HP(c_data[i], 1600), 6) # non-sensical here, since we cannot HP filter the data
                        # end

                        # # See how many points fall within the confidence intervals
                        # r_data      = log_transformation(qu_outside_est[esttag][j, cond_cex])
                        # r_data      = 100 * HP(r_data .- HP(r_data, 1600), 6)
                        # within_stat = floor(Int, (count(c_data[2] .<= r_data .<= c_data[3]) ./ length(r_data)) * 100)


                        # if ty == "normal" # For average, we currently do not have a comparison to the data
                        #     # Plotting the confidence intervals 
                        #     Plots.plot!(xaxis[cond_cex],
                        #         c_data[2],
                        #         fillrange = c_data[3],
                        #         fillalpha = 0.1,
                        #         fillcolor = obj != "top" ? palette(:glasbey_bw_n256)[j] : select_color(plot_name),
                        #         la=0.0,
                        #         lc=:white,
                        #         lw=4, dpi=500,
                        #         label=L"\textrm{Within\, \,  bounds: %$(within_stat)\%}"
                        #         )

                        #     # Plotting the CEX data 
                        #     Plots.scatter!(xaxis[cond_cex],
                        #         c_data[1],
                        #         marker=markers[j],
                        #         markercolor=:black,
                        #         markersize=5,
                        #         la=0.5,
                        #         lw=4, dpi=500,
                        #         label=L"\textrm{Corr. \,\, with\,\, data: %$(ρ)}",
                        #     )
                        # end        
                        # end
                        Plots.savefig(path * "$meas/" * "quantiles_levels/" * plot_name * "_$meas" * "_$obj" * "_loghp_quantiles_" * esttag * "_" * detrended_or_not * label * ".pdf")

                        # Do it for all groups at once
                        if obj == "top"
                            sequences = define_sequences(grid_choice_pcf)
                            dist_dict = Dict("bottom" => [sequences[1]], "middle" => [sequences[2]], "top" => [sequences[3]])

                            plot_tag = split(esttag, " ")
                            plot_tag = join(plot_tag, "\\,\\,")
                            s_data_all = qu[:, cond[1]:cond[end]]

                            for (k, g) in dist_dict
                                n_lines = length(collect(g...))
                                outside_est = log_transformation(vec(sum(s_data[g..., :], dims=1)' ./ n_lines))

                                # if the 'outside_est' has NaNs, remove them, do the HP stuff and then put them back in
                                temp_outside_est = filter(x -> !isnan(x), outside_est) # only in the beginning for 'every 4 years'
                                nn = length(outside_est) - length(temp_outside_est)
                                NaNs_to_add = [NaN for _ in 1:nn]

                                temp_outside_est = 100 .* HP(temp_outside_est .- HP(temp_outside_est, 1600), 6)
                                outside_est = [NaNs_to_add; temp_outside_est]
                                cex_all_est = log_transformation(vec(sum(s_data_all[g..., :], dims=1)' ./ n_lines))
                                cex_all_est = 100 .* HP(cex_all_est .- HP(cex_all_est, 1600), 6)

                                Plots.plot()

                                Plots.plot!(s_axis,
                                    cex_all_est,
                                    lc=:red,
                                    label=L"\textrm{Baseline}",
                                    lw=4, dpi=500, ls=:solid,
                                )

                                Plots.plot!(s_axis,
                                    outside_est,
                                    ylabel=M == "Consum" ? L"\% \Delta \textrm{\,\, Consumption}" : L"\% \Delta \textrm{\,\, %$(M)}",
                                    lc=:blue,
                                    xformatter=:latex,
                                    yformatter=:latex,
                                    lw=4,
                                    la=0.5,
                                    xtickfontsize=14,
                                    ytickfontsize=14,
                                    legendfontsize=10,
                                    guidefontsize=14,
                                    xticks=(s_axis[1:40:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(s_dts[1:40:end])]),
                                    legend=:best,
                                    label=esttag == "7 factors" ? L"\textrm{7 \,\, Factors}" : esttag == "6 factors" ? L"\textrm{6 \,\, Factors}" : esttag == "less factors" ? L"\textrm{Less \,\, Factors}" : esttag == "higher order15" ? L"\textrm{Higher \,\, Order}" : esttag == "less AF" ? L"\textrm{Less \,\,Agg.\,\, Factors}" : esttag == "less DF and AF" ? L"\textrm{Compact\,\, Model}" : esttag == "Γ estimated" ? L"\textrm{Γ\,\, estimated}" : esttag == "Γ all" ? L"\textrm{Γ-10}" : esttag == "more AF" ? L"\textrm{More\,\,Agg.\,\,Factors}" : esttag == "Γ all 85" ? L"\textrm{Γ-12}" : esttag == "PP CEX" ? L"\textrm{\Gamma_{aug}\,\, CEX}" : esttag == "PP SCF" ? L"\textrm{\Gamma_{aug}\,\, SCF}" : esttag == "PP CEX SCF" ? L"\textrm{\Gamma_{aug}\,\, CEX-SCF}" : L"\textrm{Less \,\, Data}",
                                    dpi=500, ls=:dot,
                                )



                                # Correlations 
                                other_model_tags = [
                                    "6 factors",
                                    "7 factors",
                                    "less DF and AF",
                                    "less AF",
                                    "more AF",
                                    "Γ estimated",
                                    "Γ all",
                                    "Γ all 85",
                                    "PP CEX",
                                    "PP SCF",
                                    "PP CEX SCF"
                                ]

                                local non_overlap_ids, overlap_ids
                                if occursin("every 4 years", esttag)
                                    # Idea here is we want to get the last observation the same across both estimations -> move data point in starved model to the right 1
                                    # non_overlap_ids = [j for (i, j) in enumerate(cond_cex) if i % 4 != 2]
                                    # overlap_ids = [j for (i, j) in enumerate(cond_cex) if i % 4 == 2]
                                    non_overlap_ids = [j for (i, j) in enumerate(cond_cex) if (i - 1) % 16 < 4]
                                    overlap_ids = [j for (i, j) in enumerate(cond_cex) if (i - 1) % 16 ≥ 4]
                                elseif occursin("excluding housing cycle", esttag)
                                    ids_in_question = findall(x -> x >= QuarterlyDate(2004, 1) && x <= QuarterlyDate(2009, 4), dts)
                                    non_overlap_ids = [i for i in cond_cex if i ∈ ids_in_question]
                                    overlap_ids = [i for i in cond_cex if i ∉ ids_in_question]
                                elseif occursin("excluding housing cycle short", esttag)
                                    ids_in_question = findall(x -> x >= QuarterlyDate(2007, 4) && x <= QuarterlyDate(2011, 4), dts)
                                    non_overlap_ids = [i for i in cond_cex if i ∈ ids_in_question]
                                    overlap_ids = [i for i in cond_cex if i ∉ ids_in_question]
                                elseif occursin("excluding recent 20 quarters", esttag)
                                    ids_in_question = findall(x -> x >= QuarterlyDate(2020, 1) && x <= QuarterlyDate(2024, 4), dts)
                                    non_overlap_ids = [i for i in cond_cex if i ∈ ids_in_question]
                                    overlap_ids = [i for i in cond_cex if i ∉ ids_in_question]
                                elseif esttag ∈ other_model_tags
                                    non_overlap_ids = [i for i in cond_cex if i ∈ []]
                                    overlap_ids = [i for i in cond_cex] # basically all data points
                                end

                                if occursin("excluding housing cycle", esttag) || occursin("excluding housing cycle short", esttag) || occursin("excluding recent 20 quarters", esttag)
                                    ymin = (nanminimum(outside_est) .- 0.01 * nanminimum(outside_est)) # plus because it is a negative number
                                    ymax = (nanmaximum(outside_est) .+ 0.01 * nanmaximum(outside_est))

                                    Plots.plot!(
                                        xaxis[ids_in_question],
                                        repeat([ymin], length(ids_in_question)),
                                        fillrange=ymax,
                                        fillalpha=0.3,
                                        color=:gray,
                                        label=g == 1 ? L"\textrm{Missing\,\, data}" : "",
                                    )
                                end
                                Plots.savefig(path * "$meas/" * "quantiles_levels/" * plot_name * "_$(meas)_" * "$(k)_" * "$(esttag)_loghp_" * detrended_or_not * label * ".pdf")

                                if occursin("excluding housing cycle", esttag) || occursin("excluding housing cycle short", esttag) || occursin("excluding recent 20 quarters", esttag) || occursin("every 4 years", esttag)
                                    local ids_to_use
                                    if occursin("excluding recent 20 quarters", esttag) # since there is no extrapolation, I must stop it at the last observation
                                        ids_to_use = findall(x -> x >= QuarterlyDate(2020, 1) && x <= QuarterlyDate(2024, 4), dts)
                                        filter!(x -> x >= cond[1], ids_to_use)
                                        filter!(x -> x <= cond[end], ids_to_use)
                                        ids_to_use = (ids_to_use .- cond[1]) .+ 1
                                        # ids_to_use = collect(ids_to_use[1]:ids_to_use[end]) # Since the other esttags are in a single interval. every 4 years is multiple intervals.
                                    elseif occursin("excluding housing cycle", esttag) # this one I need to just stop at the end of the housing cycle
                                        ids_to_use = findall(x -> x >= QuarterlyDate(2004, 1) && x <= QuarterlyDate(2009, 4), dts)
                                        ids_to_use = (ids_to_use .- cond[1]) .+ 1
                                    elseif occursin("excluding housing cycle short", esttag) # this one I need to just stop at the end of the housing cycle
                                        ids_to_use = findall(x -> x >= QuarterlyDate(2007, 4) && x <= QuarterlyDate(2011, 4), dts)
                                        ids_to_use = (ids_to_use .- cond[1]) .+ 1
                                    elseif occursin("every 4 years", esttag)
                                        ids_to_use = collect(1:length(dts))
                                        filter!(x -> x ∉ overlap_ids, ids_to_use)
                                        filter!(x -> x >= cond[1], ids_to_use)
                                        filter!(x -> x <= cond[end], ids_to_use)
                                        ids_to_use = (ids_to_use .- cond[1]) .+ 1
                                    end

                                    ρ1 = round(nancor(cex_all_est, outside_est), digits=2) # correlate estimates with estimates of full CEX model
                                    ρ2 = round(nancor(cex_all_est[ids_to_use], outside_est[ids_to_use]), digits=2) # correlate estimates with estimates of full CEX model

                                    # Correlation for entire series and subsection
                                    correlations_dict_cycle[meas][esttag][k] .= [ρ1, ρ2]
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # Build correlations table from dictionary 
    if compare_to_other_est && ty == "normal" && (plot_name == "SCF" || plot_name == "CEX")
        # Save the dictionary to a .jld2 file
        jldsave(path * "correlations/" * "correlations_dict_" * data_name * ".jld2"; corr=correlations_dict)
        jldsave(path * "correlations/" * "correlations_dict_cycle_" * data_name * ".jld2"; corr=correlations_dict_cycle)

        generate_correlations_table(path, data_name, obs_meas, correlations_dict)
        generate_correlations_table(path, data_name * "_cycle", obs_meas, correlations_dict_cycle)
    end

    # Convert this dictionary into a latex table 
    export_stat_dict_to_latex(within_stat_dict, ty, plot_name, path, label)

    return within_stat_dict
end

function generate_correlations_table(path, data_name, obs_meas, c_dict)
    for meas in obs_meas
        table = """
        \\begin{table}[h!]
        \\centering
        \\caption{Correlations}
        \\begin{tabular}{|l|c|c|c|}
        \\hline
        """

        # Loop through models and add rows
        for model in collect(keys(c_dict[meas]))
            table *= "\\multicolumn{4}{|c|}{" * model * "} \\\\ \\hline\n"
            table *= "& Bottom & Middle & Top \\\\ \\hline\n"
            table *= "Entire series & "
            table *= "$(c_dict[meas][model]["bottom"][1]) & $(c_dict[meas][model]["middle"][1]) & $(c_dict[meas][model]["top"][1])" * " \\\\ \\hline\n"
            table *= "Specific timeframe & "
            table *= "$(c_dict[meas][model]["bottom"][2]) & $(c_dict[meas][model]["middle"][2]) & $(c_dict[meas][model]["top"][2])" * " \\\\ \\hline\n"
        end

        # End the table
        table *= """
        \\end{tabular}
        \\end{table}
        """

        # Write the LaTeX table to a .tex file
        open(path * "correlations/" * "$(meas)_correlations_table_across_models" * data_name * ".tex", "w") do io
            write(io, table)
        end
    end
end

# function HP_mat(x, λ)
#     hp_mat = similar(x)
#     println(size(hp_mat))
#     for j in axes(x, 1)
#         hp_mat[j, :] = HP(vec(x[j, :]), λ)
#     end
#     return hp_mat
# end


function compute_copula_within_stat(data_dict, confidence_intervals, base_jump, end_jump, data_name, dimension, grid_choice_cop)
    # How many of the estimates fall within the copula data intervals  

    # Data Intervals 
    cop_ci_l = confidence_intervals[data_name]["ci_l"]["copula"][:, base_jump:end-end_jump]
    cop_ci_u = confidence_intervals[data_name]["ci_u"]["copula"][:, base_jump:end-end_jump]

    # Estimates 
    T = size(data_dict["copulas"]["data"])[end]
    cop_est = reshape(data_dict["copulas"]["data"], (grid_choice_cop^dimension, T))[:, base_jump:end-end_jump]

    # Conditioning on the periods where the data is observed <=> where intervals are observed 
    cond = data_name != "consensus" ? vec(any(.!isnan.(cop_ci_l), dims=1)) : [1]
    observed_periods = collect(1:T)[cond]


    # See how many points fall within the confidence intervals, but also accounting for the NaNs from smaller dimensional copula 
    num = 0
    den = 0

    for observed_period in observed_periods
        observed_rows = vec(.!isnan.(cop_ci_l[:, observed_period]))

        num += count(cop_ci_l[observed_rows, observed_period] .<= cop_est[observed_rows, observed_period] .<= cop_ci_u[observed_rows, observed_period])
        den += length(cop_est[observed_rows, observed_period])
    end

    # Return the ratio 
    return "$num" * "/" * "$den"
end


function filter_rows(matrix)
    filtered_rows = Int[]
    for i in 1:size(matrix, 1)
        if any(x -> !isnan(x) && !iszero(x), matrix[i, :])
            push!(filtered_rows, i)
        end
    end
    return filtered_rows
end


"""
`export_measurements_to_latex(data, filename)`

Export a LaTeX table with aggregated measurements.

# Arguments
- `data`: Dictionary with measures and nested data series.
- `filename`: Name of the output .tex file.
"""
function export_stat_dict_to_latex(within_stat_dict, ty, data_name, path, label)

    # Start LaTeX table
    latex_table = L"\begin{tabular}{|c|c|}\n\hline\nMeasure & %$(data_name) \\\ \hline\n"

    # Loop over measures
    for (measure, details) in within_stat_dict
        numerator = 0
        denominator = 0

        if measure != "copula"
            for (_, series) in details
                for (_, ratio) in series
                    # Assume ratio is a string like "2/4"
                    num, den = split(ratio, "/")
                    numerator += parse(Float64, num)
                    denominator += parse(Float64, den)
                end
            end
            aggregate_value = numerator / denominator

            # Add row to table
            row_data = "$measure & $aggregate_value \\\\ \\hline"
            latex_table *= L"%$(row_data)\n"
        else
            num, den = split(details, "/") # details = ratio here 
            numerator += parse(Float64, num)
            denominator += parse(Float64, den)

            aggregate_value = numerator / denominator

            # Add row to table
            row_data = "$measure & $aggregate_value \\\\ \\hline"
            latex_table *= L"%$(row_data)\n"
        end
    end

    # Closing the table environment
    latex_table *= L"\end{tabular}"
    detrended_or_not = ty == "normal" ? "" : "_detrended"
    filename = path * data_name * "_within_stat_" * detrended_or_not * label * ".tex"

    # Write to file
    mkpath(dirname(filename))
    open(filename, "w") do file
        write(file, string(latex_table))
    end
end



function log_transformation(x)
    # result = similar(x)  # Create a new vector with the same size and type as x
    result = sign.(x) .* log.(abs.(x) .+ 1) #log.(x[j,:])

    return result
end


"""
    compute_cross_conditional_means(copulas, new_data_pcf, measures, cop_grid; share_spec)

Compute cross-conditional share-group means from the discretized copula and
already-integrated marginal levels.

Returns a `Dict{String, Matrix{Float64}}` where keys are e.g. `"consum_by_income"`
and values are `(length(share_spec) × T)` matrices of group means.
"""
function compute_cross_conditional_means(
    copulas,          # D-dim copula density after generate_copula_densities, shape (G,...,G, T)
    new_data_pcf,     # Vector of D matrices, each (G × T) — already-integrated levels
    measures,         # sorted measure names (strings)
    cop_grid::Int;    # G
    share_spec = [0.5, 0.4, 0.1]
)
    D = length(measures)
    T = size(new_data_pcf[1], 2)
    cop_shape = ntuple(_ -> cop_grid, D)

    # Map share groups to copula grid index ranges
    cum_shares = cumsum(share_spec)
    group_ranges = Vector{UnitRange{Int}}(undef, length(share_spec))
    for s in eachindex(share_spec)
        lo_idx = s == 1 ? 1 : floor(Int, cum_shares[s-1] * cop_grid) + 1
        hi_idx = floor(Int, cum_shares[s] * cop_grid)
        group_ranges[s] = lo_idx:hi_idx
    end

    result = Dict{String, Matrix{Float64}}()
    for (oi, out_meas) in enumerate(measures)
        for (ci, cond_meas) in enumerate(measures)
            oi == ci && continue
            result["$(out_meas)_by_$(cond_meas)"] = zeros(length(share_spec), T)
        end
    end

    for t in 1:T
        cop_t = reshape(copulas[ntuple(_ -> (:), D)..., t], cop_shape)

        for (oi, out_meas) in enumerate(measures)
            for (ci, cond_meas) in enumerate(measures)
                oi == ci && continue
                key = "$(out_meas)_by_$(cond_meas)"

                # Both the outcome and conditioning PCFs must be non-NaN
                if any(isnan, @view(new_data_pcf[oi][:, t])) ||
                   any(isnan, @view(new_data_pcf[ci][:, t]))
                    result[key][:, t] .= NaN
                    continue
                end

                # Marginalize copula to the (oi, ci) plane by summing over
                # all other dimensions, skipping NaN cells.  This handles
                # data sources that only observe a subset of measures: the
                # copula density is non-NaN on the observed-measure slice
                # and NaN elsewhere (see generate_copula_densities).
                for (si, grp) in enumerate(group_ranges)
                    num = 0.0; den = 0.0
                    for idx in CartesianIndices(cop_shape)
                        idx[ci] in grp || continue
                        c_val = cop_t[idx]
                        isnan(c_val) && continue
                        num += c_val * new_data_pcf[oi][idx[oi], t]
                        den += c_val
                    end
                    result[key][si, t] = den > 0 ? num / den : NaN
                end
            end
        end
    end
    return result
end


function export_functional_data(data_vector, ty, data_name, type, obs_data, func_data, time_params, user_t, model_options, posterior_bounds=false, plot=false, model_step=:inner; max_order::Union{Int,Nothing} = nothing)
    """Takes the data, creates a dataframe and exports as csv.
    'ty' is either "normal" or "average" type
    """

    @unpack measures, case, equivalized, bottom_coded, estimator, tag = model_options
    @unpack grid_cop, grid_pcf, integral_pcf_grid = estimator
    sort!(measures) # just in case

    # Resolve max_order: default to full polynomial order from estimator
    max_order = isnothing(max_order) ? grid_pcf - 1 : min(max_order, grid_pcf - 1)

    @unpack tmin, tmax = time_params  # m = model
    @unpack gdp_series = obs_data
    @unpack confidence_intervals, func_dict = func_data

    local bot, mid, top
    if integral_pcf_grid == 10 || integral_pcf_grid == 100 || integral_pcf_grid == 20
        bot, mid, top = "bottom50", "next40", "top10"
    elseif integral_pcf_grid == 5
        bot, mid, top = "bottom40", "next40", "top20"
    end

    obs_meas = data_name == "consensus" ? measures : get_obs_meas(func_dict, data_name, measures, top=top)

    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid, integral_cop_grid = estimator
    end

    grid_data_size_pcf = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf
    grid_data_size_cop = typeof(estimator) <: SeriesEstimator ? integral_cop_grid : grid_cop

    # Making label 
    m_label = measures_folder(measures)
    folder = type == :optimization ? "from_optimization" : "from_mcmc"
    equiv = equivalized == true ? "eq" : ""
    botcod = !isempty(bottom_coded) ? "bc" : ""
    label = "_$case" * "_$equiv" * "$botcod"

    # Objects 
    copulas = copy(data_vector[1])  # It is necessary for it to be a copy!
    pcfs = copy(data_vector[2])  # gridp * d, T

    # Dimensions 
    D = length(measures)
    T = size(copulas)[end]

    # Reshape N-dimensional copula to N x T matrix
    prodt = grid_data_size_cop^D
    cop_size = tuple(vcat([grid_data_size_cop for i in 1:D], [T])...)

    CI = CartesianIndices(cop_size)
    CI_r = reshape(CI, (prodt, T))
    mat = zeros(prodt, T)

    # Ideally, I want to transpose CI_r, but I cant. So ...
    for (i, ci) in enumerate(CI_r)
        mat[i] = parse(Int64, join(ci.I[1:end-1]))  # convert string to Int, end -1 is there to remove time index 
    end

    # Getting column names for percentile functions 
    stubs = ["quantiles", "levels", "shares"]

    grid_points_pcf = select_grid_points(grid_data_size_pcf)

    # combining data 
    data_pcf = reshape(pcfs, (grid_pcf * D, T))

    # Transform pcf data to original form 
    gdp_series[!, "date"] = QuarterlyDate.(gdp_series[!, "time"])
    correction_names = [meas * "_per_hh" for meas in measures]
    filter!(row -> row.date >= QuarterlyDate(tmin["year"], tmin["quarter"]), gdp_series)
    filter!(row -> row.date <= QuarterlyDate(tmax["year"], tmax["quarter"]), gdp_series)

    select_series = select(gdp_series, correction_names)
    split_pcfs = [data_pcf[I, :] for I in Iterators.partition(axes(data_pcf, 1), grid_pcf)]  # split by measure 

    local data_cop
    if typeof(estimator) <: HistogramEstimator
        if model_step != :forecast
            for t in 1:T
                data_pcf[:, t] .= vcat([split_pcfs[m][:, t] .* select_series[t, correction_names[m]] for m in eachindex(split_pcfs)]...)
            end
        end
    elseif typeof(estimator) <: SeriesEstimator
        # Generate container to store the data of choice 
        new_data_pcf = [zeros(integral_pcf_grid, T) for _ in 1:D]
        intervals = vcat([0.0] .+ 1e-6, grid_points_pcf)

        # We need to generate new data_pcf, which are the average quantiles over the intervals (e.g., deciles)
        for m in eachindex(new_data_pcf)
            for t in 1:T
                if all(isnan.(split_pcfs[m][:, t]))
                    new_data_pcf[m][:, t] .= NaN
                else
                    for i in 1:integral_pcf_grid
                        # Using coefs, generate pcf function and then integrate pcf function over diff. intervals
                        integral, _ = quadgk(u -> reverse_inverse_hyperbolic_sine(eval_quantile_function(split_pcfs[m][:, t], max_order, u))[1] .* select_series[t, correction_names[m]], intervals[i], intervals[i+1], rtol=1e-8)

                        # Undo treatment of data => gives us average quantile within the interval
                        new_data_pcf[m][i, t] = integral / (intervals[i+1] - intervals[i]) #reverse_inverse_hyperbolic_sine(integral)[1] .* select_series[t, correction_names[m]] #./ (intervals[i+1] - intervals[i])
                    end
                end
            end
        end

        # Share-group integration: bot50/mid40/top10
        share_spec = [0.5, 0.4, 0.1]
        share_intervals = vcat([0.0 + 1e-6], cumsum(share_spec)[1:end-1], [1.0 - 1e-6])
        share_data_pcf = [zeros(length(share_spec), T) for _ in 1:D]
        corr_mat = Matrix(select_series[:, correction_names])
        integrate_quantile_functions!(share_data_pcf, split_pcfs, grid_pcf, share_intervals, corr_mat; max_order = max_order)

        data_pcf = vcat([new_data_pcf[m] for m in eachindex(new_data_pcf)]...)

        # Save copula coefficients before overwriting with densities
        copulas_coefs = copy(copulas)

        copulas = generate_copula_densities(copulas, measures, integral_cop_grid)

        # Cross-conditional share means on a finer grid (deciles → 10×10×10 = 1000 cells)
        xcond_grid = 10
        copulas_xcond = generate_copula_densities(copy(copulas_coefs), measures, xcond_grid)
        xcond_intervals = vcat([0.0 + 1e-6], select_grid_points(xcond_grid))
        xcond_pcf = [zeros(xcond_grid, T) for _ in 1:D]
        integrate_quantile_functions!(xcond_pcf, split_pcfs, grid_pcf, xcond_intervals, corr_mat; max_order = max_order)
        cross_cond_data = compute_cross_conditional_means(
            copulas_xcond, xcond_pcf, sort(measures), xcond_grid
        )

        # 'data_cop' is (x, x, x, T). Reshape to (x^d, T)
        data_cop = reshape(copulas, (integral_cop_grid^D, T))
    end


    # How to get the names of all variables 
    column_names = String[]
    abv = join([measures[i][1] for i in eachindex(measures)]) # take first letters of each measure 
    for ci in mat'[1, :]
        push!(column_names, abv * "_" * "$ci")
    end

    for stub in stubs
        for meas in measures
            for group in grid_points_pcf
                push!(column_names, stub * meas * "_$group")
            end
        end
    end

    # Generate levels and shares from the percentile functions
    levels, shares = generate_shares_levels(data_pcf, model_options, gdp_series)
    data = vcat(data_cop, data_pcf, levels, shares)
    E = DataFrame(Matrix(data'), column_names)
    E[!, "time"] = QuarterlyDate(tmin["year"], tmin["quarter"]):Quarter(1):QuarterlyDate(tmax["year"], tmax["quarter"])

    # Generate microdata 
    # Export Kendall's Tau as well 
    # if plot == true
    #     kendalls_tau(data_cop, data_pcf, data_name, folder, time_params, model_options, false, true, func_data, gdp_series)
    #     compute_tail_dependence(data_cop, data_pcf, data_name, folder, time_params, model_options)
    # end

    # Export data 
    data_tag = ty == "normal" ? "" : "_detrended"
    if model_step != :forecast
        # Export micro-data
        create_micro_df(copulas, data_pcf, data_tag, data_name, folder, time_params, model_options) #TODO: check that the order is correct for the pcfs 

        init_path = BASE_PATH
        path = init_path * "/7_Results/$m_label" * "$tag" * "/$folder/data"
        mkpath(path)

        # Export functional data. For 'detrended', the marginals itself are uninformative. It is later divided by the average, which is when its informative.
        CSV.write(path * "/" * data_name * "_functional_data" * data_tag * "$label" * ".csv", select(E, "time", :))
    end

    # Place data in a dictionary for time series plotting
    data_dict = create_time_series_dictionary([copulas, data_pcf, levels, shares], estimator, measures)

    # Inject share-group means (bot50/mid40/top10) into the dictionary
    if typeof(estimator) <: SeriesEstimator && @isdefined(share_data_pcf)
        share_labels = ["bot50", "mid40", "top10"]
        for (i, meas) in enumerate(sort(measures))
            for (si, sl) in enumerate(share_labels)
                data_dict[meas]["quantiles"]["common series"][sl] = share_data_pcf[i][si, :]
            end
        end
    end

    # Inject cross-conditional share means into data_dict
    if typeof(estimator) <: SeriesEstimator && @isdefined(cross_cond_data)
        cross_share_labels = ["bot50", "mid40", "top10"]
        for (key, mat) in cross_cond_data
            data_dict[key] = Dict{String, Any}()
            for (si, sl) in enumerate(cross_share_labels)
                data_dict[key][sl] = mat[si, :]
            end
        end
    end

    # Everything has been reconstructed now. Creating a new select_series object, which has the averages that correspond to the data estimates 
    # To do this, generate total levels and divide the number of households  
    local avg_series
    if occursin("HANK", data_name)
        @unpack gdp_series = obs_data
        avg_series = deepcopy(gdp_series)
    else
        println(data_name)
        avg_series = generate_average_series(data_dict, gdp_series, measures)
    end

    if plot == true
        within_stat_dict = generate_specific_plots(data_dict, ty, func_data, data_name, time_params, user_t, model_options, type, avg_series, gdp_series)
        return within_stat_dict, data_dict

    elseif model_step == :forecast
        return data_dict, avg_series

    else
        return data_dict
    end
end


function generate_average_series(data_dict, gdp_series, measures)
    avg_series = DataFrame()
    tot_hhs = gdp_series[!, "tot_hhs"]

    for meas in measures
        # Generate the level 
        try
            avg_series[!, meas*"_per_hh"] = vec(sum(data_dict[meas]["levels"]["data"], dims=1)' ./ tot_hhs)
        catch e
            @warn "Could not compute average series for $meas. Error: $e"
        end
    end

    # select!(avg_series, Not(:date))
    return avg_series
end


"""
`export_stat_dict_to_latex(stat_dict, path, label)`

Export a LaTeX table with aggregated statistics for multiple datasets.

# Arguments
- `stat_dict`: Dictionary with dataset names as keys, and further dictionaries of measures as values.
- `path`: Path to save the .tex file.
- `label`: Label to append to the filename.
"""

# Helper: return the dict that contains the per-measure stats (usually the "normal" one)
function _within_stats(stat_dict, dataset; regime::String="normal")
    d = stat_dict[dataset]

    # Case 1: stat_dict[dataset]["normal"]
    if isa(d, AbstractDict) && haskey(d, regime)
        return d[regime]
    end

    # Case 2: stat_dict[dataset][something]["normal"] (e.g. ["a"]["normal"])
    if isa(d, AbstractDict)
        # fast-path for your known key
        if haskey(d, "a") && isa(d["a"], AbstractDict) && haskey(d["a"], regime)
            return d["a"][regime]
        end

        # generic fallback: find the first sub-dict that has the regime key
        for (_, sub) in d
            if isa(sub, AbstractDict) && haskey(sub, regime)
                return sub[regime]
            end
        end
    end

    error("Could not find \"$regime\" stats for dataset=$(dataset). Top-level keys: $(collect(keys(d)))")
end

function export_combined_stat_dict_to_latex(stat_dict, measures, path, label)
    datasets = collect(keys(stat_dict))
    num_datasets = length(datasets)

    # Optional: escape underscores for LaTeX headers
    dataset_names = [replace(string(d), "_" => "\\_") for d in datasets]

    latex_table = L"\begin{tabular}{|" * "c"^(num_datasets + 2) * "}\n\\hline\n"
    header_row = "Measure & " * join(dataset_names, " & ") * " & Overall & " * L"\\ \hline\n"
    latex_table *= L"%$(header_row)"

    for measure in measures
        row_data = string(measure)

        overall_num = 0.0
        overall_den = 0.0

        for dataset in datasets
            within_stat_dict = _within_stats(stat_dict, dataset; regime="normal")

            if haskey(within_stat_dict, measure)
                numerator, denominator = 0.0, 0.0

                if measure != "copula"
                    for (_, series) in within_stat_dict[measure]
                        for (_, ratio) in series
                            num, den = split(ratio, "/")
                            numerator += parse(Float64, num)
                            denominator += parse(Float64, den)
                        end
                    end
                    aggregate_value = round(numerator / denominator, digits=2) * 100
                else
                    num, den = split(within_stat_dict["copula"], "/")
                    numerator += parse(Float64, num)
                    denominator += parse(Float64, den)
                    aggregate_value = denominator != 0 ? ceil(numerator * 100 / denominator) : "-"
                end

                overall_num += numerator
                overall_den += denominator
            else
                aggregate_value = "-"
            end

            row_data *= " & $aggregate_value %"
        end

        overall_ratio = overall_den != 0 ? Int(ceil(overall_num * 100 / overall_den)) : "-"
        row_data *= " & $overall_ratio %"
        row_data *= L" \\\\ \\hline"

        latex_table *= L"%$(row_data)\n"
    end

    latex_table *= L"\end{tabular}"
    filename = path * "combined_stat_" * label * ".tex"

    open(filename, "w") do file
        write(file, string(latex_table))
    end
end

# function export_combined_stat_dict_to_latex(stat_dict, measures, path, label)
#     # Determine the number of datasets
#     datasets = keys(stat_dict)
#     num_datasets = length(datasets)

#     # Start LaTeX table
#     latex_table = L"\begin{tabular}{|" * "c"^(num_datasets + 2) * "}\n\\hline\n"
#     header_row = "Measure & " * join(datasets, " & ") * " & Overall & " * L"\\ \hline\n"
#     latex_table *= L"%$(header_row)"

#     # Loop over measures
#     for measure in measures
#         row_data = measure

#         # For each measure + copula, generate overall within-stat
#         overall_num = 0
#         overall_den = 0

#         # Loop over datasets
#         for dataset in datasets
#             within_stat_dict = stat_dict[dataset]["normal"]

#             # Aggregate value calculation
#             if haskey(within_stat_dict, measure)
#                 numerator, denominator = 0, 0
#                 if measure != "copula"
#                     for (_, series) in within_stat_dict[measure]
#                         for (_, ratio) in series
#                             num, den = split(ratio, "/")
#                             numerator += parse(Float64, num)
#                             denominator += parse(Float64, den)
#                         end
#                     end
#                     aggregate_value = round(numerator / denominator, digits=2) * 100
#                     overall_num += numerator
#                     overall_den += denominator
#                 else
#                     num, den = split(within_stat_dict["copula"], "/")
#                     numerator += parse(Float64, num)
#                     denominator += parse(Float64, den)

#                     aggregate_value = denominator != 0 ? ceil(numerator * 100 / denominator) : "-"
#                     overall_num += numerator
#                     overall_den += denominator
#                 end
#             else
#                 aggregate_value = "-"  # Measure not present in this dataset
#                 overall_num += 0
#                 overall_den += 0
#             end

#             row_data *= " & $aggregate_value %"
#         end

#         local overall_ratio
#         try
#             overall_ratio = Int(ceil(overall_num * 100 / overall_den))
#         catch e
#             overall_ratio = "-"
#         end

#         row_data *= " & $overall_ratio %"

#         row_data *= L" \\\\ \\hline"
#         latex_table *= L"%$(row_data)\n"
#     end

#     # Closing the table environment
#     latex_table *= L"\end{tabular}"
#     filename = path * "combined_stat_" * label * ".tex"

#     # Write to file
#     open(filename, "w") do file
#         write(file, string(latex_table))
#     end
# end


function generate_shares_levels(data_pcf, model_options, gdp_series)
    """
    With the means of the deciles representative of the distribution, we can compute the shares and levels of the distribution.
    """
    # Each decile represents the mean of 10 percentiles so we need to multiply it by 10 to get the levels of each decile 
    @unpack measures, case, estimator, equivalized, bottom_coded, tag = model_options
    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid = estimator
    else
        @unpack grid_pcf, grid_cop = estimator
    end
    grid_choice = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf

    D = length(measures)
    T = size(data_pcf, 2)
    shares = zeros(grid_choice * D, T)
    lvls = zeros(grid_choice * D, T)

    # Housing 
    local tot_hhs
    try
        tot_hhs = gdp_series[!, "tot_hhs"]
    catch e
        tot_hhs = ones(T) .* 5000 # if not available, assume 1 household
    end

    # The grid granularity defines the percent we need    
    n_grups = tot_hhs ./ grid_choice

    for I in Iterators.partition(1:size(data_pcf, 1), grid_choice)
        for t in 1:T
            lvls[I, t] = data_pcf[I, t] .* n_grups[t]
            t_lvl = sum(lvls[I, t])
            shares[I, t] .= lvls[I, t] ./ t_lvl
        end
    end

    return lvls, shares
end



#TODO: do the rest of the labels 
function label_qs(grid)
    local plot_labels
    if grid == 10
        plot_labels = reshape([L"\textrm{1st\,\,Decile}", L"\textrm{2nd\,\,Decile}", L"\textrm{3rd\,\,Decile}", [L"\textrm{%$(i)th\,\,Decile}" for i in 4:1:10]...], (1, length(10:10:100)))
    elseif grid == 20
        plot_labels = reshape([L"\textrm{%$(i)th}" for i in 5:5:100], (1, length(5:5:100)))
    elseif grid == 5
        plot_labels = reshape([L"\textrm{%$(i)th}" for i in 20:20:100], (1, length(20:20:100)))
    elseif grid == 100
        plot_labels = reshape([L"\textrm{%$(i)th}" for i in 1:1:100], (1, length(1:1:100)))
    end
    return plot_labels
end

function generate_agg_labels(grid)
    local agg_labels
    if grid == 10 || grid == 20 || grid == 100
        agg_labels = [L"\textrm{Top\, 10}" L"\textrm{Next\, 40}" L"\textrm{Bottom\, 50}"]
    elseif grid == 5
        agg_labels = [L"\textrm{Top\, 20}" L"\textrm{Next\, 40}" L"\textrm{Bottom\, 40}"]
    end
    return agg_labels
end


function define_sequences(grid)
    # Quantiles 
    local sequences
    if grid == 10
        sequences = [1:5, 6:9, 10:10]
    elseif grid == 5
        sequences = [1:2, 3:4, 5:5]
    elseif grid == 20
        sequences = [1:10, 11:18, 19:20]
    elseif grid == 100
        sequences = [1:50, 51:90, 91:100]
    end
    return sequences
end


function define_sequences_q(grid)
    # Quantiles 
    local sequences
    if grid == 10
        sequences = [1:5, 6:9, 10]
    elseif grid == 5
        sequences = [1:2, 3:4, 5:5]
    elseif grid == 20
        sequences = [1:10, 11:18, 19:20]
    elseif grid == 100
        sequences = [1:50, 51:90, 91:100]
    end
    return sequences
end


function generate_microdata_implicates(draws, k, param_sizes, priors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources, tag)
    micro_full_df = DataFrame()
    m_label = measures_folder(measures)
    init_path = BASE_PATH
    @unpack estimator = model_options
    @unpack integral_pcf_grid = estimator
    # Generate draws 
    postd = jldopen(init_path * "/posterior_draws" * "/" * m_label * "_$tag.jld2", "r")
    itr, ch, par = size(postd["d_chains"])
    param_mat = zeros(par, draws)

    all_draws = postd["d_chains"][end, :, :][:, :] # take last draw from each chain

    # Generate 50 draws, with sufficient distance between them
    draw_ids = sample(1:ch, draws, replace=false)

    for i in 1:draws
        param_mat[:, i] = all_draws[draw_ids[i], :]
    end

    # Matrix to store the data
    # - we want to store 3 series for each measure, for each draw
    q_dict = Dict()
    T = time_params.tot_periods
    for meas in measures
        q_dict[meas] = zeros(3, T, draws)
    end

    # For each parameter vector given, generate the microdata implicates
    pb = Progress(draws, desc="Creating micro data implicates for $k")
    for p in 1:draws
        dv, _ = reconstruct_data(param_mat[:, p], param_sizes, priors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources)
        imp_df = get_implicate!(dv[k]["normal"], q_dict, p, k, kind_of_plots, obs_data, func_data, time_params, user_t, model_options)

        # Add column to imp_df that indicates the parameter vector
        imp_df[!, :implicate] .= p
        append!(micro_full_df, imp_df)

        next!(pb)
    end

    # Keep the 95% quantiles of the implicates from q_dict 
    lower_bound = Dict()
    upper_bound = Dict()

    for meas in measures
        lower_bound[meas], upper_bound[meas] = conservative_bounds(q_dict[meas], 3)
    end

    save_name = k == "CEX_pooled" ? "CEX" : k

    # Export the bounds
    mkpath(init_path * "/7_Results/$m_label" * "$tag" * "/from_mcmc/data")
    jldsave(init_path * "/7_Results/$m_label" * "$tag" * "/from_mcmc/data/" * save_name * "_bounds" * ".jld2"; lb=lower_bound, ub=upper_bound, q_dict=q_dict)

    # Export the data
    path = init_path * "/7_Results/$m_label" * "$tag" * "/from_mcmc/data"
    CSV.write(path * "/" * save_name * "_micro_data_w_imp" * ".csv", micro_full_df)
end


function get_implicate!(data_vector, q_dict, draw, data_name, type, obs_data, func_data, time_params, user_t, model_options)
    """Takes the data, creates a dataframe and exports as csv."""
    # Wrestling with time
    @unpack tmin, tmax = time_params  # m = model

    base_jump, end_jump = find_subset_frame(tmin, tmax, tmin, tmax)

    @unpack measures, case, equivalized, bottom_coded, estimator, tag = model_options
    @unpack grid_cop, grid_pcf = estimator
    sort!(measures) # just in case 

    @unpack gdp_series = obs_data
    @unpack confidence_intervals, func_dict = func_data

    obs_meas = data_name == "consensus" ? measures : get_obs_meas(func_dict, data_name, measures)

    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid, integral_cop_grid = estimator
    end
    grid_data_size_pcf = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf
    grid_data_size_cop = typeof(estimator) <: SeriesEstimator ? integral_cop_grid : grid_cop

    # Making label 
    m_label = measures_folder(measures)
    folder = type == :optimization ? "from_optimization" : "from_mcmc"
    equiv = equivalized == true ? "eq" : ""
    botcod = !isempty(bottom_coded) ? "bc" : ""
    label = "_$case" * "_$equiv" * "$botcod"

    # Objects 
    copulas = copy(data_vector[1])  # It is necessary for it to be a copy!
    pcfs = copy(data_vector[2])  # gridp * d, T

    # Dimensions 
    D = length(measures)
    T = size(copulas)[end]

    # Reshape N-dimensional copula to N x T matrix
    prodt = grid_data_size_cop^D
    cop_size = tuple(vcat([grid_data_size_cop for i in 1:D], [T])...)

    CI = CartesianIndices(cop_size)
    CI_r = reshape(CI, (prodt, T))
    mat = zeros(prodt, T)

    # Ideally, I want to transpose CI_r, but I cant. So ...
    for (i, ci) in enumerate(CI_r)
        mat[i] = parse(Int64, join(ci.I[1:end-1]))  # convert string to Int, end -1 is there to remove time index 
    end

    grid_points_pcf = select_grid_points(grid_data_size_pcf)
    grid_points_cop = select_grid_points(grid_data_size_cop)

    # combining data 
    data_pcf = reshape(pcfs, (grid_pcf * D, T))

    # Transform pcf data to original form 
    gdp_series[!, "date"] = QuarterlyDate.(gdp_series[!, "time"])
    correction_names = [meas * "_per_hh" for meas in measures]
    filter!(row -> row.date >= QuarterlyDate(tmin["year"], tmin["quarter"]), gdp_series)
    filter!(row -> row.date <= QuarterlyDate(tmax["year"], tmax["quarter"]), gdp_series)

    select_series = select(gdp_series, correction_names)
    split_pcfs = [data_pcf[I, :] for I in Iterators.partition(axes(data_pcf, 1), grid_pcf)]  # split by measure 

    local data_cop
    if typeof(estimator) <: HistogramEstimator
        if model_step != :forecast
            for t in 1:T
                data_pcf[:, t] .= vcat([split_pcfs[m][:, t] .* select_series[t, correction_names[m]] for m in eachindex(split_pcfs)]...)
            end
        end
    elseif typeof(estimator) <: SeriesEstimator
        # Generate container to store the data of choice 
        new_data_pcf = [zeros(integral_pcf_grid, T) for _ in 1:D]
        intervals = vcat([0.0] .+ 1e-6, grid_points_pcf)

        # We need to generate new data_pcf, which are the average quantiles over the intervals (e.g., deciles)
        for m in eachindex(new_data_pcf)
            for t in 1:T
                if all(isnan.(split_pcfs[m][:, t]))
                    new_data_pcf[m][:, t] .= NaN
                else
                    for i in 1:integral_pcf_grid
                        # Using coefs, generate pcf function and then integrate pcf function over diff. intervals 
                        # try
                        integral, _ = quadgk(u -> reverse_inverse_hyperbolic_sine(eval_quantile_function(split_pcfs[m][:, t], grid_pcf - 1, u))[1] .* select_series[t, correction_names[m]], intervals[i], intervals[i+1], rtol=1e-8)
                        # catch E
                        #     println(split_pcfs[m][:, t])
                        #     println(select_series[t, correction_names[m]])
                        #     turtles
                        # end

                        # Undo treatment of data => gives us average quantile within the interval 
                        new_data_pcf[m][i, t] = integral / (intervals[i+1] - intervals[i]) #reverse_inverse_hyperbolic_sine(integral)[1] .* select_series[t, correction_names[m]] #./ (intervals[i+1] - intervals[i])
                    end
                end
            end
        end

        data_pcf = vcat([new_data_pcf[m] for m in eachindex(new_data_pcf)]...)
        copulas = generate_copula_densities(copulas, measures, integral_cop_grid)
    end

    # Generate microdata 
    imp_df = create_micro_df(copulas, data_pcf, data_name, "", folder, time_params, model_options; return_df=true) #TODO: check that the order is correct for the pcfs 

    # Get data for uncertainty around the estimates
    levels, shares = generate_shares_levels(data_pcf, model_options, gdp_series)
    indices = [I for I in Iterators.partition(axes(data_pcf, 1), grid_data_size_pcf)]
    sequences = define_sequences(grid_data_size_pcf)

    avg_series = DataFrame()
    tot_hhs = gdp_series[!, "tot_hhs"]

    for (i, meas) in enumerate(measures)
        # Generate the level 
        avg_series[!, meas*"_per_hh"] = vec(sum(levels[indices[i], :], dims=1)' ./ tot_hhs)
        q_mat = data_pcf[indices[i], base_jump:end-end_jump] ./ avg_series[base_jump:end-end_jump, meas*"_per_hh"]'

        # Aggregate 'q_mat' to 3 series
        for (s, g) in enumerate(sequences)
            q_dict[meas][s, :, draw] .= vec(sum(q_mat[g, :], dims=1)' ./ length(g))
        end
    end

    return imp_df
end


function conservative_bounds(data, dim)
    # Sort along the specified dimension
    sorted_data = sort(data, dims=dim)

    # Determine the size of the draw dimension
    draws = size(data, dim)

    # Calculate indices for 2.5% and 97.5% bounds
    lower_idx = max(1, floor(Int, draws * 0.025))  # Ensure at least 1
    upper_idx = min(draws, ceil(Int, draws * 0.975))  # Ensure at most `draws`

    # Extract bounds
    lower_bound = selectdim(sorted_data, dim, lower_idx)
    upper_bound = selectdim(sorted_data, dim, upper_idx)

    return lower_bound, upper_bound
end


function export_table_to_tex_with_strings(measures, type)
    m_label = measures_folder(measures)
    folder = type == :optimization ? "from_optimization" : "from_mcmc"
    init_path = BASE_PATH
    path = init_path * "/7_Results/$m_label" * "$tag" * "/$folder/plots/"
    mkpath(path)

    # Get first letter and capitalize from 'measures'
    abv = join([measures[i][1] for i in eachindex(measures)])
    ABV = uppercase(abv)


    # Headers 
    headers = [
        "excluding housing cycle",
        "excluding recent 20 quarters",
        "excluding housing cycle short",
        "every 4 years",
        # "PP CEX excluding housing cycle",
        # "PP CEX excluding housing cycle short",
        # "PP CEX excluding recent 20 quarters",
        # "PP CEX every 4 years"
    ]

    # Import tables 
    for data_name in ["SCF", "CEX"]
        cycle_dict = jldopen(path * "correlations/" * "correlations_dict_cycle_" * data_name * ".jld2")["corr"]

        # Initialize the LaTeX table string
        table = """
        \\begin{table}[h!]
        \\centering
        \\begin{tabular}{l c c c c c c c c c}
        \\toprule\\toprule
        """

        for header in headers
            first_m = collect(keys(cycle_dict))[1]
            println(collect(keys(cycle_dict[first_m])))
            if header ∉ collect(keys(cycle_dict[first_m]))
                nothing
            else
                header_label = header == "excluding housing cycle" ? "Excluding Housing Cycle" : header == "excluding recent 20 quarters" ? "Excluding Last 4 Years" : header == "excluding housing cycle short" ? "Excluding Housing Cycle \\#2" : header == "PP CEX excluding housing cycle" ? "CEX Factor excluding housing cycle" : header == "PP CEX excluding housing cycle short" ? "CEX Factor excluding housing cycle short" : header == "PP CEX excluding recent 20 quarters" ? "CEX Factor excluding recent 20 quarters" : header == "PP CEX every 4 years" ? "CEX Factor every 4 years" : header == "excluding housing cycle wealth" ? "Excluding Housing Cycle Wealth" : "Every 4 Years"

                # Add header section
                table *= """
                &\\multicolumn{9}{c}{\\textbf{$header_label}} \\vspace{2mm}\\\\
                \\multirow{2}{*}{\\textbf{Condition}} & \\multicolumn{3}{c}{\\textbf{Bottom}} & \\multicolumn{3}{c}{\\textbf{Middle}} & \\multicolumn{3}{c}{\\textbf{Top}} \\\\
                """

                obs_meas = collect(keys(cycle_dict))
                for _ in eachindex(measures)
                    for meas in measures
                        cap_meas = uppercase(meas[1])
                        table *= "& \\textbf{$cap_meas}"
                    end
                end

                table *= " \\\\\n"

                # Add condition rows
                for (zz, condition) in enumerate(["Entire Series", "Specific Timeframe"])
                    table *= condition
                    for k in ["bottom", "middle", "top"]
                        for m in measures
                            try
                                ρ_values = cycle_dict[m][header][k]
                                table *= " & " * "$(ρ_values[zz])"
                            catch ee
                                println(ee)
                                table *= " & " * "-"
                            end
                        end
                    end
                    table *= " \\\\\n"
                end
                # Add spacing between headers
                table *= "\\vspace{2mm}\\\\\n"
            end
        end

        # End the table
        table *= """
        \\bottomrule\\bottomrule
        \\end{tabular}
        \\caption{Your Caption Here}
        \\label{tab:your_label}
        \\end{table}
        """

        # Write the LaTeX table to a .tex file
        filename = path * "correlations/" * "correlations_table_" * data_name * ".tex"
        mkpath(dirname(filename))
        open(filename, "w") do io
            write(io, table)
        end
    end
end


function generate_correlations_table_for_external_comparisons(data_choice, measures, tag, type, series_type)
    init_path = BASE_PATH
    meas_folder = measures_folder(measures)
    corr_path = init_path * "/7_Results/$(meas_folder)$(tag)/$type/plots/correlations"
    ext_corr = jldopen(corr_path * "/external_correlations.jld2", "r")["ext_corr"]
    recon_corr = CSV.read(corr_path * "/correlations.csv", DataFrame)
    # ext_corr["ext_corr"]["quantiles"]["wealth"][series_type]["top"]["WID-DFA"]
    # SCF_WID_quantiles_cycle_bottom50_income


    # Generate a table with two panels: one for income and one for wealth ... in each panel, we should have 3 correlation matries, one for bottom, middle, and top quantiles
    # For now we focus on cycles 
    # Initialize the LaTeX table string
    table = """
    \\begin{table}[h!]
    \\centering
    \\begin{tabular}{l c c c}
    \\toprule\\toprule
    """

    table *= """
    & \\multicolumn{1}{c}{Baseline-WID} & \\multicolumn{1}{c}{Baseline-DFA} & \\multicolumn{1}{c}{WID-DFA}\\\\
    \\cmidrule(lr){2-4}
    """

    seg_dict = Dict("Bottom" => "bottom50", "Middle" => "next40", "Top" => "top10")
    seg_othdict = Dict("Bottom" => "bot", "Middle" => "mid", "Top" => "top")

    for (q, quantile) in seg_dict
        table *= """
          & \\multicolumn{3}{c}{\\textbf{$(q)}} \\\\
        """
        for meas in ["income", "wealth"]
            table *= "$(uppercasefirst(meas))"
            # Access the correlation values for this measure and quantile
            corr_baseline_wid = recon_corr[!, data_choice*"_WID_quantiles_"*series_type*"_"*quantile*"_"*meas][1]

            local corr_baseline_dfa
            try
                corr_baseline_dfa = recon_corr[!, data_choice*"_DFA_quantiles_"*series_type*"_"*quantile*"_"*meas][1]
            catch ee
                corr_baseline_dfa = "-"
            end

            local corr_wid_dfa
            try
                corr_wid_dfa = ext_corr["quantiles"][meas][series_type][lowercasefirst(seg_othdict[q])]["WID-DFA"]
            catch ee
                println(ee)
                corr_wid_dfa = "-"
            end

            # Create part of the first row             
            table *= """
            & $(corr_baseline_wid) & $(corr_baseline_dfa) & $(corr_wid_dfa) \\\\
            """
        end
    end

    # End the table
    table *= """
    \\bottomrule\\bottomrule
    \\end{tabular}
    \\caption{Your Caption Here}
    \\label{tab:your_label}
    \\end{table}
    """

    # Write the LaTeX table to a .tex file
    filename = corr_path * "/external_correlations_" * data_choice * ".tex"
    open(filename, "w") do io
        write(io, table)
    end
end