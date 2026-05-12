function run_counterfactual_distributions(blind_to_dict, params, obs_data, model_options)
    # The idea is to run the same model with the original parameters, but different data 
    # This is done by removing some observations from the data and running the model 
    # The projection matrix will still be the old one, but the data will be different
    # This is done to see how the model, given posterior params, reacts to missing data
    
    func_data, time_params, model_elements = data_prep(obs_data, model_options);
    _, param_sizes, priors, meas_ind, Σ_ids       = set_params(model_elements, time_params, model_options)
    
    @unpack data_sources                          = func_data
    @unpack gdp_series, agg_data = obs_data
    @unpack estimation_object, measures, grid, reconstruction_to_show     = model_options

    # First, construct estimates with old settings, extract shocks   
    dv, _ = reconstruct_data(params, param_sizes, priors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources)

    # Using new settings = remove some data 
    c_dv       = Dict()

    # Collect MVs, one per scenario 
    @unpack MV = model_elements
    set_of_MVs = filter_measurement_vector(MV, blind_to_dict, data_sources, model_options)

    for scenario in collect(keys(blind_to_dict))
        model_elements.MV = set_of_MVs[scenario]
        c_dv[scenario], _    = reconstruct_data(params, param_sizes, priors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources, reconstruction_to_show=reconstruction_to_show)
    end

    # Plot data against one another 
    generate_missing_info_plots(dv, c_dv, gdp_series, time_params, model_options)
end

function filter_measurement_vector(MV, blind_to, data_sources, model_options)
    # Copula first, then percentiles, by dataset 
    @unpack measures, grid = model_options
    set_of_MVs = Dict()
    dims       = length(measures)
    immutable  = grid + (dims - 1) * (grid - 1)
    cop_size   = grid^(dims) - immutable 

    # For each scenario ... 
    for k in collect(keys(blind_to))
        copy_MV    = deepcopy(MV) # amazing
        # Get indices of data sources we are blind to 
        d_id       = [findall(occursin.(data_sources, source)) for source in collect(keys(blind_to[k]))]
        m_ids      = Vector{Any}(undef, length(d_id))
        for (j,d) in enumerate(d_id)
            m_ids[j]                                       = findall(occursin.(measures, blind_to[k][data_sources[j]]))
            cop_slice                                      = retrieve_cop_rows(copy_MV[d...][1:cop_size, :], grid, dims, m_ids[j])
            pcf_rows                                       = cat([(cop_size) + 1 + (x-1) .* grid:(cop_size) + x .* grid for x in m_ids[j]]..., dims=1)
            copy_MV[d...][vcat([cop_slice, pcf_rows]...), :] .= NaN
        end
        set_of_MVs[k] = copy(copy_MV)
    end
    return set_of_MVs
end

function retrieve_cop_rows(copula, grid, dims, m_ids)
    local cop_slice
    if dims == 2
        if length(m_ids) > 0
            cop_slice = 1:size(copula, 1)
        else
            cop_slice = []
        end
    # Copula on a decile grid would have 1000 rows minus 28 = 972 rows
    elseif dims == 3
        # Preliminaries, Defining immutable portion 
        f(c)         = sum((==(1)).(c.I)) >= length(c.I) - 1 
        cop_ids      = CartesianIndices((grid, grid)) # for identifying immutable of first slice in 3D case!
        imut_part    = filter(f, cop_ids)

        # We want everything outside the first slice. What is the first slice? It is observed minus immutable.
        # How can we find the immutable of the first slice? That is 

        # If there are two or more missings, then slice is everything
        if length(m_ids) >= 2
            cop_slice = 1:size(copula, 1)
        # If there is one missing, we just take the first slice as observed, everything else missing
        elseif length(m_ids) == 1
            cop_slice      = grid^(dims-1) - length(imut_part) + 1:size(copula, 1)
        elseif length(m_ids) == 0
            cop_slice = []
        end
    end
    return cop_slice
end


function run_counterfactual_aggregates(counterfactuals_dict, model_options, params, func_data, time_params, model_elements, param_sizes, priors, meas_ind, Σ_ids) 
    @unpack data_sources         = func_data
    @unpack gdp_series, agg_data = obs_data 
    @unpack tmin, tmax           = time_params
    user_t     = (deepcopy(tmin), deepcopy(tmax))

    # Using old settings 
    @info("Running counterfactuals with old settings")
    dv, dε_smoothed      = reconstruct_data(params, param_sizes, priors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources)

    # Create dictionary for each option 
    counterfactuals_data_dict = Dict()

    # Fill matrix with density, percentile function measurements 
    avg_series_cf = Dict()
     for k in collect(keys(counterfactuals_dict))
        avg_series_cf[k] = Dict() 
        @info("Generating new aggregate factors from muting $k")
        u       = generate_counterfactual_agg_factors(agg_data, counterfactuals_dict[k], time_params, model_elements, model_options);    
        @pack! model_elements = u;
        @info("Reconstructing with new aggregate factors")
        counterfactuals_data_dict[k], _    = reconstruct_data(params, param_sizes, priors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources; dε_smoothed=dε_smoothed)
        
        for (c, data_name) in enumerate(keys(counterfactuals_data_dict[k]))
            counterfactuals_data_dict[k][data_name], avg_series_cf[k][data_name]               = export_functional_data(counterfactuals_data_dict[k][data_name], data_name, :mcmc, obs_data, func_data, time_params, user_t, model_options,false,false,:forecast)
        end
    end

    avg_series = Dict()
    for (c, data_name) in enumerate(keys(dv))
        dv[data_name], avg_series[data_name]               = export_functional_data(dv[data_name], data_name, :mcmc, obs_data, func_data, time_params, user_t, model_options,false,false,:forecast)
    end

    # Plot data against one another 
    @info("Plotting counterfactuals.")
    @unpack measures, estimator, case, equivalized, bottom_coded, estimation_object, tag = model_options
    if typeof(estimator) <: SeriesEstimator @unpack integral_pcf_grid, integral_cop_grid = estimator else @unpack grid_pcf, grid_cop = estimator end
    grid_choice_pcf = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf

    equiv     = equivalized == true ? "eq" : ""
    botcod    = !isempty(bottom_coded) ? "bc" : ""
    label     = "$case" * "_$equiv" * "$botcod"

    # Time params 
    @unpack tmin, tmax     = time_params  # m = model
    smin, smax             = user_t    # what user wants 

    @assert smin["year"] >= tmin["year"]
    @assert smax["year"] <= tmax["year"]
    
    # Distinguishing by plots 
    generate_c_shares_levels_quantiles(dv, counterfactuals_data_dict, avg_series, avg_series_cf, smin["year"], smax["year"], tmin, tmax, grid_choice_pcf, label, measures, tag) 

    # compare_to_counterfactual(dv, counterfactuals_data_dict, avg_series, avg_series_cf, time_params, user_t, model_options)
    # generate_counterfactual_plots(dv, counterfactuals_data_dict, avg_series, avg_series_cf, time_params, model_options)
end

    
function generate_counterfactual_agg_factors(agg_data, counterfactuals_dict_choice, time_params, model_elements, model_options)
    @unpack u_proj, agg_count = model_elements
    @unpack agg_lags          = model_options
    @unpack tmin, tmax        = time_params

    new_aggs                 = deepcopy(agg_data)
    
    # new_aggs[:, "date"] = QuarterlyDate.(aggregates[:, "time"]) 

    # Calculate the quarter correction
    q_correction = tmin["quarter"] > agg_lags ? tmin["quarter"] - agg_lags : 4 - (agg_lags - tmin["quarter"]) % 4

    # Calculate the year correction
    y_correction = tmin["quarter"] > agg_lags ? tmin["year"] : trunc(tmin["year"] - ((agg_lags + 3 - 1) / 4) + 1)

    filter!(row -> row.date >= QuarterlyDate(y_correction, q_correction), new_aggs)

    b = filter(row -> row.date <= QuarterlyDate(tmax["year"], tmax["quarter"]), new_aggs)
    # local b
    # if until != false
    #     b = filter(row -> row.date <= QuarterlyDate(until["year"], until["quarter"]), new_aggs)
    # else
    # end
    dropped_rows_ub = nrow(new_aggs) - nrow(b)
    println(dropped_rows_ub)

    select!(new_aggs, Not(["date", "time", "year", "quarter"]))
        
    if counterfactuals_dict_choice != []
        for s in counterfactuals_dict_choice
            new_aggs[:, s] .= 0  
        end
    else
        new_aggs .= 0  
    end

    
    new_aggs                 = copy(transpose(Matrix{Float64}(new_aggs)))  # removes date column, K x T
    super_aggs               = new_aggs[:, 1:end-agg_lags]  

    for i in 1:agg_lags
        super_aggs = vcat(super_aggs, new_aggs[:, i+1:end-agg_lags+i])
    end

    super_aggs              .= (super_aggs .- mean(super_aggs, dims=2)) ./ std(super_aggs, dims=2)
    super_aggs               = super_aggs[mapslices(col -> all((!isnan).(col)), super_aggs, dims = 2)[:], :]  # Only complete data 

    # Generate new factors 
    F_f              = u_proj' * Matrix(super_aggs)
    u                = F_f[1:agg_count, 1:end-dropped_rows_ub]  

    return u
end

    

# # random data
# rd = rand(100, 10)
# M    = fit(PCA, Matrix(rd'), pratio=1; method=:svd)  
# pcs  = MultivariateStats.transform(M, Matrix(rd'))

# # Get larger proj 
# proj    = projection(M)
# inv(proj)
# proj'

# return proj, pcs, M 


function generate_missing_info_plots(dv, c_dv, gdp_series, time_params, model_options)
    @unpack estimation_object, measures, grid = model_options
    @unpack tmin, tmax = time_params  
    D                  = length(measures)
    random_key         = first(keys(dv))
    T                  = size(dv[random_key][2], 2)

    # Transform pcf data to original form
    gdp_series[!, "date"]    = QuarterlyDate.(gdp_series[!, "time"])
    filter!(row -> row.date >= QuarterlyDate(tmin["year"], tmin["quarter"]), gdp_series)
    filter!(row -> row.date <= QuarterlyDate(tmax["year"], tmax["quarter"]), gdp_series)

    for (i,v) in enumerate(values(c_dv))
        for source in keys(v)
            if i == 1
                for t in 1:T
                    # v[source][2][:, t]          .= (abs.(v[source][2][:, t]) ./ sign.(v[source][2][:, t])).^3 .* gdp_series[t, :real_gdp_pc]
                    # dv[source][2][:, t]         .= (abs.(dv[source][2][:, t]) ./ sign.(dv[source][2][:, t])).^3 .* gdp_series[t, :real_gdp_pc]
                    v[source][2][:, t]          .= sinh.(v[source][2][:, t])  .* gdp_series[t, :real_gdp_pc]
                    dv[source][2][:, t]         .= sinh.(dv[source][2][:, t]) .* gdp_series[t, :real_gdp_pc]
                end
                levels, shares     = generate_shares_levels(dv[source][2], model_options, gdp_series)
                dv[source]         = create_time_series_dictionary([dv[source]..., levels, shares], grid, measures)
            else
                for t in 1:T
                    # v[source][2][:, t]          .= (abs.(v[source][2][:, t]) ./ sign.(v[source][2][:, t])).^3 .* gdp_series[t, :real_gdp_pc]
                    v[source][2][:, t]          .= sinh.(v[source][2][:, t]) .* gdp_series[t, :real_gdp_pc]
                end
            end
            levels, shares     = generate_shares_levels(v[source][2], model_options, gdp_series)
            v[source]          = create_time_series_dictionary([v[source]..., levels, shares], grid, measures)
        end
    end

    user_t  = (Dict("year" => 1965), tmax)
    compare_to_nonmissing(dv, c_dv, time_params, user_t, model_options)
end


# function generate_counterfactual_plots(dv, counterfactuals_data_dict, gdp_series, time_params, model_options)
#     @unpack estimation_object, measures, grid = model_options
#     @unpack tmin, tmax = time_params  
#     D                  = length(measures)
#     random_key         = first(keys(dv))
#     T                  = size(dv[random_key][2], 2)

#     # Transform pcf data to original form 
#     gdp_series[!, "date"]    = QuarterlyDate.(gdp_series[!, "time"])
#     correction_names         = [meas * "_per_hh" for meas in measures]
#     filter!(row -> row.date >= QuarterlyDate(tmin["year"], tmin["quarter"]), gdp_series)
#     filter!(row -> row.date <= QuarterlyDate(tmax["year"], tmax["quarter"]), gdp_series)
#     select_series            = select(gdp_series, correction_names)
    
#     for (i,v) in enumerate(values(counterfactuals_data_dict))
#         for source in collect(keys(v))
#             if i == 1
#                 v_split_pcfs        = [v[source][2][I, :] for I in Iterators.partition(axes(v[source][2], 1), grid)]  # split by measure 
#                 dv_split_pcfs       = [dv[source][2][I, :] for I in Iterators.partition(axes(dv[source][2], 1), grid)]  # split by measure 
#                 # for t in 1:T
#                 #     v[source][2][:, t]          .= vcat([v_split_pcfs[m][:, t] .* select_series[t, correction_names[m]] for m in eachindex(v_split_pcfs)]...) #v[source][2][:, t]  .* gdp_series[t, :real_gdp_pc]
#                 #     dv[source][2][:, t]         .= vcat([dv_split_pcfs[m][:, t] .* select_series[t, correction_names[m]] for m in eachindex(dv_split_pcfs)]...) #dv[source][2][:, t] .* gdp_series[t, :real_gdp_pc]
#                 # end
#                 levels, shares     = generate_shares_levels(dv[source][2], model_options, gdp_series)
#                 dv[source]         = create_time_series_dictionary([dv[source]..., levels, shares], grid, measures)
#             else
#                 v_split_pcfs        = [v[source][2][I, :] for I in Iterators.partition(axes(v[source][2], 1), grid)]  # split by measure 
#                 # for t in 1:T
#                 #     # v[source][2][:, t]          .= (abs.(v[source][2][:, t]) ./ sign.(v[source][2][:, t])).^3 .* gdp_series[t, :real_gdp_pc]
#                 #     # v[source][2][:, t]          .= sinh.(v[source][2][:, t]) .* gdp_series[t, :real_gdp_pc]
#                 #     v[source][2][:, t]          .= vcat([v_split_pcfs[m][:, t] .* select_series[t, correction_names[m]] for m in eachindex(v_split_pcfs)]...) #v[source][2][:, t] .* gdp_series[t, :real_gdp_pc]
#                 # end
#             end
#             levels, shares     = generate_shares_levels(v[source][2], model_options, gdp_series)
#             v[source]          = create_time_series_dictionary([v[source]..., levels, shares], grid, measures)
#         end
#     end
#     user_t  = (Dict("year" => 1965), tmax)

#     compare_to_counterfactual(dv, counterfactuals_data_dict, avg_series, time_params, user_t, model_options)
# end


function compare_to_nonmissing(dv, c_dv, time_params, user_t, model_options)
    @unpack measures, grid, case, equivalized, bottom_coded, estimation_object, blind_to, tag = model_options
    equiv     = equivalized == true ? "eq" : ""
    botcod    = !isempty(bottom_coded) ? "bc" : ""
    label     = "$case" * "_$equiv" * "$botcod"
    m_label   = measures_folder(measures)

    # Time params 
    @unpack tmin, tmax     = time_params  # m = model
    smin, smax             = user_t    # what user wants 

    @assert smin["year"] >= tmin["year"]
    @assert smax["year"] <= tmax["year"]

    dts       = QuarterlyDate(tmin["year"], tmin["quarter"]) : Quarter(1) : QuarterlyDate(tmax["year"], tmax["quarter"])
    init_path = BASE_PATH
    path      = init_path * "/7_Results/$m_label" * tag * "/other_results/counterfactuals/from_missing_information/correlations/"
    mkpath(path)
    xaxis     = collect(1:length(dts))

    # Distinguishing by plots. m = missing 
    generate_m_shares_levels_quantiles(dv, c_dv, smin["year"], smax["year"], tmin, tmax, grid, label, measures, blind_to, tag)

    # Exporting Kendalls Tau 
    D              = length(measures)
    random_key     = first(keys(dv))
    prodt          = prod(size(dv[random_key]["copulas"]["data"])[1:end-1])
    T              = size(dv[random_key]["copulas"]["data"])[end]
    pcfs           = zeros(grid * D, T)

    # TODO: put all the taus together in one plot/csv, wrap these all in "try's"
    # kend_τ = Vector{Matrix{Float64}}(undef, length(keys(dv)))
    # kend_τ = Dict()
    # for (i, source) in enumerate(collect(keys(dv)))
        # c_data_cop            = reshape(c_dv["Wealth"][source]["copulas"]["data"], (prodt, T))
        # pcfs                 .= vcat([c_dv["Wealth"][source][meas]["quantiles"]["data"] for meas in measures]...)
        # kend_τ["Wealth"]      = kendalls_tau(c_data_cop, pcfs, source, "other_results/counterfactuals/from_missing_information", time_params, model_options, true, false, false)

        # c_data_cop            = reshape(c_dv["Income"][source]["copulas"]["data"], (prodt, T))
        # pcfs                 .= vcat([c_dv["Income"][source][meas]["quantiles"]["data"] for meas in measures]...)
        # kend_τ["Income"]      = kendalls_tau(c_data_cop, pcfs, source, "other_results/counterfactuals/from_missing_information", time_params, model_options, true, false, false)

        # c_data_cop            = reshape(c_dv["Both"][source]["copulas"]["data"], (prodt, T))
        # pcfs                 .= vcat([c_dv["Both"][source][meas]["quantiles"]["data"] for meas in measures]...)
        # kend_τ["Both"]        = kendalls_tau(c_data_cop, pcfs, source, "other_results/counterfactuals/from_missing_information", time_params, model_options, true, false, false)

        # Plot 
        # Plots.plot(xaxis, kend_τ[1], xlabel = L"\textrm{Year}", ylabel = L"\textrm{Kendall's\,\,Tau}", xformatter=:latex, yformatter=:latex, xticks=(xaxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_,date) in enumerate(dts[1:20:end])]), legend=:outertopright, label=L"\textrm{No Wealth\,\,%$(source)}", lw=1, dpi=500)
        # Plots.plot!(xaxis, kend_τ[2], label=L"\textrm{No\,\, Income\,\,%$(source)}", lw=1, dpi=500)
        # Plots.plot!(xaxis, kend_τ[3], label=L"\textrm{No\,\, data\,\, from\,\,%$(source)}", lw=1, dpi=500)
        # Plots.savefig(path * "/correlations/" * source * "_kendalls_tau_" * label * ".pdf")

        # combs     = join.(combinations(measures, 2), "\\,")
        # for (i, comb) in enumerate(combs)
        #     for (j, scenario) in enumerate(collect(keys(kend_τ)))
        #         if j == 1
        #             Plots.plot(
        #                 xaxis, 
        #                 kend_τ[scenario][:, i], 
        #                 xlabel = L"\textrm{Year}", 
        #                 ylabel = L"\textrm{Kendall's\,\,Tau}", 
        #                 xformatter=:latex, 
        #                 yformatter=:latex, 
        #                 xticks=(xaxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_,date) in enumerate(dts[1:20:end])]), 
        #                 legend=:outertopright, 
        #                 legendtitle=L"\textrm{%$(comb)}", 
        #                 label=L"\textrm{%$(scenario)\,\,%$(source)}", 
        #                 lw=1, dpi=500
        #                 )
        #         else
        #             Plots.plot!(xaxis, kend_τ[scenario][:, i], label=L"\textrm{%$(scenario)\,\,%$(source)}", lw=1, dpi=500)
        #         end
        #     end
        #     Plots.savefig(path * source * "_" * join(split(comb, "\\,"), "_") * "_kendalls_tau_" * label * ".pdf")
        # end
    # end
end


# function compare_to_counterfactual(data_dict, counterfactuals_data_dict, avg_series, avg_series_cf, time_params, user_time_params, model_options)
    # @unpack measures, integral_pcf_grid, integral_cop_grid, case, equivalized, bottom_coded, estimation_object = model_options
    # equiv     = equivalized == true ? "eq" : ""
    # botcod    = !isempty(bottom_coded) ? "bc" : ""
    # label     = "$case" * "_$equiv" * "$botcod"
    # m_label   = measures_folder(measures)

    # # Time params 
    # @unpack tmin, tmax     = time_params  # m = model
    # smin, smax             = user_time_params    # what user wants 

    # @assert smin["year"] >= tmin["year"]
    # @assert smax["year"] <= tmax["year"]
    
    # # Distinguishing by plots 
    # generate_c_shares_levels_quantiles(data_dict, counterfactuals_data_dict, avg_series, avg_series_cf, smin["year"], smax["year"], tmin, tmax, integral_pcf_grid, integral_cop_grid, label, measures, tag) 

    # dts       = QuarterlyDate(tmin["year"], tmin["quarter"]) : Quarter(1) : QuarterlyDate(tmax["year"], tmax["quarter"]) 
    # init_path = BASE_PATH
    # path      = init_path * "/7_Results/$m_label/other_results/counterfactuals/from_aggregates/correlations/"
    # xaxis     = collect(1:length(dts))
    # combs     = join.(combinations(measures, 2), "\\,")
    # kend_τ    = Dict([scenario => zeros(length(dts), length(combs)) for scenario in collect(keys(counterfactuals_data_dict))])
    
    # generate_c_copula_plots(data_dict, counterfactuals_data_dict, smin["year"], smax["year"], tmin, tmax, grid, label, measures)

    # # Exporting Kendalls Tau 
    # D              = length(measures)
    # random_key     = first(keys(data_dict))
    # prodt          = prod(size(data_dict[random_key]["copulas"]["data"])[1:end-1])
    # T              = size(data_dict[random_key]["copulas"]["data"])[end]
    # c_data_cop     = zeros(prodt, T)
    # pcfs           = zeros(grid * D, T)

    # TODO: put all the taus together in one plot/csv 
    # for source in collect(keys(data_dict))
    #     plots = []  # Array to store plots
        
    #     for scenario in collect(keys(counterfactuals_data_dict))
    #         c_data_cop           .= reshape(counterfactuals_data_dict[scenario][source]["copulas"]["data"], (prodt, T))
    #         pcfs                 .= vcat([counterfactuals_data_dict[scenario][source][meas]["quantiles"]["data"] for meas in measures]...)
    #         kend_τ[scenario]      = kendalls_tau(c_data_cop, pcfs, scenario * "_" * source * "_counterfactual_", "other_results/counterfactuals/from_aggregates", time_params, model_options, true, false, false)
    #     end

    #     # Include reconstruction in these counterfactual plots for comparison 
    #     c_data_cop              .= reshape(data_dict[source]["copulas"]["data"], (prodt, T))
    #     pcfs                    .= vcat([data_dict[source][meas]["quantiles"]["data"] for meas in measures]...)
    #     kend_τ["reconstruction"] = kendalls_tau(c_data_cop, pcfs, "reconstruction_" * source * "_counterfactual_", "other_results/counterfactuals/from_aggregates", time_params, model_options, true, false, false)
        
        # for (i, comb) in enumerate(combs)
        #     for (j, scenario) in enumerate(collect(keys(counterfactuals_data_dict)))
        #         if j == 1
        #             Plots.plot(xaxis, kend_τ[scenario][:, i], xlabel = L"\textrm{Year}", ylabel = L"\textrm{Kendall's\,\,Tau}", xformatter=:latex, yformatter=:latex, xticks=(xaxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_,date) in enumerate(dts[1:20:end])]), legend=:outertopright, legendtitle=L"\textrm{%$(comb)}", label=L"\textrm{%$(scenario)\,\,%$(source)}", lw=1, dpi=500)
        #         else
        #             Plots.plot!(xaxis, kend_τ[scenario][:, i], label=L"\textrm{%$(scenario)\,\,%$(source)}", lw=1, dpi=500)
        #         end
        #     end
        #     Plots.plot!(xaxis, kend_τ["reconstruction"][:, i], label=L"\textrm{reconstruction\,\,%$(source)}", lw=1, dpi=500)
            
        #     Plots.savefig(path * source * "_" * join(split(comb, "\\,"), "_") * "_kendalls_tau_" * label * ".pdf")
        # end
    # end
# end


function generate_m_shares_levels_quantiles(data_dict, counterfactuals_data_dict, smin, smax, tmin, tmax, grid, label, measures, blind_to, run_tag="") 
    # Subset data for the reconstruction 
    base_jump = (smin - tmin["year"]) * 4 + 1  
    end_jump  = (tmax["year"] - smax) * 4  
    m_label   = measures_folder(measures)

    # Dates and such
    dts       = QuarterlyDate(smin, tmin["quarter"]) : Quarter(1) : QuarterlyDate(smax, tmax["quarter"])
    init_path = BASE_PATH
    path      = init_path * "/7_Results/$m_label" * run_tag * "/other_results/counterfactuals/from_missing_information/"
    for m in measures
        mkpath(path * "/$m/quantiles_levels/")
        mkpath(path * "/$m/shares/")
    end
    xaxis     = collect(1:length(dts))

    # Import hh_gdp 
    hh_gdp = CSV.read(init_path * "/2_Data_processing/HH_GDP.csv", DataFrame, header = 1, delim = ",")
    filter!(row -> row.date >= QuarterlyDate(smin, tmin["quarter"]), hh_gdp)
    filter!(row -> row.date <= QuarterlyDate(smax, tmax["quarter"]), hh_gdp)


    # Storing for later 
    cf_dict = Dict()
    for k in keys(counterfactuals_data_dict)
        cf_dict[k] = Dict()
        cf_dict[k]["levels"]  = zeros(length(dts), 3)
        cf_dict[k]["shares"]  = zeros(length(dts), 3)
        cf_dict[k]["bottomq"] = zeros(3, length(dts))
        cf_dict[k]["topq"]    = zeros(length(dts))
    end

    # Pre-allocate vectors 
    lv10   = zeros(length(dts))
    sh10   = zeros(length(dts))
    lv40   = zeros(length(dts))
    sh40   = zeros(length(dts))
    lv50   = zeros(length(dts))
    sh50   = zeros(length(dts))
    qu     = zeros(3, length(dts))
    top_qu = zeros(length(dts)) 
    all_lv = zeros(length(dts), 3)
    all_sh = zeros(length(dts), 3)

    # The label for the plots 
    miss_label = L"\textrm{Blind\,\,to:}" # * "\n" * join([L"\textrm{%$(k)\,\,%$(v)}" for k in keys(blind_to) for v in blind_to[k]], "\n")

    # Plotting
    for source in collect(keys(data_dict))
        for meas in measures 
            println("Generating plots for $meas")
            if grid == 10
            # For the reconstruction 
                # Top 10
                lv10 .= data_dict[source][meas]["levels"]["common series"]["top10"][base_jump:end-end_jump]
                sh10 .= data_dict[source][meas]["shares"]["common series"]["top10"][base_jump:end-end_jump]
                
                # Next 40
                lv40 .= data_dict[source][meas]["levels"]["common series"]["next40"][base_jump:end-end_jump]
                sh40 .= data_dict[source][meas]["shares"]["common series"]["next40"][base_jump:end-end_jump]

                # Botttom 50
                lv50 .= data_dict[source][meas]["levels"]["common series"]["bottom50"][base_jump:end-end_jump]
                sh50 .= data_dict[source][meas]["shares"]["common series"]["bottom50"][base_jump:end-end_jump]

                # The bottom 20th, median, the 90th and the top 
                qu       .= data_dict[source][meas]["quantiles"]["data"][[2, 5, 9], base_jump:end-end_jump] ./ hh_gdp[!, "hh_contribution"]'
                top_qu   .= data_dict[source][meas]["quantiles"]["data"][grid, base_jump:end-end_jump] ./ hh_gdp[!, "hh_contribution"]

                all_lv .= hcat(lv10, lv40, lv50)
                all_sh .= hcat(sh10, sh40, sh50)

                # Loop over every scenario 
                for (k,v) in counterfactuals_data_dict
                    lv10 .= v[source][meas]["levels"]["common series"]["top10"][base_jump:end-end_jump]
                    sh10 .= v[source][meas]["shares"]["common series"]["top10"][base_jump:end-end_jump]

                    lv40 .= v[source][meas]["levels"]["common series"]["next40"][base_jump:end-end_jump]
                    sh40 .= v[source][meas]["shares"]["common series"]["next40"][base_jump:end-end_jump]

                    lv50 .= v[source][meas]["levels"]["common series"]["bottom50"][base_jump:end-end_jump]
                    sh50 .= v[source][meas]["shares"]["common series"]["bottom50"][base_jump:end-end_jump]

                    cf_dict[k]["bottomq"]  .= v[source][meas]["quantiles"]["data"][[2, 5, 9], base_jump:end-end_jump] ./ hh_gdp[!, "hh_contribution"]'
                    cf_dict[k]["topq"]     .= v[source][meas]["quantiles"]["data"][grid, base_jump:end-end_jump] ./ hh_gdp[!, "hh_contribution"]

                    cf_dict[k]["levels"] .= hcat(lv10, lv40, lv50)
                    cf_dict[k]["shares"] .= hcat(sh10, sh40, sh50)
                end

            elseif grid == 5
                error("Plot is not possible yet!")
            elseif grid == 20
                # TODO:
            end
            # opacity    = [0.8, 0.5, 0.2] 
            line_style     = [:solid, :dash, :dot, :dashdot]
            M              = uppercasefirst(meas)
            series(source) = [L"\textrm{%$(source)\,\,Top\,10}" L"\textrm{%$(source)\,\,Next\,40}" L"\textrm{%$(source)\,\,Bottom\,50}"]
            cf_series(k)   = [L"\textrm{%$(k)\,\,Top\,10}" L"\textrm{%$(k)\,\,Next\,40}" L"\textrm{%$(k)\,\,Bottom\,50}"]
            tag            = ["top10_", "next40_", "bottom50_"]

            # Plot the shares against one another 
            for i in axes(all_sh, 2)
                Plots.plot(xaxis, 
                            all_sh[:, i], 
                            xlabel = L"\textrm{Year}",
                            ylabel = L"\textrm{%$(M)}\,\,\textrm{shares}",
                            xformatter=:latex, 
                            yformatter=:latex, 
                            xticks=(xaxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_,date) in enumerate(dts[1:20:end])]),
                            legend=:outertopright,
                            legendtitle=miss_label,
                            label=series(source)[i],
                            lw=1, dpi=500
                        )
                for (j,k) in enumerate(collect(keys(cf_dict)))
                    Plots.plot!(xaxis, cf_dict[k]["shares"][:, i], lc=:black, ls=line_style[j], label=cf_series(k)[i], lw=1, dpi=500)
                end
                Plots.savefig(path * "/$meas/" * "/shares/" * source * "_" * tag[i] * label * ".pdf")
            end

            # Levels 
            for i in axes(all_lv, 2)
                Plots.plot(xaxis, 
                            all_lv[:, i], 
                            xlabel = L"\textrm{Year}",
                            ylabel = L"\textrm{%$(M)}\,\,\textrm{levels}",
                            xformatter=:latex, 
                            yformatter=:latex, 
                            xticks=(xaxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_,date) in enumerate(dts[1:20:end])]),
                            legend=:outertopright,
                            legendtitle=miss_label,
                            label=series(source)[i],
                            lw=1, dpi=500
                        )
                for (j,k) in enumerate(collect(keys(cf_dict)))
                    Plots.plot!(xaxis, cf_dict[k]["levels"][:, i], lc=:black, ls=line_style[j], label=cf_series(k)[i], lw=1, dpi=500)
                end
                Plots.savefig(path * "/$meas/" * "/quantiles_levels/" * source * "_" * "levels_" * tag[i] * label * ".pdf")
            end
            # Quantiles 
            qseries(source) = [L"\textrm{%$(source)\,\,\,20th\,pct.}" L"\textrm{%$(source)\,\,\,Median}" L"\textrm{%$(source)\,\,\,90th\,pct.}"]
            qcfseries(k)    = [L"\textrm{%$(k)\,\,20th\,pct.}" L"\textrm{%$(k)\,\,Median}" L"\textrm{%$(k)\,\,90th\,pct.}"]
            qtag            = ["20pct_", "50pct_", "90pct_"]
            for i in axes(all_lv, 2)
                Plots.plot(xaxis, 
                            qu'[:, i],
                            xlabel = L"\textrm{Year}",
                            ylabel = L"\textrm{%$(M)}\,\,\textrm{quantiles}",
                            xformatter=:latex, 
                            yformatter=:latex, 
                            xticks=(xaxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_,date) in enumerate(dts[1:20:end])]),
                            legendtitle=miss_label,
                            legend=:outertopright,
                            label=qseries(source)[i], 
                            lw=1, dpi=500
                        )

                for (j,k) in enumerate(collect(keys(cf_dict)))
                    Plots.plot!(xaxis, cf_dict[k]["bottomq"]'[:, i], lc=:black, ls=line_style[j], label=qcfseries(k)[i], lw=1, dpi=500)
                end
                Plots.savefig(path * "/$meas/" * "/quantiles_levels/" * source * "_" * qtag[i] * label * ".pdf")
            end

            # Top quantile  
            Plots.plot(xaxis, 
                        top_qu,
                        xlabel = L"\textrm{Year}",
                        ylabel = L"\textrm{%$(M)}\,\,\textrm{quantiles}",
                        xformatter=:latex, 
                        yformatter=:latex, 
                        xticks=(xaxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_,date) in enumerate(dts[1:20:end])]),
                        legend=:outertopright,
                        legendtitle=miss_label,
                        lc=:blue,
                        label=L"\textrm{%$(source)\,\,\,Top\,10}", 
                        lw=1, dpi=500
                    )
            for (j,k) in enumerate(collect(keys(cf_dict)))
                Plots.plot!(xaxis, cf_dict[k]["topq"], ls=line_style[j], lc=:black, label=L"\textrm{%$(k)\,\,Top\,10}", lw=1, dpi=500)
            end
            Plots.savefig(path * "/$meas/" * "/quantiles_levels/" * source * "_" * "top10pct_" * label * ".pdf")
        end
    end
end


function generate_c_shares_levels_quantiles(data_dict, counterfactuals_data_dict, avg_series, avg_series_cf, smin, smax, tmin, tmax, grid, label, measures, tag)
    # Subset data for the reconstruction 
    base_jump = (smin - tmin["year"]) * 4 + 1  
    end_jump  = (tmax["year"] - smax) * 4  
    m_label   = measures_folder(measures)

    # Dates and such
    dts       = QuarterlyDate(smin, tmin["quarter"]) : Quarter(1) : QuarterlyDate(smax, tmax["quarter"])
    init_path = BASE_PATH
    path      = init_path * "/7_Results/$m_label" * tag * "/other_results/counterfactuals/from_aggregates/"
    for m in measures
        mkpath(path * "/$m/quantiles_levels/")
    end
    xaxis     = collect(1:length(dts))

    # Importing HH GDP 
    # hh_gdp = CSV.read(, DataFrame, header=true, delim=',')
    # hh_gdp                 = DataFrame(XLSX.readtable(init_path * "/2_Data_processing/HH_GDP.xlsx", "HH_GDP", header=true,))
    # hh_gdp[!, "date"]      = QuarterlyDate.(hh_gdp[!, "date"])
    # filter!(row -> row.date >= QuarterlyDate(smin, tmin["quarter"]), hh_gdp)
    # filter!(row -> row.date <= QuarterlyDate(smax, tmax["quarter"]), hh_gdp)


    # Storing for later 
    cf_dict = Dict()
    for k in keys(counterfactuals_dict)
        cf_dict[k] = Dict()
        # cf_dict[k]["levels"]  = zeros(length(dts), 3)
        # cf_dict[k]["shares"]  = zeros(length(dts), 3)
        cf_dict[k]["bottomq"] = zeros(3, length(dts))
        cf_dict[k]["topq"]    = zeros(length(dts))
    end

    # Pre-allocate vectors 
    # lv10   = zeros(length(dts))
    # sh10   = zeros(length(dts))
    # lv40   = zeros(length(dts))
    # sh40   = zeros(length(dts))
    # lv50   = zeros(length(dts))
    # sh50   = zeros(length(dts))
    qu     = zeros(3, length(dts))
    top_qu = zeros(length(dts)) 
    # all_lv = zeros(length(dts), 3)
    # all_sh = zeros(length(dts), 3)

    for source in collect(keys(data_dict))
        for meas in measures 
            println("Generating plots for $meas")
            # if grid == 10
            # For the reconstruction 
                # # Top 10
                # lv10 = data_dict[source][meas]["levels"]["common series"]["top10"][base_jump:end-end_jump]
                # sh10 = data_dict[source][meas]["shares"]["common series"]["top10"][base_jump:end-end_jump]
                
                # # Next 40
                # lv40 = data_dict[source][meas]["levels"]["common series"]["next40"][base_jump:end-end_jump]
                # sh40 = data_dict[source][meas]["shares"]["common series"]["next40"][base_jump:end-end_jump]

                # # Botttom 50
                # lv50 = data_dict[source][meas]["levels"]["common series"]["bottom50"][base_jump:end-end_jump]
                # sh50 = data_dict[source][meas]["shares"]["common series"]["bottom50"][base_jump:end-end_jump]

            # The bottom 20th, median, the 90th and the top 
            qu       = data_dict[source][meas]["quantiles"]["data"][[2, 5, 9], base_jump:end-end_jump] ./ avg_series[source][base_jump:end-end_jump, meas * "_per_hh"]'
            top_qu   = vec(data_dict[source][meas]["quantiles"]["data"][grid, base_jump:end-end_jump]) ./ vec(avg_series[source][base_jump:end-end_jump, meas * "_per_hh"]') # issues with Julia reading both as matrices

            # all_lv = hcat(lv10, lv40, lv50)
            # all_sh = hcat(sh10, sh40, sh50)

                # Loop over every scenario 
                for (k,v) in counterfactuals_data_dict
                    # lv10 = v[source][meas]["levels"]["common series"]["top10"][base_jump:end-end_jump]
                    # sh10 = v[source][meas]["shares"]["common series"]["top10"][base_jump:end-end_jump]

                    # lv40 = v[source][meas]["levels"]["common series"]["next40"][base_jump:end-end_jump]
                    # sh40 = v[source][meas]["shares"]["common series"]["next40"][base_jump:end-end_jump]

                    # lv50 = v[source][meas]["levels"]["common series"]["bottom50"][base_jump:end-end_jump]
                    # sh50 = v[source][meas]["shares"]["common series"]["bottom50"][base_jump:end-end_jump]

                    cf_dict[k]["bottomq"]  = v[source][meas]["quantiles"]["data"][[2, 5, 9], base_jump:end-end_jump] ./ avg_series_cf[k][source][base_jump:end-end_jump, meas * "_per_hh"]'
                    cf_dict[k]["topq"]     = vec(v[source][meas]["quantiles"]["data"][grid, base_jump:end-end_jump]) ./ vec(avg_series_cf[k][source][base_jump:end-end_jump, meas * "_per_hh"]')

                    # cf_dict[k]["levels"] = hcat(lv10, lv40, lv50)
                    # cf_dict[k]["shares"] = hcat(sh10, sh40, sh50)
                end

            # opacity    = [0.8, 0.5, 0.2] 
            line_style = [:dash, :dot, :dashdot, :dashdotdot, :solid, :dash]
            line_color = [:black, :black, :black, :black, :black, :red, :red]
            M          = uppercasefirst(meas)
            
            # Plot the shares against one another 
            series(source) = [L"\textrm{%$(source)\,\,Top\,10}" L"\textrm{%$(source)\,\,Next\,40}" L"\textrm{%$(source)\,\,Bottom\,50}"]
            cf_series(k)   = [L"\textrm{%$(k)\,\,Top\,10}" L"\textrm{%$(k)\,\,Next\,40}" L"\textrm{%$(k)\,\,Bottom\,50}"]

            tag    = ["top10_", "next40_", "bottom50_"]

            # Quantiles 
            qseries = ["20th", "50th", "90th"]
            qcfseries(k)    = k == "AP" ? L"\textrm{Asset\,\, prices \,\, muted}" : k == "RE" ? L"\textrm{Real\,\, economy\,\, muted}" : k == "Nothing" ? L"\textrm{All\,\, muted}" : k == "MP" ? L"\textrm{Monetary\,\, muted}" : k == "UN" ? L"\textrm{Employment\,\, muted}" : L"\textrm{No\,\,data\,\,from\,\,%$(k)}"
            qtag    = ["20pct_", "50pct_", "90pct_"]

            plot_name        = occursin("CEX", source) ? "CEX" : source 

            for i in axes(qu, 1)
                cond     = source != "consensus" ? findall(!isnan, qu[i, :]) : [1] # .!isnan.(all_lv_o[1][:, j])

                if isempty(cond)
                    println(source)
                    println(qu'[:, i])
                    break
                end

                s_axis   = xaxis[cond[1]:cond[end]] # start at the first observation 
                s_data   = qu[i, cond[1]:cond[end]]
                s_dts    = dts[cond[1]:cond[end]]

                for (j,k) in enumerate(collect(keys(cf_dict)))
                    println(size(cf_dict[k]["bottomq"]))
                    println(size(cf_dict[k]["topq"]))
                    Plots.plot(s_axis, 
                                s_data,
                                xlabel = L"\textrm{Year}",
                                ylabel = L"\textrm{%$(M)\,\, rel.\,\, to\,\, average}",
                                xformatter=:latex, 
                                yformatter=:latex, 
                                xtickfontsize=10,
                                ytickfontsize=10,
                                legendfontsize=10,
                                guidefontsize=14,
                                lc = select_color(plot_name),
                                xticks=(s_axis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_,date) in enumerate(s_dts[1:20:end])]),
                                legend=:best,
                                label=L"\textrm{Model}", 
                                lw=3, dpi=500
                            )

                    Plots.plot!(s_axis, cf_dict[k]["bottomq"][i, cond[1]:cond[end]], lc=:lightblue, la=0.7, ls=line_style[j], label=qcfseries(k), lw=2, dpi=500)
                    Plots.savefig(path * "/$meas/" * "/quantiles_levels/" * k * "_" * plot_name * "_" * qtag[i] * label * ".pdf")
                end
            end

            # Top quantile  
            println(size(top_qu))
            cond     = source != "consensus" ? findall(!isnan, top_qu) : [1] # .!isnan.(all_lv_o[1][:, j])
            if isempty(cond)
                break
            end
            s_axis   = xaxis[cond[1]:cond[end]] # start at the first observation 
            s_dts    = dts[cond[1]:cond[end]]

            for (j,k) in enumerate(collect(keys(cf_dict)))
                Plots.plot(s_axis, 
                            top_qu[cond[1]:cond[end]],
                            xlabel = L"\textrm{Year}",
                            ylabel = L"\textrm{%$(M)\,\, rel.\,\, to\,\, average}",
                            xformatter=:latex, 
                            yformatter=:latex, 
                            xtickfontsize=10,
                            ytickfontsize=10,
                            legendfontsize=10,
                            guidefontsize=14,
                            lc = select_color(plot_name),
                            xticks=(s_axis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_,date) in enumerate(s_dts[1:20:end])]),
                            legend=:best,
                            label=L"\textrm{Model}", 
                            lw=3, dpi=500
                        )
                Plots.plot!(s_axis, cf_dict[k]["topq"][cond[1]:cond[end]], lc=:lightblue, la=0.7, ls=line_style[j], label=qcfseries(k), lw=2, dpi=500)
                Plots.savefig(path * "/$meas/" * "/quantiles_levels/" * k * "_" * plot_name * "_" * "top10pct_" * label * ".pdf")
            end
        end
    end
end
