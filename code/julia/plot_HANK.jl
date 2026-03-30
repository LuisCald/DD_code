
function generate_quantiles_shares_levels_HANK(data_dict, ty, func_data, data_name, smin, smax, tmin, tmax, estimator, label, type, measures, time_params, select_series, gdp_series, posterior_bounds, compare_to_other_est, tag, freq)

    # --- helper: linear interpolation onto a target grid ---
    # We use Interpolations.jl (already used elsewhere in the repo) instead of an
    # undefined `interpolate(...)` helper.
    #
    # Behavior:
    # - linear interpolation between observed points
    # - flat extrapolation at the ends
    function interp1_linear(x_obs::AbstractVector, y_obs::AbstractVector, x_tgt::AbstractVector)
        @assert length(x_obs) == length(y_obs)
        itp = linear_interpolation(x_obs, y_obs, extrapolation_bc=Line())
        return itp.(x_tgt)
    end

    @unpack func_dict, year_vec, data_sources, confidence_intervals = func_data
    @unpack time_dict = time_params
    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid, integral_cop_grid = estimator
    else
        @unpack grid_pcf, grid_cop = estimator
    end

    grid_choice_pcf = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf
    grid_choice_cop = typeof(estimator) <: SeriesEstimator ? integral_cop_grid : grid_cop

    # Subset data for the reconstruction
    base_jump, end_jump = find_subset_frame(smin, smax, tmin, tmax)
    m_label = measures_folder(measures)
    folder = type == :optimization ? "from_optimization" : "from_mcmc"
    init_path = BASE_PATH
    path = init_path * "/7_Results/$m_label" * "$tag" * "/$folder/plots/"

    # Estimation dates
    dts = QuarterlyDate(smin["year"], smin["quarter"]):Quarter(1):QuarterlyDate(smax["year"], smax["quarter"])

    # Import the confidence intervals from the full HANK model 

    # Find the integer from tag 
    economy_number = split(tag, " ")[end]
    # file_name = init_path * "/2_Data_processing/confidence_intervals/ci_draws_illiqd_and_income_and_liquid_quintiles_series_HANK full $(economy_number)_all.jld2"
    # raw_ci = jldopen(file_name, "r")["ci"]

    # file_name_raw_data = "/home/luisc/Distributional_Dynamics/7_Results/illiqd_and_income_and_liquid HANK full $(economy_number)/from_mcmc/data/HANK_full_$(economy_number).jld2"
    # func_dict_full_hank = jldopen(file_name_raw_data, "r")["data"]
    # Import truth
    truth_file = DATA_PROCESSING * "/HANK_truth_$(economy_number).csv"
    truth_data = CSV.read(truth_file, DataFrame)

    # # Find the bounds of each 
    # hank_year_vec = repeat([i for i in collect(tmin["year"]:tmax["year"])], inner=4)
    # confidence_intervals = Dict("HANK" => Dict())
    # confidence_intervals["HANK"]["ci_u"], confidence_intervals["HANK"]["ci_l"] = construct_confidence_intervals(raw_ci, 0.025, 0.975, measures, hank_year_vec, estimator)

    # Correct the confidence intervals to be in the same time frame as the estimation
    hank_full_dts = QuarterlyDate(tmin["year"], 1):Quarter(1):QuarterlyDate(tmax["year"], 4)

    # Find the indices of the estimation time frame in the HANK full time frame
    hank_indices = findall(x -> x in dts, hank_full_dts)

    # Sometimes it may not have a year-quarter column 
    if :year ∉ names(truth_data) && :quarter ∉ names(truth_data)
        truth_data.time = QuarterlyDate.(truth_data.time)
        truth_data[:, :year] = year.(truth_data.time)
        truth_data[:, :quarter] = quarter.(truth_data.time)
    end

    truth_data_filtered = filter(row -> QuarterlyDate(row.year, row.quarter) in dts, truth_data)
    # for meas in measures
    #     for o in ["quantiles"]
    #         # println("HANK ci_u size: ", size(confidence_intervals["HANK"]["ci_u"][meas][o]))
    #         # println("HANK full size: ", size(func_dict_full_hank["HANK full $(economy_number)"][meas][o]["data"]))
    #         # confidence_intervals["HANK"]["ci_u"][meas][o] = confidence_intervals["HANK"]["ci_u"][meas][o][:, hank_indices]
    #         # confidence_intervals["HANK"]["ci_l"][meas][o] = confidence_intervals["HANK"]["ci_l"][meas][o][:, hank_indices]
    #         # func_dict_full_hank["HANK full $(economy_number)"][meas][o]["data"] = func_dict_full_hank["HANK full $(economy_number)"][meas][o]["data"][:, hank_indices]
    #     end
    # end


    local bot, mid, top
    if grid_choice_pcf == 10 || grid_choice_pcf == 100 || grid_choice_pcf == 20
        bot, mid, top = "bottom50", "next40", "top10"
    elseif grid_choice_pcf == 5
        bot, mid, top = "bottom40", "next40", "top20"
    end

    # Get the observed measures
    obs_meas = get_obs_meas(func_dict, data_name, measures; top=top)

    within_stat_dict = Dict()

    # Create a dictionary to hold the correlations
    corr_dict = Dict()
    plot_name = data_name

    # Generate the within statistic for the copulas 
    dimension = length(measures)

    # if data_name != "consensus" && ty == "normal" # ty == "average" doesn't work because it has no comparison to the data
    #     within_stat_dict["copula"] = compute_copula_within_stat(data_dict, confidence_intervals, base_jump, end_jump, data_name, dimension, grid_choice_cop)
    # end

    # println("select_series size: ", size(select_series))
    # println("select_series content: ", select_series[base_jump:end-end_jump, :])
    # println("base_jump: ", base_jump, " end_jump: ", end_jump)


    for meas in obs_meas # TODO: this issue here is that not all measures are observed ... ofc, we can use the reconstructed data but not the confidence intervals

        # All quantiles 
        qu = data_dict[meas]["quantiles"]["data"][:, base_jump:end-end_jump] ./ select_series[base_jump:end-end_jump, meas*"_per_hh"]'

        # For the (o)bservations 
        qu_o = Vector{Any}(undef, 4)

        # All quantiles 
        qu_o[1] = func_dict[data_name][meas]["quantiles"]["data"][:, base_jump:end-end_jump] ./ select_series[base_jump:end-end_jump, meas*"_per_hh"]'
        # qu_o[2] = confidence_intervals["HANK"]["ci_l"][meas]["quantiles"][:, base_jump:end-end_jump] ./ select_series[base_jump:end-end_jump, meas*"_per_hh"]'
        # qu_o[3] = confidence_intervals["HANK"]["ci_u"][meas]["quantiles"][:, base_jump:end-end_jump] ./ select_series[base_jump:end-end_jump, meas*"_per_hh"]'

        # Quantiles from full HANK model
        # smin_full = Dict("year" => 2001, "quarter" => 1)
        # smax_full = Dict("year" => 2024, "quarter" => 4)
        # base_jump_full, end_jump_full = find_subset_frame(smin, smax, smin_full, smax_full)
        # println("HANK full jumps: ", base_jump_full, " ", end_jump_full)
        # qu_o[4] = func_dict_full_hank["HANK full $(economy_number)"][meas]["quantiles"]["data"] ./ select_series[base_jump:end-end_jump, meas*"_per_hh"]'

        # Plots 
        markers = [:diamond, :utriangle, :circle, :star5, :cross, :diamond, :utriangle, :circle, :star5, :cross, :diamond, :utriangle, :circle, :star5, :cross, :diamond, :utriangle, :circle, :star5, :cross] #collect(range(0.3, 1, length=length(data_sources)))

        # x-axis
        xaxis = collect(1:length(dts))
        M = uppercasefirst(meas)

        # Important tag 
        detrended_or_not = ty == "average" ? "detrended_" : ""

        # All possible line-styles
        line_styles = [:solid :dash :dot :dashdot :dashdotdot :solid :dash :dot :dashdot :dashdotdot] # has to be a matrix

        within_stat_dict[meas] = Dict()
        corr_dict[meas] = Dict()

        # Doing quantiles now 
        label_quantiles = label_qs(grid_choice_pcf) # e.g., L"\textrm{%$(i)th}"
        sequences = define_sequences(grid_choice_pcf) # gets indices of each group
        dist_dict = Dict("bottom" => [sequences[1]], "middle" => [sequences[2]], "top" => [sequences[3]])
        log_qu = log_transformation(deepcopy(qu))


        for (obj, dist) in dist_dict
            corr_dict[meas][obj] = Dict()
            within_stat_dict[meas][obj] = Dict()

            cond = data_name != "consensus" ? findall(!isnan, qu_o[1][dist[1][1], :]) : [1] # .!isnan.(all_lv_o[1][:, j])

            est_ids = cond[1]:cond[end]
            s_axis = xaxis[est_ids] # start at the first observation
            s_data = qu[:, est_ids]
            s_dts = dts[est_ids]
            intd = 40

            # Plotting model estimates
            for (lsⱼ, j) in enumerate(dist[1])
                Plots.plot()
                Plots.plot!(s_axis,
                    s_data[j, :],
                    ylabel=L"\textrm{%$(M)\, \, rel.\,  to\,\,  average}",
                    lc=:red,
                    xformatter=:latex,
                    yformatter=:latex,
                    xtickfontsize=14,
                    ytickfontsize=14,
                    legendfontsize=10,
                    guidefontsize=14,
                    # xticks=(s_axis[1:intd:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(s_dts[1:intd:end])]),
                    legend=:best,
                    label="",
                    lw=4, dpi=500, ls=:solid,
                )

                # Plots.plot!(s_axis,
                #     qu_o[2][j, :][est_ids],
                #     fillrange=qu_o[3][j, :][est_ids],
                #     fillalpha=0.3,
                #     fillcolor=:red,
                #     la=0.0,
                #     lc=:white,
                #     lw=4, dpi=500,
                #     label="",
                # )

                # Find the variable from truth_df, and then scale it
                q_tag = integral_pcf_grid == 5 ? "q" : integral_pcf_grid == 10 ? "" : integral_pcf_grid == 20 ? "v" : "unknown"
                truth_var_name = Symbol(meas * "$(j)" * q_tag)
                # meas_name = meas == "consum" ? "consumption" : meas
                avg_name = Symbol(meas * "_per_hh")
                truth_values = truth_data_filtered[:, truth_var_name] ./ truth_data_filtered[base_jump:end-end_jump, meas*"_per_hh"] #select_series[base_jump:end-end_jump, meas*"_per_hh"] # truth_data_filtered[:, avg_name] #

                Plots.plot!(s_axis,
                    truth_values[est_ids],
                    # qu_o[4][j, :][est_ids],
                    lc=:blue,
                    lw=4,
                    la=0.5,
                    xtickfontsize=14,
                    ytickfontsize=14,
                    legendfontsize=10,
                    guidefontsize=14,
                    dpi=500, ls=:dot,
                    label=""
                )

                # Compute correlation between model and truth
                corr_n = round(cor(s_data[j, :], truth_values[est_ids]) * 100, digits=0)

                # R² to truth (demeaned): 1 - MSE / Var(truth) on demeaned series
                model_dm = s_data[j, :] .- mean(s_data[j, :])
                truth_dm = truth_values[est_ids] .- mean(truth_values[est_ids])
                mse_val = sum((model_dm .- truth_dm).^2) / length(est_ids)
                var_truth = var(truth_values[est_ids])
                r2_val = var_truth > 0 ? round((1.0 - mse_val / var_truth) * 100, digits=0) : NaN

                # Keep data points for scatter plot
                c_data = Vector{Any}(undef, 3)
                c_data[1] = qu_o[1][j, :][cond]

                # c_data[2] = qu_o[2][j, est_ids]
                # c_data[3] = qu_o[3][j, est_ids]

                # See how many points fall within the confidence intervals
                # r_data = qu[j, est_ids] # estimates that correspond to the indices of the data points 
                # num = count(c_data[2] .<= r_data .<= c_data[3])
                # den = length(r_data)
                # within_stat = floor(Int, (num ./ den) * 100)

                # within_stat_dict[meas][obj][label_quantiles[:, j][1]] = "$num" * "/" * "$den"

                Plots.scatter!(xaxis[cond],
                    c_data[1],
                    marker=:utriangle, #markers[j],
                    markercolor=:black,
                    ms=5,
                    la=0.5,
                    lw=2, dpi=500,
                    label=""
                    # label=L"\textrm{HANK \,\,Data\,\, Used}",
                )

                # Interpolation line of the data points
                Plots.plot!(xaxis[cond],
                    c_data[1],
                    la=0.5,
                    lw=2, dpi=500,
                    label="",
                    ls=:dash
                    # label=L"\textrm{HANK \,\,Data\,\, Used}",
                )

                # Store the correlation in the dictionary
                try
                    # Remove the .0 from the number
                    corr_dict[meas][obj][label_quantiles[:, j][1]] = Int(corr_n)
                catch e
                    corr_dict[meas][obj][label_quantiles[:, j][1]] = "-"
                end

                # Plots.plot!([], [], ls=:dash, lc=:black, la=0.0, label=L"\textrm{Within\, \,  bounds: %$(within_stat)\%}",)
                Plots.plot!([], [], ls=:dot, lc=:black, la=0.0, label=L"\textrm{Corr.\, to\,\, Truth\,: %$(corr_n)\%}",)
                Plots.plot!([], [], ls=:dot, lc=:black, la=0.0, label=L"R^2\textrm{\, to\,\, Truth\,: %$(r2_val)\%}",)
                Plots.savefig(path * "$meas/" * "quantiles_levels/" * plot_name * "_$meas" * "_$obj" * "_quantiles_$(j)" * detrended_or_not * label * ".pdf")
            end
        end

        # --- Share-group comparison: bot50, mid40, top10 (relative to mean) ---
        share_labels = ["bot50", "mid40", "top10"]
        share_plot_titles = [L"\textrm{Bottom\,\, 50\%}", L"\textrm{Middle\,\, 40\%}", L"\textrm{Top\,\, 10\%}"]

        for (si, sl) in enumerate(share_labels)
            # Model reconstruction
            if !haskey(data_dict[meas]["quantiles"]["common series"], sl)
                continue
            end
            model_share = data_dict[meas]["quantiles"]["common series"][sl][base_jump:end-end_jump] ./ select_series[base_jump:end-end_jump, meas*"_per_hh"]

            # Truth
            truth_col = Symbol(meas * "_" * sl)
            if !hasproperty(truth_data_filtered, truth_col)
                @warn "Truth column $truth_col not found, skipping share-group plot for $meas $sl"
                continue
            end
            truth_share = truth_data_filtered[:, truth_col] ./ select_series[base_jump:end-end_jump, meas*"_per_hh"]

            # Observed data (from func_dict)
            obs_share = nothing
            if haskey(func_dict[data_name][meas]["quantiles"]["common series"], sl)
                raw = func_dict[data_name][meas]["quantiles"]["common series"][sl]
                obs_share = raw[base_jump:end-end_jump] ./ select_series[base_jump:end-end_jump, meas*"_per_hh"]
            end

            # Non-NaN range
            s_axis = xaxis
            est_ids = 1:length(dts)

            # Plot: red = model reconstruction
            Plots.plot()
            Plots.plot!(s_axis,
                model_share[est_ids],
                ylabel=L"\textrm{%$(M)\, \, rel.\,  to\,\,  average}",
                # title=share_plot_titles[si],
                lc=:red,
                xformatter=:latex,
                yformatter=:latex,
                xtickfontsize=14,
                ytickfontsize=14,
                legendfontsize=10,
                guidefontsize=14,
                legend=:best,
                label="",
                lw=4, dpi=500, ls=:solid,
            )

            # Blue dotted = truth
            Plots.plot!(s_axis,
                truth_share[est_ids],
                lc=:blue,
                lw=4,
                la=0.5,
                xtickfontsize=14,
                ytickfontsize=14,
                legendfontsize=10,
                guidefontsize=14,
                dpi=500, ls=:dot,
                label=""
            )

            # Correlation to truth
            valid = findall(i -> !isnan(model_share[i]) && !isnan(truth_share[i]), est_ids)
            corr_n = length(valid) > 2 ? round(cor(model_share[valid], truth_share[valid]) * 100, digits=0) : NaN

            # R² to truth (demeaned)
            if length(valid) > 2
                model_dm = model_share[valid] .- mean(model_share[valid])
                truth_dm = truth_share[valid] .- mean(truth_share[valid])
                mse_s = sum((model_dm .- truth_dm).^2) / length(valid)
                var_s = var(truth_share[valid])
                r2_s = var_s > 0 ? round((1.0 - mse_s / var_s) * 100, digits=0) : NaN
            else
                r2_s = NaN
            end

            # Observed data: black scatter + dashed interpolation
            if obs_share !== nothing
                obs_cond = findall(!isnan, obs_share)
                if !isempty(obs_cond)
                    Plots.scatter!(xaxis[obs_cond],
                        obs_share[obs_cond],
                        marker=:utriangle,
                        markercolor=:black,
                        ms=5,
                        la=0.5,
                        lw=2, dpi=500,
                        label=""
                    )
                    Plots.plot!(xaxis[obs_cond],
                        obs_share[obs_cond],
                        la=0.5,
                        lw=2, dpi=500,
                        label="",
                        ls=:dash
                    )
                end
            end

            Plots.plot!([], [], ls=:dot, lc=:black, la=0.0, label=L"\textrm{Corr.\, to\,\, Truth\,: %$(corr_n)\%}",)
            Plots.plot!([], [], ls=:dot, lc=:black, la=0.0, label=L"R^2\textrm{\, to\,\, Truth\,: %$(r2_s)\%}",)
            Plots.savefig(path * "$meas/" * "quantiles_levels/" * plot_name * "_$meas" * "_share_$(sl)" * detrended_or_not * label * ".pdf")
        end
    end

    # --- Cross-conditional comparison: E[outcome | conditioning ∈ group] ---
    cross_labels = ["bot50", "mid40", "top10"]
    cross_plot_titles = [L"\textrm{Bottom\,\, 50\%}", L"\textrm{Middle\,\, 40\%}", L"\textrm{Top\,\, 10\%}"]
    group_names = ["bottom", "middle", "top"]
    meas_display = Dict("consum" => "Consumption", "income" => "Income", "wealth" => "Wealth")
    detrended_or_not = ty == "average" ? "detrended_" : ""

    for out_meas in measures
        for cond_meas in measures
            key = "$(out_meas)_by_$(cond_meas)"
            is_self = (out_meas == cond_meas)

            # For self-conditionals, use the share-group data under the measure's own name
            if is_self
                haskey(data_dict, out_meas) || continue
                haskey(data_dict[out_meas], "quantiles") || continue
                haskey(data_dict[out_meas]["quantiles"], "common series") || continue
            else
                haskey(data_dict, key) || continue
            end

            # Ensure output directory exists
            mkpath(path * "cross_conditional/")

            corr_dict[key] = Dict()
            if !haskey(corr_dict, "$(key)_r2")
                corr_dict["$(key)_r2"] = Dict()
            end
            for (si, sl) in enumerate(cross_labels)
                # DDFM reconstruction
                if is_self
                    haskey(data_dict[out_meas]["quantiles"]["common series"], sl) || continue
                    model_series = data_dict[out_meas]["quantiles"]["common series"][sl][base_jump:end-end_jump]
                else
                    haskey(data_dict[key], sl) || continue
                    model_series = data_dict[key][sl][base_jump:end-end_jump]
                end

                # Truth: self-conditionals use "consum_bot50", cross use "consum_by_income_bot50"
                truth_col = is_self ? Symbol("$(out_meas)_$(sl)") : Symbol("$(key)_$(sl)")
                if !hasproperty(truth_data_filtered, truth_col)
                    @warn "Truth column $truth_col not found, skipping cross-conditional for $key $sl"
                    continue
                end
                mean_col = Symbol("$(out_meas)_per_hh")
                truth_series = truth_data_filtered[:, truth_col] ./
                               truth_data_filtered[:, mean_col]
                model_normed = model_series ./ select_series[base_jump:end-end_jump, string(mean_col)]

                # Correlation
                valid = findall(i -> !isnan(model_normed[i]) && !isnan(truth_series[i]),
                                1:length(truth_series))
                corr_n = length(valid) > 2 ?
                    round(Int, cor(model_normed[valid], truth_series[valid]) * 100) : NaN
                corr_dict[key][group_names[si]] = corr_n

                # R² (detrended): 1 - MSE / Var(truth)
                if length(valid) > 2
                    model_dm = model_normed[valid] .- mean(model_normed[valid])
                    truth_dm = truth_series[valid] .- mean(truth_series[valid])
                    mse_val = sum((model_dm .- truth_dm).^2) / length(valid)
                    var_truth = var(truth_series[valid])
                    corr_dict["$(key)_r2"][group_names[si]] = var_truth > 0 ?
                        round(Int, (1.0 - mse_val / var_truth) * 100) : NaN
                else
                    corr_dict["$(key)_r2"][group_names[si]] = NaN
                end

                # Plot: DDFM (red) vs Truth (blue dotted)
                s_axis = collect(1:length(dts))
                Out_M = get(meas_display, out_meas, uppercasefirst(out_meas))
                Cond_M = get(meas_display, cond_meas, uppercasefirst(cond_meas))

                Plots.plot()
                Plots.plot!(s_axis,
                    model_normed,
                    ylabel=L"\textrm{%$(Out_M)\, \, rel.\,  to\,\,  average}",
                    lc=:red,
                    xformatter=:latex, yformatter=:latex,
                    xtickfontsize=14, ytickfontsize=14,
                    legendfontsize=10, guidefontsize=14,
                    legend=:best, label="",
                    lw=4, dpi=500, ls=:solid,
                )
                Plots.plot!(s_axis,
                    truth_series,
                    lc=:blue, lw=4, la=0.5,
                    xtickfontsize=14, ytickfontsize=14,
                    legendfontsize=10, guidefontsize=14,
                    dpi=500, ls=:dot, label=""
                )
                Plots.plot!([], [], ls=:dot, lc=:black, la=0.0,
                    label=L"\textrm{Corr.\, to\,\, Truth\,: %$(corr_n)\%}",)
                Plots.savefig(path * "cross_conditional/" * plot_name *
                    "_$(key)_$(sl)_" * detrended_or_not * label * ".pdf")
            end
        end
    end

    # --- Sample cross-conditional correlations ---
    sample_corr_dict = Dict{String, Dict}()
    sample_map = Dict(
        "a" => ("PSID", ["consum", "income", "wealth"]),
        "c" => ("CEX",  ["consum", "income"]),
        "d" => ("SCF",  ["income", "wealth"]),
    )
    for (sample_letter, (sample_name, sample_vars)) in sample_map
        sample_key = "HANK $(sample_letter) $(economy_number)"
        sample_corr_dict[sample_name] = Dict()
        sample_path = DATA_PROCESSING * "/HANK_$(sample_name)_$(economy_number).csv"
        if !isfile(sample_path)
            @warn "Sample file not found: $sample_path"
            continue
        end
        for out_meas in measures
            for cond_meas in measures
                key = "$(out_meas)_by_$(cond_meas)"

                # Check if sample has both variables
                if !(out_meas in sample_vars) || !(cond_meas in sample_vars)
                    sample_corr_dict[sample_name][key] = Dict(
                        gn => "--" for gn in group_names
                    )
                    continue
                end

                sample_df = compute_sample_cross_conditional(
                    sample_path, out_meas, cond_meas; shares = [0.5, 0.4, 0.1]
                )
                if sample_df === nothing
                    sample_corr_dict[sample_name][key] = Dict(
                        gn => "--" for gn in group_names
                    )
                    continue
                end

                # Align sample with truth via linear interpolation
                truth_times = QuarterlyDate.(truth_data_filtered.year, truth_data_filtered.quarter)
                sample_df_filt = filter(row -> row.time in Set(truth_times), sample_df)
                sample_corr_dict[sample_name][key] = Dict()
                if !haskey(sample_corr_dict[sample_name], "$(key)_r2")
                    sample_corr_dict[sample_name]["$(key)_r2"] = Dict()
                end

                is_self_sample = (out_meas == cond_meas)
                for (si, sl) in enumerate(cross_labels)
                    col_name = "$(key)_$(sl)"
                    truth_col = is_self_sample ? Symbol("$(out_meas)_$(sl)") : Symbol(col_name)
                    if !hasproperty(truth_data_filtered, truth_col) ||
                       !(col_name in names(sample_df_filt))
                        sample_corr_dict[sample_name][key][group_names[si]] = "--"
                        sample_corr_dict[sample_name]["$(key)_r2"][group_names[si]] = "--"
                        continue
                    end

                    mean_col = Symbol("$(out_meas)_per_hh")

                    # Sample observed points (sparse)
                    sample_times_raw = sample_df_filt[:, :time]
                    sample_vals_raw  = sample_df_filt[:, Symbol(col_name)]

                    # Remove NaN from sample before interpolation
                    valid_sample = findall(i -> !isnan(sample_vals_raw[i]), 1:length(sample_vals_raw))
                    if length(valid_sample) < 2
                        sample_corr_dict[sample_name][key][group_names[si]] = "--"
                        sample_corr_dict[sample_name]["$(key)_r2"][group_names[si]] = "--"
                        continue
                    end

                    # Interpolate sample onto all truth quarters
                    x_obs = Float64.(Dates.value.(sample_times_raw[valid_sample]))
                    y_obs = Float64.(sample_vals_raw[valid_sample])
                    x_tgt = Float64.(Dates.value.(truth_times))
                    itp = linear_interpolation(x_obs, y_obs, extrapolation_bc=Line())
                    sample_interp = itp.(x_tgt)

                    # Normalize both by mean
                    truth_vals  = truth_data_filtered[:, truth_col] ./ truth_data_filtered[:, mean_col]
                    sample_vals = sample_interp ./ truth_data_filtered[:, mean_col]

                    valid_s = findall(i -> !isnan(truth_vals[i]) && !isnan(sample_vals[i]),
                                     1:length(truth_vals))
                    if length(valid_s) > 2
                        sample_corr_dict[sample_name][key][group_names[si]] =
                            round(Int, cor(truth_vals[valid_s], sample_vals[valid_s]) * 100)

                        # R² (detrended)
                        model_dm = sample_vals[valid_s] .- mean(sample_vals[valid_s])
                        truth_dm = truth_vals[valid_s] .- mean(truth_vals[valid_s])
                        mse_val = sum((model_dm .- truth_dm).^2) / length(valid_s)
                        var_truth = var(truth_vals[valid_s])
                        sample_corr_dict[sample_name]["$(key)_r2"][group_names[si]] = var_truth > 0 ?
                            round(Int, (1.0 - mse_val / var_truth) * 100) : NaN
                    else
                        sample_corr_dict[sample_name][key][group_names[si]] = "--"
                        sample_corr_dict[sample_name]["$(key)_r2"][group_names[si]] = "--"
                    end
                end
            end
        end
    end

    # Export correlations to a JLD2 file (includes cross-conditional + sample correlations)
    file_name = path * "correlations/correlations_$(data_name)"
    jldsave(file_name * ".jld2"; correlations=corr_dict, sample_correlations=sample_corr_dict)

    return within_stat_dict
end


function generate_quantiles_shares_levels_HANK_full(data_dict, ty, func_data, data_name, smin, smax, tmin, tmax, estimator, label, type, measures, time_params, select_series, gdp_series, posterior_bounds, compare_to_other_est, tag, freq)

    @unpack func_dict, year_vec, data_sources, confidence_intervals = func_data
    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid, integral_cop_grid = estimator
    else
        @unpack grid_pcf, grid_cop = estimator
    end

    grid_choice_pcf = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf
    grid_choice_cop = typeof(estimator) <: SeriesEstimator ? integral_cop_grid : grid_cop

    # Subset data for the reconstruction
    base_jump, end_jump = find_subset_frame(smin, smax, tmin, tmax)
    m_label = measures_folder(measures)
    folder = type == :optimization ? "from_optimization" : "from_mcmc"
    init_path = BASE_PATH
    path = init_path * "/7_Results/$m_label" * "$tag" * "/$folder/plots/"

    # Estimation dates
    dts = QuarterlyDate(smin["year"], smin["quarter"]):Quarter(1):QuarterlyDate(smax["year"], smax["quarter"])
    println("dates of user: ", dts)

    # Import the model estimates from the incomplete-HANK model
    estimates_dict = get_estimates_for_comparison_HANK(data_name, ty, time_params, measures, estimator)

    # # Find the integer from tag 
    # economy_number = split(tag, " ")[end]
    # file_name = init_path * "/2_Data_processing/confidence_intervals/ci_draws_illiqd_and_income_and_liquid_quintiles_series_HANK $(economy_number)_all.jld2"
    # raw_ci = jldopen(file_name, "r")["ci"]

    # file_name_raw_data = "/home/luisc/Distributional_Dynamics/7_Results/illiqd_and_income_and_liquid HANK full/other_results/func_dict_HANK_full_$(economy_number).jld2"
    # func_dict_full_hank = jldopen(file_name_raw_data, "r")["data_obs"]

    # # Find the bounds of each 
    # hank_year_vec = repeat([i for i in collect(2001:2024)], inner=4)
    # confidence_intervals = Dict("HANK" => Dict())
    # confidence_intervals["HANK"]["ci_u"], confidence_intervals["HANK"]["ci_l"] = construct_confidence_intervals(raw_ci, 0.025, 0.975, measures, hank_year_vec, estimator)

    # # Correct the confidence intervals to be in the same time frame as the estimation
    # hank_full_dts = QuarterlyDate(2001, 1):Quarter(1):QuarterlyDate(2024, 4)

    # # Find the indices of the estimation time frame in the HANK full time frame
    # hank_indices = findall(x -> x in dts, hank_full_dts)
    # println("HANK indices: ", hank_indices)
    # for meas in measures
    #     for o in ["quantiles"]
    #         println("HANK ci_u size: ", size(confidence_intervals["HANK"]["ci_u"][meas][o]))
    #         println("HANK full size: ", size(func_dict_full_hank["HANK full"][meas][o]["data"]))
    #         confidence_intervals["HANK"]["ci_u"][meas][o] = confidence_intervals["HANK"]["ci_u"][meas][o][:, hank_indices]
    #         confidence_intervals["HANK"]["ci_l"][meas][o] = confidence_intervals["HANK"]["ci_l"][meas][o][:, hank_indices]
    #         func_dict_full_hank["HANK full"][meas][o]["data"] = func_dict_full_hank["HANK full"][meas][o]["data"][:, hank_indices]
    #     end
    # end


    local bot, mid, top
    if grid_choice_pcf == 10 || grid_choice_pcf == 100 || grid_choice_pcf == 20
        bot, mid, top = "bottom50", "next40", "top10"
    elseif grid_choice_pcf == 5
        bot, mid, top = "bottom40", "next40", "top20"
    end

    # Get the observed measures
    obs_meas = get_obs_meas(func_dict, data_name, measures; top=top)

    within_stat_dict = Dict()
    plot_name = data_name

    # Generate the within statistic for the copulas 
    dimension = length(measures)

    # if data_name != "consensus" && ty == "normal" # ty == "average" doesn't work because it has no comparison to the data
    #     within_stat_dict["copula"] = compute_copula_within_stat(data_dict, confidence_intervals, base_jump, end_jump, data_name, dimension, grid_choice_cop)
    # end

    println("select_series size: ", size(select_series))
    println("select_series content: ", select_series[base_jump:end-end_jump, :])
    println("base_jump: ", base_jump, " end_jump: ", end_jump)

    # Create a dictionary to hold the correlations
    corr_dict = Dict()
    within_stat_dict = Dict()

    for df_k in keys(estimates_dict)
        corr_dict[df_k] = Dict()
        within_stat_dict[df_k] = Dict()
        for meas in obs_meas # TODO: this issue here is that not all measures are observed ... ofc, we can use the reconstructed data but not the confidence intervals

            # All quantiles 
            qu = data_dict[meas]["quantiles"]["data"][:, base_jump:end-end_jump] ./ select_series[base_jump:end-end_jump, meas*"_per_hh"]'

            # For the (o)bservations 
            qu_o = Vector{Any}(undef, 1)

            # # All quantiles 
            # qu_o[1] = func_dict[data_name][meas]["quantiles"]["data"][:, base_jump:end-end_jump] ./ select_series[base_jump:end-end_jump, meas*"_per_hh"]'
            # qu_o[2] = confidence_intervals[data_name]["ci_l"][meas]["quantiles"][:, base_jump:end-end_jump] ./ select_series[base_jump:end-end_jump, meas*"_per_hh"]'
            # qu_o[3] = confidence_intervals[data_name]["ci_u"][meas]["quantiles"][:, base_jump:end-end_jump] ./ select_series[base_jump:end-end_jump, meas*"_per_hh"]'

            # Quantiles from incomplete HANK model
            # smin_full = Dict("year" => 2001, "quarter" => 1)
            # smax_full = Dict("year" => 2024, "quarter" => 4)
            # base_jump_full, end_jump_full = find_subset_frame(smin, smax, smin_full, smax_full)
            # println("HANK full jumps: ", base_jump_full, " ", end_jump_full)
            qu_o[1] = estimates_dict[df_k][meas]["quantiles"]["data"] ./ select_series[base_jump:end-end_jump, meas*"_per_hh"]'
            if all(isnan.(qu_o[1]))
                println("All NaN for df_k: ", df_k, " meas: ", meas)
                continue
            end

            # Plots 
            markers = [:diamond, :utriangle, :circle, :star5, :cross, :diamond, :utriangle, :circle, :star5, :cross, :diamond, :utriangle, :circle, :star5, :cross, :diamond, :utriangle, :circle, :star5, :cross] #collect(range(0.3, 1, length=length(data_sources)))

            # x-axis
            xaxis = collect(1:length(dts))
            M = uppercasefirst(meas)

            # Important tag 
            detrended_or_not = ty == "average" ? "detrended_" : ""

            # All possible line-styles
            line_styles = [:solid :dash :dot :dashdot :dashdotdot :solid :dash :dot :dashdot :dashdotdot] # has to be a matrix

            within_stat_dict[df_k][meas] = Dict()

            # Doing quantiles now 
            label_quantiles = label_qs(grid_choice_pcf) # e.g., L"\textrm{%$(i)th}"
            sequences = define_sequences(grid_choice_pcf) # gets indices of each group
            dist_dict = Dict("bottom" => [sequences[1]], "middle" => [sequences[2]], "top" => [sequences[3]])


            for (obj, dist) in dist_dict
                within_stat_dict[df_k][meas][obj] = Dict()

                cond = data_name != "consensus" ? findall(!isnan, qu_o[1][dist[1][1], :]) : [1] # .!isnan.(all_lv_o[1][:, j])

                est_ids = cond[1]:cond[end]
                s_axis = xaxis[est_ids] # start at the first observation
                s_data = qu[:, est_ids]
                s_dts = dts[est_ids]
                intd = 40

                # Plotting model estimates
                for (lsⱼ, j) in enumerate(dist[1])
                    Plots.plot()
                    Plots.plot!(s_axis,
                        s_data[j, :],
                        ylabel=L"\textrm{%$(M)\, \, rel.\,  to\,\,  average}",
                        lc=:red,
                        xformatter=:latex,
                        yformatter=:latex,
                        xtickfontsize=14,
                        ytickfontsize=14,
                        legendfontsize=10,
                        guidefontsize=14,
                        # xticks=(s_axis[1:intd:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(s_dts[1:intd:end])]),
                        legend=:best,
                        label="",
                        lw=4, dpi=500, ls=:solid,
                    )

                    Plots.plot!(s_axis,
                        qu_o[1][j, :][est_ids],
                        lc=:blue,
                        lw=4,
                        la=0.5,
                        xtickfontsize=14,
                        ytickfontsize=14,
                        legendfontsize=10,
                        guidefontsize=14,
                        dpi=500, ls=:dot,
                        label=""
                    )

                    # Compute correlation between full model and incomplete HANK model
                    corr_n = round(cor(s_data[j, :], qu_o[1][j, :][est_ids]), digits=2) * 100

                    # c_data = Vector{Any}(undef, 3)
                    # c_data[1] = qu_o[1][j, :][cond]
                    # c_data[2] = qu_o[2][j, est_ids]
                    # c_data[3] = qu_o[3][j, est_ids]

                    # # See how many points fall within the confidence intervals
                    # r_data = qu[j, est_ids] # estimates that correspond to the indices of the data points 
                    # num = count(c_data[2] .<= r_data .<= c_data[3])
                    # den = length(r_data)
                    # within_stat = floor(Int, (num ./ den) * 100)

                    # within_stat_dict[meas][obj][label_quantiles[:, j][1]] = "$num" * "/" * "$den"

                    # Plots.scatter!(xaxis[cond],
                    #     c_data[1],
                    #     marker=:utriangle, #markers[j],
                    #     markercolor=:black,
                    #     ms=5,
                    #     la=0.5,
                    #     lw=2, dpi=500,
                    #     label=""
                    #     # label=L"\textrm{HANK \,\,Data\,\, Used}",
                    # )
                    # Store the correlation in the dictionary
                    corr_dict[df_k]["$meas"*"_"*"$obj"*"_"*"$(label_quantiles[:, j][1])"] = corr_n

                    # Plots.plot!([], [], ls=:dash, lc=:black, la=0.0, label=L"\textrm{Within\, \,  bounds: %$(within_stat)\%}",)
                    Plots.plot!([], [], ls=:dash, lc=:black, la=0.0, label=L"\textrm{Corr.: %$(corr_n)\%}",)
                    Plots.savefig(path * "$meas/" * "quantiles_levels/" * df_k * "_$meas" * "_$obj" * "_quantiles_$(j)" * detrended_or_not * label * ".pdf")
                end
            end
        end
    end

    # Export correlations to a JLD2 file
    economy_number = split(tag, " ")[end]
    file_name = path * "correlations/correlations_HANK_$(economy_number)"
    jldsave(file_name * ".jld2"; correlations=corr_dict)

    return within_stat_dict
end


function get_estimates_for_comparison_HANK(data_name, ty, time_p, measures, estimator)
    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid, integral_cop_grid = estimator
    else
        @unpack grid_pcf, grid_cop = estimator
    end

    grid_choice_pcf = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf
    grid_choice_cop = typeof(estimator) <: SeriesEstimator ? integral_cop_grid : grid_cop

    plot_name = data_name

    @unpack year_vec = time_p
    init_path = BASE_PATH
    objects = sort(["quantiles", "levels", "shares"])
    meas_folder = measures_folder(measures)

    # Download estimates from incomplete HANK models 
    economy_number = split(tag, " ")[end]
    data_tag = ty == "normal" ? "" : "_detrended"

    # Actual estimates
    inc_tag = "HANK $(economy_number)" # incomplete HANK model

    estimates_dict = Dict()

    for df in ["a", "b", "c", "d"]
        local estimates
        if economy_number == 1
            estimates = CSV.read(init_path * "/7_Results/$meas_folder " * inc_tag * "/from_mcmc/data/HANK $(df)_functional_data" * data_tag * "_A non-diag_.csv", DataFrame)
        else
            estimates = CSV.read(init_path * "/7_Results/$meas_folder " * inc_tag * "/from_mcmc/data/HANK $(df) $(economy_number)_functional_data" * data_tag * "_A non-diag_.csv", DataFrame)
        end

        # Decompose this into copulas, pcfs, levels and shares
        T = nrow(estimates)
        D = length(measures)
        colons = ntuple(_ -> (:), D) # TODO: needs to be corrected for 4 dimensional case 
        cop_size = tuple([grid_choice_cop for i in 1:D]...)

        cop_mat_size = tuple([grid_choice_cop for i in 1:D]..., T)

        copulas = fill!(Array{Float64}(undef, cop_mat_size...), NaN)

        # Reshape objects
        data_pcf = Matrix(transpose(Matrix(select(estimates, filter(x -> occursin("quantiles", x), names(estimates))))))
        levels = Matrix(transpose(Matrix(select(estimates, filter(x -> occursin("levels", x), names(estimates))))))
        shares = Matrix(transpose(Matrix(select(estimates, filter(x -> occursin("shares", x), names(estimates))))))

        # Putting it all together 
        estimates_dict2 = create_time_series_dictionary([copulas, data_pcf, levels, shares], estimator, sort(measures))

        @unpack tmin, tmax = time_p

        # Filter cex_all_data to match the overall estimation time frame. cex_all_data is already in quarters and is from 1962Q3 to 2021Q4. cex_all_estimates as well.
        min_time = minimum(estimates[!, "time"]) # the first period of the estimation of the outside estimates 
        max_time = maximum(estimates[!, "time"]) # the last period of the estimation of the outside estimates

        pushup = QuarterlyDate(tmin["year"], tmin["quarter"]) - QuarterlyDate(min_time) # Tricky part is if min_time is later than tmin ... I guess i could add NaNs for missing periods. Makes sense!
        pushdown = QuarterlyDate(max_time) - QuarterlyDate(tmax["year"], tmax["quarter"]) # max_time solved similarly

        estimates_dict["HANK $df"] = fill_data_dicts(estimates_dict2, measures, pushup, pushdown, estimator)
    end

    return estimates_dict
end


"""
    generate_averaged_hank_avg_corr_table(economy_ids, measures, tag, type; output_path)

Load cross-conditional correlations from multiple HANK economy runs, average them,
and generate LaTeX table `hank_avg_corr.tex`.

`economy_ids` is e.g. `2:10` — the economies to average over.
"""
function generate_averaged_hank_avg_corr_table(
    economy_ids,
    measures::Vector{String},
    type::String;          # "from_mcmc" or "from_optimization"
    data_sources::Union{String, Vector{String}} = "a",   # single or multiple: ["a", "c", "d"]
    output_path::String = ""
)
    # Normalize to vector
    data_sources = data_sources isa String ? [data_sources] : data_sources
    m_label = measures_folder(measures)
    group_names = ["bottom", "middle", "top"]
    group_labels = ["Bottom \$<50\$", "Middle \$50\\text{-}90\$", "Top \$>90\$"]
    cond_vars = ["consum", "income", "wealth"]
    out_vars = ["consum", "income", "wealth"]
    section_titles = ["By Consumption Groups", "By Income Groups", "By Wealth Groups"]

    # Collect all correlations
    all_ddfm = Dict{String, Dict{String, Vector{Float64}}}()
    all_sample = Dict{String, Dict{String, Dict{String, Vector{Float64}}}}()

    for econ_id in economy_ids
        econ_tag = " HANK $(econ_id)"
        init_path = BASE_PATH
        corr_path = init_path * "/7_Results/$(m_label)$(econ_tag)/$type/plots/correlations/"

        for ds in data_sources
            data_name = "HANK $(ds) $(econ_id)"
            corr_file = corr_path * "correlations_$(data_name).jld2"
            if !isfile(corr_file)
                @warn "Correlations file not found: $corr_file"
                continue
            end

            jld = jldopen(corr_file, "r")
            corr_dict = jld["correlations"]
            sample_dict = haskey(jld, "sample_correlations") ? jld["sample_correlations"] : Dict()
            close(jld)

            # Accumulate DDFM correlations (same across sources, but averaging handles duplicates)
            for out_var in out_vars
                for cond_var in cond_vars
                    key = "$(out_var)_by_$(cond_var)"
                    haskey(corr_dict, key) || continue
                    if !haskey(all_ddfm, key)
                        all_ddfm[key] = Dict(gn => Float64[] for gn in group_names)
                    end
                    for gn in group_names
                        if haskey(corr_dict[key], gn)
                            v = corr_dict[key][gn]
                            if v isa Number && !isnan(v)
                                push!(all_ddfm[key][gn], Float64(v))
                            end
                        end
                    end
                end
            end

            # Accumulate sample correlations
            for (sname, sdict) in sample_dict
                if !haskey(all_sample, sname)
                    all_sample[sname] = Dict{String, Dict{String, Vector{Float64}}}()
                end
                for (key, gdict) in sdict
                    if !haskey(all_sample[sname], key)
                        all_sample[sname][key] = Dict(gn => Float64[] for gn in group_names)
                    end
                    for gn in group_names
                        if haskey(gdict, gn)
                            v = gdict[gn]
                            if v isa Number && !isnan(v)
                                push!(all_sample[sname][key][gn], Float64(v))
                            end
                        end
                    end
                end
            end
        end
    end

    # Build table with averaged correlations
    _fmt(vals) = isempty(vals) ? "--" : string(round(Int, mean(vals)))

    # Map data_sources to column headers: data_sources[i] → measures[i] column
    source_label = Dict("a" => "\\textit{A}", "b" => "\\textit{B}",
                        "c" => "\\textit{C}", "d" => "\\textit{D}")
    # Build per-measure sample label and survey name lookup
    measure_sample_label = Dict{String,String}()
    measure_sample_name  = Dict{String,String}()
    survey_name = Dict("a" => "PSID", "b" => "CPS", "c" => "CEX", "d" => "SCF")
    for (i, m) in enumerate(sort(measures))
        ds = i <= length(data_sources) ? data_sources[i] : data_sources[end]
        measure_sample_label[m] = get(source_label, ds, "Sample")
        measure_sample_name[m]  = get(survey_name, ds, "Sample")
    end

    sm = sort(measures)
    table = """\\begin{tabular}{l cc cc cc}
\\toprule
 & \\multicolumn{2}{c}{Consumption} & \\multicolumn{2}{c}{Income} & \\multicolumn{2}{c}{Wealth} \\\\
\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}
 & DDFM & Sample $(measure_sample_label[sm[1]]) & DDFM & Sample $(measure_sample_label[sm[2]]) & DDFM & Sample $(measure_sample_label[sm[3]]) \\\\
\\midrule
"""

    for (sec_idx, cond_var) in enumerate(cond_vars)
        table *= "\\multicolumn{7}{l}{\\textit{$(section_titles[sec_idx])}} \\\\\n"

        for (gi, gn) in enumerate(group_names)
            row = group_labels[gi]
            for out_var in out_vars
                key = "$(out_var)_by_$(cond_var)"

                # Averaged DDFM correlation
                ddfm_val = haskey(all_ddfm, key) ? _fmt(all_ddfm[key][gn]) : "--"

                # Averaged sample correlation — use the survey assigned to the outcome column
                sample_val = "--"
                sname = measure_sample_name[out_var]
                if haskey(all_sample, sname)
                    sdata = all_sample[sname]
                    if haskey(sdata, key) && haskey(sdata[key], gn)
                        sample_val = _fmt(sdata[key][gn])
                    end
                end

                row *= " & $ddfm_val & $sample_val"
            end
            table *= "$row \\\\\n"
        end
        if sec_idx < length(cond_vars)
            table *= "\\addlinespace\n"
        end
    end

    table *= """\\bottomrule
\\end{tabular}
"""

    # Determine output path
    if isempty(output_path)
        tag = " HANK $(first(economy_ids))-$(last(economy_ids))"
        output_path = BASE_PATH * "/7_Results/$(m_label)$(tag)/plots/correlations/"
    end
    mkpath(output_path)
    outfile = joinpath(output_path, "hank_avg_corr.tex")
    open(outfile, "w") do io
        write(io, table)
    end
    @info "Wrote averaged hank_avg_corr table ($(length(economy_ids)) economies) to $outfile"
end


"""
Pick the most informative sample for a given (conditioning_var, outcome_var) cell.

Rules:
- CEX has consum + income (no wealth)
- SCF has income + wealth (no consumption)
- PSID has all three, but is less frequent → use as fallback

For "By Consumption Groups" (cond=consum): outcome=income → CEX; outcome=wealth → PSID
For "By Income Groups"     (cond=income): outcome=consum → CEX; outcome=wealth → SCF
For "By Wealth Groups"     (cond=wealth): outcome=consum → PSID; outcome=income → SCF
"""
function _pick_sample_for_cell(cond_var::String, out_var::String)
    if cond_var == "consum"
        out_var == "income" && return "CEX"
        out_var == "wealth" && return "PSID"
    elseif cond_var == "income"
        out_var == "consum" && return "CEX"
        out_var == "wealth" && return "SCF"
    elseif cond_var == "wealth"
        out_var == "consum" && return "PSID"
        out_var == "income" && return "SCF"
    end
    return nothing
end
