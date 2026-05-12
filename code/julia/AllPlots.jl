function logspace(start, stop, num; base=10.0)
    logstart, logstop = log(base, start), log(base, stop)
    logspace_vals = exp.(LinRange(logstart, logstop, num))
    return logspace_vals
end

function scale_transformation(series_type, series)
    if series_type == "quantiles"
        return series ./ 10000
    elseif series_type == "levels"
        return series ./ 1_000_000_000_000
    end
end

function symlog(x, n=-3)
    result = zeros(size(x, 1))
    for i in axes(x, 1)
        result[i] = sign(x[i]) * (log10(1 + abs(x[i]) / (10^n)))
    end
    result
end

function symlogformatter(x, n)
    if sign(x) == 0
        L"10^{%$(Int(n))}"
    else
        s = sign(x) == 1 ? "+" : "-"
        nexp = sign(x) * (abs(x) + n)

        if sign(x) == -1
            nexp = -nexp
        end

        if mod(nexp, 1) == 0
            nexp = Int(nexp)
        end

        L"%$(s)10^{%$(nexp)}"
    end
end

function gen_proof_of_concept_copulas(d_data_dict::Dict, r_data_dict::Dict, source_of_id, measures, model_options, time_p, fully_observed_id, estimation_id)
    """Plotting copula visuals"""
    @unpack equivalized, bottom_coded, estimator = model_options
    @unpack time_dict, tmin, tmax = time_p
    @unpack grid_cop, integral_cop_grid, integral_pcf_grid = estimator

    # Get some observations 
    random_measure = first(measures)
    all_periods = axes(d_data_dict[random_measure]["quantiles"]["data"], 2)
    fully_observed_periods = all_periods[fully_observed_id]
    interval = Int(floor(length(fully_observed_periods) / 3))
    the_periods_id = fully_observed_periods[1:interval:end]

    # println(all_periods)
    # println(fully_observed_id)
    # println(estimation_id)

    # Get the IDs 'the_years' we want 
    dts = QuarterlyDate(tmin["year"], tmin["quarter"]):Quarter(1):QuarterlyDate(tmax["year"], tmax["quarter"])
    dts_to_keep = dts[estimation_id[1:interval:end]]

    # Preliminaries 
    D = length(measures)
    axis = select_grid_points(integral_cop_grid)
    colons = ntuple(_ -> (:), D)

    # Path stuff
    init_path = BASE_PATH
    m_label = measures_folder(measures)
    path = init_path * "/7_Results/proof_of_concept/$m_label/"
    mkpath(path)

    # Generate side by side copulas
    grid_tup = tuple([integral_cop_grid for _ in 1:D]...)
    count = D * length(collect(the_periods_id))
    legendfl = [false for _ in 1:count]
    legendfl[1] = true
    j = 0

    # Find indices of the observations from d_data_dict["copulas"]["data"][colons..., :] where the column is not all NaN
    # Find columns that do not consist entirely of NaNs
    non_nan_columns = [col for col in 1:size(d_data_dict["copulas"]["data"], D + 1) if all(!isnan, d_data_dict["copulas"]["data"][colons..., col])]

    # println(d_data_dict["copulas"]["data"][colons..., non_nan_columns])
    # println(non_nan_columns)

    # First generate time series of KL divergence 
    # Identify indices that are before covid and then 2010s 
    b4_cov_T_ids = [est_id for est_id in estimation_id if dts[est_id] <= QuarterlyDate(2019, 2)]
    the2010s_T_ids = [est_id for est_id in estimation_id if dts[est_id] >= QuarterlyDate(2010, 1)]

    b4_cov_ids = [id for (i, id) in enumerate(fully_observed_periods) if dts[estimation_id[i]] <= QuarterlyDate(2019, 2)]
    the2010s_ids = [id for (i, id) in enumerate(fully_observed_periods) if dts[estimation_id[i]] >= QuarterlyDate(2010, 1)]

    # b4_cov_ids   = [i for (i, id) in enumerate(fully_observed_periods) if dts[estimation_id[i]] <= QuarterlyDate(2019, 2)]
    # the2010s_ids = [i for (i, id) in enumerate(fully_observed_periods) if dts[estimation_id[i]] >= QuarterlyDate(2010, 1)]


    C̄_nocov = reshape(nanmean(d_data_dict["copulas"]["data"][colons..., b4_cov_ids], dims=D + 1), grid_tup)
    C̄_2010s = reshape(nanmean(d_data_dict["copulas"]["data"][colons..., the2010s_ids], dims=D + 1), grid_tup)

    KL_nocov = Vector{Float64}(undef, length(b4_cov_ids))
    KL_nocov_o = Vector{Float64}(undef, length(b4_cov_ids))
    KL_2010s = Vector{Float64}(undef, length(the2010s_ids))
    KL_2010s_o = Vector{Float64}(undef, length(the2010s_ids))

    #TODO: COMMENT -> d_data_dict has only fully observed observations, so 'y' should be a number from 1:length(observed)

    for (i, y) in enumerate(b4_cov_ids)
        KL_nocov[i] = compute_Kullback_Leibler_divergence(d_data_dict["copulas"]["data"][colons..., y], r_data_dict["copulas"]["data"][colons..., y])
        KL_nocov_o[i] = compute_Kullback_Leibler_divergence(d_data_dict["copulas"]["data"][colons..., y], C̄_nocov; compare_to_average=true)
    end

    for (i, y) in enumerate(the2010s_ids)
        KL_2010s[i] = compute_Kullback_Leibler_divergence(d_data_dict["copulas"]["data"][colons..., y], r_data_dict["copulas"]["data"][colons..., y])
        KL_2010s_o[i] = compute_Kullback_Leibler_divergence(d_data_dict["copulas"]["data"][colons..., y], C̄_2010s; compare_to_average=true)
    end

    # Plot this as a time series 
    select_annual = source_of_id == "CEX" ? 4 : 1
    select_dts = source_of_id == "CEX" ? 12 : 1
    intval = source_of_id == "CEX" ? select_annual + select_dts : source_of_id == "SCF" ? 2 : 1
    xaxis = collect(1:1:length(b4_cov_ids))
    sxaxis = collect(1:1:length(the2010s_ids))

    Plots.plot(
        xaxis[1:select_annual:end],
        KL_nocov[1:select_annual:end],
        # xlabel = L"\textrm{Years}",
        ylabel=L"\textrm{Kullback–Leibler\,\, divergence}",
        xformatter=:latex,
        yformatter=:latex,
        xticks=(xaxis[1:intval:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(dts[b4_cov_T_ids][1:intval:end])]),
        linewidth=4,
        xtickfontsize=14,
        ytickfontsize=14,
        legendfontsize=12,
        guidefontsize=14,
        legend=length(measures) == 2 ? false : :right,
        lc=:red,
        ls=:solid,
        label=L"\textrm{Information\,\,Loss}",
    )

    # Save plot

    Plots.plot!(
        xaxis[1:select_annual:end],
        KL_nocov_o[1:select_annual:end],
        # xlabel = L"\textrm{Years}",
        linewidth=4,
        lc=:black,
        ls=:dash,
        label=L"\textrm{Data\,\,Variation}",
    )
    # println(KL_nocov_o)
    # println(KL_nocov)

    Plots.savefig(path * "KL_divergence_no_covid" * ".pdf")


    # Plot this as a time series 
    Plots.plot(
        sxaxis[1:select_annual:end],
        KL_2010s[1:select_annual:end],
        # xlabel = L"\textrm{Years}",
        ylabel=L"\textrm{Kullback–Leibler\,\, divergence}",
        xformatter=:latex,
        yformatter=:latex,
        xticks=(sxaxis[1:intval:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(dts[the2010s_T_ids][1:intval:end])]),
        linewidth=4,
        xtickfontsize=14,
        ytickfontsize=14,
        legendfontsize=12,
        guidefontsize=14,
        legend=false,
        lc=:red,
        ls=:solid,
        label=L"\textrm{Model\,\,Variation}",
    )

    # Save plot

    Plots.plot!(
        sxaxis[1:select_annual:end],
        KL_2010s_o[1:select_annual:end],
        # xlabel = L"\textrm{Years}",
        linewidth=4,
        lc=:black,
        ls=:dash,
        label=L"\textrm{Data\,\,Variation}",
    )

    Plots.savefig(path * "KL_divergence_2010s" * ".pdf")

    measure1 = uppercasefirst(measures[1])
    measure2 = uppercasefirst(measures[2])
    m1 = uppercasefirst(measures[1])[1]
    m2 = uppercasefirst(measures[2])[1]

    for (i, y) in enumerate(the_periods_id)
        if length(measures) == 2
            Plots.surface(
                axis,
                axis,
                d_data_dict["copulas"]["data"][colons..., y],  # copula 
                xlabel=L"\textrm{%$(measure1)}",
                ylabel=L"\textrm{%$(measure2)}",
                zlabel=L"dC(%$(m1), %$(m2))",
                xformatter=:latex,
                yformatter=:latex,
                zformatter=:latex,
                xtickfontsize=10,
                ytickfontsize=10,
                legendfontsize=10,
                guidefontsize=10,
                legend=false,
                camera=(30, 10),
                size=(400, 400),
                color=:winter,
                display_option=Plots.GR.OPTION_SHADED_MESH)
            Plots.pdf(path * "copula_proof_" * "$(dts_to_keep[i])" * ".pdf")

            Plots.surface(
                axis,
                axis,
                r_data_dict["copulas"]["data"][colons..., y],  # copula 
                xlabel=L"\textrm{%$(measure1)}",
                ylabel=L"\textrm{%$(measure2)}",
                zlabel=L"dC(%$(m1), %$(m2))",
                xformatter=:latex,
                yformatter=:latex,
                zformatter=:latex,
                legend=false,
                xtickfontsize=10,
                ytickfontsize=10,
                legendfontsize=10,
                guidefontsize=10,
                camera=(30, 10),
                color=:winter,
                size=(400, 400),
                display_option=Plots.GR.OPTION_SHADED_MESH)

            Plots.pdf(path * "copula_proof_" * "$(dts_to_keep[i])" * "_Approximated" * ".pdf")
        end
        degree = -2.0 #TODO: interpretation -> 
        for m in measures
            j += 1
            M = uppercasefirst(m)
            Plots.plot(
                axis,
                symlog(reshape(scale_transformation("quantiles", d_data_dict[m]["quantiles"]["data"][:, y]), (integral_pcf_grid, 1)), degree),
                xlabel=L"\textrm{Percentile\,\, Grid}",
                ylabel=L"\textrm{%$(M)\,\,\,\, (100k\, USD,\, 2019\, base\,\,  year)}",
                xformatter=:latex,
                # yformatter=:latex,
                yformatter=x -> symlogformatter(x, degree),
                xticks=(1:integral_pcf_grid, [L"\textrm{%$(i)th}" for i in (100/integral_pcf_grid):(100/integral_pcf_grid):100]),
                linewidth=3,
                xtickfontsize=14,
                ytickfontsize=14,
                legendfontsize=12,
                guidefontsize=14,
                legend=legendfl[j],
                ylims=define_ylims(symlog(scale_transformation("quantiles", d_data_dict[m]["quantiles"]["data"][:, end]), degree)),
                label=L"\textrm{Observed}",
                # yscale=:log10,
                size=(400, 400),
            )

            Plots.plot!(
                axis,
                symlog(reshape(scale_transformation("quantiles", r_data_dict[m]["quantiles"]["data"][:, y]), (integral_pcf_grid, 1)), degree),
                yformatter=x -> symlogformatter(x, degree),
                linewidth=3,
                linestyle=:dash,
                linecolor=:orange,
                label=L"\textrm{Approximated}",
                # yscale=:log10,
                size=(400, 400)
            )

            Plots.pdf(path * "pcfs_proof_" * m * "_$(dts_to_keep[i])" * ".pdf")
        end
    end
end

function define_ylims(series)
    yMax = maximum(series) .+ 0.1 * maximum(series)
    yMin = minimum(series) .- 0.1 * minimum(series)
    return (yMin, yMax)
end

function find_plot_ybounds(d_data_dict, measures)

    yMax = Dict()
    yMin = Dict()
    for m in measures
        yMax[m] = Dict()
        yMin[m] = Dict()
        for k in collect(keys(d_data_dict[m]))
            yMax[m][k] = maximum(d_data_dict[m][k]["data"]) .+ 0.1 * maximum(d_data_dict[m][k]["data"])
            yMin[m][k] = minimum(d_data_dict[m][k]["data"]) .- 0.1 * minimum(d_data_dict[m][k]["data"])
        end
    end

    return yMax, yMin
end

# Bootstrap of the reconstructions vs. observed 
function gen_proof_of_concept_figure(d_data_dict::Dict, r_data_dict::Dict, model_options, time_p, estimation_id, years, source)
    @unpack measures, equivalized, bottom_coded, blind_to, estimator = model_options
    @unpack tmin, tmax = time_p

    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid = estimator
    else
        @unpack grid_pcf, grid_cop = estimator
    end
    grid_choice = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf

    # ci_u    = confidence_intervals["ci_u"]
    # ci_l    = confidence_intervals["ci_l"]

    # Plots reflect same scale (ex-post)
    @unpack time_dict, tmin, tmax = time_p

    # Get some observations 
    random_measure = first(measures)
    all_periods = axes(d_data_dict[random_measure]["quantiles"]["data"], 2)

    xaxis = collect(1:1:length(all_periods))
    objects = sort(["quantiles", "levels", "shares"])
    dts = QuarterlyDate(tmin["year"], tmin["quarter"]):Quarter(1):QuarterlyDate(tmax["year"], tmax["quarter"])
    dts_to_keep = dts[estimation_id]
    fully_obs_id = estimation_id #.- estimation_id[1] .+ 1

    # Define path to save figures
    init_path = BASE_PATH
    m_label = measures_folder(measures)
    path = init_path * "/7_Results/proof_of_concept/$m_label/"
    mkpath(path)

    # Generate functions so their out of the inner scope defined by the if-else control flow
    local series, series_l
    if grid_choice == 10 || grid_choice == 20 || grid_choice == 100
        series = ["bottom50", "next40", "top10"]
        series_l = [L"\textrm{Bottom \,50}", L"\textrm{Next\,40}", L"\textrm{Top\,10}"]
    elseif grid_choice == 5
        series = ["bottom40", "next40", "top20"]
        series_l = [L"\textrm{Bottom \,40}", L"\textrm{Next\,40}", L"\textrm{Top\,20}"]
    end

    color_tag = [:green, :sienna, :blue]

    # Generate vector of objects depending on the grid size 
    interval_dts = measures == ["consum", "income"] ? 7 : 4
    for (n, m) in enumerate(measures)
        M = uppercasefirst(m)
        for (c, o) in enumerate(objects)
            if o == "shares"
                Plots.plot()
                for i in eachindex(series)
                    dataₒ = d_data_dict[m][o]["common series"][series[i]]
                    cond = vec(.!any(isnan.(dataₒ), dims=2))
                    select_annual = source == "CEX" && M == "Consum" ? 4 : 1

                    # Condition
                    cond = findall(cond)
                    cond = cond[1:select_annual:end]        # same downsampling as sxaxis


                    sxaxis = xaxis[cond]
                    sdata = dataₒ[cond]
                    # sci_l  = ci_l[m][o][i, :][cond]
                    # sci_u  = ci_u[m][o][i, :][cond]

                    Plots.plot!(
                        sxaxis,
                        sdata, # d_lines[id[1], :, n], 
                        xlabel=L"\textrm{Years}",
                        ylabel=L"\textrm{%$(M)}\,\,\textrm{%$(o)}",
                        xformatter=:latex,
                        yformatter=:latex,
                        xticks=(sxaxis[1:interval_dts:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(dts_to_keep[cond][1:interval_dts:end])]),
                        # xticks=(xaxis[1:8:end], [L"%$(y)" for y in years[1:8:end]]),
                        legend=:best,
                        label=series_l[i],
                        lw=2,
                        xtickfontsize=14,
                        ytickfontsize=14,
                        legendfontsize=12,
                        guidefontsize=14,
                        fillalpha=0.3,
                        lc=color_tag[i],
                        dpi=500
                    )

                    Plots.plot!(
                        sxaxis,
                        sdata,
                        label="",
                        lw=2,
                        lc=color_tag[i],
                        dpi=500,
                        seriesalpha=0.3,
                        # yerror=(sdata - sci_l, sci_u - sdata)
                    )



                    Plots.plot!(
                        sxaxis,
                        hcat(r_data_dict[m][o]["common series"][series[i]][cond], r_data_dict[m][o]["common series"][series[i]][cond], r_data_dict[m][o]["common series"][series[i]][cond]), #r_lines[id, :, n]', 
                        label=["" "" ""],
                        linestyle=:dash,
                        linecolor=:goldenrod1,
                        dpi=500,
                        lw=2
                    )
                end

                Plots.pdf(path * m * "_" * o * "_proof" * ".pdf")

            elseif o == "levels"
                Plots.plot()
                for i in eachindex(series)
                    dataₒ = scale_transformation(o, d_data_dict[m][o]["common series"][series[i]])
                    cond = vec(.!any(isnan.(dataₒ), dims=2))
                    select_annual = source == "CEX" && M == "Consum" ? 4 : 1

                    # Condition
                    cond = findall(cond)
                    cond = cond[1:select_annual:end]        # same downsampling as sxaxis

                    sxaxis = xaxis[cond]
                    sdata = dataₒ[cond]
                    # sci_l  = scale_transformation(o, ci_l[m][o][i, :][cond])
                    # sci_u  = scale_transformation(o, ci_u[m][o][i, :][cond])


                    Plots.plot!(
                        sxaxis,
                        sdata, # d_lines[id[1], :, n], 
                        xlabel=L"\textrm{Years}",
                        ylabel=L"\textrm{%$(M)}\,\,\textrm{%$(o)\,\,\,\, (1T\, USD,\, 2019\, base\,\,  year)}",
                        xformatter=:latex,
                        yformatter=:latex,
                        xticks=(sxaxis[1:interval_dts:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(dts_to_keep[cond][1:interval_dts:end])]),
                        legend=:outertopright,
                        lc=color_tag[i],
                        label=series_l[i],
                        lw=2,
                        xtickfontsize=14,
                        ytickfontsize=14,
                        legendfontsize=12,
                        guidefontsize=14,
                    )

                    Plots.plot!(
                        sxaxis,
                        sdata,
                        lc=color_tag[i],
                        label="",
                        lw=2,
                        seriesalpha=0.3,
                        dpi=500,
                        # yerror=(sdata - sci_l, sci_u - sdata)
                    )
                    # The reconstruction 
                    Plots.plot!(
                        sxaxis,
                        scale_transformation(o, r_data_dict[m][o]["common series"][series[i]])[cond],
                        linestyle=:dash,
                        linecolor=:goldenrod1,
                        label=L"\textrm{Approximated}",
                        lw=2,
                        dpi=500,
                    )

                end
                Plots.pdf(path * m * "_" * o * ".pdf")
            else
                sequences = define_sequences(grid_choice) # sequences = [1, 2:5, 6, 7:9, 10]
                color_choices = [:blue, :green, :sienna, :black, :red]
                file_tags = ["bottom", "middle", "top"]
                circlemaker(x, y, r) = Plots.Shape(r * sind.(0:10:360) .+ x, r * cosd.(0:10:360) .+ y)

                for (s, sequence) in enumerate(sequences)
                    Plots.plot()
                    dataₒ = scale_transformation(o, d_data_dict[m][o]["data"][sequence, :])
                    # cond   = length(sequence) == 1 ? vec(.!any(isnan.(dataₒ), dims=2)) : vec(.!any(isnan.(dataₒ), dims=1))
                    cond = length(sequence) == 1 ? vec(.!any(isnan.(r_data_dict[m][o]["data"][sequence..., :]), dims=2)) : vec(.!any(isnan.(r_data_dict[m][o]["data"][sequence, :]), dims=1))
                    select_annual = source == "CEX" && M == "Consum" ? 4 : 1

                    # Condition
                    idx_cols = findall(cond)
                    idx_cols = idx_cols[1:select_annual:end]        # same downsampling as sxaxis

                    sxaxis = xaxis[idx_cols]
                    sdata = length(sequence) == 1 ? dataₒ[idx_cols] : dataₒ[:, idx_cols]
                    l_seq = length(sequence)

                    for (i, j) in enumerate(sequence)
                        # sci_l  = scale_transformation(o, ci_l[m][o][j, :][cond])
                        # sci_u  = scale_transformation(o, ci_u[m][o][j, :][cond])

                        indices = l_seq == 1 ? (:, :) : (i, :)
                        Plots.plot!(
                            sxaxis,
                            sdata[indices...],
                            # xlabel = L"\textrm{Years}",
                            ylabel=L"\textrm{\,\,in\,\, 10k\, USD}",
                            xformatter=:latex,
                            yformatter=:latex,
                            xticks=(sxaxis[1:interval_dts:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(dts_to_keep[idx_cols][1:interval_dts:end])]),
                            label="", #label_q(j), #L"\textrm{Bottom\,50}",, 
                            linewidth=3,
                            ls=:dot,
                            xtickfontsize=16,
                            ytickfontsize=16,
                            guidefontsize=16,
                            lc=:gray, #color_choices[i],
                            dpi=500,
                            legend=false
                        )
                        fs = 12
                        if M == "Consum"
                            if j <= 5
                                if j == 1
                                    Plots.annotate!(sxaxis[end-10], sdata[indices...][end-6] .* 1.10, Plots.text(label_q(j), fs, :bottom))
                                else
                                    Plots.annotate!(sxaxis[end-10], sdata[indices...][end-4], Plots.text(label_q(j), fs, :bottom))
                                end
                            elseif j >= 7 && j <= 9
                                Plots.annotate!(sxaxis[end-10], sdata[indices...][end-4], Plots.text(label_q(j), fs, :bottom))
                            elseif j == 6
                                Plots.annotate!(sxaxis[end-10], sdata[indices...][end-10] .* 0.85, Plots.text(label_q(j), fs, :bottom))
                            else
                                Plots.annotate!(sxaxis[end-10], sdata[indices...][end-5] .* 1.015, Plots.text(label_q(j), fs, :bottom))
                            end
                        elseif M == "Income"
                            if j <= 5
                                Plots.annotate!(sxaxis[end-4], sdata[indices...][end-3] * 1.08, Plots.text(label_q(j), fs, :bottom))
                            elseif j >= 6 && j <= 9
                                Plots.annotate!(sxaxis[end-4], sdata[indices...][end-3] .* 1.08, Plots.text(label_q(j), fs, :bottom))
                            else
                                Plots.annotate!(sxaxis[end-5], sdata[indices...][end-4] .* 1.06, Plots.text(label_q(j), fs, :bottom))
                            end
                        elseif M == "Wealth"
                            if j <= 5
                                if j ∈ [4, 5]
                                    Plots.annotate!(sxaxis[end-6], sdata[indices...][end-5] .* 1.06, Plots.text(label_q(j), fs, :bottom))
                                elseif j == 2
                                    Plots.annotate!(sxaxis[end-2], [-2], Plots.text(label_q(j), fs, :bottom)) # sdata[indices...][end-5] .* 0.01
                                elseif j == 3
                                    Plots.annotate!(sxaxis[end-9], sdata[indices...][5], Plots.text(label_q(j), fs, :bottom))
                                else
                                    Plots.annotate!(sxaxis[end-2], sdata[indices...][end-2] .* 0.8, Plots.text(label_q(j), fs, :bottom))
                                end
                            elseif j >= 6 && j <= 9
                                if j == 6
                                    Plots.annotate!(sxaxis[end-7], sdata[indices...][end-7] .* 0.4, Plots.text(label_q(j), fs, :bottom))
                                else
                                    Plots.annotate!(sxaxis[end-3], sdata[indices...][end] .* 0.72, Plots.text(label_q(j), fs, :bottom))
                                end
                            else
                                Plots.annotate!(sxaxis[end-5], sdata[indices...][end-5] .* 1.05, Plots.text(label_q(j), fs, :bottom))
                            end
                        end
                        # Separated since i dont want a label for the intervals 
                        # M = rand(12,2)
                        # circles = circlemaker.([M[:,1], M[:,2]]..., (0.02,))
                        # p = Plots.plot(M[:,1], M[:,2], ls=:dot, c=:red, framestyle=:box, grid=false)
                        # Plots.plot!(circles, ratio=1, fc=:transparent, lc=:darkred, lw=1, label=false)
                        # Plots.scatter!([-1e16], [circles[1].y[1]], mc=:white, msc=:darkred, lw=1, label="I am transparent!", xlims=Plots.xlims(p), ylims=Plots.ylims(p))

                        if i != l_seq
                            Plots.scatter!(sxaxis, sdata[indices...], marker=:square, mc=:white, msc=:black, msw=4, ms=6, label="", dpi=500)
                        else
                            Plots.scatter!(sxaxis, sdata[indices...], marker=:square, mc=:white, msc=:black, msw=4, ms=6, label="", dpi=500)
                            Plots.scatter!([], [], marker=:square, mc=:white, msc=:black, msw=4, dpi=500)
                        end
                    end
                    println(size(scale_transformation(o, r_data_dict[m][o]["data"][sequence, :])))

                    Plots.scatter!(
                        sxaxis,
                        length(sequence) == 1 ? scale_transformation(o, r_data_dict[m][o]["data"][sequence, :])[:, idx_cols]' : scale_transformation(o, r_data_dict[m][o]["data"][sequence, idx_cols])',
                        label=["" "" "" ""],
                        dpi=500,
                        msc=:orange,
                        mc=:white,
                        ms=3,
                        msw=2
                    )

                    Plots.scatter!(
                        sxaxis,
                        scale_transformation(o, r_data_dict[m][o]["data"][sequence[end], idx_cols]),
                        # label=L"\textrm{Approximation}",
                        msc=:orange,
                        mc=:white,
                        ms=3,
                        msw=2,
                        dpi=500,
                    )

                    Plots.pdf(path * m * "_" * o * "_$(file_tags[s])_quantiles" * ".pdf")
                end
            end
        end
    end
end



function f_b(o, grid)
    local plot_label
    if grid == 10 || grid == 20 || grid == 100
        plot_label = o == "quantiles" ? L"\textrm{Median}" : L"\textrm{Bottom\,50}"
    elseif grid == 5
        plot_label = o == "quantiles" ? L"40\textrm{th}" : L"\textrm{Bottom\,40}"
    end
    return plot_label
end


function f_n(o, grid)
    local plot_label
    if grid == 10 || grid == 20 || grid == 100
        plot_label = o == "quantiles" ? L"90\textrm{th}" : L"\textrm{Next\,40}"
    elseif grid == 5
        plot_label = o == "quantiles" ? L"80\textrm{th}" : L"\textrm{Next\,40}"
    end
    return plot_label
end


function f_t(o, grid)
    local plot_label
    if grid == 10 || grid == 20 || grid == 100
        plot_label = L"\textrm{Top\,10}"
    elseif grid == 5
        plot_label = L"\textrm{Top\,20}"
    end
    return plot_label
end

label_q(i) = i == 1 ? L"\textbf{%$(i)st \,\,Decile}" : i == 2 ? L"\textbf{%$(i)nd \,\,Decile}" : i == 3 ? L"\textbf{%$(3)rd \,\,Decile}" : L"\textbf{%$(i)th \,\,Decile}"