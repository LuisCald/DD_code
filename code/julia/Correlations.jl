function create_micro_df(copula_obj, pcfs, data_tag, k, folder, time_params, method_options; return_df=false)
    @unpack measures, equivalized, bottom_coded, case, estimator, tag = method_options
    @unpack tmin, tmax = time_params
    @unpack confidence_intervals = func_data
    @unpack grid_cop = estimator

    # if data_tag == "_detrended" 
    #     return 
    # end

    if typeof(estimator) <: SeriesEstimator
        @unpack integral_cop_grid = estimator
    end
    grid_size_data_cop = typeof(estimator) <: SeriesEstimator ? integral_cop_grid : grid_cop

    # Size of objects 
    dim = length(measures)
    # cop_size = tuple([grid_size_data_cop for i in 1:dim]...)
    m_grid = measures .* "grid"


    # Dates of the series 
    q_dates = QuarterlyDate(tmin["year"], tmin["quarter"]):Quarter(1):QuarterlyDate(tmax["year"], tmax["quarter"])
    T = length(q_dates)

    micro_full_df = DataFrame() #micro_df, Symbol.(["cop_share", "grid_point", measures..., m_grid...]))

    cop_ind = copula_case(measures, measures)

    # For each year ...
    # p = Progress(T, desc="Creating micro data for $k")
    # Threads.@threads 

    println("Creating micro data for $k")
    for y in eachindex(q_dates)
        if any(isfinite, copula_obj[cop_ind..., y])
            # Make dataframe 
            micro_df = construct_micro_dataset(copula_obj[cop_ind..., y], pcfs[:, y], grid_size_data_cop) #construct_micro_dataset_weighted(copula_y, pcfs_y, grid)
            micro_DF = DataFrame(micro_df, Symbol.(["cop_share", "grid_point", measures..., m_grid...])) # 
            micro_DF[!, "time"] .= q_dates[y]

            append!(micro_full_df, micro_DF)
            # next!(p)
        else
            nothing
        end
    end
    # finish!(p)

    # Sort dataframe by time 
    try
        sort!(micro_full_df, :time)
    catch ee
        println(ee)
    end

    # Export the data
    m_label = measures_folder(measures)
    init_path = BASE_PATH
    path = init_path * "/7_Results/$m_label" * "$tag" * "/$folder/data"
    mkpath(path)
    equiv = equivalized == true ? "eq" : ""
    botcod = !isempty(bottom_coded) ? "bc" : ""
    label = "_$case" * "_$equiv" * "$botcod"
    CSV.write(path * "/" * k * "_micro_data" * data_tag * "$label" * ".csv", micro_full_df)

    if return_df
        return micro_full_df
    end

end


function kendalls_tau(copula_obj, pcfs, k, folder, time_params, method_options, return_option, plot_tau, plot_data, gdp_series)
    @unpack measures, equivalized, bottom_coded, case, grid, tag = method_options
    @unpack tmin, tmax = time_params
    @unpack confidence_intervals = func_data

    # Size of objects 
    dim = length(measures)
    T = size(copula_obj, 2)
    cop_size = tuple([grid for i in 1:dim]...)
    m_grid = measures .* "grid"

    # Dates of the series 
    q_dates = QuarterlyDate(tmin["year"], tmin["quarter"]):Quarter(1):QuarterlyDate(tmax["year"], tmax["quarter"])

    # For the different combinations of variables
    combs = join.(combinations(measures, 2), "\\,")
    df_combs = join.(combinations(measures, 2), "_")

    # Containers for Kendalls Tau and Pearson Correlation
    kend_tau = zeros(T, length(combs))
    pwcorr = zeros(T, length(combs))

    # Containers for the granular correlations
    local g_kend_tau, g_pwcorr
    if grid >= 10 && dim >= 2  #TODO: why grid >= 10? 
        # Generate containers for the granular correlations 
        g_kend_tau = zeros(T, length(combs) * 4)
        g_pwcorr = zeros(T, length(combs) * 4)

        # For empircal correlations, later
        g_ekend_tau = zeros(T, length(combs) * 4)
        g_epwcorr = zeros(T, length(combs) * 4)
    end


    micro_full_df = DataFrame() #micro_df, Symbol.(["cop_share", "grid_point", measures..., m_grid...]))


    # For each year ...
    for y in eachindex(q_dates)
        if any(isfinite, copula_obj[:, y])
            # Get series for households 
            # tot_hhs            = filter(row -> row.date == q_dates[y], gdp_series)[:, "tot_hhs"][1] # gdp_series is already subset to the time frame of the estimation 
            copula_y = reshape(copula_obj[:, y], cop_size)
            pcfs_y = pcfs[:, y]
            # condition          = (!isnan).(pcfs[:, y])

            # Make dataframe 
            micro_df = construct_micro_dataset(copula_y, pcfs_y, grid) #construct_micro_dataset_weighted(copula_y, pcfs_y, grid)
            micro_DF = DataFrame(micro_df, Symbol.(["cop_share", "grid_point", measures..., m_grid...])) # 
            micro_DF[!, "time"] .= q_dates[y]

            # Compute correlations from the micro df 
            kend_tau[y, :], pwcorr[y, :] = compute_correlations(micro_df, measures)
            # kend_tau[y, scomb+1:end] .= NaN
            # pwcorr[y, scomb+1:end]   .= NaN

            # Compute granular correlations from the micro df 
            if grid >= 10 && dim == 2
                g_kend_tau[y, :], g_pwcorr[y, :] = compute_granular_correlations(micro_df, measures, grid)
                # g_kend_tau[y, scomb+1:end] .= NaN
                # g_pwcorr[y, scomb+1:end]   .= NaN
            end

            # if y == 1
            #     micro_full_df     = micro_DF
            # else
            append!(micro_full_df, micro_DF)
            # end


        else
            kend_tau[y, :] .= NaN
            pwcorr[y, :] .= NaN
            # micro_full_df        = DataFrame(zeros(10, 2+dim), Symbol.(["time", "grid_point", measures...])) #TODO: something random

            if grid >= 10 && dim == 2
                g_kend_tau[y, :] .= NaN
                g_pwcorr[y, :] .= NaN
            end

        end
    end

    # Generate name combinations
    kend_tau = DataFrame(kend_tau, Symbol.([df_combs...]))
    kend_tau[!, "time"] = q_dates

    pwcorr = DataFrame(pwcorr, Symbol.([df_combs...]))
    pwcorr[!, "time"] = q_dates

    # micro_full_df       = DataFrame(micro_full_df, Symbol.(["time", "grid_point", measures...]))

    # Export and plot Kendall's Tau 
    xaxis = collect(1:length(q_dates))

    m_label = measures_folder(measures)
    equiv = equivalized == true ? "eq" : ""
    botcod = !isempty(bottom_coded) ? "bc" : ""
    label = "_$case" * "_$equiv" * "$botcod"
    line_sty = [:solid :dash :dot :dashdot :dashdotdot]

    # Export
    init_path = BASE_PATH
    path = init_path * "/7_Results/$m_label" * tag * "/$folder/data"
    mkpath(path)
    CSV.write(path * "/" * k * "_kendalls_tau" * "$label" * ".csv", select(kend_tau, "time", :))
    CSV.write(path * "/" * k * "_pearson" * "$label" * ".csv", select(pwcorr, "time", :))
    CSV.write(path * "/" * k * "_micro_data" * "$label" * ".csv", micro_full_df)

    # Treatment of the granular correlations
    if grid >= 10 && dim == 2
        col_names = []
        df_names = []

        # Idea: Sort everything with the first columns being the combinations of the lowest possible dimension (2) then everything after 

        # Generate the combinations, where each combination is a pair 
        # comb           = obs_dims == dim ? combinations(1:dim, 2) : obs_dims < dim && obs_dims >= 2 ? obs_ids : println("Is this a 4 dimensional case? Check this")

        # Generate combinations for the granular correlations
        comb = combinations(1:dim, 2)
        for c in comb
            for p in ["low", "high"]
                for q in ["low", "high"]
                    push!(df_names, "$p" * "_$(measures[c[1]])" * "," * "_$q" * "_$(measures[c[2]])")
                    push!(col_names, "$p" * "\\,$(measures[c[1]])" * "\\," * "\\,$q" * "\\,$(measures[c[2]])")
                end
            end
        end

        sort!(col_names) #TODO: why sort here but not df_names? Check that this makes sense 

        # Prepare for export 
        if dim == 2
            g_kend_tau = DataFrame(g_kend_tau, Symbol.(df_names))
            g_kend_tau[!, "time"] = q_dates

            g_pwcorr = DataFrame(g_pwcorr, Symbol.(df_names))
            g_pwcorr[!, "time"] = q_dates

            CSV.write(path * "/" * k * "_granular_kendalls_tau" * "$label" * ".csv", select(g_kend_tau, "time", :))
            CSV.write(path * "/" * k * "_granular_pearson" * "$label" * ".csv", select(g_pwcorr, "time", :))
        end
    end

    # Computing Kendalls tau from the OBSERVED data 
    if plot_data != false && k != "consensus"
        @unpack func_dict = func_data
        T = size(func_dict[k]["copulas"]["data"])[end]
        copula_y = zeros(cop_size...)

        # Containers for (e)mpirical correlations
        ekend_tau = zeros(T, length(combs))
        epwcorr = zeros(T, length(combs))

        for y in 1:T
            cop = func_dict[k]["copulas"]["data"][:, y]
            if any(isfinite, cop) == true
                tot_hhs = filter(row -> row.date == q_dates[y], gdp_series)[:, "tot_hhs"][1] # gdp_series is already subset to the time frame of the estimation 

                # # In case not all measures are observed, but just enough
                # scomb      = length(combinations(obs_meas, 2))
                # comb_bound = scomb * 4 

                pcfs_y = zeros(grid * dim)

                # Fill containers 
                copula_y .= reshape(cop, cop_size)
                pcfs_y .= vcat([func_dict[k][m]["quantiles"]["data"][:, y] for m in measures]...)

                # Create micro df for the respective year 
                micro_df = construct_micro_dataset(copula_y, pcfs_y, grid, tot_hhs)

                # Compute empirical correlations for the year 
                ekend_tau[y, :], epwcorr[y, :] = compute_correlations(micro_df, measures)
                # ekend_tau[y, scomb+1:end]                       .= NaN
                # epwcorr[y, scomb+1:end]                         .= NaN 

                if grid >= 10 && dim == 2 #TODO: dimension equal to 2 because otherwise, we have 12 combinations -> too many lines in the plot, but will comeback to this.  
                    # Compute empirical granular correlations for the year 
                    # The order in which these are spit out is like on line 131
                    g_ekend_tau[y, :], g_epwcorr[y, :] = compute_granular_correlations(micro_df, measures, grid)
                    # g_ekend_tau[y, comb_bound+1:end]   .= NaN
                    # g_epwcorr[y, comb_bound+1:end]     .= NaN
                end
            else
                ekend_tau[y, :] .= NaN
                epwcorr[y, :] .= NaN
                if grid >= 10 && dim == 2
                    g_ekend_tau[y, :] .= NaN
                    g_epwcorr[y, :] .= NaN
                end
            end
        end
    end

    if plot_tau == true
        local path
        if length(folder) > 17
            path = init_path * "/7_Results/$m_label" * tag * "/$folder/correlations"
        else
            path = init_path * "/7_Results/$m_label" * tag * "/$folder/plots/correlations"
        end
        mkpath(path)

        # Plotting Kendalls Tau 
        if plot_data == false # so no empirical data in the plot 
            # Plotting correlations
            Plots.plot(
                xaxis,
                convert.(Float64, Matrix(select(kend_tau, Not("time")))),
                xlabel=L"\textrm{Year}",
                ylabel=L"\textrm{%$(k)\,\,Kendall's\,\,Tau}",
                xformatter=:latex,
                yformatter=:latex,
                xticks=(xaxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(q_dates[1:20:end])]),
                legend=:outertopright,
                label=reshape([L"\textrm{%$(c)}" for c in combs], (1, length(combs))),
                lw=2, dpi=500
            )

            Plots.savefig(path * "/" * k * "_kendalls_tau" * label * ".pdf")

            # Plotting granular correlations
            if grid >= 10 && dim == 2
                Plots.plot(
                    xaxis,
                    convert.(Float64, Matrix(select(g_kend_tau, Not("time")))),
                    xlabel=L"\textrm{Year}",
                    ylabel=L"\textrm{%$(k)\,\,Kendall's\,\,Tau}",
                    xformatter=:latex,
                    yformatter=:latex,
                    xticks=(xaxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(q_dates[1:20:end])]),
                    legend=:outertopright,
                    label=reshape([L"\textrm{%$(c)}" for c in col_names], (1, length(col_names))),
                    lw=2, dpi=500
                )

                Plots.savefig(path * "/" * k * "_granular_kendalls_tau" * label * ".pdf")
            end

            # Plotting normal correlation 
            Plots.plot(
                xaxis,
                convert.(Float64, Matrix(pwcorr)),
                xlabel=L"\textrm{Year}",
                ylabel=L"\textrm{%$(k)\,\,Pearson\,\,Correlation}",
                xformatter=:latex,
                yformatter=:latex,
                xticks=(xaxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(q_dates[1:20:end])]),
                legend=:outertopright,
                label=reshape([L"\textrm{%$(c)}" for c in combs], (1, length(combs))),
                lw=2, dpi=500
            )

            Plots.savefig(path * "/" * k * "_pearson_correlation" * label * ".pdf")

            # Plotting granular normal correlation
            if grid >= 10 && dim == 2
                Plots.plot(
                    xaxis,
                    convert.(Float64, Matrix(g_pwcorr)),
                    xlabel=L"\textrm{Year}",
                    ylabel=L"\textrm{%$(k)\,\,Pearson\,\,Correlation}",
                    xformatter=:latex,
                    yformatter=:latex,
                    xticks=(xaxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(q_dates[1:20:end])]),
                    legend=:outertopright,
                    label=reshape([L"\textrm{%$(c)}" for c in col_names], (1, length(col_names))),
                    lw=2, dpi=500
                )

                Plots.savefig(path * "/" * k * "_granular_pearson_correlation" * label * ".pdf")
            end

        elseif plot_data != false && k != "consensus"
            # Start the series at the first observation 
            cond = findall(!isnan, Matrix(ekend_tau)[:, 1])
            c_data = ekend_tau[cond, :]
            ci_l = confidence_intervals[k]["ci_l"]["correlations"]["ktau"][:, cond]
            ci_u = confidence_intervals[k]["ci_u"]["correlations"]["ktau"][:, cond]
            start = length(cond) == 0 ? 1 : cond[1]
            sxaxis = xaxis[start:end]
            sq_dates = q_dates[start:end]
            sdata = convert.(Float64, Matrix(select(kend_tau, Not("time"))))[start:end, :]
            csize = size(sdata, 2)

            Plots.plot(
                sxaxis,
                sdata,
                xlabel=L"\textrm{Year}",
                ylabel=L"\textrm{%$(k)\,\,Kendall's\,\,Tau}",
                xformatter=:latex,
                yformatter=:latex,
                ls=line_sty[:, 1:csize],
                xticks=(sxaxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(sq_dates[1:20:end])]),
                legend=:outertopright,
                label=reshape([L"\textrm{%$(c)}" for c in combs], (1, length(combs))),
                lw=2, dpi=500
            )

            Plots.scatter!(
                xaxis[cond],
                c_data,
                label=[L"\textrm{%$(k)\,\, obs}" ["" for i in 1:length(combs)-1]...],
                marker=:diamond,
                markercolor=:black,
                lw=2, dpi=500,
                fillcolor=:black,
                fillalpha=0.2,
                yerror=(c_data - ci_l', ci_u' - c_data),
            )

            Plots.savefig(path * "/" * k * "_kendalls_tau" * label * ".pdf")

            # Plotting granular correlations
            if grid >= 10 && dim == 2
                Plots.plot(
                    xaxis,
                    convert.(Float64, Matrix(select(g_kend_tau, Not("time")))),
                    xlabel=L"\textrm{Year}",
                    ylabel=L"\textrm{%$(k)\,\,Kendall's\,\,Tau}",
                    xformatter=:latex,
                    yformatter=:latex,
                    xticks=(xaxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(q_dates[1:20:end])]),
                    legend=:outertopright,
                    label=reshape([L"\textrm{%$(c)}" for c in col_names], (1, length(col_names))),
                    lw=2, dpi=500
                )

                # Plotting the data
                c_data = g_ekend_tau[.!isnan.(Matrix(g_ekend_tau)[:, 1]), :]
                ci_l = confidence_intervals[k]["ci_l"]["correlations"]["g_ktau"][:, findall(!isnan, Matrix(g_ekend_tau)[:, 1])]
                ci_u = confidence_intervals[k]["ci_u"]["correlations"]["g_ktau"][:, findall(!isnan, Matrix(g_ekend_tau)[:, 1])]

                Plots.scatter!(
                    xaxis[findall(!isnan, Matrix(g_ekend_tau)[:, 1])],
                    c_data,
                    label=[L"\textrm{%$(k)\,\, obs}" ["" for i in 1:length(col_names)-1]...],
                    marker=:diamond,
                    markercolor=:black,
                    lc=:black,
                    lw=2, dpi=500,
                    yerror=(c_data - ci_l', ci_u' - c_data),
                )

                Plots.savefig(path * "/" * k * "_granular_kendalls_tau" * label * ".pdf")
            end

            # Plotting normal correlation
            # Start the series at the first observation 
            cond = findall(!isnan, Matrix(epwcorr)[:, 1])
            c_data = epwcorr[cond, :]
            ci_l = confidence_intervals[k]["ci_l"]["correlations"]["pcorr"][:, cond]
            ci_u = confidence_intervals[k]["ci_u"]["correlations"]["pcorr"][:, cond]
            start = length(cond) == 0 ? 1 : cond[1]
            sxaxis = xaxis[start:end]
            sq_dates = q_dates[start:end]
            sdata = convert.(Float64, Matrix(select(pwcorr, Not("time"))))[start:end, :]

            Plots.plot(
                sxaxis,
                sdata,
                xlabel=L"\textrm{Year}",
                ylabel=L"\textrm{%$(k)\,\,Pearson\,\,Correlation}",
                xformatter=:latex,
                yformatter=:latex,
                ls=line_sty[:, 1:csize],
                xticks=(sxaxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(sq_dates[1:20:end])]),
                legend=:outertopright,
                label=reshape([L"\textrm{%$(c)}" for c in combs], (1, length(combs))),
                lw=2, dpi=500
            )



            Plots.scatter!(
                xaxis[cond],
                c_data,
                label=[L"\textrm{%$(k)\,\, obs}" "" ""],
                marker=:diamond,
                markercolor=:black,
                # lc=:black,
                lw=2, dpi=500,
                fillalpha=0.2,
                yerror=(c_data - ci_l', ci_u' - c_data),
            )

            Plots.savefig(path * "/" * k * "_pearson_correlation" * label * ".pdf")

            # Plotting granular correlations
            if grid >= 10 && dim == 2
                Plots.plot(
                    xaxis,
                    convert.(Float64, Matrix(select(g_pwcorr, Not("time")))),
                    xlabel=L"\textrm{Year}",
                    ylabel=L"\textrm{%$(k)\,\,Pearson\,\,Correlation}",
                    xformatter=:latex,
                    yformatter=:latex,
                    xticks=(xaxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(q_dates[1:20:end])]),
                    legend=:outertopright,
                    label=reshape([L"\textrm{%$(c)}" for c in col_names], (1, length(col_names))),
                    lw=2, dpi=500
                )

                # Plotting the data
                c_data = g_epwcorr[.!isnan.(Matrix(g_epwcorr)[:, 1]), :]
                ci_l = confidence_intervals[k]["ci_l"]["correlations"]["g_pcorr"][:, findall(!isnan, Matrix(g_epwcorr)[:, 1])]
                ci_u = confidence_intervals[k]["ci_u"]["correlations"]["g_pcorr"][:, findall(!isnan, Matrix(g_epwcorr)[:, 1])]

                Plots.scatter!(
                    xaxis[findall(!isnan, Matrix(g_epwcorr)[:, 1])],
                    c_data,
                    label=[L"\textrm{%$(k)\,\, obs}" ["" for i in 1:length(col_names)-1]...],
                    marker=:diamond,
                    markercolor=:black,
                    lc=:black,
                    lw=2, dpi=500,
                    fillcolor=:black,
                    fillalpha=0.1,
                    yerror=(c_data - ci_l', ci_u' - c_data),
                )

                Plots.savefig(path * "/" * k * "_granular_pearson_correlation" * label * ".pdf")
            end

        end
    end

    if return_option == true
        return convert.(Float64, Matrix(select(kend_tau, Not("time"))))
    end #TODO: should be all columns?
end



# function compute_kendalls_tau_2d(empirical_copula::Matrix{T}) where T <: Real
#     n = size(empirical_copula, 1)  # number of rows/columns in the matrix

#     # Compute the marginal probabilities for each variable
#     p = sum(empirical_copula, dims=2) ./ n
#     q = sum(empirical_copula, dims=1) ./ n

#     # Compute the mean ranks for each variable
#     mean_ranks = (1:n) .+ (n+1)/2

#     # Compute the observed and expected number of concordant and discordant pairs
#     n_conc = 0
#     n_disc = 0
#     for i in 1:n-1
#         for j in i+1:n
#             conc = 0
#             disc = 0
#             for k in 1:n
#                 if empirical_copula[i,k] <= empirical_copula[j,k]
#                     conc += 1
#                 else
#                     disc += 1
#                 end
#             end
#             n_conc += conc * p[i] * p[j]
#             n_disc += disc * p[i] * p[j]
#         end
#     end

#     # Compute Kendall's tau and its variance
#     tau = (n_conc - n_disc) / (n * (n-1) / 2)
#     var_tau = 2 * (2*n + 5) / (9*n * (n-1))

#     # Compute the standard error and p-value
#     se = sqrt(var_tau)
#     z = tau / se
#     pvalue = 2 * (1 - cdf(Normal(0, 1), abs(z)))

#     return tau, pvalue
# end


function compute_kendalls_tau(cop)
    """
    Computes Kendall's Tau correlation coefficient from a given copula density matrix.
    From: https://www.uni-muenster.de/Physik.TP/~lemm/seminarSS08/JHE-2007.pdf
    From: http://www.scielo.org.co/pdf/rce/v43n1/0120-1751-rce-43-01-3.pdf
    """
    n = size(cop, 1)
    edf = copula_cdf(cop)
    c_edf = copula_integral(edf, n)
    tau = 4 * c_edf - 1

    return tau
end


# function compute_correlations(micro_df, measures)
#     """
#     Computes Kendall's Tau correlation coefficient from a given copula density matrix.
#     From: https://www.uni-muenster.de/Physik.TP/~lemm/seminarSS08/JHE-2007.pdf
#     From: http://www.scielo.org.co/pdf/rce/v43n1/0120-1751-rce-43-01-3.pdf
#     """

#     # Kendalls Tau for all pairs of variables 
#     dim = length(measures)
#     rvs = convert.(Float64, Matrix(micro_df[:, 3:end-dim]))

#     R"""
#     library(wdm)
#     pwco <- wdm($rvs, method = "pearson", weights = $(micro_df[:,1]))  # matrices
#     taus <- wdm($rvs, method = "kendall", weights = $(micro_df[:,1]))  # matrices
#     """

#     @rget pwco
#     @rget taus

#     # Replace missing with NaN 
#     pwco = coalesce.(pwco, NaN)
#     taus = coalesce.(taus, NaN)

#     # println(pwco)
#     # println(taus)
#     # taus = corkendall(convert.(Float64, Matrix(micro_df[:, 3:end])))  # col1 = cop shares, col2 = cop group  #TODO: could fail in some cases because of the combination of numbers 
#     # pwco = cor(micro_df[:, 3:end]) #, Weights(dataset[:, 1]))
#     comb = combinations(1:dim, 2)

#     return [taus[c...] for c in comb], [pwco[c...] for c in comb]
# end

# function compute_granular_correlations(micro_df, measures, grid)
#     dim = length(measures)

#     # Define groups based on measures 
#     group_dict = find_groups(grid, measures)

#     # Divide the dataframe into these groups based on the value of a row 
#     split_df = split_micro_data(micro_df, group_dict)

#     # Container for the correlations 
#     comb = combinations(1:dim, 2)
#     tau_vec = Vector{Vector{Float64}}(undef, length(split_df))
#     pwco_vec = Vector{Vector{Float64}}(undef, length(split_df))

#     # compute the correlations for the different groups 
#     for (i, k) in enumerate(sort(collect(keys(split_df))))

#         # Find columns associated with the group
#         sub_meas = split(k, "_")[2:2:end] #TODO: for 2D case, this is equal to the measures. For 3D, there are 3 pairs so ..
#         m_ind = [findfirst(x -> x == elem, measures) for elem in sub_meas]

#         # Compute the correlations
#         dff = split_df[k][:, 3:end-dim]

#         R"""
#         library(wdm)
#         pwco <- wdm($dff, method = "pearson", weights = $(split_df[k][:, 1]))  # matrices
#         taus <- wdm($dff, method = "kendall", weights = $(split_df[k][:, 1]))  # matrices
#         """

#         @rget pwco
#         @rget taus

#         # Replace missing with NaN 
#         pwco = coalesce.(pwco, NaN)
#         taus = coalesce.(taus, NaN)

#         # taus        = corkendall(convert.(Float64, Matrix(dff)))  # col1 = cop shares, col2 = cop group  #TODO: could fail in some cases because of the combination of numbers 
#         # pwco        = cor(dff) 

#         # From these correlations, extract the relevant ones
#         tau_vec[i] = [taus[m_ind...]]
#         pwco_vec[i] = [pwco[m_ind...]]
#     end

#     return reduce(vcat, tau_vec), reduce(vcat, pwco_vec)
# end

# a = rand(100,3)
# comb        = combinations(1:3, 2)
# taus        = corkendall(a)
# pwco        = cor(a)
# tau_vec     = [taus[meas_ind...]]
# pwco_vec    = [pwco[c...] for c in comb]

# tau_l = [tau_vec, tau_vec, tau_vec]
# reduce(vcat, tau_l)

# a = rand(100,2)
# comb        = combinations(1:2, 2)
# taus        = corkendall(a)
# pwco        = cor(a)
# tau_vec     = [taus[c...] for c in comb]
# pwco_vec    = [pwco[c...] for c in comb]

# tau_l = [tau_vec, tau_vec, tau_vec]
# reduce(vcat, tau_l)

# string_value = "high_income_low_debt"
# sub_meas = split(string_value, "_")[2:2:end]
# meas_ind = [findfirst(x -> x == elem, ["income", "wealth", "debt"]) for elem in sub_meas]



function split_micro_data(micro_df, group_dict)
    # Column of numbers
    group_vec = Dict()
    groups = collect(keys(group_dict))
    for k in groups
        group_vec[k] = []
    end

    # Iterate over each element in the column
    for (index, value) in enumerate(micro_df[:, 2])
        for i in eachindex(groups)
            if value ∈ group_dict[groups[i]]
                push!(group_vec[groups[i]], index)
            end
        end
    end

    # Split the matrix based on the groups
    split_df = Dict()
    for i in eachindex(groups)
        split_df[groups[i]] = micro_df[group_vec[groups[i]], :]
    end

    return split_df
end

function find_groups(grid, measures)
    # Assuming the variables are already defined
    half_point = floor(Int, grid / 2)
    group_dict = Dict()
    dim = length(measures)

    # Generate the combinations
    comb = combinations(1:dim, 2)

    for c in comb
        for p in ["low", "high"]
            for q in ["low", "high"]
                group_dict["$p"*"_$(measures[c[1]])"*"_$q"*"_$(measures[c[2]])"] = []
            end
        end
    end

    if dim == 2
        # Check the length of number_part
        for first_number in 1:grid
            for second_number in 1:grid
                # Cases
                if first_number <= half_point && second_number <= half_point
                    push!(group_dict["low"*"_$(measures[1])"*"_low"*"_$(measures[2])"], parse(Int, string(first_number) * string(second_number)))

                elseif first_number > half_point && second_number > half_point
                    push!(group_dict["high"*"_$(measures[1])"*"_high"*"_$(measures[2])"], parse(Int, string(first_number) * string(second_number)))

                elseif first_number > half_point && second_number <= half_point
                    push!(group_dict["high"*"_$(measures[1])"*"_low"*"_$(measures[2])"], parse(Int, string(first_number) * string(second_number)))

                elseif first_number <= half_point && second_number > half_point
                    push!(group_dict["low"*"_$(measures[1])"*"_high"*"_$(measures[2])"], parse(Int, string(first_number) * string(second_number)))
                end
            end
        end

    elseif dim == 3
        # Check the length of number_part
        for c in comb
            for first_number in 1:grid
                for second_number in 1:grid
                    for third_number in 1:grid
                        gp = parse(Int, string(first_number) * string(second_number) * string(third_number))  # grid point 

                        first_gp = parse(Int, string(gp)[c[1]])
                        second_gp = parse(Int, string(gp)[c[2]])

                        # Cases
                        if first_gp <= half_point && second_gp <= half_point
                            push!(group_dict["low"*"_$(measures[c[1]])"*"_low"*"_$(measures[c[2]])"], gp)

                        elseif first_gp > half_point && second_gp > half_point
                            push!(group_dict["high"*"_$(measures[c[1]])"*"_high"*"_$(measures[c[2]])"], gp)

                        elseif first_gp > half_point && second_gp <= half_point
                            push!(group_dict["high"*"_$(measures[c[1]])"*"_low"*"_$(measures[c[2]])"], gp)

                        elseif first_gp <= half_point && second_gp > half_point
                            push!(group_dict["low"*"_$(measures[c[1]])"*"_high"*"_$(measures[c[2]])"], gp)
                        end
                    end
                end
            end
        end
    end

    for k in collect(keys(group_dict))
        group_dict[k] = unique(group_dict[k])
    end

    return group_dict
end


# TODO: only works if pcf_grid and cop_grid are the same
function construct_micro_dataset(cop, pcfs, grid)
    dim = length(size(cop))

    dataset = zeros(prod(size(cop)), 2 * dim + 2)
    if dim == 3
        for i in 1:grid
            for j in 1:grid
                for k in 1:grid
                    combined_element = parse(Int, string(i) * string(j) * string(k))
                    dataset[(i-1)*(grid*grid)+(j-1)*grid+k, :] = [cop[i, j, k], combined_element, pcfs[i], pcfs[grid+j], pcfs[grid+grid+k], i, j, k]
                end
            end
        end

    elseif dim == 2
        for i in 1:grid
            for j in 1:grid
                combined_element = parse(Int, string(i) * string(j))
                dataset[(i-1)*grid+j, :] = [cop[i, j], combined_element, pcfs[i], pcfs[grid+j], i, j]
            end
        end
    end

    # Scale the predictions s.t. it all fits within [0,1] and sums to 1
    # NaN_cond   = (!isnan).(dataset[:,1])
    # min_weight = minimum(dataset[NaN_cond, 1])
    # max_weight = maximum(dataset[NaN_cond, 1]) 
    # sum_weight = sum(dataset[NaN_cond, 1]) # basically 1

    # # min-max scale 
    # dataset[:, 1] .= (dataset[:, 1] .- min_weight) ./ (max_weight .- min_weight) 
    # dataset[:, 1] .= (dataset[:, 1] ./ sum_weight) .* tot_hhs


    # dataset        = dataset[dataset[:,1] .>= 0, :]

    # Drop rows which have negative weight 
    # dataset = dataset[dataset[:,1] .<= 1, :]
    # mult    = size(dataset, 2) - 2 == 3 ? 1000000 : 10000

    # # Generate dataset     
    # tau_dataset = vcat([repeat(dataset[i,:]', inner=(max.(Int_NaN.(ceil.(dataset[i,1] * mult)), 0),1)) for i in axes(dataset,1)]...)

    return dataset
end

function construct_micro_dataset_weighted(cop, pcfs, grid)
    dim = length(size(cop))
    dataset = zeros(prod(size(cop)), 2 * dim + 2)
    if dim == 3
        for i in 1:grid
            for j in 1:grid
                for k in 1:grid
                    combined_element = parse(Int, string(i) * string(j) * string(k))
                    dataset[(i-1)*(grid*grid)+(j-1)*grid+k, :] = [cop[i, j, k], combined_element, pcfs[i], pcfs[grid+j], pcfs[grid+grid+k], i, j, k]
                end
            end
        end

    elseif dim == 2
        for i in 1:grid
            for j in 1:grid
                combined_element = parse(Int, string(i) * string(j))
                dataset[(i-1)*grid+j, :] = [cop[i, j], combined_element, pcfs[i], pcfs[grid+j], i, j]
            end
        end
    end

    # Drop rows where they are less than .0001. Why? x 10000 # TODO: why not use tot_hhs
    dataset = dataset[dataset[:, 1].>=0.000001, :]
    mult = size(dataset, 2) - 2 == 3 ? 1000000 : 10000

    # Generate dataset     
    tau_dataset = vcat([repeat(dataset[i, :]', inner=(max.(Int_NaN.(ceil.(dataset[i, 1] * mult)), 0), 1)) for i in axes(dataset, 1)]...)

    return tau_dataset
end

function Int_NaN(x)
    if isnan(x) == true
        return 0
    else
        return Int(x)
    end
end

# a = rand(2,2)
# vcat([repeat(a[i,:]', inner=(1,1)) for i in 1:2]...)

# dataset = rand(1000,2) 
# dataset = (dataset ./ sum(dataset, dims=2)) 
# dataset = hcat(dataset, 1:size(dataset,1)) 
# tau_dataset = vcat([repeat(dataset[i,:]', inner=(max.(Int.(floor.(round.(dataset[i,1] * 10000, digits=2))), 0),1)) for i in axes(dataset,1)]...)
# unique(tau_dataset[:,3])
# combine(groupby(tau, :grid_point, nrow => :my_desired_name))

function rep(x::AbstractVector, lengths::AbstractVector{T}) where T<:Integer
    if length(x) != length(lengths)
        throw(DimensionMismatch("vector lengths must match"))
    end
    res = similar(x, sum(lengths))
    i = 1
    for idx in axes(x)
        tmp = x[idx]
        for _ in 1:lengths[idx]
            res[i] = tmp
            i += 1
        end
    end
    return res
end

function copula_integral(cop, grid)
    local c_edf
    if length(size(cop)) == 2
        vx = range(0, 1, length=grid)
        vy = range(0, 1, length=grid)
        c_edf = trapz((vx, vy), cop)
    elseif length(size(cop)) == 3
        vx = range(0, 1, length=grid)
        vy = range(0, 1, length=grid)
        vz = range(0, 1, length=grid)
        c_edf = trapz((vx, vy, vz), cop)
    end

    return c_edf
end


function compute_tail_dependence(copula_obj, pcfs, data_name, folder, time_params, method_options)
    @unpack measures, equivalized, bottom_coded, case, grid, tag = method_options
    @unpack tmin, tmax = time_params

    # Some params and initialization of the tail dependence vectors
    T = size(copula_obj)[end]
    q_dates = QuarterlyDate(tmin["year"], tmin["quarter"]):Quarter(1):QuarterlyDate(tmax["year"], tmax["quarter"])
    pairs = 1:grid-1
    dimension = length(measures)

    local tail_dep_l, tail_dep_u
    if dimension == 2
        copula_y = zeros(grid, grid)
        tail_dep_l = zeros(T, length(pairs))
        tail_dep_u = zeros(T, length(pairs))

        for (q, p) in enumerate(pairs)
            for i in eachindex(q_dates)
                copula_y .= reshape(copula_obj[:, i], (grid, grid))
                result = compute_tail_dependence_2d(copula_y, p)
                tail_dep_u[i, q] = result[1]
                tail_dep_l[i, q] = result[2]
            end
        end

        # Export tail_dep_l and tail_dep_u to excel 
        half_way_p = floor(Int, grid / 2)
        tail_dep_u = DataFrame(tail_dep_u[:, half_way_p+1:grid-1], Symbol.(["tdu_$(i)" for i in half_way_p+1:grid-1]))
        tail_dep_l = DataFrame(tail_dep_l[:, 1:half_way_p], Symbol.(["tdl_$(i)" for i in 1:half_way_p]))

        # Merge the tail dependence dataframes into a single DataFrame
        tail_dep = hcat(tail_dep_l, tail_dep_u)
        tail_dep[!, "time"] = q_dates

        # Export the tail dependence dataframes to excel
        m_label = measures_folder(measures)
        equiv = equivalized == true ? "eq" : ""
        botcod = !isempty(bottom_coded) ? "bc" : ""
        label = "_$case" * "_$equiv" * "$botcod"

        # Export
        init_path = BASE_PATH
        path = init_path * "/7_Results/$m_label" * tag * "/$folder/data"
        mkpath(path)
        CSV.write(path * "/" * data_name * "_tail_dependence" * "$label" * ".csv", select(tail_dep, "time", :))

        # Plot Tail Dependence 
        tail_mat = Matrix(select(tail_dep, Not(:time)))
        seq = 1:50:size(tail_mat, 1)

        Plots.plot()
        for j in eachindex(seq)
            Plots.plot!(
                1:grid-1,
                tail_mat[seq[j], :],
                xlabel=L"\textrm{Lower\,\,Tail\,\,and\,\,Upper\,\,Tail}",
                ylabel=L"\textrm{Tail\,\,Dependence}",
                xformatter=:latex,
                yformatter=:latex,
                label=L"%$(q_dates[seq[j]])"[1:5] * "\$",
                xticks=(1:grid-1, [L"0.%$(i)" for i in 1:grid-1]),
                legend=:outertopright,
                lw=2, dpi=500
            )
        end
        # for i in 2:length(seq)
        #     Plots.plot!(1:grid-1, tail_mat[seq[i], :], label="", lw=2, dpi=500)
        # end

        # Plots.plot!(1:grid-1, tail_mat[seq[end], :], xformatter=:latex, yformatter=:latex,  label = L"%$(q_dates[seq[end]])"[1:5] * "\$", lw=2, dpi=500)
        mid_point = floor(Int, grid / 2)
        Plots.vline!([mid_point], yformatter=:latex, linestyle=:dash, linecolor=:grey, label="")

        path = init_path * "/7_Results/$m_label" * tag * "/$folder/plots/correlations"
        mkpath(path)
        Plots.savefig(path * "/" * data_name * "_tail_dependence" * "$label" * ".pdf")


    elseif dimension == 3

        copula_y = zeros(grid, grid, grid)
        itr_combs = combinations(measures, 2) # For all combinations: generate_combinations(measures)
        combs = join.(itr_combs, "_")
        m_cols = [m * "grid" for m in measures]

        # In this case, we need to do the tail dependence for each pair of variables
        tail_dep_l = Dict()
        tail_dep_u = Dict()

        for c in combs
            tail_dep_l[c] = zeros(T, length(pairs)) # 9 pairs := (1,1) (2,2) ... 
            tail_dep_u[c] = zeros(T, length(pairs))
        end

        cop_2d = zeros(grid, grid)
        for (q, p) in enumerate(pairs)
            for i in eachindex(q_dates)
                tot_hhs = filter(row -> row.date == q_dates[i], gdp_series)[:, "tot_hhs"][1] # gdp_series is already subset to the time frame of the estimation 
                copula_y .= reshape(copula_obj[:, i], (grid, grid, grid))
                pcfs_y = pcfs[:, i]
                artificial_df = construct_micro_dataset(copula_y, pcfs_y, grid, tot_hhs)
                artificial_df = DataFrame(artificial_df, ["copula_weights", "grid_point", measures..., m_cols...])

                for (cc, c) in enumerate(itr_combs)
                    # Create new artificial dataset that groups over the variables in c, creates weighted mean
                    grouped_df = combine(groupby(artificial_df, [Symbol.(m * "grid") for m in c]), :copula_weights => nan_sum => :total_copula_weight)

                    # Construct the copula from the grouped data, using the "grid" columns as indices 
                    for fid in axes(cop_2d, 1)
                        for sid in axes(cop_2d, 2)
                            cond = (grouped_df[!, Symbol.(c[1] * "grid")] .== fid) .&& (grouped_df[!, Symbol.(c[2] * "grid")] .== sid)
                            cop_2d[fid, sid] = grouped_df[cond, :total_copula_weight][1]
                        end
                    end

                    # Compute now the tail dependence for the 2D copula
                    result = compute_tail_dependence_2d(cop_2d, p)
                    tail_dep_u[combs[cc]][i, q] = result[1]
                    tail_dep_l[combs[cc]][i, q] = result[2]
                end
            end
        end

        # Export the tail dependence dataframes to excel
        m_label = measures_folder(measures)
        equiv = equivalized == true ? "eq" : ""
        botcod = !isempty(bottom_coded) ? "bc" : ""
        label = "_$case" * "_$equiv" * "$botcod"
        init_path = BASE_PATH
        path = init_path * "/7_Results/$m_label" * tag * "/$folder/data"
        mkpath(path)

        # Export tail_dep_l and tail_dep_u to excel # TODO: to check indices
        half_way_p = floor(Int, grid / 2)
        for (cc, c) in enumerate(combs)
            println(c)
            tail_dep_u[c*"_df"] = DataFrame(tail_dep_u[c][:, half_way_p+1:grid-1], Symbol.(["tdu_$(i)" for i in half_way_p+1:grid-1]))
            tail_dep_l[c*"_df"] = DataFrame(tail_dep_l[c][:, 1:half_way_p], Symbol.(["tdl_$(i)" for i in 1:half_way_p]))

            # Merge the tail dependence dataframes into a single DataFrame
            tail_dep = hcat(tail_dep_l[c*"_df"], tail_dep_u[c*"_df"])
            tail_dep[!, "time"] = q_dates

            # Export 
            CSV.write(path * "/" * data_name * "_" * c * "_tail_dependence" * "$label" * ".csv", select(tail_dep, "time", :))

            # Plot Tail Dependence 
            tail_mat = Matrix(select(tail_dep, Not(:time)))
            seq = collect(45:45:size(tail_mat, 1))
            line_sty = [:solid :dash :dot :dashdot :dashdotdot :solid :dash :dot :dashdot :dashdotdot]
            m_comb = collect(itr_combs)[cc]

            Plots.plot()
            for i in eachindex(seq)
                println(L"%$(q_dates[seq[i]])"[1:5] * "\$")
                # println(q_dates[seq[i]])

                Plots.plot!(
                    pairs,
                    tail_mat[seq[i], :],
                    xlabel=L"\textrm{Lower\,\,Tail\,\,and\,\,Upper\,\,Tail}",
                    ylabel=L"\textrm{Tail\,\,Dependence,\,\, λ_{(%$(m_comb[1]),\, %$(m_comb[2]))}}",
                    xformatter=:latex,
                    yformatter=:latex,
                    label=L"%$(q_dates[seq[i]])"[1:5] * "\$",
                    ls=line_sty[:, i],
                    xticks=(pairs, [L"0.%$(i)" for i in pairs]),
                    legend=:outertopright,
                    lw=2, dpi=500
                )
            end

            # Plots.plot!(1:grid-1, tail_mat[seq[end], :], xformatter=:latex, yformatter=:latex,  label = L"%$(q_dates[seq[end]])"[1:5] * "\$", lw=2, dpi=500)
            Plots.vline!([half_way_p], yformatter=:latex, linestyle=:dash, linecolor=:grey, label="")

            path = init_path * "/7_Results/$m_label" * tag * "/$folder/plots/correlations"
            mkpath(path)
            Plots.savefig(path * "/" * data_name * "_" * c * "_tail_dependence" * "$label" * ".pdf")
        end
    end
end

# Generate a sum over NaNs function 
function nan_sum(x)
    if all(isnan, x)
        return NaN
    else
        return nansum(x)
    end
end

function compute_tail_dependence_2d(copula_hist, u)
    """
    From https://wisostat.uni-koeln.de/fileadmin/sites/statistik/pdf_publikationen/TDCSchmidt.pdf,
    equation 1.4
    and 
    https://escholarship.org/content/qt07x6p3bk/qt07x6p3bk_noSplash_c6faf62d3de3c34d81607b2465a48c15.pdf?t=q8hduy
    equation 6.6
    https://wisostat.uni-koeln.de/fileadmin/sites/statistik/pdf_publikationen/FrahmJunkerSchmidt.pdf
    """
    n = size(copula_hist, 1)
    m = length(size(copula_hist))
    edf = copula_cdf(copula_hist)

    # Compute the tail dependence coefficient
    k = u / 10
    p = 10 - u
    upper_tail = 2 - (1 / k) + (1 / k) * edf[p, p] # (Huang, 1992)
    lower_tail = edf[u, u] / k

    # upper_tail = (1 - m * (u/n) + edf[u, u]) / (1 - (m-1)*(u/n)) # (1/(1000 * z)) * sum(edf[u:end,u:end]) #
    # lower_tail = edf[u, u] / ((m-1)*(u/n))  # (1/z) * edf[u, u] #    
    # Return the tail dependence coefficients
    return upper_tail, lower_tail
end

function copula_cdf(copula_hist)
    local copula_cdf
    # Compute the empirical copula CDF
    if length(size(copula_hist)) == 2
        n = size(copula_hist, 1)
        copula_cdf = zeros(n, n)
        for i in 1:n
            for j in 1:n
                copula_cdf[i, j] = sum(copula_hist[1:i, 1:j])
            end
        end
    elseif length(size(copula_hist)) == 3
        n = size(copula_hist, 1)
        copula_cdf = zeros(n, n, n)
        for i in 1:n
            for j in 1:n
                for k in 1:n
                    copula_cdf[i, j, k] = sum(copula_hist[1:i, 1:j, 1:k])
                end
            end
        end
    end
    return copula_cdf
end



# function compute_tail_dependence_3d(copula_hist, measure_combo, measures, u)
#     """Tail dependence is only for 2 variables, but we take a vector of 3 variables, where estimates are unconditional on the declared 3rd variable."""

#     n   = size(copula_hist, 1)
#     m   = length(size(copula_hist))
#     edf = copula_cdf(copula_hist)

#     # See who is third, measures are always in order 
#     third_var       = setdiff(measures, measure_combo)
#     id_of_third_var = findfirst(measures .== third_var)

#     # Create a vector where all elements are equal to u except for the id of the third variable
#     U = Vector{Int64}(undef, m)
#     U[id_of_third_var] = n
#     U[setdiff(1:m, id_of_third_var)] .= u


#     # Compute the tail dependence coefficient
#     upper_tail = (1 - m * (u/n) + edf[U[1], U[2], U[3]]) / (1 - (m-1)*(u/n)) # (1/(1000 * z)) * sum(edf[u:end,u:end]) #
#     lower_tail = edf[U[1], U[2], U[3]] / ((m-1)*(u/n))  # (1/z) * edf[u, u] #    

#     return upper_tail, lower_tail
# end

function compute_tail_dependence_3d(artificial_df, measure_combo, measures, u)
    """Tail dependence is only for 2 variables, but we take a vector of 3 variables, where estimates are unconditional on the declared 3rd variable."""

    # find the IDs of the measures in measure_combo, where the IDs are the location of the measures in 'measures'
    measure_combo_id = [findfirst(measures .== m) for m in measure_combo]

    # The first column of the df is the copula weights, so we shift indices by 1
    pcfs = DataFrame(artificial_df[:, measure_combo_id.+1], Symbol.(measure_combo))

    # Sort each column, extract unique values (should only be 10), fill quantiles with 0 and assign quantiles based on order of unique values 
    for c in Symbol.(measure_combo)
        sort!(pcfs[:, c])
        unique_vals = sort(unique(pcfs[!, c]))
        pcfs[!, String(c)*"_quantiles"] .= 0

        # assign quantiles 
        for (i, v) in enumerate(unique_vals)
            pcfs[pcfs[!, c].==v, String(c)*"_quantiles"] .= i
        end
    end

    # Conditional probabilities: upper  # TODO: this is conditional on the second measure. Would be good to do reverse while we are here. Results will be similar.  
    numerator_u = sum((pcfs[:, measures[measure_combo_id[1]]*"_quantiles"] .>= u) .& (pcfs[:, measures[measure_combo_id[2]]*"_quantiles"] .>= u))
    denominator_u = sum(pcfs[:, measures[measure_combo_id[2]]*"_quantiles"] .>= u)

    # Conditional Probabilities: lower 
    numerator_l = sum((pcfs[:, measures[measure_combo_id[1]]*"_quantiles"] .< u) .& (pcfs[:, measures[measure_combo_id[2]]*"_quantiles"] .< u))
    denominator_l = sum(pcfs[:, measures[measure_combo_id[2]]*"_quantiles"] .< u)


    # Compute the tail dependence coefficient
    upper_tail = numerator_u / denominator_u
    lower_tail = numerator_l / denominator_l

    return upper_tail, lower_tail
end



function generate_combinations(measures)
    combinations = []

    # Generate combinations of two variables
    for i in eachindex(measures)
        for j in (i+1):length(measures)
            v1 = measures[i]
            v2 = measures[j]
            push!(combinations, [v1, v2])
            push!(combinations, [v2, v1])
        end
    end

    # # Generate combinations of three variables
    # if length(measures) == 3
    #     for i in eachindex(measures)
    #         for j in (i + 1):length(measures)
    #             for k in (j + 1):length(measures)
    #                 v1 = measures[i]
    #                 v2 = measures[j]
    #                 v3 = measures[k]
    #                 push!(combinations, (v1, v2, v3))
    #                 push!(combinations, (v1, v3, v2))
    #                 push!(combinations, (v2, v1, v3))
    #                 push!(combinations, (v2, v3, v1))
    #                 push!(combinations, (v3, v1, v2))
    #                 push!(combinations, (v3, v2, v1))
    #             end
    #         end
    #     end
    # end

    return combinations
end
