# # using CategoricalArrays
# # using CSV 
# # using Coverage
# # using Random
# # using GalacticOptim
# # using DualNumbers
# # using DataFramesMeta 
# # using PyPlot
# # using DataSets
# # using MatrixEquations
# # using SparseArrays
# # using RecursiveArrayTools

# # using StaticArrays
# using Pkg
# Pkg.activate("env_dd_v19/")
# # using BAT 
# # using ValueShapes
# using Optimization 
# using OptimizationOptimJL
# using ReverseDiff
# using ForwardDiff
# using Optim
# using DataFrames
# using FreqTables
# using DualNumbers
# using AdvancedHMC
# using CSV
# using LineSearches
# using Zygote 
# # using DIMESampler
# using Distributed 
# using Printf
# using ProgressBars
# using StatsFuns
# using Interpolations
# using GR
# using Random
# GR.inline("pdf")

# # using GalacticOptim: AutoReverseDiff
# using GenericLinearAlgebra: svd
# using Distributions
# using StatsPlots
# using Plots  
# using MCMCChains
# using HypothesisTests
# using GLM
# using HPFilter

# using ProgressBars
# using LaTeXTabulars
# using Dates 
# using PeriodicalDates
# using KernelDensity
# using Measures
# using MultivariateStats
# using PlotlyJS
# using Parameters 
# using FFTW
# using XLSX
# using StatsBase
# using LaTeXStrings
# using LinearAlgebra
# using BenchmarkTools
# using TimerOutputs
# using Test
# using DelimitedFiles
# using PDMats
# # using MCMCDiagnosticTools
# using JLD2
# using Trapz
# using Combinatorics
# using MatrixEquations
# using FiniteDiff
# using Latexify 
# using BlackBoxOptim
# using KissSmoothing
# using RCall 
# using NaNStatistics
# using JuMP
# using HiGHS

# generate_relative_to_peak_plots()

function generate_relative_to_peak_plots(; tag::AbstractString="")

    # Import recession dates
    init_path  = DATA_PROCESSING
    nber_dates = CSV.read(init_path * "/nber_dates.csv", DataFrame)

    # just the housing cycle 
    filter!(x -> x.peak == nber_dates[!, :peak][7], nber_dates)

    # overwrite peak date 
    nber_dates[!, :peak] .= "August 2008 (2008Q3)"

    # Drop row associated with 1980 recession 
    # filter!(x -> x.peak != "January 1980 (1980Q1)", nber_dates)
    # a = nber_dates[:, :peak][7] # Covid
    # filter!(x -> x.peak != a, nber_dates)

    # From the first column, the current format is like this: December 1969 (1969Q4) --- we only want whats in the parentheses 
    for point in ["peak", "trough"]
        nber_dates[!, point] = [x[end-6:end-1] for x in nber_dates[!, point]]

        # Add a hyphen in character 5
        nber_dates[!, point] = [x[1:4] * "-" * x[5:end] for x in nber_dates[!, point]]

        # Convert to dates, but just year quarter
        nber_dates[!, point] = QuarterlyDate.(nber_dates[!, point])
    end

    # Now I want to create a dictionary of dataframes, where each dataframe is a different recession 
    nber_dates_dict  = Dict()
    
    # Make a key for each group 
    for ext in ["peak"]
        nber_dates_dict[ext] = Dict()
        for (i, d) in enumerate(nber_dates[!, ext])
            nber_dates_dict[ext][d] = DataFrame()
            
            # create a column of dates that is 8 quarters before and 20 quarters after
            nber_dates_dict[ext][d][!, "date"] = d - Quarter(20) : Quarter(1) : d + Quarter(20)
        end
    end

    # pop!(nber_dates_dict["peak"], QuarterlyDate(1981, 3))
    # pop!(nber_dates_dict["peak"], QuarterlyDate(1969, 4))
    # pop!(nber_dates_dict["peak"], QuarterlyDate(1990, 3))
    # pop!(nber_dates_dict["peak"], QuarterlyDate(1973, 4))



    # # Identifying the different kinds of recessions 
    # demand_dates = [QuarterlyDate(1969, 4), QuarterlyDate(1981, 3), QuarterlyDate(2001, 1), QuarterlyDate(2007, 4)]
    # supply_dates = [QuarterlyDate(1973, 4), QuarterlyDate(1990, 3)] # QuarterlyDate(2019, 4)
    # both_dates   = [] #[QuarterlyDate(1980, 1)]
    # recession_types = Dict()

    # for d in nber_dates[!, "peak"]
    #     if d ∈ demand_dates
    #         recession_types[d] = "D"
    #     elseif d ∈ supply_dates
    #         recession_types[d] = "S"
    #     elseif d ∈ both_dates
    #         recession_types[d] = "B"
    #     end
    # end
    sets_of_measures = [["consum", "income", "wealth"]]
    opttag = "from_mcmc"

    for m_set in sets_of_measures
        measures_folder       = m_set[1] * "_and_" * m_set[2] * "_and_" * m_set[3]
        file_path             = BASE_PATH * "/7_Results/$(measures_folder)" * tag * "/$opttag/data/PSID_micro_data_A non-diag_.csv"
        micro_data            = CSV.read(file_path, DataFrame)
        micro_data[!, "time"] = QuarterlyDate.(micro_data[!, "time"])
        plot_path             = BASE_PATH * "/7_Results/consumption_cyclicality" * tag
        mkpath(plot_path)

        # For cop_share, zero out the values that are less than 0 
        # micro_data[!, "cop_share"] = [x < 0 ? 0 : x for x in micro_data[!, "cop_share"]]
        filter!(x -> x.cop_share > 0, micro_data)


        # Define groups, which are described by their grid points
        groups  = Dict("low_$(m_set[2])_low_$(m_set[3])" => [1:10, 1:5, 1:5], "low_$(m_set[2])_high_$(m_set[3])" => [1:10, 1:5, 6:10], "high_$(m_set[2])_low_$(m_set[3])" => [1:10, 6:10, 1:5], "high_$(m_set[2])_high_$(m_set[3])" => [1:10, 6:10, 6:10]) # "low_$(m_set[2])_top_$(m_set[3])" => [1:10, 1:5, 10], "top_$(m_set[2])_low_$(m_set[3])" => [1:10, 10, 1:5], "top_$(m_set[2])_top_$(m_set[3])" => [1:10, 10, 10])
        gr_labs = Dict("low_$(m_set[2])_low_$(m_set[3])" => L"\textrm{Low\,\, %$(m_set[2]),\,\, Low\,\, %$(m_set[3])}", "low_$(m_set[2])_high_$(m_set[3])" => L"\textrm{Low\,\, %$(m_set[2]),\,\, High\,\, %$(m_set[3])}", "high_$(m_set[2])_low_$(m_set[3])" => L"\textrm{High\,\, %$(m_set[2]),\,\, Low\,\, %$(m_set[3])}", "high_$(m_set[2])_high_$(m_set[3])" => L"\textrm{High\,\, %$(m_set[2]),\,\, High\,\, %$(m_set[3])}")
        l_styles = [:dash :solid :dot :dashdot :dashdotdot]
        logocolors = Colors.JULIA_LOGO_COLORS
        l_colors   = [logocolors.blue, logocolors.red, logocolors.green, logocolors.purple]
         

        # Now, merge my data with the nber_dates_dict
        consum_df = Dict()

        sorted_dict = sort(collect(groups), by=x->x[2])

        for ext in collect(keys(nber_dates_dict))
            # Plots.plot()
            p1 = Plots.plot()
            p2 = Plots.plot()
            p3 = Plots.plot()

            qq = 0
            rec_min = []
            for (n, g) in sorted_dict
                qq += 1
                consum_df[ext] = DataFrame()
                println(n)

                # filter data based on the grid points
                s_micro_data = micro_data[findall(x -> x ∈ g[1], micro_data[!, "consumgrid"]), :]
                s_micro_data = s_micro_data[findall(x -> x ∈ g[2], s_micro_data[!, "$(m_set[2])grid"]), :]
                s_micro_data = s_micro_data[findall(x -> x ∈ g[3], s_micro_data[!, "$(m_set[3])grid"]), :]
            
                for (i, d) in nber_dates_dict[ext]
                    # Subset to the dates of interest, but keeping all columns
                    ss_micro_data = filter(row -> row.time .<= i + Quarter(20), s_micro_data) # 
                    ss_micro_data = filter(row -> row.time .>= i - Quarter(20), ss_micro_data)

                    # correct dates for plotting, making the recession date zero
                    ss_micro_data[!, "actual_time"] = ss_micro_data[!, "time"] 

                    append!(consum_df[ext], ss_micro_data)
                end
                
                # Create 1 df for consumption, 1 for weights 
                n_recessions = length(nber_dates_dict[ext])
                final_df_w   = combine(groupby(consum_df[ext], [:actual_time]), [:cop_share] => (x) -> sum(x) ./ n_recessions)
                final_df     = combine(groupby(consum_df[ext], [:actual_time]), [:consum, :cop_share] => (x,y) -> mean(x, weights(y)))

                # Consumption first  
                psid_obs_dates           = QuarterlyDate(2003, 2) : Quarter(8) : QuarterlyDate(2021, 2)
                recession_dates          = QuarterlyDate(2007, 4) : Quarter(1) :  QuarterlyDate(2009, 2)
                cond_psid                = findall(x -> x ∈ collect(psid_obs_dates), final_df[!, :actual_time])
                cond_recession           = findall(x -> x ∈ collect(recession_dates), final_df[!, :actual_time])
                final_df_w[!, "time"]    = final_df_w[!, "actual_time"] .- QuarterlyDate(2008, 3) #TODO: hard-coded 
                final_df[!, "time"]      = final_df[!, "actual_time"] .- QuarterlyDate(2008, 3) #TODO: hard-coded 
                peak_cond                = Dates.value.(final_df[:, "time"]) .== 0
                consum_peak              = final_df[:, "consum_cop_share_function"][peak_cond]
                consum_peak_rounded      = Int.(round.(consum_peak, digits=0))
                final_df[:, :consum_rel] = final_df[:, "consum_cop_share_function"] ./ consum_peak

                # Detrended version stuff 
                trend          = HP(final_df[:, "consum_cop_share_function"], 1600)
                detrended_data = final_df[:, "consum_cop_share_function"] .- trend

                # Make relative to trough 
                detrended_data .+= mean(final_df[:, "consum_cop_share_function"])
                consum_peak_detrended = detrended_data[peak_cond]
                final_df[!, "dt_consum_rel"] = (detrended_data ./ consum_peak_detrended) # in % terms

                # Now, weights
                peak_cond_w               = Dates.value.(final_df_w[:, "time"]) .== 0
                share_peak_w              = final_df_w[:, "cop_share_function"][peak_cond_w]
                share_peak_rounded_w      = round.(share_peak_w, digits=2)
                final_df_w[:, :share_rel] = final_df_w[:, :cop_share_function] ./ share_peak_w

                # Consumption First  
                Plots.plot!(p1, axes(final_df[:, :actual_time]), 
                final_df[:, :consum_rel], 
                label="", #gr_labs[n], 
                lw=4, 
                ls=l_styles[qq],
                lc=l_colors[qq],
                legend=:outerbottom,
                legend_columns=2,
                legendfontsize=10,
                xtickfontsize=10,
                ytickfontsize=10,
                guidefontsize=14,
                xformatter=:latex, 
                yformatter=:latex, 
                ylabel=L"\textrm{Consumption\,\, rel.\,\, to\,\, %$(ext)}", 
                xlabel=L"\textrm{PSID\,\, observation\,\,date}", 
                xticks=(cond_psid, [L"%$(date)"[1:5] * "\$" for (_,date) in enumerate(final_df[cond_psid, :actual_time])]),
                )
                Plots.plot!([], [], ls=l_styles[qq], lc=l_colors[qq], label=gr_labs[n], lw=1) 

                # Plots.plot!(p1, [], [], gridlinewidth=4, gridstyle=:dash, grid=:x, gridalpha=0.2, label="")
                Plots.plot!(p1, [], [], gridlinewidth=2, gridstyle=:dash, grid=:x, gridalpha=0.2, label="")

                append!(rec_min, [nanminimum(final_df[:, :consum_rel])])
                # # Plotting the confidence intervals 
                # Plots.plot!(p1, collect(axes(final_df[:, :actual_time])...)[cond_recession],
                # [minimum(final_df[:, :consum_rel])  for _ in eachindex(cond_recession)],
                # fillrange = ones(length(cond_recession)),
                # fillalpha = 0.2,
                # fillcolor = :gray,
                # la=0.0,
                # lc=:gray,
                # dpi=500,
                # label="",
                # )
                
                # Plots.bar!(p1, collect(axes(final_df[:, :actual_time])...)[cond_recession], ones(length(cond_recession)), color = "grey", alpha = 0.3, label = "")
                

                # Consumption Detrended   
                Plots.plot!(p3, axes(final_df[:, :actual_time]), 
                final_df[!, "dt_consum_rel"], 
                label=gr_labs[n], 
                lw=4, 
                ls=l_styles[qq],
                legend=:outerbottom,
                legend_columns=2,
                legendfontsize=10,
                xtickfontsize=10,
                ytickfontsize=10,
                guidefontsize=14,
                xformatter=:latex, 
                yformatter=:latex, 
                ylabel=L"\textrm{Consumption\,\, rel.\,\, to\,\, %$(ext)}", 
                xlabel=L"\textrm{PSID\,\, observation\,\,date}", 
                xticks=(cond_psid, [L"%$(date)"[1:5] * "\$" for (_,date) in enumerate(final_df[cond_psid, :actual_time])]),
                )
                Plots.plot!(p3, [], [], gridlinewidth=4, gridstyle=:dot, grid=:x, gridalpha=0.1, label="")

                Plots.plot!(p3, collect(axes(final_df[:, :actual_time])...)[cond_recession],
                [minimum(final_df[:, :dt_consum_rel]) for _ in eachindex(cond_recession)],
                fillrange = ones(length(cond_recession)),
                fillalpha = 0.1,
                fillcolor = :gray,
                la=0.0,
                lc=:white,
                lw=4, dpi=500,
                label="",
                )

                # Weights 
                if n != "overall"
                    Plots.plot!(p2, final_df_w[:, :time], 
                    final_df_w[:, :share_rel], 
                    label=gr_labs[n], 
                    lw=2, 
                    ls=l_styles[qq],
                    xformatter=:latex, 
                    yformatter=:latex, 
                    legend=:outertopright,
                    legendfontsize=7,
                    ylabel=L"\textrm{Share\,\, rel.\,\, to\,\, %$(ext)}", 
                    xlabel=L"\textrm{quarters\,\, from \,\, %$(ext)}", 
                    xticks=(-7:1:7, [L"%$(f)" for f in -7:1:7])
                    )

                    Plots.plot!(p2, final_df_w[:, :time], 
                    final_df_w[:, :share_rel], 
                    la=0.0,
                    label=L"\textrm{At\,\, peak: \, %$(share_peak_rounded_w[1])}",
                    )
                end

                if qq == length(groups)
                    println(rec_min)
                    # Plotting the confidence intervals 
                    Plots.plot!(p1, collect(axes(final_df[:, :actual_time])...)[cond_recession],
                    [minimum(rec_min) for _ in eachindex(cond_recession)],
                    fillrange = ones(length(cond_recession)),
                    fillalpha = 0.4,
                    fillcolor = :gray,
                    ylimits=(minimum(rec_min), 1.12),
                    la=0.0,
                    lc=:gray,
                    dpi=500,
                    label="",
                    )

                    Plots.savefig(p1, plot_path * "/$measures_folder" * "_consumption_rel_$ext.pdf")
                    Plots.savefig(p2, plot_path * "/$measures_folder" * "_weight_rel_$ext.pdf")
                    Plots.savefig(p3, plot_path * "/$measures_folder" * "_detrended_consumption_rel_$ext.pdf")
                end
            end
        end
    end
end


# # Now, merge my data with the nber_dates_dict
# param_dict = Dict()

# p1 = Plots.plot()
# p2 = Plots.plot()

# # Now, merge my data with the nber_dates_dict
# consum_df = Dict()

# for m_set in sets_of_measures 
#     measures_folder = m_set[1] * "_and_" * m_set[2] * "_and_" * m_set[3]
#     file_path = dirname(pwd())[end-7:end] == "Dynamics" ? "/Users/lc/Dropbox/Distributional_Dynamics/7_Results/$(measures_folder)/from_optimization/data/PSID_micro_data_A non-diag_.csv" : init_path * "/PSID_micro_data_A non-diag_.csv"
#     micro_data = CSV.read(file_path, DataFrame)
#     micro_data[!, "time"] = QuarterlyDate.(micro_data[!, "time"])

#     # For cop_share, zero out the values that are less than 0 
#     micro_data[!, "cop_share"] = [x < 0 ? 0 : x for x in micro_data[!, "cop_share"]]


#     # Define groups, which are described by their grid points
#     groups  = Dict("overall" => [1:10, 1:10, 1:10], "low_$(m_set[2])_low_$(m_set[3])" => [1:10, 1:5, 1:5], "low_$(m_set[2])_high_$(m_set[3])" => [1:10, 1:5, 6:10], "high_$(m_set[2])_low_$(m_set[3])" => [1:10, 6:10, 1:5], "high_$(m_set[2])_high_$(m_set[3])" => [1:10, 6:10, 6:10]) # "low_$(m_set[2])_top_$(m_set[3])" => [1:10, 1:5, 10], "top_$(m_set[2])_low_$(m_set[3])" => [1:10, 10, 1:5], "top_$(m_set[2])_top_$(m_set[3])" => [1:10, 10, 10])
#     gr_labs = Dict("overall" => L"\textrm{Overall}", "low_$(m_set[2])_low_$(m_set[3])" => L"\textrm{Low\,\, %$(m_set[2]),\,\, Low\,\, %$(m_set[3])}", "low_$(m_set[2])_high_$(m_set[3])" => L"\textrm{Low\,\, %$(m_set[2]),\,\, High\,\, %$(m_set[3])}", "high_$(m_set[2])_low_$(m_set[3])" => L"\textrm{High\,\, %$(m_set[2]),\,\, Low\,\, %$(m_set[3])}", "high_$(m_set[2])_high_$(m_set[3])" => L"\textrm{High\,\, %$(m_set[2]),\,\, High\,\, %$(m_set[3])}")
#     l_styles = [:dash :solid :dot :dashdot :dashdotdot]

#     # Now, merge my data with the nber_dates_dict
#     consum_df = Dict()

#     sorted_dict = sort(collect(groups), by=x->x[2])
#     for _ in 1:1
#         qq = 0 
#         for (n, g) in sorted_dict
#             qq += 1
#             if n == "overall"
#                 nothing 
#             else
#                 consum_df        = DataFrame()

#                 println(n)

#                 # filter data based on the grid points
#                 s_micro_data = micro_data[findall(x -> x ∈ g[1], micro_data[!, "consumgrid"]), :]
#                 s_micro_data = s_micro_data[findall(x -> x ∈ g[2], s_micro_data[!, "$(m_set[2])grid"]), :]
#                 s_micro_data = s_micro_data[findall(x -> x ∈ g[3], s_micro_data[!, "$(m_set[3])grid"]), :]
                
#                 # Now we want all observations around the recessions/booms 
#                 for (i, d) in nber_dates_dict["peak"]
#                     # Subset to the dates of interest, but keeping all columns
#                     ss_micro_data = filter(row -> row.time .<= i + Quarter(8), s_micro_data) # 
#                     ss_micro_data = filter(row -> row.time .>= i - Quarter(8), ss_micro_data)        

#                     # correct dates for plotting, making the recession date zero
#                     ss_micro_data[!, "time"] = ss_micro_data[!, "time"] .- i

#                     append!(consum_df, ss_micro_data)
#                 end

#                 # Now we perform the optimization 
#                 peak_cond = Dates.value.(consum_df[:, "time"]) .== 0
#                 X     = [Matrix(consum_df[Dates.value.(consum_df[:, "time"]) .== i, [:income, :wealth]]) for i in -8:8] # 0 is filtered out later
#                 ω̄     = consum_df[peak_cond, :cop_share][:,:]
#                 ω̄     = reshape(ω̄, (1, size(ω̄, 1)))

#                 # Define model parameters 
#                 model       = Model(HiGHS.Optimizer)
#                 m           = length(ω̄)
#                 T           = size(unique(consum_df[:, :time]), 1)

#                 set_attribute(model, "time_limit", 1800.0)

#                 # Define omega_t variables for each t in T, excluding the peak
#                 @variable(model, 0 <= ω[t=filter(t -> t != 9, 1:T), j=1:m] <= 1)
#                 # @variable(model, 0 .<= ω[t=filter(t -> t != 9, 1:T)] .<= 1)

#                 # Objective function
#                 @objective(model, Min, 0.5 * sum((ω[t, j] - ω̄[j])^2 for t in filter(t -> t != 9, 1:T), j in 1:m))
#                 # @objective(model, Min, 0.5 * sum((ω[t] - ω̄[j])^2 for t in filter(t -> t != 9, 1:T)))

#                 # Vector equality constraints for omega_t' X = omega_peak' X
#                 # for i in axes(X, 2)  # For each measure in X
#                     for t in filter(t -> t != 9, 1:T)
#                         @constraint(model, sum(ω[t, j] * X[t][j,:] for j in 1:m) == sum(ω̄[j] * X[9][j,:] for j in 1:m))
#                     end
#                 # end
                
#                 # for t in filter(t -> t != 9, 1:T)
#                 #     @constraint(model, ω[t] * X[t] == ω̄ * X[9])
#                 # end
                

#                 # Weights in period must sum to weights in peak 
#                 for t in filter(t -> t != 9, 1:T)
#                     @constraint(model, sum(ω[t, j] for j in 1:m) == sum(ω̄[j] for j in 1:m))
#                 end

#                 optimize!(model)
#                 new_weights = Matrix(value.(ω))
#                 param_dict[n] = new_weights
#                 # weights_cf = jldopen("/Users/lc/Dropbox/Distributional_Dynamics/5_Code/param_dict_counterfactual.jld2", "r")
#                 # weights_cf = weights_cf["cf"]
#                 # new_weights = Matrix(weights_cf[n])

#                 # Allocate weights back into df 
#                 sort!(consum_df, :time)
#                 t, k                     = size(new_weights)
#                 half_t                   = Int(t/2)
#                 consum_df[!, :cop_share] = vcat(vec(reshape(new_weights[1:half_t, :], (1, k * half_t))), vec(ω̄), vec(reshape(new_weights[half_t+1:end, :], (1, k * half_t))))

#                 # Generate weighted averages 
#                 E_C     = combine(groupby(consum_df, [:time]), [:consum, :cop_share] => (x,y) -> mean(x, weights(y)))

#                 # Make relative to peak/trough 
#                 peak_cond            = Dates.value.(E_C[:, "time"]) .== 0
#                 consum_peak          = E_C[:, "consum_cop_share_function"][peak_cond]
#                 consum_peak_rounded          = round.(consum_peak, digits=2)
#                 E_C[!, "consum_rel"] = round.(E_C[:, "consum_cop_share_function"] ./ consum_peak, digits=6)

#                 # Normal version first
#                 Plots.plot!(p1, E_C[:, :time], 
#                 E_C[:, :consum_rel], 
#                 label=gr_labs[n], 
#                 lw=2, 
#                 xformatter=:latex, 
#                 yformatter=:latex, 
#                 legend=:outertopright,
#                 ylabel=L"\textrm{Consumption\,\, rel.\,\, to\,\, peak}", 
#                 xlabel=L"\textrm{quarters\,\, from \,\, peak}", 
#                 xticks=(-7:1:7, [L"%$(f)" for f in -7:1:7])
#                 )

#                 Plots.plot!(p1, E_C[:, :time], 
#                 E_C[:, :consum_rel], 
#                 la=0.0,
#                 label=L"\textrm{At\,\, peak: \, %$(consum_peak_rounded[1])}",
#                 )

#                 if qq == length(groups)
#                     Plots.savefig(p1, "Counterfactual.pdf")
#                 end

#                 # detrend 
#                 # df = DataFrame(Y = E_C[:, "consum_cop_share_function"], X = 1:length(unique(E_C[:, "time"])))
                
#                 # Fit linear model
#                 # model = lm(@formula(Y ~ X), df)

#                 # Detrend by subtracting fitted values from the original data
#                 # detrended_data = residuals(model)
#                 # detrended_data = E_C[:, "consum_cop_share_function"] .- HP(E_C[:, "consum_cop_share_function"], 100)
#                 # detrended_data .+= mean(E_C[:, "consum_cop_share_function"])

#                 # consum_peak_detrended = round.(detrended_data[peak_cond], digits=6)
#                 # E_C[!, "consum_rel"]  = (detrended_data ./ consum_peak_detrended) # in % terms

#                 # Plots.plot!(p2, E_C[:, :time], 
#                 # E_C[:, :consum_rel], 
#                 # label=gr_labs[n], 
#                 # lw=2, 
#                 # xformatter=:latex, 
#                 # yformatter=:latex, 
#                 # legend=:outertopright,
#                 # ylabel=L"\textrm{Consumption\,\, rel.\,\, to\,\, peak}", 
#                 # xlabel=L"\textrm{quarters\,\, from \,\, peak}", 
#                 # xticks=(-7:1:7, [L"%$(f)" for f in -7:1:7])
#                 # )

#                 # Plots.plot!(p2, E_C[:, :time], 
#                 # E_C[:, :consum_rel], 
#                 # la=0.0,
#                 # label=L"\textrm{At\,\, peak: \, %$(consum_peak[1])}",
#                 # )
#                 # if qq == length(groups)
#                 #     Plots.savefig(p2, "CounterfactualDetrended.pdf")
#                 # end
#             end
#         end
#     end
# end
