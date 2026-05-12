#TODO: check that this is correct. I think the grid_choice_pcf condition is for 20, not 5
# Next 40
# local qq 
# if grid_choice_pcf == 5
#     qq = "q"
# else
#     qq = ""
# end


function compare_to_external_sources(dv, ty, func_data, obs_data, user_params, time_params, model_options, type, label, data_bounds=false)
    # Comparing our estimates with the observed data, both survey anfunction compare_to_external_sources(dv, func_data, user_params, time_params, method_options, type, label, data_bounds=false)
    @unpack measures, lags, case, estimator, tag = model_options
    @unpack func_dict, year_vec = func_data
    @unpack gdp_series = obs_data
    smin, smax = user_params

    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid, integral_cop_grid = estimator
    else
        @unpack grid_pcf, grid_cop = estimator
    end
    grid_choice_pcf = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf
    grid_choice_cop = typeof(estimator) <: SeriesEstimator ? integral_cop_grid : grid_cop


    # Path situation + other
    data_sources = setdiff(sort(collect(keys(dv))), ["consensus"])
    m_label = measures_folder(measures)
    folder = type == :optimization ? "from_optimization" : "from_mcmc"
    init_path = BASE_PATH
    path = init_path * "/7_Results/$m_label" * "$tag" * "/$folder/plots/"
    mkpath(path)
    series_t = ["top10", "next40", "bottom50"]

    # Creating and filling dictionaries to loop over 
    top = Dict("consensus" => Dict())
    bot = Dict("consensus" => Dict())
    mid = Dict("consensus" => Dict())

    # Defining the series 
    series = Dict()
    if grid_choice_pcf == 10 || grid_choice_pcf == 20 || grid_choice_pcf == 100
        series["bot"] = "bottom50"
        series["mid"] = "next40"
        series["top"] = "top10"
    elseif grid_choice_pcf == 5
        series["bot"] = "bottom40"
        series["mid"] = "next40"
        series["top"] = "top20"
    end

    # For xaxis 
    dates = Dict()

    # External sources
    external_sources = retrieve_external_sources(measures)

    q_ext_s = ["DFA", "WIDq"]
    agr_dates = Dict()
    dist_dict = Dict("top" => top, "bot" => bot, "mid" => mid)
    dist_label = Dict("top" => "top 10", "bot" => "bottom 50", "mid" => "next 40")
    all_recons = [data_sources..., "consensus"]


    # Everything will relative to the average 
    gdp_series[!, "date"] = QuarterlyDate.(gdp_series[!, "time"])
    correction_names = [meas * "_per_hh" for meas in measures]
    select_series = select(gdp_series, [correction_names..., "date"])

    compare_to_DFA(deepcopy(dv), ty, time_params, type, case, measures, func_dict, gdp_series, tag)

    # Second, the other sources 
    for obj in ["levels", "quantiles", "shares"]
        # Create aggreable dates dict per object 
        agr_dates[obj] = Dict()

        # Create Approximation dictionary 
        for dt in [top, bot, mid]
            dt[obj] = Dict()
        end

        for meas in measures
            for dt in [top, bot, mid]
                dt[obj][meas] = Dict()
                dt[obj][meas]["lb"] = Dict()
                dt[obj][meas]["ub"] = Dict()
            end

            agr_dates[obj][meas] = Dict()

            for (ext_s, df) in external_sources[obj][meas]
                # Get dates from external source 
                dates = Dict()
                dropmissing!(df)

                # read in CSV 
                # a = CSV.read("/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/validation/WIDq/income_quantiles_WIDq.csv", DataFrame)
                # Date.(a[!, :time], "yyyyq")
                year_min = ext_s ∉ q_ext_s ? minimum(df[!, "year"]) : minimum(Dates.year.(df[:, "time"]))
                year_max = ext_s ∉ q_ext_s ? maximum(df[!, "year"]) : maximum(Dates.year.(df[:, "time"]))

                if ext_s ∉ q_ext_s
                    # Dates of the external sources  
                    # dates["tmin"]    = Dict("year" => year_min, "quarter" => minimum(df[df[!, "year"] .== year_min, "quarter"]))
                    # dates["tmax"]    = Dict("year" => year_max, "quarter" => maximum(df[df[!, "year"] .== year_max, "quarter"])) 
                    df[!, "quarter"] .= 4
                    dates["tmin"] = Dict("year" => year_min, "quarter" => 4)
                    dates["tmax"] = Dict("year" => year_max, "quarter" => 4)
                else
                    df[!, "year"] = Dates.year.(df[!, "time"])
                    df[!, "quarter"] = Dates.quarterofyear.(df[!, "time"])
                    dates["tmin"] = Dict("year" => year_min, "quarter" => minimum(df[df[!, "year"].==year_min, "quarter"]))
                    dates["tmax"] = Dict("year" => year_max, "quarter" => maximum(df[df[!, "year"].==year_max, "quarter"]))
                end

                df[!, :date] = QuarterlyDate.(df[!, "year"], df[!, "quarter"])

                # Figure out agreeable dates between user specified time and external source dates 
                ext_freqtype = ext_s ∉ q_ext_s ? :annual : :quarterly
                agr_dates[obj][meas][ext_s] = find_agreeable_dates(dates, ext_freqtype, smin, smax)

                # Filter by year and quarter 
                filter!(row -> row.date >= QuarterlyDate(agr_dates[obj][meas][ext_s]["tmin"]["year"], agr_dates[obj][meas][ext_s]["tmin"]["quarter"]), df)
                filter!(row -> row.date <= QuarterlyDate(agr_dates[obj][meas][ext_s]["tmax"]["year"], agr_dates[obj][meas][ext_s]["tmax"]["quarter"]), df)
                col_names = setdiff(names(df), ["year", "quarter", "date", "time"])

                # Subset and log-transform of the data first
                dates_dict = Dict("tmin" => smin, "tmax" => smax)

                # Filter the correction series to dates dict for data 
                sub_select_series = filter(row -> row.date >= QuarterlyDate(dates_dict["tmin"]["year"], dates_dict["tmin"]["quarter"]), select_series)
                filter!(row -> row.date <= QuarterlyDate(dates_dict["tmax"]["year"], dates_dict["tmax"]["quarter"]), sub_select_series)

                for r in all_recons
                    for (k, v) in dist_dict
                        # Only subsets to smin and smax 
                        v[obj][meas][r] = subset_to_cutoff(dv[r][ty][meas][obj]["common series"][series[k]], dates_dict, tmin, tmax, obj)
                        # if ty == "average"
                        #     v[obj][meas][r] = log_transformation(v[obj][meas][r])
                        # end
                    end
                end

                # Filter the correction series to dates dict for external source
                sub_select_series = filter(row -> row.date >= QuarterlyDate(agr_dates[obj][meas][ext_s]["tmin"]["year"], agr_dates[obj][meas][ext_s]["tmin"]["quarter"]), select_series)
                filter!(row -> row.date <= QuarterlyDate(agr_dates[obj][meas][ext_s]["tmax"]["year"], agr_dates[obj][meas][ext_s]["tmax"]["quarter"]), sub_select_series)

                # Since the external sources are mostly annual, we need to average the data over the quarters
                if ext_freqtype == :annual
                    sub_select_series = average_out_quarters(sub_select_series, ext_s, meas)
                end

                # Log transform of the External Source 
                if obj == "levels"
                    for c in col_names
                        raw_data = ext_s ∉ q_ext_s ? df[!, c] : df[!, c] .* 1000000
                        df[!, c] = log_transformation(raw_data)
                    end

                elseif obj == "quantiles"
                    for c in col_names
                        if ext_s == "DFA" && c == "bottom50"
                            df[!, c] .+= abs(minimum(df[!, c]) * 2) # Since it goes negative very quickly
                        end
                        df[!, c] = log_transformation(df[!, c])
                    end

                elseif obj == "shares"
                    for c in col_names
                        raw_data = ext_s ∉ q_ext_s ? df[!, c] : df[!, c] ./ 100
                        df[!, c] = raw_data
                    end
                end
            end
        end
    end

    # Correlation matrix
    corr_df = DataFrame(rand(1, 1), :auto)
    cycle_dict = Dict("cycle" => :cycle, "raw" => :raw)

    # Estimation dates (user specified)
    dts = QuarterlyDate(smin["year"], smin["quarter"]):Quarter(1):QuarterlyDate(smax["year"], smax["quarter"])
    xa = collect(1:length(dts))

    plot_tag = ty == "normal" ? "" : "_detrended_"

    # Get correlations between external sources
    ext_corr = Dict()

    # External source with Approximation 
    objects = collect(keys(external_sources))

    for obj in objects
        ext_corr[obj] = Dict()
        for meas in measures
            M = uppercasefirst(meas)
            ext_corr[obj][meas] = Dict()

            for (series_type, choice) in cycle_dict
                z = 1
                ext_corr[obj][meas][series_type] = Dict()
                for (seg, dt) in dist_dict
                    ext_corr[obj][meas][series_type][seg] = Dict()

                    for (j, source) in enumerate(all_recons) #data_sources
                        # First, check whether the plot is worth making 
                        possible_ext_sources = collect(keys(external_sources[obj][meas]))

                        if length(possible_ext_sources) == 0
                            nothing
                        else
                            # First, get the data, which corresponds to (smin, smax) -> use obj, meas to get agr_dates
                            not_nan = .!isnan.(dt[obj][meas][source])
                            cycle_data = fill(NaN, length(not_nan))

                            if sum(not_nan) > 0
                                cycle_data[not_nan] = get_raw_or_cycle(dt[obj][meas][source][not_nan], ty, false, :quarterly, choice, obj, meas)
                            end

                            # Generate vector of opacities: one 100% and the rest 0%
                            opacities = [1.0, 0.0]

                            for xx in 1:2
                                # The approximation  
                                Plots.plot(xa,
                                    cycle_data,
                                    xformatter=:latex,
                                    yformatter=:latex,
                                    ylabel=L"\textrm{Cyclical\,\,Component}",
                                    xticks=(xa[1:40:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(dts[1:40:end])]),
                                    label="",#label=L"\textrm{Model}",
                                    lc=:red,
                                    la=opacities[xx],
                                    xtickfontsize=14,
                                    ytickfontsize=14,
                                    legendfontsize=10,
                                    guidefontsize=14,
                                    legend=false,
                                    lw=4, dpi=500
                                )

                                # Create data vector which stores the ext. source data so we can create correlations for the respective plot 
                                ext_source_data = Dict()

                                # Second, plot the external sources 
                                for (i, ext_s) in enumerate(possible_ext_sources)
                                    ext_freqtype = ext_s ∉ q_ext_s ? :annual : :quarterly
                                    dat_freqtype = ext_freqtype == :annual ? 4 : 1
                                    ext_data = get_raw_or_cycle(external_sources[obj][meas][ext_s][:, series[seg]], ty, true, ext_freqtype, choice, obj, meas)

                                    # Subset the x_axis and data to the same dates as the external source
                                    x_axis = subset_to_cutoff(xa, agr_dates[obj][meas][ext_s], smin, smax, "dates")
                                    data_sub = subset_to_cutoff(dt[obj][meas][source], agr_dates[obj][meas][ext_s], smin, smax, "no log") #TODO: check if it overwrites 

                                    # Average over the quarters to get an annual estimate, if the external source is annual
                                    if ext_freqtype == :annual
                                        # Convert Matrix to Dataframe and add dates column
                                        data_sub = DataFrame(data_sub[:, :], :auto)
                                        data_sub[!, :dates] = QuarterlyDate(agr_dates[obj][meas][ext_s]["tmin"]["year"], agr_dates[obj][meas][ext_s]["tmin"]["quarter"]):Quarter(1):QuarterlyDate.(agr_dates[obj][meas][ext_s]["tmax"]["year"], agr_dates[obj][meas][ext_s]["tmax"]["quarter"])
                                        data_sub[!, :year] = Dates.year.(data_sub[!, :dates])
                                        data_sub[!, :quarter] = Dates.quarterofyear.(data_sub[!, :dates])

                                        # Take the average by year 
                                        data_sub = combine(groupby(data_sub, [:year]), ["x1"] => mean) # :auto generates :x1

                                        # Convert to vector 
                                        data_sub = vec(data_sub[:, :x1_mean])

                                    end

                                    # Extract correlation from the Approximation and the external source
                                    cond = .!isnan.(data_sub)
                                    get_correlation = sum(cond)
                                    local ρ
                                    if get_correlation > 0
                                        cc = get_raw_or_cycle(data_sub[cond], ty, false, ext_freqtype, choice, obj, meas)
                                        ρ = round(cor(cc, ext_data[cond]), digits=2)
                                    else
                                        ρ = NaN
                                    end
                                    col = source * "_" * ext_s * "_" * obj * "_" * series_type * "_" * series_t[z] * "_$meas"
                                    corr_df[!, col] = [ρ]

                                    # Store the external data 
                                    dates_of_ext_s = ext_freqtype == :annual ? [QuarterlyDate(agr_dates[obj][meas][ext_s]["tmin"]["year"], agr_dates[obj][meas][ext_s]["tmin"]["quarter"]):Year(1):QuarterlyDate(agr_dates[obj][meas][ext_s]["tmax"]["year"], agr_dates[obj][meas][ext_s]["tmax"]["quarter"])] : [QuarterlyDate(agr_dates[obj][meas][ext_s]["tmin"]["year"], agr_dates[obj][meas][ext_s]["tmin"]["quarter"]):Quarter(1):QuarterlyDate.(agr_dates[obj][meas][ext_s]["tmax"]["year"], agr_dates[obj][meas][ext_s]["tmax"]["quarter"])]
                                    ext_source_data[ext_s] = DataFrame(hcat(dates_of_ext_s..., ext_data), [:dates, Symbol(ext_s)])

                                    # Plot the external sources first and then the data ... 
                                    # if ext_s == "DFA" && seg == "bot"
                                    #     # Replace negative values with NaN
                                    #     ext_data = [x < 0 ? NaN : x for x in ext_data]
                                    #     ext_data = [x > 0.5 ? NaN : x for x in ext_data]
                                    #     println(ext_data)
                                    # end
                                    Plots.plot!(x_axis[1:dat_freqtype:end],
                                        ext_data,
                                        lc=select_color(ext_s),
                                        ls=select_linestyle(ext_s),
                                        label="",
                                        lw=4, dpi=500
                                    )


                                    Plots.plot!([], [], ls=select_linestyle(ext_s), lc=select_color(ext_s), lw=2, label=L"\textrm{%$(ext_s)}")

                                    # Once we've gone through all the datatsets
                                    if i == length(possible_ext_sources) && xx == 2 && j == length(all_recons)
                                        # Define Combinations 
                                        ext_combs = collect(combinations(possible_ext_sources, 2))

                                        if length(ext_combs) != 0
                                            # Merge the datasets in 'ext_source_data' together based on dates 
                                            merged_df = DataFrame()

                                            # Iterate through the dictionary and merge data frames
                                            for (key, df) in ext_source_data
                                                if isempty(merged_df)
                                                    merged_df = df
                                                else
                                                    merged_df = outerjoin(merged_df, df, on=:dates)
                                                end
                                            end

                                            # Now run correlations 
                                            for (i, comb) in enumerate(ext_combs)
                                                merged_df_sub = deepcopy(merged_df)
                                                # First, drop rows which are NaN for each column 
                                                for c in comb
                                                    # Convert missings to NaN 
                                                    replace!(merged_df_sub[!, c], missing => NaN)
                                                    merged_df_sub = merged_df_sub[.!isnan.(merged_df_sub[!, c]), :]
                                                end
                                                cor_label = join(comb, "-")
                                                ext_corr[obj][meas][series_type][seg][cor_label] = round(cor(merged_df_sub[:, Symbol(comb[1])], merged_df_sub[:, Symbol(comb[2])]), digits=2)
                                            end
                                        end
                                    end
                                end
                                Plots.savefig(path * "/$meas/" * "/external_comparisons/" * meas * "_" * obj * "_$source" * "_" * series_type * "_" * series[seg] * "_" * plot_tag * label * "_$(opacities[xx])" * ".pdf")
                            end
                        end
                    end
                    z += 1
                end
            end
        end
    end

    #Save correlations related to external sources
    jldsave(path * "/correlations/" * "external_correlations.jld2"; ext_corr=ext_corr)

    # Regarding correlation table 
    select!(corr_df, Not(:x1))
    select!(corr_df, sort(names(corr_df)))
    CSV.write(path * "/correlations/correlations" * plot_tag * ".csv", corr_df)


    # x_objects = join.(split.(names(corr_df), "_"), "\n ")
    # Plots.plot(
    #     axes(Matrix(corr_df), 2),
    #     Matrix(corr_df)[1, :],
    #     xlabel = L"\textrm{Series}",
    #     ylabel = L"\textrm{Correlation\,\,with\,\,%$(ext_s)}",
    #     fontfamily="Computer Modern",
    #     lw=4,
    #     xticks=(axes(Matrix(corr_df), 2), [object for object in x_objects]),
    #     legend=false
    #     )

    # Plots.savefig(path * "/correlations/" * "$ext_s" * "_correlations.pdf")
end


function average_out_quarters(series, ext_s, meas)
    series[!, :year] = Dates.year.(series[!, :date])

    # Take the average by year 
    local data_sub
    if ext_s == "WID" && meas == "income"
        # Only save quarter 4
        series[!, :quarter] = Dates.quarterofyear.(series[!, :date])
        data_sub = filter(row -> row.quarter == 4, series)
    else
        data_sub = combine(groupby(series, [:year]), [meas * "_per_hh"] => mean) # :auto generates :x1

        # rename mean column 
        rename!(data_sub, Symbol(meas * "_per_hh_mean") => meas * "_per_hh")

        # data_sub = vec(data_sub[:, meas * "_per_hh" * "_mean"])
    end
    return data_sub
end

function generate_wealth_by_income_df(source, ty, measures, opttag, gdp_series, tag)
    file_tag = ty == "normal" ? "" : "_detrended"
    # Import data 
    meas_folder = measures_folder(measures)
    init_path = BASE_PATH
    micro_df = CSV.read(init_path * "/7_Results/$(meas_folder)" * "$tag" * "/$(opttag)/data/" * "$(source)_micro_data" * file_tag * "_A non-diag_.csv", DataFrame)
    micro_df[!, "time"] = QuarterlyDate.(micro_df[!, "time"])

    # Create wealth by income. But first, define groups 
    micro_df[!, "income_groups"] = [micro_df[i, "incomegrid"] <= 4 ? "low" : micro_df[i, "incomegrid"] <= 8 .&& micro_df[i, "incomegrid"] >= 5 ? "middle" : "high" for i in 1:nrow(micro_df)]

    # Create wealth by income, by taking the weighted sum of wealth over the groups, per quarter 
    replace!(micro_df[!, :cop_share], NaN => 0.0)

    # Merge gdp_series on time
    micro_df = innerjoin(micro_df, gdp_series, on=:time)
    micro_df[!, :tot_hhs_in_group] = micro_df[!, :cop_share] .* micro_df[!, :tot_hhs]

    divide_by_agg = ty == "normal" ? 1 : gdp_series[!, m*"_per_hh"]
    for m in measures
        micro_df[!, m*"levels"] = micro_df[!, m] .* micro_df[!, "tot_hhs_in_group"] ./ divide_by_agg
    end

    final_df = combine(groupby(micro_df, [:time, :income_groups]), [:wealthlevels, :tot_hhs_in_group] .=> x -> sum(x))

    # Create a wide dataframe where we have 1 column for each wealth to income group 
    wide_df1 = unstack(final_df, :time, :income_groups, :wealthlevels_function)
    wide_df2 = unstack(final_df, :time, :income_groups, :tot_hhs_in_group_function)

    # Rename columns 
    rename!(wide_df1, :low => :wealthlevels_bottom40, :middle => :wealthlevels_next40, :high => :wealthlevels_top20)
    rename!(wide_df2, :low => :tot_hhs_bottom40, :middle => :tot_hhs_next40, :high => :tot_hhs_top20)

    # Combine wide_df1 and wide_df2
    wide_df = innerjoin(wide_df1, wide_df2, on=:time)

    # Rename time to dates 
    rename!(wide_df, :time => :dates)

    return wide_df
end


"""
    weighted_conditional_bin_means_simple(outcome, conditioning, w, shares)

Compute weighted mean of `outcome` for observations in share-defined groups of
`conditioning`.  Simplified version for sample micro data (plain vectors).
"""
function weighted_conditional_bin_means_simple(
    outcome::AbstractVector,
    conditioning::AbstractVector,
    w::AbstractVector,
    shares::AbstractVector{<:Real},
)
    @assert length(outcome) == length(conditioning) == length(w)
    n = length(shares)
    totw = sum(w)
    totw > 0 || return fill(NaN, n)

    ix = sortperm(conditioning)
    ws = w[ix]
    os = outcome[ix]

    targets = Float64.(shares) .* totw
    num = zeros(Float64, n)
    den = zeros(Float64, n)

    bin = 1
    filled = 0.0
    tol = 1e-14 * totw

    @inbounds for i in eachindex(os)
        wi = ws[i]
        oi = os[i]
        wi <= 0 && continue
        while wi > 0 && bin <= n
            rem = targets[bin] - filled
            if rem <= tol
                bin += 1
                filled = 0.0
                continue
            end
            take = min(wi, rem)
            num[bin] += oi * take
            den[bin] += take
            wi -= take
            filled += take
            if filled >= targets[bin] - tol
                bin += 1
                filled = 0.0
            end
        end
    end

    means = fill(NaN, n)
    @inbounds for b = 1:n
        means[b] = den[b] > 0 ? num[b] / den[b] : NaN
    end
    return means
end


"""
    compute_sample_cross_conditional(sample_csv_path, outcome_var, conditioning_var; shares)

Compute cross-conditional share-group means from a sample micro CSV.
Returns a DataFrame with time and columns like `<outcome>_by_<conditioning>_bot50`, etc.,
or `nothing` if the required variables are not present in the sample.
"""
function compute_sample_cross_conditional(
    sample_csv_path::String,
    outcome_var::String,
    conditioning_var::String;
    shares = [0.5, 0.4, 0.1]
)
    df = CSV.read(sample_csv_path, DataFrame)
    # Check both variables exist
    if !(Symbol(outcome_var) in propertynames(df)) ||
       !(Symbol(conditioning_var) in propertynames(df))
        return nothing
    end

    df[!, :time] = QuarterlyDate.(df[!, :time])
    dates = sort(unique(df.time))

    # Share labels
    cum = 0.0
    labels = String[]
    for (i, s) in enumerate(shares)
        pct = round(Int, s * 100)
        if i == 1
            push!(labels, "bot$(pct)")
        elseif i == length(shares)
            push!(labels, "top$(pct)")
        else
            push!(labels, "mid$(pct)")
        end
        cum += s
    end

    result = DataFrame(time = dates)
    for sl in labels
        result[!, "$(outcome_var)_by_$(conditioning_var)_$(sl)"] = zeros(length(dates))
    end

    for (di, d) in enumerate(dates)
        sub = filter(row -> row.time == d, df)
        cond_x = Float64.(sub[!, Symbol(conditioning_var)])
        out_x = Float64.(sub[!, Symbol(outcome_var)])
        w = Float64.(sub[!, :weight])
        means = weighted_conditional_bin_means_simple(out_x, cond_x, w, shares)
        for (si, sl) in enumerate(labels)
            result[di, "$(outcome_var)_by_$(conditioning_var)_$(sl)"] = means[si]
        end
    end
    return result
end


function find_tmin_tmax(time_dict, k)

    year_max = maximum(keys(time_dict[k]))
    qtr_max = maximum(time_dict[k][year_max])

    year_min = minimum(keys(time_dict[k]))
    qtr_min = minimum(time_dict[k][year_min])

    # Create two dictionaries
    tminₖ = Dict("year" => year_min, "quarter" => qtr_min)
    tmaxₖ = Dict("year" => year_max, "quarter" => qtr_max)

    return tminₖ, tmaxₖ
end

make_range(ymin, qmin, ymax, qmax) = Dict(
    "tmin" => Dict("year" => Int(ymin), "quarter" => Int(qmin)),
    "tmax" => Dict("year" => Int(ymax), "quarter" => Int(qmax))
)

function df_range(df)
    @assert !isempty(df) "DataFrame is empty"
    dmin = minimum(skipmissing(df.dates))
    dmax = maximum(skipmissing(df.dates))
    return make_range(year(dmin), quarter(dmin), year(dmax), quarter(dmax))
end


function intersect_ranges(ranges)
    to_tuple(d) = (Int(d["year"]), Int(d["quarter"]))
    starts = [to_tuple(r["tmin"]) for r in ranges]
    ends = [to_tuple(r["tmax"]) for r in ranges]
    tmin = maximum(starts)
    tmax = minimum(ends)
    if tmin > tmax
        return Dict("tmin" => Dict("year" => missing, "quarter" => missing),
            "tmax" => Dict("year" => missing, "quarter" => missing))
    end
    return make_range(tmin[1], tmin[2], tmax[1], tmax[2])
end

function compare_to_DFA(dv, ty, time_params, type, case, measures, func_dict, gdp_series, tag)
    @unpack tmin, tmax, time_dict = time_params # 'time_dict' gives us the dates of the data that entered the estimation

    """Comparing the cyclicality generated from the DFA to ours."""

    # First check that any measures we've just estimated are in here 
    if "wealth" in measures && "income" in measures
        @info("Wealth and income are in measures, performing wealth by income DFA comparison")
    else
        return "no DFA wealth by income comparison since income or wealth are not in measures"
    end

    dimension = length(measures)
    label = "$dimension" * "D" * "_$case"

    # Load data
    init_path = BASE_PATH
    file_path = init_path * "/2_Data_processing/validation/DFA/DFA.xlsx"
    m_label = measures_folder(measures)
    folder = type == :optimization ? "from_optimization" : "from_mcmc"

    # Get DFA data 
    wealth_by_income = DataFrame(XLSX.readtable(file_path, "wealth_by_income", header=true,))

    # Extract dates
    dfa_dates = Dict()
    wealth_by_income[:, "time"] = QuarterlyDate.(wealth_by_income[:, "time"])

    # Find earliest and latest quarter
    earliest_period = minimum(wealth_by_income[:, "time"])
    latest_period = maximum(wealth_by_income[:, "time"])

    wealth_by_income[!, "year"] = Dates.year.(wealth_by_income[:, "time"])
    wealth_by_income[!, "quarter"] = Dates.quarterofyear.(wealth_by_income[!, "time"])
    dfa_dates["tmin"] = Dict(
        "year" => year(earliest_period),
        "quarter" => quarter(earliest_period)
    )
    dfa_dates["tmax"] = Dict(
        "year" => year(latest_period),
        "quarter" => quarter(latest_period)
    )
    # dfa_dates["tmax"] = Dict("year" => maximum(wealth_by_income[!, "year"]), "quarter" => maximum(wealth_by_income[wealth_by_income[!, "year"].==maximum(wealth_by_income[!, "year"]), "quarter"]))
    wealth_by_income[!, "dates"] = QuarterlyDate.(dfa_dates["tmin"]["year"], dfa_dates["tmin"]["quarter"]):Quarter(1):QuarterlyDate(dfa_dates["tmax"]["year"], dfa_dates["tmax"]["quarter"])


    # Create dictionaries 
    top = Dict{String,Dict}()
    bot = Dict{String,Dict}()
    mid = Dict{String,Dict}()

    top_dfa = Dict{String,Dict}()
    bot_dfa = Dict{String,Dict}()
    mid_dfa = Dict{String,Dict}()

    for dt in [top, bot, mid]
        for k in [collect(keys(dv))...]
            dt[k] = Dict{String,Dict}()
            for j in ["shares", "levels", "quantiles"]
                dt[k][j] = Dict()
            end
        end
    end

    dsources = sort(setdiff(collect(keys(dv)), ["consensus"]))

    # Importing our data for wealth by income
    obs_k = []
    agr_dates = Dict()
    wealth_by_incomeₖ = Dict()

    for k in collect(dsources)
        # Find the observed measures 
        obs_meas = k == "consensus" ? measures : get_obs_meas(func_dict, k, measures)

        if "income" in obs_meas && "wealth" in obs_meas
            # Find tmin and tmax for 'k'
            idₖ = findall(x -> x == k, dsources)[1]

            # Bounds set by user
            tminₖ, tmaxₖ = find_tmin_tmax(time_dict, idₖ)

            # Data from our estimates
            micro_df = generate_wealth_by_income_df(k, ty, measures, folder, gdp_series, tag) #CSV.read(init_path * "/2_Data_processing/" * k * "_wealth_by_income" * ".csv", DataFrame)

            # Find the agreeable dates between the series 
            agr_dates[k] = find_agreeable_dates(dfa_dates, :quarterly, tminₖ, tmaxₖ)

            # 2) micro's own window
            micro_rng = df_range(micro_df)

            # 3) final intersection across all three windows
            final_rng = intersect_ranges([agr_dates[k], micro_rng])

            # Materialize bounds once
            tminQ = QuarterlyDate(final_rng["tmin"]["year"], final_rng["tmin"]["quarter"])
            tmaxQ = QuarterlyDate(final_rng["tmax"]["year"], final_rng["tmax"]["quarter"])

            # Subset the DFA data and micro data based on the agreeable dates
            wealth_by_incomeₖ[k] = filter(x -> tminQ <= x.dates <= tmaxQ, wealth_by_income)
            filter!(x -> tminQ <= x.dates <= tmaxQ, micro_df)

            push!(obs_k, k)

            # filter!(x -> x.dates >= QuarterlyDate(agr_dates[k]["tmin"]["year"], agr_dates[k]["tmin"]["quarter"]), micro_df)
            # filter!(x -> x.dates <= QuarterlyDate(agr_dates[k]["tmax"]["year"], agr_dates[k]["tmax"]["quarter"]), micro_df)

            # Levels and quantiles are transformed
            top[k]["levels"]["wealth_by_income"] = micro_df[!, "wealthlevels_top20"]
            bot[k]["levels"]["wealth_by_income"] = micro_df[!, "wealthlevels_bottom40"]
            mid[k]["levels"]["wealth_by_income"] = micro_df[!, "wealthlevels_next40"]

            top[k]["quantiles"]["wealth_by_income"] = ty == "normal" ? log_transformation(top[k]["levels"]["wealth_by_income"] ./ micro_df[!, "tot_hhs_top20"]) : top[k]["levels"]["wealth_by_income"] ./ micro_df[!, "tot_hhs_top20"]
            bot[k]["quantiles"]["wealth_by_income"] = ty == "normal" ? log_transformation(bot[k]["levels"]["wealth_by_income"] ./ micro_df[!, "tot_hhs_bottom40"]) : bot[k]["levels"]["wealth_by_income"] ./ micro_df[!, "tot_hhs_bottom40"]
            mid[k]["quantiles"]["wealth_by_income"] = ty == "normal" ? log_transformation(mid[k]["levels"]["wealth_by_income"] ./ micro_df[!, "tot_hhs_next40"]) : mid[k]["levels"]["wealth_by_income"] ./ micro_df[!, "tot_hhs_next40"]
        end
    end

    # Remove the datasets that don't have wealth by income 
    dsources = copy(obs_k)

    for dt in [top_dfa, bot_dfa, mid_dfa]
        for k in dsources
            dt[k] = Dict{String,Dict}()
            for j in ["quantiles"]
                dt[k][j] = Dict()
            end
        end
    end

    # Wealth by income 
    # Shares
    for k in dsources
        # top_dfa[k]["shares"]["wealth"]          = wealth_by_incomeₖ[k][!, "sharetop20"] ./ 100
        # bot_dfa[k]["shares"]["wealth"]          = wealth_by_incomeₖ[k][!, "sharebottom40"] ./ 100
        # mid_dfa[k]["shares"]["wealth"]          = wealth_by_incomeₖ[k][!, "sharenext40"] ./ 100

        # # Levels
        # top_dfa[k]["levels"]["wealth"]          = log_transformation(wealth_by_incomeₖ[k][!, "wealthlevels_top20"] .* 1000000)
        # bot_dfa[k]["levels"]["wealth"]          = log_transformation(wealth_by_incomeₖ[k][!, "wealthlevels_bottom40"] .* 1000000)
        # mid_dfa[k]["levels"]["wealth"]          = log_transformation(wealth_by_incomeₖ[k][!, "wealthlevels_next40"] .* 1000000)

        # Quantiles 
        top_dfa[k]["quantiles"]["wealth"] = log_transformation(wealth_by_incomeₖ[k][!, "quantile_top20"])
        bot_dfa[k]["quantiles"]["wealth"] = log_transformation(wealth_by_incomeₖ[k][!, "quantile_bottom40"])
        mid_dfa[k]["quantiles"]["wealth"] = log_transformation(wealth_by_incomeₖ[k][!, "quantile_next40"])
    end

    # Plot the series against each other, one plot for each dictionary
    path = init_path * "/7_Results/$m_label" * "$tag" * "/$folder/plots/wealth/external_comparisons/"
    mkpath(path)

    series_t = ["top20", "bottom40", "next40"]
    series_y = [["top", "20"], ["bottom", "40"], ["next", "40"]]

    corr_df = DataFrame(rand(1, 1), :auto)

    plot_tag = ty == "average" ? "_detrended_" : ""

    foo = zip([top, bot, mid], [top_dfa, bot_dfa, mid_dfa])

    for (q, dt) in enumerate(foo)
        for j in ["quantiles"]
            for (i, k) in enumerate(dsources)
                # Wealth by income 
                deviations_DFA = dt[2][k][j]["wealth"] .- HP(dt[2][k][j]["wealth"], 1600)
                deviations_DFA .= HP(deviations_DFA, 6)

                Plots.plot()

                xaxis = QuarterlyDate(agr_dates[k]["tmin"]["year"], agr_dates[k]["tmin"]["quarter"]):Quarter(1):QuarterlyDate(agr_dates[k]["tmax"]["year"], agr_dates[k]["tmax"]["quarter"])

                plot_int = k == "SIPP2" ? 4 : 40

                # Object of interest 
                local deviations_r
                if ty == "normal"
                    deviations_r = dt[1][k][j]["wealth_by_income"] .- HP(dt[1][k][j]["wealth_by_income"], 1600)
                    deviations_r .= HP(deviations_r, 6)
                elseif ty == "average"
                    deviations_r = HP(dt[1][k][j]["wealth_by_income"], 6)
                end

                # println(k)
                # println(size(deviations_DFA))
                # println(size(deviations_r))
                # println(deviations_r)
                # println(deviations_DFA)

                # Extract Correlation between DFA and Approximation
                col = k * "_" * j * "_cycle_" * series_y[q][1] * series_y[q][2] * "_wealthbyincome"
                ρ = round(cor(deviations_DFA, deviations_r), digits=2)
                corr_df[!, col] = [ρ]

                Plots.plot!(
                    axes(deviations_r),
                    deviations_r,
                    lc=:red, #select_color(k),
                    # xlabel = L"\textrm{Year}",  
                    ylabel=L"\textrm{Cyclical\,\, Component}",#L"\textrm{Wealth\,\, of\,\, %$(series_y[q][1])\,\,%$(series_y[q][2])\,\,Income\,\,%$(j[1:end-1])}",
                    xformatter=:latex,
                    yformatter=:latex,
                    xtickfontsize=14,
                    ytickfontsize=14,
                    legendfontsize=10,
                    guidefontsize=14,
                    xticks=(collect(axes(deviations_r))[1][1:plot_int:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(xaxis[1:plot_int:end])]),
                    legend=false, #series_y[q][1] == "top" ? :best : false, #:outertopright,
                    label="", #L"\textrm{Model}",#L"\textrm{%$(k)\,\,Approximation}",
                    # label = L"\textrm{%$(k)\,-\,DFA\,-\, ρ: %$(ρ)}",
                    lw=4, dpi=500
                )

                # Wealth by income  
                Plots.plot!(
                    axes(deviations_DFA),
                    deviations_DFA,
                    lc=select_color("DFA"),
                    ls=select_linestyle("DFA"),
                    lw=4,
                    legendmarkerstroke=2,
                    label="", #L"\textrm{DFA}",
                )

                # Plots.plot!(
                #     axes(deviations_r), 
                #     deviations_r, 
                #     la=0.0,
                #     label = L"\textrm{%$(k)\,-\,DFA\,-\, corr: %$(ρ)}",
                #     lw=4, dpi=500 
                # )                        

                Plots.savefig(path * k * "_wealth_byincome_" * j * "_DFA_cycle_" * series_t[q] * "_" * plot_tag * label * ".pdf")
            end

            # Plotting the raw version, Wealth by income 
            if ty == "normal"
                local ρ
                for (i, k) in enumerate(dsources)
                    # deviations_DFA = dt[2][k][j]["wealth"] .- HP(dt[2][k][j]["wealth"], 1600)
                    # deviations_DFA .= HP(deviations_DFA, 6)

                    xaxis = QuarterlyDate(agr_dates[k]["tmin"]["year"], agr_dates[k]["tmin"]["quarter"]):Quarter(1):QuarterlyDate(agr_dates[k]["tmax"]["year"], agr_dates[k]["tmax"]["quarter"])

                    Plots.plot()
                    # Extract Correlation between DFA and Approximation
                    col = k * "_" * j * "_" * series_y[q][1] * series_y[q][2] * "_wealthbyincome"
                    ρ = round(cor(dt[2][k][j]["wealth"], dt[1][k][j]["wealth_by_income"]), digits=2)
                    corr_df[!, col] = [ρ]

                    plot_int = k == "SIPP2" ? 4 : 20

                    Plots.plot!(
                        axes(dt[1][k][j]["wealth_by_income"]),
                        dt[1][k][j]["wealth_by_income"],
                        # xlabel = L"\textrm{Year}",
                        ylabel=L"\textrm{Cyclical\,\, Component}",
                        xformatter=:latex,
                        yformatter=:latex,
                        xticks=(collect(axes(dt[2][k][j]["wealth"]))[1][1:plot_int:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(xaxis[1:plot_int:end])]),
                        legend=false,
                        xtickfontsize=10,
                        ytickfontsize=10,
                        legendfontsize=10,
                        guidefontsize=14,
                        label="", #L"\textrm{Model}",
                        # label = L"\textrm{%$(k)\,-\,DFA\,-\, corr.: %$(ρ)}",
                        lw=4, dpi=500
                    )

                    Plots.plot!(
                        axes(dt[2][k][j]["wealth"]),
                        dt[2][k][j]["wealth"],
                        label="",
                        # label = L"\textrm{DFA}",
                        # legendmarkerstroke=2,
                        lc=select_color("DFA"),
                        ls=select_linestyle("DFA"),
                        lw=4
                    )
                    Plots.plot!([], [], ls=select_linestyle("DFA"), lc=select_color("DFA"), lw=2, label=L"\textrm{DFA}")

                    # Plots.plot!(
                    # axes(dt[k][j]["wealth_by_income"]), 
                    # dt[k][j]["wealth_by_income"], 
                    # la=0.0,
                    # label = L"\textrm{%$(k)\,-\,DFA\,-\, corr.: %$(ρ)}",
                    # lw=4, dpi=500 
                    # )    

                    Plots.savefig(path * k * "_wealth_byincome_" * j * "_DFA_" * series_t[q] * "_" * plot_tag * label * ".pdf")
                end
            end
        end
    end

    select!(corr_df, Not(:x1))
    select!(corr_df, sort(names(corr_df)))
    # x_objects = join.(split.(names(corr_df), "_"), "\n ")


    # Plots.plot(
    #     axes(Matrix(corr_df), 2),
    #     Matrix(corr_df)[1, :],
    #     xlabel = L"\textrm{Series}",
    #     ylabel = L"\textrm{Correlation\,\,with\,\,DFA}",
    #     fontfamily="Computer Modern",
    #     lw=4,
    #     xticks=(axes(Matrix(corr_df), 2), [object for object in x_objects]),
    #     legend=false
    #     )

    # Plots.savefig(init_path * "/7_Results/$m_label/$folder/plots/correlations/" * "DFA_correlations" * ".pdf")

    # Export DFA correlations 
    CSV.write(init_path * "/7_Results/$m_label" * "$tag" * "/$folder/plots/correlations/" * "DFA_correlations" * plot_tag * ".csv", corr_df)
end


function get_raw_or_cycle(series, ty, ext_source, type, choice, object, meas, select_series=false)
    λ = type == :quarterly ? 1600 : 6


    if choice == :raw
        # if ty == "normal"
        return series
        # for the rest, if "detrended/average", we don't detrend -- only smooth -- if external source, we detrend and smooth 
        # elseif ty == "average" && ext_source 
        #     return series .- HP(series, λ)
        # elseif ty == "average" && ext_source == false
        #     # find linear trend 
        #     # n    = length(series)
        #     # time = 1:n
        #     # X    = [ones(n) time]
        #     # coeffs = X \ series
        #     # trend = X * coeffs

        #     return series .- nanmean(series) 
        # end
    else
        # Take the average of the series by year and return the output
        if object == "quantiles" && select_series != false
            corr_series = series .* select_series[:, meas*"_per_hh"]
            new_series = corr_series .- HP(corr_series, λ)
        else
            # if ty == "normal"
            new_series = series .- HP(series, λ)
            # elseif ty == "average" && ext_source
            #     new_series  = series .- HP(series, λ)
            # elseif ty == "average" && ext_source == false
            #     # find linear trend 
            #     # n    = length(series)
            #     # time = 1:n
            #     # X    = [ones(n) time]
            #     # coeffs = X \ series
            #     # trend = X * coeffs
            #     new_series  = series .- nanmean(series)
            # end
        end

        # Additional smoothing
        if type == :quarterly
            new_series .= HP(new_series, 6)
        end

        return new_series
    end
end



function compare_to_data(dv, ty, func_data, obs_data, user_t, time_params, model_options, type, label, data_bounds=false)
    @unpack measures, lags, estimator, tag = model_options
    @unpack tmin, tmax = time_params
    @unpack gdp_series = obs_data
    @unpack func_dict, year_vec = func_data

    # Getting the average series of the data
    # avg_series = generate_average_series(dv, gdp_series, measures)

    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid, integral_cop_grid = estimator
    else
        @unpack grid_pcf, grid_cop = estimator
    end
    grid_choice_pcf = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf
    grid_choice_cop = typeof(estimator) <: SeriesEstimator ? integral_cop_grid : grid_cop

    data_sources = setdiff(sort(collect(keys(dv))), ["consensus"])
    smin, smax = user_t
    m_label = measures_folder(measures)
    dates = Dict()

    # Common annual series 
    top = Dict()
    bot = Dict()
    mid = Dict()

    for dt in [top, bot, mid]
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
                # dt["consensus"][obj] = Dict()
            end
        end
        # for obj in ["shares", "quantiles"]
        #     dt["consensus"][obj] = Dict()
        # end
    end

    # Defining the series 
    series = Dict()
    if grid_choice_pcf == 10 || grid_choice_pcf == 20 || grid_choice_pcf == 100
        series["bot"] = "bottom50"
        series["mid"] = "next40"
        series["top"] = "top10"
    elseif grid_choice_pcf == 5
        series["bot"] = "bottom40"
        series["mid"] = "next40"
        series["top"] = "top20"
    end

    # Path situation
    folder = type == :optimization ? "from_optimization" : "from_mcmc"
    init_path = BASE_PATH
    path = init_path * "/7_Results/$m_label" * "$tag" * "/$folder/plots/"
    mkpath(path)

    dates = Dict()

    # Everything will be relative to the average 
    gdp_series[!, "date"] = QuarterlyDate.(gdp_series[!, "time"])
    correction_names = [meas * "_per_hh" for meas in measures]
    select_series = select(gdp_series, [correction_names..., "date"])

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

        # Filter the aggregate correction 
        sub_select_series = filter(x -> x.date >= QuarterlyDate(dates[meas]["tmin"]["year"], dates[meas]["tmin"]["quarter"]), select_series)
        agg_corr = sub_select_series[!, meas*"_per_hh"]

        # @info(dates[meas])
        for obj in ["shares", "quantiles"] #TODO: levels?

            # # Adjust frequencies, Confine to dates
            # top["consensus"][obj][meas] = subset_to_cutoff(dv["consensus"][ty][meas][obj]["common series"][series["top"]], dates[meas], tmin, tmax, obj, agg_corr)
            # bot["consensus"][obj][meas] = subset_to_cutoff(dv["consensus"][ty][meas][obj]["common series"][series["bot"]], dates[meas], tmin, tmax, obj, agg_corr)
            # mid["consensus"][obj][meas] = subset_to_cutoff(dv["consensus"][ty][meas][obj]["common series"][series["mid"]], dates[meas], tmin, tmax, obj, agg_corr)


            # Adjust frequencies, Confine to dates
            for source in data_sources
                if data_bounds != false
                    top[source]["lb"][obj][meas] = subset_to_cutoff(data_bounds[source]["lb"][meas][obj]["common series"][series["top"]], dates[meas], tmin, tmax, obj, agg_corr)
                    bot[source]["lb"][obj][meas] = subset_to_cutoff(data_bounds[source]["lb"][meas][obj]["common series"][series["bot"]], dates[meas], tmin, tmax, obj, agg_corr)
                    mid[source]["lb"][obj][meas] = subset_to_cutoff(data_bounds[source]["lb"][meas][obj]["common series"][series["mid"]], dates[meas], tmin, tmax, obj, agg_corr)

                    top[source]["ub"][obj][meas] = subset_to_cutoff(data_bounds[source]["ub"][meas][obj]["common series"][series["top"]], dates[meas], tmin, tmax, obj, agg_corr)
                    bot[source]["ub"][obj][meas] = subset_to_cutoff(data_bounds[source]["ub"][meas][obj]["common series"][series["bot"]], dates[meas], tmin, tmax, obj, agg_corr)
                    mid[source]["ub"][obj][meas] = subset_to_cutoff(data_bounds[source]["ub"][meas][obj]["common series"][series["mid"]], dates[meas], tmin, tmax, obj, agg_corr)
                end

                # For the Approximations 
                top[source][obj][meas] = subset_to_cutoff(dv[source][ty][meas][obj]["common series"][series["top"]], dates[meas], tmin, tmax, obj, agg_corr)
                bot[source][obj][meas] = subset_to_cutoff(dv[source][ty][meas][obj]["common series"][series["bot"]], dates[meas], tmin, tmax, obj, agg_corr)
                mid[source][obj][meas] = subset_to_cutoff(dv[source][ty][meas][obj]["common series"][series["mid"]], dates[meas], tmin, tmax, obj, agg_corr)

                # For the observed 
                top[source*"_obs"][obj][meas] = subset_to_cutoff(func_dict[source][meas][obj]["common series"][series["top"]], dates[meas], tmin, tmax, obj, agg_corr)
                bot[source*"_obs"][obj][meas] = subset_to_cutoff(func_dict[source][meas][obj]["common series"][series["bot"]], dates[meas], tmin, tmax, obj, agg_corr)
                mid[source*"_obs"][obj][meas] = subset_to_cutoff(func_dict[source][meas][obj]["common series"][series["mid"]], dates[meas], tmin, tmax, obj, agg_corr)
            end

            # Plots of both data, 2 lines per plot (we take the first quarter for the year)
            dts = QuarterlyDate(dates[meas]["tmin"]["year"], dates[meas]["tmin"]["quarter"]):Quarter(1):QuarterlyDate(dates[meas]["tmax"]["year"], dates[meas]["tmax"]["quarter"])
            xaxis = collect(1:length(dts))

            series_dict = Dict("top" => top, "mid" => mid, "bot" => bot)
            series_lab = Dict("top" => "top 10", "mid" => "next 40", "bot" => "bottom 50")
            @info("Building Approximation plots, comparing Approximations to data: $obj")
            # Top
            plot_tag = ty == "average" ? "_detrended_" : ""
            for (dist_label, dist) in series_dict
                Plots.plot(xaxis, [NaN for _ in eachindex(xaxis)], label="", xticks=(xaxis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(dts[1:20:end])]))

                for (j, source) in enumerate(data_sources)
                    plot_name = occursin("CEX", source) ? "CEX" : source
                    # Not all sources have observe the measure, so we need to check that 
                    if all(isnan.(dist[source][obj][meas]))
                        nothing
                    else
                        cond = source != "consensus" ? findall(!isnan, dist[source*"_obs"][obj][meas]) : [1] # .!isnan.(all_lv_o[1][:, j])
                        println(cond)
                        s_axis = xaxis[cond[1]:cond[end]] # start at the first observation 
                        s_data = dist[source][obj][meas][cond[1]:cond[end]]

                        Plots.plot!(s_axis,
                            s_data,
                            xlabel="",
                            ylabel=obj == "quantiles" ? L"\textrm{%$(M)\,\,rel.\,  to\,\, average}" : L"\textrm{%$(M)}\,\,\textrm{%$(obj)}",
                            xformatter=:latex,
                            yformatter=:latex,
                            # xticks=(s_axis[1:20:end], [L"%$(date)"[1:5] * "\$" for (_,date) in enumerate(s_dts[1:20:end])]),
                            legend=:best,
                            lw=4, dpi=500,
                            xtickfontsize=14,
                            ytickfontsize=14,
                            guidefontsize=14,
                            lc=select_color(plot_name),
                            yerror=data_bounds != false ? (s_data - dist[source]["lb"][obj][meas][cond[1]:end], dist[source]["ub"][obj][meas][cond[1]:end] - s_data) : nothing,
                            label=L"\textrm{%$(plot_name)}"
                        )
                        # Plot the observed data
                        c_data = dist[source*"_obs"][obj][meas][cond]

                        Plots.plot!(xaxis[cond],
                            c_data,
                            marker=:diamond,
                            markercolor=:black,
                            lc=:black,
                            la=0.3,
                            lw=3, dpi=500,
                            label=j == length(data_sources) ? L"\textrm{Linear \,\, Interpolation}" : ""
                        )
                    end
                end
                Plots.savefig(path * "/$meas/" * "/data_comparisons/" * meas * "_" * obj * "_" * series[dist_label] * plot_tag * "_" * label * ".pdf")
            end
        end
    end
end


function arsinh(x)
    return log(x + sqrt(x^2 + 1))  # the big is too handle large numbers 
end

function subset_to_cutoff(q_series, agr_dates, smin, smax, obj, agg_corr=false)

    # println(length(q_series))
    dts = QuarterlyDate(smin["year"], smin["quarter"]):Quarter(1):QuarterlyDate(smax["year"], smax["quarter"])

    a = length(dts)
    b = length(dts[dts.>=QuarterlyDate(agr_dates["tmin"]["year"], agr_dates["tmin"]["quarter"])])
    c = length(dts[dts.<=QuarterlyDate(agr_dates["tmax"]["year"], agr_dates["tmax"]["quarter"])])

    subset_series = q_series[a-b+1:end-(a-c)]

    # if ty == "average"
    #     return subset_series
    # else
    if obj == "quantiles" && agg_corr != false
        subset_series = subset_series ./ agg_corr

    elseif obj == "quantiles" && agg_corr == false
        # subset_series = log_transformation(subset_series)
        subset_series = inverse_hyperbolic_sine(subset_series)

    elseif obj == "levels"
        subset_series = inverse_hyperbolic_sine(subset_series)

    elseif obj == "shares"
        nothing
    end

    return subset_series
    # end
end


function generate_time_bounds(series, tmin, tmax)
    agr_dates = Dict()
    agr_dates["tmin"] = Dict()
    agr_dates["tmax"] = Dict()

    agr_dates["tmin"]["quarter"] = 4
    agr_dates["tmax"]["quarter"] = 4

    if minimum(series) > tmin["year"]
        agr_dates["tmin"]["year"] = minimum(series)

    elseif minimum(series) < tmin["year"]
        agr_dates["tmin"]["year"] = tmin["year"]

    elseif minimum(series) == tmin["year"]
        agr_dates["tmin"]["year"] = tmin["year"]
    end

    if maximum(series) > tmax["year"]
        agr_dates["tmax"]["year"] = tmax["year"]

    elseif maximum(series) < tmax["year"]
        agr_dates["tmax"]["year"] = maximum(series)

    elseif maximum(series) == tmax["year"]
        agr_dates["tmax"]["year"] = maximum(series) - 1
    end

    return agr_dates
end

# function generate_time_bounds(series, tmin, tmax)
#     min_diff    = minimum(series) - tmin["year"]
#     max_diff    = maximum(series) - tmax["year"]

#     local aligned_min, aligned_max
#     if min_diff > 0 
#         aligned_min = minimum(series)
#     elseif min_diff < 0 
#         aligned_min = tmin["year"]
#     elseif min_diff == 0 
#         aligned_min = tmin["year"] 
#     end

#     if max_diff > 0 
#         aligned_max = tmax["year"]
#     elseif max_diff < 0 
#         aligned_max = maximum(series)
#     elseif max_diff == 0 
#         aligned_max = tmax["year"] - 1 
#     end

#     return (amin = aligned_min, amax = aligned_max)
# end


function retrieve_external_sources(measures)

    # Path 
    init_path = BASE_PATH
    fold_path = init_path * "/2_Data_processing/validation"

    external_sources = Dict()
    help_dict = Dict()
    folders = [name for name in readdir(fold_path) if isdir(joinpath(fold_path, name))]

    # We create the entire structure and fill in if it exists for that external source 
    for obj in ["levels", "quantiles", "shares"]
        external_sources[obj] = Dict()

        for meas in measures
            external_sources[obj][meas] = Dict()

        end
    end

    for ext_s in folders
        for meas in measures
            for obj in ["shares", "quantiles", "levels"]
                try
                    external_sources[obj][meas][ext_s] = CSV.read(fold_path * "/$ext_s/" * meas * "_$obj" * "_$ext_s.csv", DataFrame)
                catch e
                    # external_sources[obj][meas][ext_s]    = NaN
                    # println("No data for $ext_s, $meas, $obj")
                    nothing
                end
            end
        end
    end

    return external_sources
end

# function find_agreeable_dates(ext_dates, ext_freqtype, smin, smax)
#     # Find the agreeable dates between the series 
#     agr_dates = Dict()
#     agr_dates["tmin"] = Dict()
#     agr_dates["tmax"] = Dict()

#     # Find the agreeable dates between the series
#     if ext_freqtype == :quarterly
#         if smin["year"] < ext_dates["tmin"]["year"]
#             agr_dates["tmin"]["year"] = ext_dates["tmin"]["year"]
#             agr_dates["tmin"]["quarter"] = ext_dates["tmin"]["quarter"]

#         elseif smin["year"] > ext_dates["tmin"]["year"]
#             agr_dates["tmin"]["year"] = smin["year"]
#             agr_dates["tmin"]["quarter"] = smin["quarter"]

#         else
#             agr_dates["tmin"]["year"] = smin["year"]
#             if smin["quarter"] < ext_dates["tmin"]["quarter"]
#                 agr_dates["tmin"]["quarter"] = ext_dates["tmin"]["quarter"]
#             else
#                 agr_dates["tmin"]["quarter"] = smin["quarter"]
#             end
#         end

#         if smax["year"] > ext_dates["tmax"]["year"]
#             agr_dates["tmax"]["year"] = ext_dates["tmax"]["year"]
#             agr_dates["tmax"]["quarter"] = ext_dates["tmax"]["quarter"]
#         elseif smax["year"] < ext_dates["tmax"]["year"]
#             agr_dates["tmax"]["year"] = smax["year"]
#             agr_dates["tmax"]["quarter"] = smax["quarter"]
#         else
#             agr_dates["tmax"]["year"] = smax["year"]
#             if smax["quarter"] > ext_dates["tmax"]["quarter"]
#                 agr_dates["tmax"]["quarter"] = ext_dates["tmax"]["quarter"]
#             else
#                 agr_dates["tmax"]["quarter"] = smax["quarter"]
#             end
#         end
#     else
#         if smin["year"] < ext_dates["tmin"]["year"]
#             agr_dates["tmin"]["year"] = ext_dates["tmin"]["year"]
#             agr_dates["tmin"]["quarter"] = 4
#         elseif smin["year"] > ext_dates["tmin"]["year"]
#             agr_dates["tmin"]["year"] = smin["year"]

#             if smin["quarter"] < ext_dates["tmin"]["quarter"]
#                 agr_dates["tmin"]["quarter"] = 4
#             elseif smin["quarter"] == ext_dates["tmin"]["quarter"]
#                 agr_dates["tmin"]["quarter"] = 4
#             else # will never happen by assumption
#                 agr_dates["tmin"]["quarter"] = 4
#             end
#         else
#             agr_dates["tmin"]["year"] = smin["year"]
#             if smin["quarter"] < ext_dates["tmin"]["quarter"] #TODO: assumption: annual data is Q4 
#                 agr_dates["tmin"]["quarter"] = 4
#             elseif smin["quarter"] == ext_dates["tmin"]["quarter"]
#                 agr_dates["tmin"]["quarter"] = 4
#             else # this case will never happen by assumption 
#                 agr_dates["tmin"]["year"] = smin["year"] + 1
#                 agr_dates["tmin"]["quarter"] = 4
#             end
#         end

#         if smax["year"] > ext_dates["tmax"]["year"]
#             agr_dates["tmax"]["year"] = ext_dates["tmax"]["year"]
#             agr_dates["tmax"]["quarter"] = 4

#         elseif smax["year"] < ext_dates["tmax"]["year"]
#             agr_dates["tmax"]["year"] = smax["year"]

#             if smin["quarter"] == ext_dates["tmin"]["quarter"]
#                 agr_dates["tmax"]["quarter"] = 4
#             elseif smin["quarter"] < ext_dates["tmin"]["quarter"]
#                 agr_dates["tmax"]["year"] = agr_dates["tmax"]["year"] - 1
#                 agr_dates["tmax"]["quarter"] = 4
#             else # will never happen by assumption
#                 agr_dates["tmax"]["quarter"] = 4
#             end
#         else
#             agr_dates["tmax"]["year"] = smax["year"]
#             if smax["quarter"] > ext_dates["tmax"]["quarter"]
#                 agr_dates["tmax"]["quarter"] = 4
#             elseif smax["quarter"] == ext_dates["tmax"]["quarter"]
#                 agr_dates["tmax"]["quarter"] = 4
#             else
#                 agr_dates["tmax"]["year"] = smax["year"] - 1
#                 agr_dates["tmax"]["quarter"] = 4
#             end
#         end
#     end
#     return agr_dates
# end
"""
find_agreeable_dates(ext_dates, ext_freqtype, smin, smax) -> Dict
Return the overlapping date window between:
- a quarterly series with bounds smin..smax (Dicts with "year","quarter"),
- an external series ext_dates with "tmin"/"tmax".

If `ext_freqtype == :quarterly`, use the quarters as-is.
Otherwise (annual), align the external series to Q4 only.

Returns:
Dict("tmin"=>Dict("year"=>Y, "quarter"=>Q),
     "tmax"=>Dict("year"=>Y, "quarter"=>Q))
If there is no overlap, both entries are set to `missing`.
"""
function find_agreeable_dates(ext_dates::Dict, ext_freqtype::Symbol,
    smin::Dict, smax::Dict)
    # helpers
    yq(d::Dict) = (Int(d["year"]), Int(d["quarter"]))
    dictyq(t::Tuple{Int,Int}) = Dict("year" => t[1], "quarter" => t[2])

    if ext_freqtype == :quarterly
        tmin = max(yq(smin), yq(ext_dates["tmin"]))
        tmax = min(yq(smax), yq(ext_dates["tmax"]))
    else
        # Annual: valid points are Q4 of each year
        ext_min_year = Int(ext_dates["tmin"]["year"])
        ext_max_year = Int(ext_dates["tmax"]["year"])

        # earliest Q4 on/after smin, and latest Q4 on/before smax
        first_q4_year = Int(smin["year"]) + (Int(smin["quarter"]) > 4 ? 1 : 0)   # (>4 won't happen, but robust)
        last_q4_year = Int(smax["year"]) - (Int(smax["quarter"]) < 4 ? 1 : 0)

        tmin_year = max(first_q4_year, ext_min_year)
        tmax_year = min(last_q4_year, ext_max_year)

        tmin = (tmin_year, 4)
        tmax = (tmax_year, 4)
    end

    # no overlap
    if tmin > tmax
        return Dict(
            "tmin" => Dict("year" => missing, "quarter" => missing),
            "tmax" => Dict("year" => missing, "quarter" => missing),
        )
    end

    return Dict("tmin" => dictyq(tmin), "tmax" => dictyq(tmax))
end



# function subset_to_cutoff(time_series, dts, agr_dates)


#     a = length(dts)
#     b = length(dts[dts .>= QuarterlyDate(agr_dates["tmin"]["year"], agr_dates["tmin"]["quarter"])])
#     c = length(dts[dts .<= QuarterlyDate(agr_dates["tmax"]["year"], agr_dates["tmax"]["quarter"])])

#     return time_series[a-b+1:end-(a-c)]    
# end

function nan_floor(type, corr)

    try
        a = floor(type, corr)
        return a
    catch e
        println(corr)
        return NaN
    end
end