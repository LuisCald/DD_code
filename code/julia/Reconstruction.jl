function reconstruct_data(par_final, param_sizes, hyperpriors, meas_ind, Σ_ids, model_elements, obs_data, model_options, time_params, data_sources; reconstruction_to_show=false, dε_smoothed=false)
    """Reverse all steps taken in the data preparation."""

    @unpack lags, blind_to, case, pre_multiply, estimator, tag = model_options
    @unpack factor_count, n_less_than_one, u, pcs, βs, trend, Gⱼ = model_elements

    @unpack gdp_series, df_vec = obs_data
    @unpack tot_periods, tmin, tmax = time_params

    number_of_dfs = length(df_vec.data)
    data_names = df_vec[2]

    local x_smoothed
    if dε_smoothed != false
        A, B, _, _ = matrisize(par_final, param_sizes, case)
        T = size(u, 2)
        x_smoothed = zeros(factor_count, T)
        x̂ = zeros(factor_count)
        for t in 1:T
            x_smoothed[:, t] = state_transition(x̂, dε_smoothed[:, t], u[:, t], A, B)
            copyto!(x̂, x_smoothed[:, t])
        end
    else
        smoother_output, _, _ = likeli(model_elements, par_final, param_sizes, hyperpriors, Σ_ids, model_options; smooth=true)
        @unpack x_smoothed, dε_smoothed = smoother_output
    end

    @unpack proj, means, stds, agg_count = model_elements
    @unpack measures = model_options
    @unpack grid_pcf, grid_cop = estimator
    dimension = length(measures)

    # Plot the smoothed data
    Plots.plot()
    dts = QuarterlyDate(tmin["year"], tmin["quarter"]):Quarter(1):QuarterlyDate(tmax["year"], tmax["quarter"])
    xaxis = 1:tot_periods
    init_path = BASE_PATH
    n_factors = size(proj, 2)

    for i in 1:n_factors
        Plots.plot!(
            xaxis[5:end],
            x_smoothed[i, 5:end],
            xformatter=:latex,
            yformatter=:latex,
            # xticks=(xaxis[5:20:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(dts[5:20:end])]),
            xticks=(xaxis[5:20:end], [i for i in 5:20:tot_periods]),
            ylabel=L"\textrm{Factor\,\, Value}",
            label=L"\textrm{Factor \,\, %$(i)}",
            legend=:best) #, xticks=(collect(1:20:tot_periods), [L"%$(date)"[1:5] * "\$" for (_,date) in enumerate(gdp_series[1:20:end, :time])]))
    end
    Plots.savefig(init_path * "/7_Results/factor_analysis/distributional_factors$tag.pdf")

    aggregate_rep = hasproperty(model_options, :aggregate_rep) ? model_options.aggregate_rep : :as_states

    X = Vector{Matrix{Float64}}(undef, number_of_dfs)
    T = size(x_smoothed, 2)
    for i in eachindex(Gⱼ)
        if aggregate_rep == :as_inputs
            # :as_inputs — x_smoothed is already factor-only (4r rows)
            X[i] = (Gⱼ[i] * x_smoothed)
        else
            # :as_states — slice off aggregate rows from x_smoothed
            X[i] = (Gⱼ[i] * x_smoothed[1:end-agg_count, :])
        end
    end

    # Split by object, multiply by stds, reform object 
    add_variance!(estimator, X, stds, measures)
    add_mean!(X, means) # completely filled matrices

    # Generate two X's: one with the trend and one with average trend in

    # What kind of trend to add back?
    new_trend = select_trend(trend, "average")

    X_dict = Dict()
    X_dict["normal"] = deepcopy(X)
    X_dict["average"] = deepcopy(X)

    # Add trend back for each 
    t = collect(1:tot_periods)
    for j in eachindex(X)
        for t in 1:tot_periods
            # For copula, weights do not seem to be affected
            X_dict["normal"][j][:, t] .+= trend[j][:, t] #
            X_dict["average"][j][:, t] .+= new_trend[j] # By adding the nanmean, we can generate estimates for all periods ... to check is how similar are the copula estimates across both
        end
    end

    # Save smoothed factors and reconstructed coefficients as wide DataFrames with dates
    m_label = measures_folder(measures)
    save_dir = init_path * "/7_Results/$m_label" * "$tag" * "/from_mcmc/data"
    mkpath(save_dir)

    # 1. Smoothed factors (T × n_factors)
    factor_df = DataFrame(Matrix(x_smoothed'), ["x$i" for i in 1:size(x_smoothed, 1)])
    factor_df[!, "time"] = collect(dts)
    CSV.write(save_dir * "/smoothed_factors.csv", select(factor_df, "time", :))

    # 2. Reconstructed coefficients per data source
    for (i, ds) in enumerate(data_sources)
        n_coeff = size(X_dict["normal"][i], 1)
        col_names = ["x$k" for k in 1:n_coeff]

        df_normal = DataFrame(Matrix(X_dict["normal"][i]'), col_names)
        df_normal[!, "time"] = collect(dts)
        CSV.write(save_dir * "/$(ds)_coefficients_normal.csv", select(df_normal, "time", :))

        df_avg = DataFrame(Matrix(X_dict["average"][i]'), col_names)
        df_avg[!, "time"] = collect(dts)
        CSV.write(save_dir * "/$(ds)_coefficients_average.csv", select(df_avg, "time", :))
    end

    # Create consensus average.
    X_array_norm = cat(X_dict["normal"]..., dims=3)
    X_array_avg = cat(X_dict["average"]..., dims=3)

    X̄_norm = similar(X_dict["normal"][1])
    X̄_avg = similar(X_dict["average"][1])

    # Take the average of the 3 dimensions, skipping over NaN, element by element
    for i in axes(X_dict["normal"][1], 1)
        for j in axes(X_dict["normal"][1], 2)
            X̄_norm[i, j] = mean(filter(!isnan, X_array_norm[i, j, :]))
            X̄_avg[i, j] = mean(filter(!isnan, X_array_avg[i, j, :]))
        end
    end

    # Create Dictionary to keep track 
    X̄̄ = Dict()

    # Create keys via a loop, data sources is already sorted as well as the reconstructed data
    for i in 1:number_of_dfs
        X̄̄[data_sources[i]] = [X_dict["normal"][i], X_dict["average"][i]]
    end

    X̄̄["consensus"] = [X̄_norm, X̄_avg]

    # Add the immutable part
    for (k, v) in X̄̄
        X̄̄[k][1] = add_multidimensional_immutable(estimator, v[1], grid_cop, measures)
        X̄̄[k][2] = add_multidimensional_immutable(estimator, v[2], grid_cop, measures)
    end

    X̄̄_new = Dict()

    for (k, v) in X̄̄
        X̄̄_new[k] = Dict()
        X̄̄_new[k]["normal"] = undo_functional_treatment(estimator, v[1], measures)
        X̄̄_new[k]["average"] = undo_functional_treatment(estimator, v[2], measures)
    end

    # n_out_of_bounds  = data_diagnostics(copulas, measures)
    if reconstruction_to_show == false
        return X̄̄_new, dε_smoothed
    else
        sub_X̄̄ = Dict()

        for k in unique(["consensus", reconstruction_to_show])
            sub_X̄̄[k] = X̄̄_new[k]
        end

        return sub_X̄̄, dε_smoothed
    end
end


function select_trend(trend, type)
    new_trend = Vector{Vector{Float64}}(undef, length(trend))
    if type == "average"
        for i in eachindex(trend)
            new_trend[i] = vec(nanmean(trend[i], dims=2))
        end
    elseif type == "normal"
        new_trend = deepcopy(trend)
    end
    return new_trend
end

function state_transition(x_t, w_t, u_t, A, B)
    return A * x_t + B * u_t + w_t
end


function add_mean!(X, means)  # Add mean by dataset 
    for j in eachindex(X)
        # name = data_names[j]
        # l_X = length(X)

        # if name ∈ keys(blind_to)
        #     condition = true
        #     i = rand(1:l_X)
        #     # Test next dataset to see if it's also blind
        #     while condition
        #         next_name = data_names[i]
        #         if next_name ∈ keys(blind_to)
        #             condition = true
        #             i = rand(1:l_X)
        #         else
        #             condition = false
        #             v = copy(means[j])
        #             idx = findall(isnan, v)
        #             v[idx] = means[i][idx]
        #             X[j] = ΓF_σ .+ v
        #         end
        #     end
        # else
        X[j] .+= means[j]
        # end
    end

    return X
end

function add_variance!(estimator, X, stds, measures)
    @unpack grid_cop, grid_pcf = estimator

    dimension = length(measures)

    # Copula stuff 
    n_objs = 1 + dimension
    # rows_v             = get_rows_of_copulas(grid_cop, measures)

    cop_part, imm_part = retrieve_cop_and_imm_part(estimator, dimension)
    cop_rows = cop_part - imm_part

    for j in eachindex(X)
        copula_data = X[j][1:cop_rows, :]

        # pcfs stuff 
        pcfs_data = X[j][cop_rows+1:end, :]
        pcf_partition = size(pcfs_data, 1) ÷ dimension
        pcfs = [pcfs_data[I, :] for I in Iterators.partition(axes(pcfs_data, 1), pcf_partition)]

        # standard deviations 
        copula_data .= copula_data .* stds[1]

        start = 2

        for o in start:n_objs
            pcfs[o-1] .= pcfs[o-1] .* stds[o]
        end

        X[j] = vcat(copula_data, vcat(pcfs...))
    end
end


function find_observed_measures_count(pcfs, dimension, gridp; indices=false)
    # split the pcfs into the number of dimensions 
    pcfs_part = [pcfs[I, :] for I in Iterators.partition(axes(pcfs, 1), gridp)]
    ids = []
    # Check each partition to see if it's NaN or not
    count = 0
    for i in eachindex(pcfs_part)
        if any(isnan, pcfs_part[i])
            count += 1
        else
            push!(ids, i)
        end
    end

    obs_dims = dimension - count

    if indices == true
        return obs_dims, ids
    else
        return obs_dims
    end
end


function find_immutable_constant(topology, dimension, gridp, object=:pcf)

    # Based on how many measures are observed 
    local obs_dims

    if object == :pcf
        obs_dims = find_observed_measures_count(topology, dimension, gridp)

        if obs_dims == dimension
            return (1 / (gridp * sqrt(gridp)^(dimension - 2)), 1)
        elseif obs_dims >= 2 && obs_dims < dimension
            return (1 / (gridp * sqrt(gridp)^((dimension - 2) - (dimension - obs_dims))), sqrt(grid)^(dimension - obs_dim))
        elseif obs_dims == 1
            return (NaN, NaN)
        end

    elseif object == :copula
        non_missing_rows = sum((!isnan).(topology))
        IL = Vector{Int64}(undef, dimension)
        n = Vector{Int64}(undef, dimension)
        IL[1], n[1] = (1, 1) # No copula 

        for i in 2:dimension
            IL[i] = gridp + (i - 1) * (gridp - 1)  # immutable length
            n[i] = gridp^i - IL[i]
        end

        obs_dims = 1
        for i in 2:dimension
            if non_missing_rows == n[i]
                obs_dims = i
            end
        end

        if obs_dims == dimension
            return (1 / (gridp * sqrt(gridp)^(dimension - 2)), 1)
        elseif obs_dims >= 2 && obs_dims < dimension
            return (1 / (gridp * sqrt(gridp)^((dimension - 2) - (dimension - obs_dims))), sqrt(gridp)^(dimension - obs_dims))
        elseif obs_dims == 1
            return (NaN, NaN)
        end

    end

    # Return immutable constant based on this criteria 
end

# The issue is we have to reallocate the data to the larger container. The larger container is a 1030 x T. 1000 for the copula, 30 for the marginals. 
# Since the marginals are just vcat underneath, we remove the last the 30.  
function add_multidimensional_immutable(estimator, X, grid_cop, measures)
    # Immutable part only from copula part, so we split 
    D = length(measures)
    IL = grid_cop + (D - 1) * (grid_cop - 1)  # immutable length
    n = grid_cop^D - IL
    T = size(X, 2)

    # Break X into Copulas and Percentile Functions 
    copulas = @view(X[1:n, :])
    pcfs = @view(X[n+1:end, :])

    # Generate containers that includes both the mutable and immutable parts
    cop_cont = zeros(n + IL, T)
    cop_size = tuple([grid_cop for i in 1:D]...)
    temp = zeros(cop_size...)

    f(c) = sum((==(1)).(c.I)) >= length(c.I) - 1
    condition1 = filter(!f, CartesianIndices(size(temp)))
    condition2 = filter(f, CartesianIndices(size(temp)))

    for j in 1:T
        # Find the immutable constant 
        local imm_c, correction, immutable
        if typeof(estimator) <: HistogramEstimator
            imm_c, correction = find_immutable_constant(copulas[:, j], D, grid_cop, :copula) #1 / (gridp * sqrt(gridp)^(D-2))  # this changes empirically since some measures are unobserved, to which this equals => dimension - obs_dims in the exponent for that time t 
            immutable = vcat(imm_c, zeros(IL - 1))
        else
            imm_c = 1
            correction = 1
            immutable = vcat(imm_c, zeros(IL - 1))
        end

        # Fill 
        temp[condition1] .= copulas[:, j] .* correction # correction is for the 2D copula in the 3D setting, it's 1 otherwise 
        temp[condition2] .= immutable
        cop_cont[:, j] .= temp[:]
    end
    return vcat(cop_cont, pcfs)
end


function single_immutable(X, gridp, measures)
    """Add the immutable part to the levels and percentile functions. Only shares has an immutable part."""
    # Define large container of all coefs, mutable and immutable 
    dimension = length(measures)
    immutable_length = 1 * dimension
    mutable_length = size(X, 1)
    T = size(X, 2)
    container = zeros(mutable_length + immutable_length, T)

    # Define locations of mutable and immutable 
    # shares come first, there are 2 measures, so the condition is like for every 10 starting at s =1, fill the row with the constant. until we reach the lenght of measures - 1
    f(c) = c.I[1] ∈ [i for i in [1:gridp:gridp*(dimension-1)+1][1]]
    condition1 = filter(!f, CartesianIndices(size(container)))
    condition2 = filter(f, CartesianIndices(size(container)))
    immutable = 1 / (gridp * sqrt(gridp)^(dimension - 2))

    container[condition1] .= X[:]
    container[condition2] .= immutable
    return container
end

function data_diagnostics(copulas, measures)
    count = 0
    d = length(measures)
    for copula in eachslice(copulas; dims=d)
        count += sum(any(x -> x < 0 || x > 1, copula))
    end
    return count
end


function univariate_idct(X, gridp, measures, time_params)
    """Take X and split into levels and shares. Take the idct for each measure."""

    @unpack tmin, tmax = time_params
    # Divide matrix into measures 
    split_data = [X[I, :] for I in Iterators.partition(axes(X, 1), gridp)]
    T = size(X, 2)
    D = length(measures)

    # Subsetting aggregate data 
    # gdp_series[!, "year"]      = Dates.year.(gdp_series[!, "time"])

    # gdp_series[!, "unit_time"] = Dates.quarterofyear.(gdp_series[!, "time"])
    # filter!(row -> row.year >= tmin["year"], gdp_series)
    # filter!(row -> row.year <= tmax["year"], gdp_series)

    # Now we want the aggregate data to start at q3
    # year_data  = filter(row -> row.year == tmin["year"], gdp_series)
    # first_qtrs = length(year_data[:, :unit_time])
    # q3         = first_qtrs - 1
    # gdp_series = gdp_series[q3:end, :]  # because of the lag 

    # Take idct and undo transformations 
    for i in axes(split_data, 1)
        split_data[i] .= idct(split_data[i], 1)
    end
    #     if i <= D #|| (i > 2*D && i <= 3*D) || (i > 4*D && i <= 5*D)
    #         for t in 1:T
    #             split_data[i][:, t] .= split_data[i][:, t] .* gdp_series[t, :gdp]  # TODO: check that this series date matches split data date 
    #         end
    #     elseif i > D # (i > D && i <= 2*D) || (i > 3*D && i <= 4*D) || (i > 5*D && i <= 6*D)
    #         for t in 1:T
    #             split_data[i][:, t] .= (abs.(split_data[i][:, t]) ./ sign.(split_data[i][:, t])).^3 .* gdp_series[t, :real_gdp_pc]   
    #         end
    #     end
    # end
    return split_data
end


function undo_functional_treatment(estimator, X, measures)
    """Break the data into copula and percentile functions. Then take the Inverse DCT of each them, at each point in time."""
    @unpack grid_pcf, grid_cop = estimator


    # Some parameters 
    d = length(measures)
    T = size(X, 2)
    n = grid_cop^d

    @views copulas = X[1:n, :]
    @views pcfs = X[n+1:end, :]

    cop_size = tuple([grid_cop for i in 1:d]...)
    colons = ntuple(_ -> (:), d)

    # Convert vectors to matrices  
    cop_array = zeros(cop_size..., T)
    fill!(cop_array, NaN)
    pcf_array = zeros(grid_pcf * d, T)
    split_pcfs = [pcfs[I, :] for I in Iterators.partition(axes(pcfs, 1), grid_pcf)]

    # Perform idct for each object 
    if typeof(estimator) <: HistogramEstimator
        for m in eachindex(split_pcfs)
            split_pcfs[m] = reverse_inverse_hyperbolic_sine(idct(split_pcfs[m], 1))
        end
    else
        nothing
        # for m in eachindex(split_pcfs)
        #     split_pcfs[m]          = reverse_inverse_hyperbolic_sine(split_pcfs[m])
        # end
    end


    # The CartesianIndices for the copula
    cop_ci = CartesianIndices(tuple([grid_cop for _ in 1:d]...))

    # Perform Inverse DCT time period by time period 
    for t in 1:T
        # Percentiles first 
        pcf_array[:, t] .= vcat([split_pcfs[m][:, t] for m in eachindex(split_pcfs)]...) #TODO: redundant?
        obs_dims, id_of_observed = find_observed_measures_count(pcfs[:, t], d, grid_pcf; indices=true)

        if obs_dims <= 1
            nothing
        else
            # Find the cartesian indices for the observed measures
            sn = copula_case_validation(id_of_observed, measures) # subset of n 

            # Subset to the observed measures
            vec_of_indices = cop_ci[sn...][:]
            id_of_cop = findall(x -> x in vec_of_indices, cop_ci[:])

            scop_size = tuple([grid_cop for i in 1:obs_dims]...)
            cop_array[sn..., t] = undo_copula_treatment(reshape(copulas[id_of_cop, t], scop_size), estimator) # TODO: issue, immutable contains 
        end
    end

    return [cop_array, pcf_array]
end


function undo_functional_treatment_validation(estimator, X, grid_pcf, grid_cop, measures, gdp_series, time_dict, freq_type)
    """Break the data into copula and percentile functions. Then take the Inverse DCT of each them, at each point in time."""

    # some params 
    D = length(measures)
    T = size(X, 2)
    n = grid_cop^D
    cop_size = tuple([grid_cop for i in 1:D]...)
    colons = ntuple(_ -> (:), D)

    # Get pcfs 
    @views copulas = X[1:n, :]
    @views pcfs = X[n+1:end, :]

    # Agg data 
    gdp_series[!, "date"] = QuarterlyDate.(gdp_series[!, "time"])
    correction_vec = zeros(T, D)

    # Storing real gdp pc data
    col_names = sort([meas * "_per_hh" for meas in measures])

    if freq_type == "year"
        # Collects real gdp pc for all years, not just observed     
        for (i, y) in enumerate(sort(collect(keys(time_dict))))
            p_data = filter(row -> row.date == QuarterlyDate(y, 4), gdp_series)
            correction_vec[i, :] = Matrix(p_data[:, col_names])[:] #float.()[1]
        end

    elseif freq_type == "quarter"
        count = 1
        for y in sort(collect(keys(time_dict))) # 
            for p in sort(time_dict[y])
                p_data = filter(row -> row.date == QuarterlyDate(y, p), gdp_series)
                correction_vec[count, :] = Matrix(p_data[:, col_names])[:]
                count += 1
            end
        end
    end

    # Convert vectors to matrices  
    cop_array = zeros(cop_size..., T)
    fill!(cop_array, NaN)
    pcf_array = zeros(grid_pcf * D, T)
    split_pcfs = [pcfs[I, :] for I in Iterators.partition(axes(pcfs, 1), grid_pcf)] # split by measure 
    cop_ci = CartesianIndices(tuple([grid_cop for _ in 1:D]...))

    # Perform inverse DCT on each percentile function separately 
    if typeof(estimator) <: HistogramEstimator
        for m in eachindex(split_pcfs)
            split_pcfs[m] = reverse_inverse_hyperbolic_sine(idct(split_pcfs[m], 1))
        end

        # Perform Inverse DCT time period by time period 
        for t in 1:T
            pcf_array[:, t] .= vcat([split_pcfs[m][:, t] .* correction_vec[t, m] for m in eachindex(split_pcfs)]...)

            # Find how many observed measures there are 
            obs_dims, id_of_observed = find_observed_measures_count(pcfs[:, t], D, grid_pcf; indices=true)
            if obs_dims <= 1
                nothing
            else
                sn = copula_case_validation(id_of_observed, measures) # subset of n 
                vec_of_indices = cop_ci[sn...][:]
                id_of_cop = findall(x -> x in vec_of_indices, cop_ci[:])

                scop_size = tuple([grid_cop for i in 1:obs_dims]...)

                cop_array[sn..., t] = undo_copula_treatment(reshape(copulas[id_of_cop, t], scop_size), estimator) # TODO: issue, immutable contains 
            end
        end

        return [cop_array, pcf_array]
    else
        # Perform Inverse DCT time period by time period 
        for t in 1:T
            pcf_array[:, t] .= vcat([split_pcfs[m][:, t] for m in eachindex(split_pcfs)]...)

            # Find how many observed measures there are 
            obs_dims, id_of_observed = find_observed_measures_count(pcfs[:, t], D, grid_pcf; indices=true)

            if obs_dims <= 1
                nothing
            else
                # Find the cartesian indices for the observed measures
                sn = copula_case_validation(id_of_observed, measures) # subset of n 

                # Subset to the observed measures
                vec_of_indices = cop_ci[sn...][:]
                id_of_cop = findall(x -> x in vec_of_indices, cop_ci[:])

                scop_size = tuple([grid_cop for i in 1:obs_dims]...)
                cop_array[sn..., t] = undo_copula_treatment(reshape(copulas[id_of_cop, t], scop_size), estimator) # TODO: issue, immutable contains 
            end
        end

        return [cop_array, pcf_array]
    end
end


function perform_proof_of_concept_reconstruction(df, source, years, gdp_series, time_p, freq_type, time_dict_k, model_options; max_order::Union{Int,Nothing} = nothing)

    @unpack lags, measures, freq, plot_proof, pca_perspective, estimator, tag = model_options
    @unpack grid_pcf, grid_cop = estimator
    @unpack tmin, tmax, tot_periods = time_p

    # Resolve max_order: default to full polynomial order from estimator
    max_order = isnothing(max_order) ? grid_pcf - 1 : min(max_order, grid_pcf - 1)

    dimension = length(measures)

    # So, actual ID corresponds to the index of the estimation 1:T
    estimation_id = generate_correct_indices(years, freq_type, freq, tmin, time_dict_k)  # for adding back trend to the correct index


    # Mask of actual ID that corresponds to fully observed data --- needed for adding trend back later 
    fully_obs_id = []
    for i in axes(df, 2)
        if sum(isnan.(df[:, i])) == 0
            push!(fully_obs_id, i)
        end
    end

    estimation_id_full = estimation_id[fully_obs_id]

    if plot_proof == true
        # This is the original data 
        RD = deepcopy(df)

        # Treating 'df', getting important factors and reconstructing #

        # Perform detrending on the data 
        βs, df, trend = perform_detrending(df, time_p, years, freq, freq_type, time_dict_k)

        # Perform PCA 
        pool, means, stds = perform_standardization([df], estimator, dimension, measures)
        proj, pcs, M, _ = perform_pca(pool, measures, :functional_data, tag)

        # dealing with observed data and reconstructed, # Split by object, multiply by stds, reform object 
        X̃ = [proj * pcs]
        add_variance!(estimator, X̃, stds, measures)
        add_mean!(X̃, means)

        # Add trend back to only fully observed data, since reconstruction only depends on the fully observed data
        X = X̃[1]
        for (i, t) in enumerate(estimation_id_full)
            X[:, i] .= X[:, i] .+ trend[:, t] #.+ βs[:,1] .+ βs[:,2] .* t .+ βs[:,3] .* (2 .* t.^2 ./ tot_periods) # 
        end

        d_container = add_multidimensional_immutable(estimator, RD, grid_cop, measures)
        r_container = add_multidimensional_immutable(estimator, X, grid_cop, measures)

        # The reconstruction needs to be the same size of the data. How do we do this?
        r_container_all = fill(NaN, size(d_container))
        r_container_all[:, fully_obs_id] .= r_container # reconstruction only exists for the fully observed periods

        d_data_vector = undo_functional_treatment_validation(estimator, d_container, grid_pcf, grid_cop, measures, gdp_series, time_dict_k, freq_type)
        r_data_vector = undo_functional_treatment_validation(estimator, r_container_all, grid_pcf, grid_cop, measures, gdp_series, time_dict_k, freq_type)

        # subsetting gdp series 
        gdp_series[!, "date"] = QuarterlyDate.(gdp_series[!, "time"])
        filter!(row -> row.date >= QuarterlyDate(tmin["year"], tmin["quarter"]), gdp_series)
        filter!(row -> row.date <= QuarterlyDate(tmax["year"], tmax["quarter"]), gdp_series)


        # From weights, generate data 
        if typeof(estimator) <: SeriesEstimator
            for dv in [d_data_vector, r_data_vector]
                @unpack integral_pcf_grid, integral_cop_grid = estimator

                grid_size_data_pcf = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf
                grid_size_data_cop = typeof(estimator) <: SeriesEstimator ? integral_cop_grid : grid_cop

                # Generate container to store the data of choice 
                T = size(dv[2], 2)
                new_data_pcf = [zeros(grid_size_data_pcf, T) for _ in 1:dimension]

                # Fill 'new_data_pcf' with NaN
                for m in eachindex(new_data_pcf)
                    new_data_pcf[m] .= NaN
                end

                grid_points_pcf = select_grid_points(grid_size_data_pcf)
                intervals = vcat([0.0] .+ 1e-6, grid_points_pcf)

                # split the pcfs by measure
                split_pcfs = [dv[2][I, :] for I in Iterators.partition(axes(dv[2], 1), grid_pcf)]

                # Agg data 
                correction_vec = zeros(T, dimension)

                # Storing real gdp pc data
                col_names = sort([meas * "_per_hh" for meas in measures])

                if freq_type == "year"
                    # Collects real gdp pc for all years, not just observed     
                    for (i, y) in enumerate(sort(collect(keys(time_dict_k))))
                        p_data = filter(row -> row.date == QuarterlyDate(y, 4), gdp_series)
                        correction_vec[i, :] = Matrix(p_data[:, col_names])[:] #float.()[1]
                    end

                elseif freq_type == "quarter"
                    count = 1
                    for y in sort(collect(keys(time_dict_k))) # 
                        for p in sort(time_dict_k[y])
                            p_data = filter(row -> row.date == QuarterlyDate(y, p), gdp_series)
                            correction_vec[count, :] = Matrix(p_data[:, col_names])[:]
                            count += 1
                        end
                    end
                end

                # Integration (max_order truncates Legendre series for smoothing)
                integrate_quantile_functions!(new_data_pcf, split_pcfs, grid_pcf, intervals, correction_vec; max_order = max_order)
                dv[2] = vcat([new_data_pcf[m] for m in eachindex(new_data_pcf)]...)

                # Integrate densities
                dv[1] = generate_copula_densities(deepcopy(dv[1]), measures, grid_size_data_cop) #FIXME:
            end
        end

        # Generate shares and levels
        d_levels, d_shares = generate_shares_levels(d_data_vector[2], model_options, gdp_series)
        r_levels, r_shares = generate_shares_levels(r_data_vector[2], model_options, gdp_series)

        # Create dictionary of all data 
        d_data_dict = create_time_series_dictionary([d_data_vector..., d_levels, d_shares], estimator, measures)
        r_data_dict = create_time_series_dictionary([r_data_vector..., r_levels, r_shares], estimator, measures)

        # For copula, must be greater than 2 dimensions
        if dimension >= 2
            gen_proof_of_concept_copulas(d_data_dict, r_data_dict, source, measures, model_options, time_p, fully_obs_id, estimation_id_full)
        end

        # Generate proof of concept figures for percentile functions 
        gen_proof_of_concept_figure(d_data_dict, r_data_dict, model_options, time_p, estimation_id, sort(collect(keys(time_dict_k))), source)
    else
        nothing
    end
end

# Idea now is to compute the integrals for each decile interval, then multiply the coefficients thereafter

function check_observed_measures(A, measures)
    observed_measures = []

    dims = size(A)
    for i in 1:length(dims)
        slice_indices = ntuple(k -> k == i ? 1 : :, length(dims))
        meas_indices = findall(x -> x == Colon(), slice_indices)
        if all(isfinite.(A[slice_indices...]))
            push!(observed_measures, measures[meas_indices]...)
        end
    end

    return unique(observed_measures)
end


function generate_copula_densities(X, measures, grid_size_data_cop; given_integrals=false)
    # In the case of the reconstructions, if the dataset ONCE had 3 measures observed, then predictions will be made for all 3 

    # 'd' has to be based on the coefficients that are obsered <==> the observed measures 
    T = size(X)[end]
    D = length(size(X)[1:end-1])

    # This has to be of the same dimensionality of the estimation 
    obs_cop_size = tuple([grid_size_data_cop for i in 1:D]...) # doesnt affect SeriesEstimator
    new_dv = fill!(Array{Float64}(undef, obs_cop_size..., T), NaN)  #zeros(cop_size..., T)
    # grid_cop      = size(X, 1)
    # cop_size      = tuple([grid_cop for i in 1:D]...) 
    colons = copula_case(measures, measures)

    # Grid points 
    x = select_grid_points(grid_size_data_cop)
    x[end] = x[end] - 1e-6  # to avoid numerical issues with the last point

    # Compute the integrals first
    N = size(X, 1) - 1
    integrals = given_integrals == false ? precompute_integrals(N, x) : given_integrals

    # Threads.@threads 
    for t in 1:T
        # At each period, check what is observed 
        obs_meas = check_observed_measures(X[colons..., t], measures)
        cop_ind = copula_case(obs_meas, measures)
        obs_d = length(obs_meas)

        if obs_d <= 1
            nothing
        else
            # XX no longer needed — copula_cdf_estimator now uses array indices directly
            # XX = obs_d == 2 ? [[x[i], x[j]] for i in eachindex(x), j in eachindex(x)] : [[x[i], x[j], x[k]] for i in eachindex(x), j in eachindex(x), k in eachindex(x)]

            cop_w = X[cop_ind..., t]
            cop_cdf = obs_d == 2 ? [copula_cdf_estimator(cop_w, integrals, [i, j]) for i in eachindex(x), j in eachindex(x)] : [copula_cdf_estimator(cop_w, integrals, [i, j, k]) for i in eachindex(x), j in eachindex(x), k in eachindex(x)]
            new_dv[cop_ind..., t] .= cdf_to_pdf(cop_cdf)
        end
    end

    return new_dv
end


# SINGLE ITERATION ... I KNOW, DUMB, BUT FOR EFFICIENCY REASONS 
function generate_copula_density(X, XX, id_x, integrals, measures, obs_meas, grid_size_data_cop)
    # Note: XX is kept in signature for compatibility but no longer used (integrals are array-indexed)
    D = length(size(X))

    # This has to be of the same dimensionality of the estimation
    obs_cop_size = tuple([grid_size_data_cop for i in 1:D]...) # doesnt affect SeriesEstimator
    new_dv = fill!(Array{Float64}(undef, obs_cop_size...), NaN)

    cop_ind = copula_case(obs_meas, measures)
    obs_d = length(obs_meas)

    if obs_d <= 1
        return NaN
    else
        cop_w = X[cop_ind...]
        # Old: passed float grid points via XX
        # cop_cdf = obs_d == 2 ? [copula_cdf_estimator(cop_w, integrals, [XX[i, j][1], XX[i, j][2]]) for i in id_x, j in id_x] : [copula_cdf_estimator(cop_w, integrals, [XX[i, j, k][1], XX[i, j, k][2], XX[i, j, k][3]]) for i in id_x, j in id_x, k in id_x]
        # New: pass integer indices directly (integrals is now an array, not a Dict)
        cop_cdf = obs_d == 2 ? [copula_cdf_estimator(cop_w, integrals, [i, j]) for i in id_x, j in id_x] : [copula_cdf_estimator(cop_w, integrals, [i, j, k]) for i in id_x, j in id_x, k in id_x]
        new_dv[cop_ind...] .= cdf_to_pdf(cop_cdf)

        return new_dv
    end
end


function precompute_integrals(N, x)
    # Old dict-based version:
    # integrals = Dict{Tuple{Int,Float64},Float64}()
    # for m in 0:N
    #     for u in x
    #         integrals[(m, u)] = integrate_legendre_polynomial(m, u)
    #     end
    # end
    # return integrals

    # Array-based version: integral_array[m+1, u_idx] for fast indexing
    integral_array = zeros(N + 1, length(x))
    for (ui, u) in enumerate(x)
        for m in 0:N
            integral_array[m + 1, ui] = integrate_legendre_polynomial(m, u)
        end
    end
    return integral_array
end


function copula_cdf_estimator(X, integrals, u_indices)
    C_N = 0.0

    # Dimension of copula
    d = length(size(X))

    # Order of the object
    N = size(X, 1) - 1

    # Ranges for the object
    ranges = [(0:N) for _ in 1:d]

    # All possible orders of the object
    m_combos = collect(Iterators.product(ranges...))

    # Look over each weight <==> looping over each m_combos
    # Now uses array indexing: integrals[m+1, u_idx] instead of Dict lookup integrals[(m, u)]
    for ci in CartesianIndices(m_combos)
        m = Tuple(m_combos[ci])
        rho_m = X[ci]
        product = 1.0

        @inbounds for j in 1:d
            product *= integrals[m[j] + 1, u_indices[j]]
        end

        C_N += rho_m * product
    end

    return C_N
end




function I_m(m, u)
    return integrate_legendre_polynomial(m, u)
end

function integrate_legendre_polynomial(m, u)
    if m == 0
        return u
    else
        integral_cop, _ = quadgk(u -> Q_m(m, u), 0, u, rtol=1e-8)

        return integral_cop
    end
end


function integrate_quantile_functions!(new_data_pcf, split_pcfs, grid_pcf, intervals, agg_corr; max_order::Int = grid_pcf - 1)
    use_order = min(max_order, grid_pcf - 1)
    for m in eachindex(new_data_pcf)
        for t in axes(new_data_pcf[m], 2)
            if all(isnan.(split_pcfs[m][:, t]))
                new_data_pcf[m][:, t] .= NaN
            else
                for i in axes(new_data_pcf[m], 1)
                    # Using coefs, generate pcf function and then integrate pcf function over diff. intervals
                    integral, _ = quadgk(u -> reverse_inverse_hyperbolic_sine(eval_quantile_function(split_pcfs[m][:, t], use_order, u))[1] * agg_corr[t, m], intervals[i], intervals[i+1], rtol=1e-8)

                    # Undo treatment of data => gives us average quantile within the interval
                    new_data_pcf[m][i, t] = integral / (intervals[i+1] - intervals[i])
                end
            end
        end
    end
end


function interval_time_correction(confidence_intervals, periods, time_p, k, freq)
    """
     The goal here is to subset the confidence intervals to the dates of estimation. 
    """

    @unpack tmin, tmax, tot_years, tot_periods, time_dict, freq_type = time_p

    start = 1
    confidence_interval_dict = Dict()

    # Intervals of the dataset 
    lower_or_upper = collect(keys(confidence_intervals))

    for bound in lower_or_upper
        confidence_interval_dict[bound] = Dict()
        measuresⱼ = setdiff(collect(keys(confidence_intervals[bound])), ["copula"])

        for m in measuresⱼ
            confidence_interval_dict[bound][m] = Dict()
            objs = collect(keys(confidence_intervals[bound][m]))

            for o in objs
                confidence_interval_dict[bound][m][o] = Dict()
                confidence_interval_dict[bound][m][o] = order_measures(confidence_intervals[bound][m][o], periods, tot_periods, freq, freq_type[k], tmin, time_dict[k], start)
            end
        end

        try
            confidence_interval_dict[bound]["copula"] = order_measures(confidence_intervals[bound]["copula"], periods, tot_periods, freq, freq_type[k], tmin, time_dict[k], start)
        catch ee
            println("don't forget copula intervals")
        end
    end

    return confidence_interval_dict
end


function find_observed_measures_using_chebycoefs(pcf_coefficients, t, measures, estimator)
    """From the coefficients, we know which measures are observed. We can find the observed measures by checking if the coefficients are not NaN."""
    @unpack granularity_pcf = estimator

    # Container
    observed_measures = []

    for x in eachindex(pcf_coefficients)
        if all(!isnan, pcf_coefficients[x][:, t])
            push!(observed_measures, measures[x])
        end
    end

    return observed_measures
end


function return_subcopula_sizes(d, combos)
    subcop_amounts = []

    for d in 2:d
        append!(subcop_amounts, length(filter(x -> length(x) == d, combos)))
    end

    return subcop_amounts
end


function undo_series_estimator(coeffs, estimator, obs_meas)
    @unpack grid_cop, granularity_cop, grid_type_cop = estimator

    obs_dim = length(obs_meas)
    cop_size = tuple([grid_cop for i in 1:obs_dim]...)

    # Weights in matrix form 
    w = reshape(coeffs, cop_size)

    # For the evaluation of the chebyshev polynomial 
    dom = [ones(obs_dim) zeros(obs_dim)]' # domain 
    O = size(w) .- 1 # order
    # g        = nodes(granularity_cop, :chebyshev_nodes, [1.0, 0.0]).points 
    g = select_grid_points(granularity_cop, grid_type_cop)
    g2 = Iterators.product(fill(g, obs_dim)...)
    X = [collect(point) for point in g2] # A matrix of vectors, where each vector is a N-Dimensional point on the copula grid

    function cheby_interp(x)
        yhat = chebyshev_evaluate(w, x, O, dom)
        return yhat
    end

    return cheby_interp.(X)
end



function undo_series_estimator_setup(coefficients, estimator, dimension, measures)
    @unpack grid_cop, grid_pcf, grid_type_cop, grid_type_pcf = estimator


    # First separate the copula coefficients from the percentile function coefficients
    cop_part, imm_part = retrieve_cop_and_imm_part(estimator, dimension)
    cop_rows = cop_part - imm_part
    cop_coefs = coefficients[1:cop_rows, :]

    pcf_coefs = coefficients[cop_rows+1:end, :]
    pcf_partition = size(pcf_coefs, 1) ÷ dimension
    pcfs_coefs2 = [pcf_coefs[I, :] for I in Iterators.partition(axes(pcf_coefs, 1), pcf_partition)]

    T = size(cop_coefs, 2)
    pcfs = [zeros(granularity_pcf, T) for _ in 1:dimension]

    # Multiply coefficients with the basis functions and undo the 'inverse_hyperbolic_sine' transformation
    for i in eachindex(pcfs)
        for t in axes(pcfs[i], 2)
            poly_int = ChebyshevT([pcfs_coefs2[i][:, t]...])
            pcfs[i][:, t] = reverse_inverse_hyperbolic_sine(poly_int.(x_values_cheb)) # TODO: this is not the original data. 
        end
    end

    # For copula, find the copula to reconstruct and then just multiply with the basis functions
    cop̂ = fill(NaN, tuple([granularity_cop for i in 1:dimension]...)..., T)


    # To get the number of sub-copulas 
    all_combos = generate_unique_combinations((1:dimension))

    # find all combos of length D, D-1, D-2, ... 2
    subcopula_amounts = return_subcopula_sizes(dimension, all_combos)

    for t in 1:T
        if all(isnan, cop_coefs[:, t])
            nothing

        else
            obs_meas = find_observed_measures_using_chebycoefs(pcfs, t, measures, estimator)
            if length(obs_meas) <= 1
                nothing

            elseif sort(obs_meas) == sort(measures) # all(!isnan, cop_coefs[:, t])
                # Extract only the largest copula
                largest_cop = cop_coefs[end-grid_cop^dimension+1:end, t]
                cop_ind = copula_case(obs_meas, measures)
                cop̂[cop_ind..., t] .= undo_series_estimator(largest_cop, estimator, obs_meas)

            else # they are a subcopula 
                # Get indices that correspond to measures that are observed e.g., [2,3]
                id_wrt_measures = findall(x -> x ∈ obs_meas, measures)
                cop_ind = copula_case(obs_meas, measures)

                # In the list of combinations, find where [2,3] lies, lets say 3rd position
                id_of_combination = findall(x -> x == id_wrt_measures, all_combos)[1]

                # See how many sub-copulas there are before [2,3]
                cumsum_sub_amounts = cumsum(subcopula_amounts)
                id_of_min_subcopula_level = findfirst(x -> id_of_combination > x, cumsum_sub_amounts)

                if typeof(id_of_min_subcopula_level) == Nothing
                    # Then you know it is in the first set of sub-copulas, which are 2D
                    start = (id_of_combination - 1) * grid_cop^2 + 1
                    term = id_of_combination * grid_cop^2
                    cop̂[cop_ind..., t] .= undo_series_estimator(cop_coefs[start:term, t], estimator, obs_meas)

                else
                    # Find 
                    r = 0
                    o = 2

                    for j in 1:id_of_min_subcopula_level # 1,2,3
                        r += grid_cop^(o) * subcopula_amounts[j]
                        o += 1
                    end

                    # This tells me how much above last level 
                    diff = id_of_combination - cumsum_sub_amounts[id_of_min_subcopula_level]

                    start = r + grid_cop^o * (diff - 1) + 1
                    term = r + grid_cop^o * diff

                    cop̂[cop_ind..., t] .= undo_series_estimator(cop_coefs[start:term, t], estimator, obs_meas)
                end
            end
        end
    end
    # For each point in time, divie the cop by the maximum element 
    cop_ind = copula_case(measures, measures)

    for t in 1:T
        cop̂[cop_ind..., t] = cop̂[cop_ind..., t] ./ nanmaximum(cop̂[cop_ind..., t])
    end

    return [cop̂, vcat(pcfs...)]
end


function move_data_to_dict(df, periods, gdp_series, model_options, time_p, data_source, k; max_order::Union{Int,Nothing} = nothing)
    @unpack estimator, lags, measures, freq, tag = model_options
    @unpack tmin, tmax, tot_years, tot_periods, time_dict, freq_type = time_p
    @unpack grid_pcf, grid_cop = estimator

    # Resolve max_order: default to full polynomial order from estimator
    max_order = isnothing(max_order) ? grid_pcf - 1 : min(max_order, grid_pcf - 1)

    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid, integral_cop_grid = estimator
    end
    grid_size_data_pcf = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf
    grid_size_data_cop = typeof(estimator) <: SeriesEstimator ? integral_cop_grid : grid_cop

    # Undo Transformations 
    D = length(measures)
    RD = copy(df)

    d_container = add_multidimensional_immutable(estimator, RD, grid_cop, measures)
    d_data_vector = undo_functional_treatment_validation(estimator, d_container, grid_pcf, grid_cop, measures, gdp_series, time_dict[k], freq_type[k])

    if typeof(estimator) <: SeriesEstimator
        # Generate container to store the data of choice 
        T = size(d_data_vector[2], 2)
        new_data_pcf = [zeros(grid_size_data_pcf, T) for _ in 1:D]
        grid_points_pcf = select_grid_points(grid_size_data_pcf)
        intervals = vcat([0.0] .+ 1e-6, grid_points_pcf)
        # intervals = vcat([0.0], grid_points_pcf)

        # split the pcfs by measure
        split_pcfs = [d_data_vector[2][I, :] for I in Iterators.partition(axes(d_data_vector[2], 1), grid_pcf)]

        # Agg data 
        gdp_series[!, "date"] = QuarterlyDate.(gdp_series[!, "time"])
        correction_vec = zeros(T, D)

        # Storing real gdp pc data
        col_names = sort([meas * "_per_hh" for meas in measures])

        if freq_type[k] == "year"
            # Collects real gdp pc for all years, not just observed     
            for (i, y) in enumerate(sort(collect(keys(time_dict[k]))))
                p_data = filter(row -> row.date == QuarterlyDate(y, 4), gdp_series)
                correction_vec[i, :] = Matrix(p_data[:, col_names])[:] #float.()[1]
            end

        elseif freq_type[k] == "quarter"
            count = 1
            for y in sort(collect(keys(time_dict[k]))) # 
                for p in sort(time_dict[k][y])
                    p_data = filter(row -> row.date == QuarterlyDate(y, p), gdp_series)
                    correction_vec[count, :] = Matrix(p_data[:, col_names])[:]
                    count += 1
                end
            end
        end

        integrate_quantile_functions!(new_data_pcf, split_pcfs, grid_pcf, intervals, correction_vec; max_order = max_order)

        # Share-group integration: bot50/mid40/top10
        share_spec = [0.5, 0.4, 0.1]
        share_intervals = vcat([0.0 + 1e-6], cumsum(share_spec)[1:end-1], [1.0 - 1e-6])
        share_data_pcf = [zeros(length(share_spec), T) for _ in 1:D]
        integrate_quantile_functions!(share_data_pcf, split_pcfs, grid_pcf, share_intervals, correction_vec; max_order = max_order)

        # We need to generate new data_pcf, which are the average quantiles over the intervals
        # for m in eachindex(new_data_pcf)
        #     for t in axes(new_data_pcf[m], 2)
        #         if all(isnan.(split_pcfs[m][:, t]))
        #             new_data_pcf[m][:, t] .= NaN
        #         else
        #             for i in axes(new_data_pcf[m], 1)
        #                 # Using coefs, generate pcf function and then integrate pcf function over diff. intervals 
        #                 integral, _            = quadgk(u -> eval_quantile_function(split_pcfs[m][:, t], grid_pcf-1, u), intervals[i], intervals[i+1], rtol=1e-8)

        #                 # Undo treatment of data => gives us average quantile within the interval 
        #                 new_data_pcf[m][i, t]  = reverse_inverse_hyperbolic_sine(integral / (intervals[i+1] - intervals[i]))[1] .* correction_vec[t, m]
        #             end
        #         end
        #     end
        # end

        d_data_vector[2] = vcat([new_data_pcf[m] for m in eachindex(new_data_pcf)]...)
    end

    # Break apart the data and assign the correct indices based on quarterly data
    cop_dim = size(d_data_vector[1])
    Πcop_dim = prod(cop_dim[1:end-1])
    T = cop_dim[end]

    # Reshape data to prepare ordering 
    copulas = reshape(d_data_vector[1], (Πcop_dim, T))

    start = 1
    copulas = order_measures(copulas, periods, tot_periods, freq, freq_type[k], tmin, time_dict[k], start)
    pcfs = order_measures(d_data_vector[2], periods, tot_periods, freq, freq_type[k], tmin, time_dict[k], start)

    gdp_series[!, "date"] = QuarterlyDate.(gdp_series[!, "time"])
    filter!(row -> row.date >= QuarterlyDate(tmin["year"], tmin["quarter"]), gdp_series)
    filter!(row -> row.date <= QuarterlyDate(tmax["year"], tmax["quarter"]), gdp_series)

    levels, shares = generate_shares_levels(pcfs, model_options, gdp_series)
    d_data_dict = create_time_series_dictionary([copulas, pcfs, levels, shares], estimator, measures)

    # Inject share-group means (bot50/mid40/top10) into the dictionary
    if typeof(estimator) <: SeriesEstimator && @isdefined(share_data_pcf)
        # Order the share data onto the full quarterly timeline (same as pcfs)
        share_pcf_cat = vcat([share_data_pcf[m] for m in eachindex(share_data_pcf)]...)
        share_pcf_ordered = order_measures(share_pcf_cat, periods, tot_periods, freq, freq_type[k], tmin, time_dict[k], 1)
        n_share = length(share_spec)
        share_indices = [I for I in Iterators.partition(axes(share_pcf_ordered, 1), n_share)]

        share_labels = ["bot50", "mid40", "top10"]
        for (i, meas) in enumerate(sort(measures))
            for (si, sl) in enumerate(share_labels)
                d_data_dict[meas]["quantiles"]["common series"][sl] = share_pcf_ordered[share_indices[i], :][si, :]
            end
        end
    end


    # Exporting raw data
    export_raw_data(d_data_dict, estimator, data_source, measures, time_p, tag)

    return d_data_dict
end


function correct_pcfs_with_aggregates(pcfs, grid_pcf, measures, gdp_series, time_dict, freq_type)
    # some params 
    D = length(measures)
    T = size(pcfs, 2)

    # Agg data 
    gdp_series[!, "date"] = QuarterlyDate.(gdp_series[!, "time"])
    correction_vec = zeros(T, D)

    # Storing real gdp pc data
    col_names = sort([meas * "_per_hh" for meas in measures])

    if freq_type == "year"
        # Collects real gdp pc for all years, not just observed     
        for (i, y) in enumerate(sort(collect(keys(time_dict))))
            p_data = filter(row -> row.date == QuarterlyDate(y, 4), gdp_series)
            correction_vec[i, :] = Matrix(p_data[:, col_names])[:] #float.()[1]
        end

    elseif freq_type == "quarter"
        count = 1
        for y in sort(collect(keys(time_dict))) # 
            for p in sort(time_dict[y])
                p_data = filter(row -> row.date == QuarterlyDate(y, p), gdp_series)
                correction_vec[count, :] = Matrix(p_data[:, col_names])[:]
                count += 1
            end
        end
    end

    pcf_array = zeros(grid_pcf * D, T)
    split_pcfs = [pcfs[I, :] for I in Iterators.partition(axes(pcfs, 1), grid_pcf)] # split by measure 

    # Perform Inverse DCT time period by time period 
    for t in 1:T
        pcf_array[:, t] .= vcat([split_pcfs[m][:, t] .* correction_vec[t, m] for m in eachindex(split_pcfs)]...)
    end

    return pcf_array
end



function store_functional_data(files, dfs, gdp_series, model_options, time_p; max_order::Union{Int,Nothing} = nothing)
    @unpack year_vec, time_dict = time_p
    func_dict = Dict()
    data_sources = sort(collect(keys(files)))

    # Stores the data
    for (k, df) in enumerate(dfs)
        func_dict[data_sources[k]] = move_data_to_dict(df, year_vec[k], gdp_series, model_options, time_p, data_sources[k], k; max_order = max_order)
    end

    return func_dict
end

function store_confidence_intervals(files, confidence_intervals, freq, time_p)
    @unpack year_vec, time_dict = time_p
    data_sources = sort(collect(keys(files)))

    # Stores the confidence intervals
    for j in collect(keys(confidence_intervals))
        k = findall(x -> x == j, data_sources)[1]
        confidence_intervals[j] = interval_time_correction(confidence_intervals[j], year_vec[k], time_p, k, freq)
    end

    return confidence_intervals
end



function collect_indices_for_quantile_groups(grid_choice_pcf)
    if grid_choice_pcf == 10
        return [1:5, 6:9, 10]
    elseif grid_choice_pcf == 5
        return [1:2, 3:4, 5]
    elseif grid_choice_pcf == 20
        return [1:10, 11:18, 19:20]
    elseif grid_choice_pcf == 100
        return [1:50, 51:90, 91:100]

    end
end



function pick_grid_for_confidence_intervals(estimator)

    @unpack grid_pcf, grid_cop = estimator

    grid_choice_pcf = copy(grid_pcf)
    grid_choice_cop = copy(grid_cop)

    return grid_choice_pcf, grid_choice_cop
end



function specify_bootstrap_container(measures, objects, estimator, year_vec, series, source, draws)
    @unpack grid_pcf, grid_cop = estimator
    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid, integral_cop_grid = estimator
    end
    grid_choice_pcf = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf
    grid_choice_cop = typeof(estimator) <: SeriesEstimator ? integral_cop_grid : grid_cop

    dimension = length(measures)
    T = length(year_vec)

    # Containers 
    sub_boot_dict = Dict()

    for m in measures
        sub_boot_dict[m] = Dict()
        for o in objects
            if o == "quantiles"
                sub_boot_dict[m][o] = zeros(grid_choice_pcf, draws, T)  #For quarterly data, this is still fine because I repeat the year for each quarter
                fill!(sub_boot_dict[m][o], NaN)
            else
                sub_boot_dict[m][o] = zeros(length(series), draws, T)
                fill!(sub_boot_dict[m][o], NaN)
            end
        end
    end

    # Copula 
    cop_n = grid_choice_cop^dimension
    sub_boot_dict["copula"] = fill(NaN, (cop_n, draws, T))

    return sub_boot_dict
end






function undo_copula_treatment(x, estimator; slice=false)
    if typeof(estimator) <: SeriesEstimator
        return x
    else
        if slice
            return idct(idct(x, 1), 2)
        else
            return idct(x)
        end
    end
end


"""
    gauss_legendre_nodes_weights(n)

Return `(nodes, weights)` for n-point Gauss-Legendre quadrature on [-1, 1].
Uses the eigenvalue method (symmetric tridiagonal companion matrix).
"""
function gauss_legendre_nodes_weights(n::Int)
    β = [i / sqrt(4i^2 - 1) for i in 1:n-1]
    J = SymTridiagonal(zeros(n), β)
    eigen_result = eigen(J)
    nodes = eigen_result.values
    weights = 2.0 .* eigen_result.vectors[1, :] .^ 2
    return nodes, weights
end

function estimate_confidence_intervals!(data, objects, series, years, time_dict, freq_type, estimator, measures, draws, gdp_series, source; noise_draws=nothing, rng=nothing, sd_scale::Real=1.0, sd_by_period::Bool=true)
    """The problem right now is I wish to generate data, but with a higher granularity for the pcfs than the copula, I basically have to generate the data twice.
    I have generate the pcfs for the higher granularity (to compare the reconstructions to) and the pcfs of the lower granularity, so when i export the data, 
    the copula weights correspond to the appropriate pcfs ... For now, i estimate both at the same granularity and comback to this later."""

    init_path = BASE_PATH

    # Find the choice of grid points for the confidence intervals
    @unpack grid_pcf, grid_cop = estimator

    # Identify years of estimation and create a dictionary with the data for each year
    year_data_dict = Dict()
    u_years = unique(years)

    for yr in u_years
        year_data_dict[yr] = filter(row -> row.year == yr, data)

        if "strata" in names(year_data_dict[yr])
            sort!(year_data_dict[yr], :income) # TODO: only works for CPS data, data generally may not have income 

            # Define the cluster, which are each block of 5 observations 
            year_data_dict[yr][!, "id_temp"] = 1:nrow(year_data_dict[yr])
            year_data_dict[yr][!, "cluster"] = ceil.(Int, year_data_dict[yr][!, "id_temp"] ./ 5)
        end
    end


    # Sampling from the distribution object. For each draw, compute the confidence intervals for quantiles, levels and shares 
    dim = length(measures)

    # Define the size of the sub_boot_dict
    sub_boot_dict = specify_bootstrap_container(measures, objects, estimator, years, series, "other", draws)

    # Sizes of different topologies + immutable
    cop_part = grid_cop^dim  # for the series estimator, this size here does NOT correspond to the order of the polynomial, but the initial granularity of the object that was approximated
    imm_part = grid_cop + (dim - 1) * (grid_cop - 1)
    pcf_part = grid_pcf * dim

    # Container for the objects to estimate the measurement and process noise
    n = cop_part + pcf_part - imm_part
    coef_boot = zeros(n, length(years), draws)
    fill!(coef_boot, NaN)

    # Define rows 
    cop_rows = 1:(cop_part-imm_part)
    cop_size = tuple([grid_cop for _ in 1:dim]...)
    cop_ci = CartesianIndices(cop_size)
    pcf_rows = [I for I in Iterators.partition(length(cop_rows)+1:length(cop_rows)+pcf_part, grid_pcf)]

    # Precompute pcf integration grid for series estimators (used in series_estimator)
    @unpack integral_pcf_grid = estimator
    grid_points_pcf = select_grid_points(integral_pcf_grid)
    intervals = vcat([0.0] .+ 1e-6, grid_points_pcf)

    # Precompute Gauss-Legendre quadrature data for each sub-interval.
    # This replaces adaptive quadgk calls with fixed-node evaluation in the bootstrap loop.
    n_gl = 16  # 16-point GL is exact for polynomials up to degree 31; more than sufficient for rtol=1e-8 on these smooth integrands
    gl_ref_nodes, gl_ref_weights = gauss_legendre_nodes_weights(n_gl)

    # For each sub-interval [a, b], map reference nodes from [-1,1] to [a,b] and precompute basis values
    gl_nodes   = Vector{Vector{Float64}}(undef, integral_pcf_grid)
    gl_weights = Vector{Vector{Float64}}(undef, integral_pcf_grid)
    gl_basis   = Vector{Matrix{Float64}}(undef, integral_pcf_grid)  # each is n_gl × grid_pcf

    for i in 1:integral_pcf_grid
        a, b = intervals[i], intervals[i+1]
        half_len = (b - a) / 2
        mid = (a + b) / 2

        # Map nodes to [a, b] and scale weights by Jacobian
        gl_nodes[i]   = mid .+ half_len .* gl_ref_nodes
        gl_weights[i] = half_len .* gl_ref_weights

        # Precompute Q_m(j, node) basis matrix: row = node, col = basis function j+1
        basis = Matrix{Float64}(undef, n_gl, grid_pcf)
        for k in 1:n_gl
            for j in 0:(grid_pcf - 1)
                basis[k, j+1] = Q_m(j, gl_nodes[i][k])
            end
        end
        gl_basis[i] = basis
    end

    # HANK special case:
    # If the sole goal is to construct DCT_boot for sigma estimation, we can avoid
    # bootstrapping the simulated microdata. Instead:
    # 1) compute the base coefficient path once from the simulated microdata;
    # 2) (optional) estimate per-coefficient SDs from empirical noise draws;
    # 3) generate draws by perturbing the base with i.i.d. N(0, sd^2) shocks.
    #
    # NOTE: Empirical (data-based) noise can be disabled by setting `sd_scale = 0.0`.
    # if occursin("HANK", source) && (noise_draws !== nothing || iszero(sd_scale))

    #     # NaN-robust SD for a vector.
    #     nanstd_vec(v) = begin
    #         s = 0.0
    #         s2 = 0.0
    #         nobs = 0
    #         @inbounds for x in v
    #             if !isnan(x)
    #                 nobs += 1
    #                 s += x
    #                 s2 += x * x
    #             end
    #         end
    #         nobs <= 1 && return NaN
    #         μ = s / nobs
    #         var = max(0.0, (s2 - nobs * μ * μ) / (nobs - 1))
    #         return sqrt(var)
    #     end

    #     # Fill coefficient vector for a given (bootstrapped) sample.
    #     function fill_coefs_for_sample!(out, sample, yr, actual_period)
    #         fill!(out, NaN)

    #         sample = coalesce.(sample, NaN)
    #         filter!("weight" => !isnan, sample)

    #         obs_measures = String[]

    #         # PCFs
    #         for (m, meas) in enumerate(measures)
    #             if meas ∉ names(sample)
    #                 continue
    #             end

    #             non_missing = filter(meas => !isnan, sample)
    #             nrow(non_missing) == 0 && continue

    #             push!(obs_measures, meas)

    #             avg_aggr = filter(row -> row.date >= QuarterlyDate(yr, actual_period), gdp_series)[!, meas*"_per_hh"][1]
    #             avg_data = mean(non_missing[:, meas], weights(non_missing[:, :weight]))
    #             multiplier = abs.(avg_aggr / avg_data)
    #             multiplier[1] > 20 && continue

    #             non_missing[!, meas] = non_missing[!, meas] .* multiplier
    #             tot_scale = avg_aggr .- mean(non_missing[:, meas], weights(non_missing[:, :weight]))
    #             non_missing[!, meas] .= non_missing[!, meas] .+ tot_scale
    #             sort!(non_missing, meas)

    #             t_rv = inverse_hyperbolic_sine(non_missing[:, meas] ./ avg_aggr)
    #             out[pcf_rows[m]] .= series_estimator(t_rv, non_missing[:, :weight], grid_pcf - 1)
    #         end

    #         # Copula coefficients
    #         unique!(obs_measures)
    #         obs_dims = length(obs_measures)
    #         if obs_dims >= 2
    #             try
    #                 cop_weights = get_copulas(sample, measures, obs_measures, estimator; with_immutable=true)
    #                 cop_weights = reshape(cop_weights, cop_size)
    #                 out[cop_rows] .= remove_immutable(cop_weights)
    #             catch
    #                 out[cop_rows] .= NaN
    #             end
    #         end

    #         return out
    #     end

    #     # Base coefficient path (n×T)
    #     coef_base = Matrix{Float64}(undef, n, length(years))
    #     fill!(coef_base, NaN)

    #     # Optional MC-noise SDs estimated by bootstrapping the simulated microdata within each period.
    #     mc_draws = draws
    #     sds_mc = mc_draws > 0 ? Matrix{Float64}(undef, n, length(years)) : nothing
    #     mc_draws > 0 && fill!(sds_mc, NaN)

    #     tmp_coef = Vector{Float64}(undef, n)
    #     mc_coef_draws = mc_draws > 0 ? Matrix{Float64}(undef, n, mc_draws) : nothing

    #     count_base = 1
    #     for yr in u_years
    #         for actual_period in time_dict[yr]
    #             period_data = filter(row -> row[freq_type] == actual_period, year_data_dict[yr])
    #             fill_coefs_for_sample!(@view(coef_base[:, count_base]), period_data, yr, actual_period)

    #             if mc_draws > 0
    #                 bN = nrow(period_data)
    #                 b_size = bN

    #                 for r in 1:mc_draws
    #                     boot_idx = rng === nothing ? rand(1:bN, b_size) : rand(rng, 1:bN, b_size)
    #                     b_sample = period_data[boot_idx, :]
    #                     fill_coefs_for_sample!(tmp_coef, b_sample, yr, actual_period)
    #                     @views mc_coef_draws[:, r] .= tmp_coef
    #                 end

    #                 @inbounds for k in 1:n
    #                     σ = nanstd_vec(@view mc_coef_draws[k, :])
    #                     sds_mc[k, count_base] = (isfinite(σ) && σ >= 0) ? σ : NaN
    #                 end
    #             end

    #             count_base += 1
    #         end
    #     end

    #     # Empirical noise (from data) is optional. If `sd_scale == 0` (or `noise_draws`
    #     # are not provided), we do NOT add any data-based noise.
    #     Tcoef = length(years)
    #     sds_data_use = zeros(Float64, n, Tcoef)

    #     if noise_draws !== nothing && !iszero(sd_scale)
    #         @assert ndims(noise_draws) == 3 "noise_draws must be a K×T×D array"
    #         @assert size(noise_draws, 1) == n "noise_draws has incompatible K dimension"

    #         # Compute a per-coefficient SD from empirical noise draws.
    #         #
    #         # Important: depending on how `noise_draws` was saved, it may contain either
    #         #   (a) residual draws (already ≈ mean-zero), or
    #         #   (b) coefficient draws whose mean varies over time.
    #         #
    #         # To avoid accidentally treating time-variation in the mean as “noise”, we
    #         # estimate dispersion *within each time period* (across draws).
    #         Tn = size(noise_draws, 2)
    #         Dn = size(noise_draws, 3)

    #         sd_by_period || error("sd_by_period=false is no longer supported; use period-by-period SDs")
    #         @assert Tcoef <= Tn "noise_draws has fewer time periods than coef_base"

    #         # Compute SDs by period; when the empirical object is not observed (e.g. PSID
    #         # only has all measures late in the sample), SDs are left as NaN and filled.
    #         sds = Matrix{Float64}(undef, n, Tn)
    #         fill!(sds, NaN)

    #         @inbounds for k in 1:n
    #             for t in 1:Tn
    #                 s = 0.0
    #                 nobs = 0
    #                 for d in 1:Dn
    #                     x = noise_draws[k, t, d]
    #                     if !isnan(x)
    #                         nobs += 1
    #                         s += x
    #                     end
    #                 end
    #                 nobs <= 1 && continue
    #                 μ = s / nobs

    #                 sse_t = 0.0
    #                 for d in 1:Dn
    #                     x = noise_draws[k, t, d]
    #                     if !isnan(x)
    #                         dx = x - μ
    #                         sse_t += dx * dx
    #                     end
    #                 end

    #                 var = max(0.0, sse_t / (nobs - 1))
    #                 σ = sqrt(var)
    #                 sds[k, t] = (isfinite(σ) && σ >= 0) ? σ : NaN
    #             end
    #         end

    #         # Carry/back-fill missing SDs across time per coefficient.
    #         n_filled = 0
    #         n_all_missing = 0
    #         @inbounds for k in 1:n
    #             row = @view sds[k, :]
    #             first_idx = nothing
    #             for t in 1:Tn
    #                 if !isnan(row[t])
    #                     first_idx = t
    #                     break
    #                 end
    #             end

    #             if first_idx === nothing
    #                 n_all_missing += 1
    #                 continue
    #             end

    #             first_val = row[first_idx]
    #             for t in 1:(first_idx-1)
    #                 if isnan(row[t])
    #                     row[t] = first_val
    #                     n_filled += 1
    #                 end
    #             end
    #             for t in (first_idx+1):Tn
    #                 if isnan(row[t])
    #                     row[t] = row[t-1]
    #                     n_filled += 1
    #                 end
    #             end
    #         end

    #         n_filled > 0 && @info("Filled $n_filled missing SD entries by carry/back-fill (HANK sigma fast-path)")
    #         n_all_missing > 0 && @info("$n_all_missing coefficients had no SD information at any time (HANK sigma fast-path)")

    #         # Any SDs still missing after fill imply no information anywhere; treat as zero noise.
    #         replace!(sds, NaN => 0.0)

    #         @views sds_data_use .= sds[:, 1:Tcoef]
    #     end

    #     total_sds = sd_scale .* sds_data_use

    #     if mc_draws > 0
    #         # Fill missing MC SDs across time per coefficient (option 2)
    #         n_filled_mc = 0
    #         n_all_missing_mc = 0
    #         @inbounds for k in 1:n
    #             row = @view sds_mc[k, :]
    #             first_idx = nothing
    #             for t in 1:Tcoef
    #                 if !isnan(row[t])
    #                     first_idx = t
    #                     break
    #                 end
    #             end
    #             if first_idx === nothing
    #                 n_all_missing_mc += 1
    #                 continue
    #             end
    #             first_val = row[first_idx]
    #             for t in 1:(first_idx-1)
    #                 if isnan(row[t])
    #                     row[t] = first_val
    #                     n_filled_mc += 1
    #                 end
    #             end
    #             for t in (first_idx+1):Tcoef
    #                 if isnan(row[t])
    #                     row[t] = row[t-1]
    #                     n_filled_mc += 1
    #                 end
    #             end
    #         end
    #         n_filled_mc > 0 && @info("Filled $n_filled_mc missing MC SD entries by carry/back-fill (HANK sigma fast-path)")
    #         n_all_missing_mc > 0 && @info("$n_all_missing_mc coefficients had no MC SD information at any time (HANK sigma fast-path)")

    #         replace!(sds_mc, NaN => 0.0)
    #         @views sds_mc_use = sds_mc[:, 1:Tcoef]
    #         total_sds = sqrt.((sd_scale .* sds_data_use) .^ 2 .+ (sds_mc_use) .^ 2)
    #     end

    #     # Generate coefficient draws around the simulated base path.
    #     for s in 1:draws
    #         Z = rng === nothing ? randn(n, Tcoef) : randn(rng, n, Tcoef)
    #         @views coef_boot[:, :, s] .= coef_base .+ total_sds .* Z
    #     end

    #     clean_sub_boot_dict!(sub_boot_dict)
    #     return sub_boot_dict, coef_boot
    # end

    # Other 
    qvec = zeros(grid_pcf)
    count = 1

    # For the copulas later 
    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid, integral_cop_grid = estimator
    end
    grid_size_data_cop = typeof(estimator) <: SeriesEstimator ? integral_cop_grid : grid_cop
    int_points = select_grid_points(grid_size_data_cop)
    int_points[end] = int_points[end] - 1e-6  # to avoid numerical issues with the last point
    id_x = eachindex(int_points)
    integrals_pre = precompute_integrals(grid_cop - 1, int_points)

    # Precompute avg_aggr values to avoid repeated DataFrame filtering inside the bootstrap loop
    avg_aggr_cache = Dict{Any, Float64}()
    for yr in u_years
        for actual_period in time_dict[yr]
            for meas in measures
                col = meas * "_per_hh"
                if col in names(gdp_series)
                    avg_aggr_cache[(yr, actual_period, meas)] = filter(row -> row.date >= QuarterlyDate(yr, actual_period), gdp_series)[!, col][1]
                end
            end
        end
    end

    for y in eachindex(u_years)
        yr = u_years[y]
        println("$y of $(length(u_years))")
        local rep_w
        if occursin("SCF", source)
            rep_w = CSV.read(init_path * "/1_Data/SCF+/replicate_weights/replicate_weights_$yr.csv", DataFrame)
        end

        for actual_period in time_dict[yr]
            # Import data for this period
            period_data = filter(row -> row[freq_type] == actual_period, year_data_dict[yr])
            b_size = nrow(period_data) # source == "HANK" ? 5000 :
            boot_indices = rand(1:nrow(period_data), b_size, draws)

            # Bootstrap
            Threads.@threads for s in 1:draws  # when you bootstrap, it should be the length of the sample
                local b_sample
                if "strata" in names(period_data) && !occursin("SCF", source)
                    # Generate empty df 
                    b_sample = DataFrame()
                    strata = unique(period_data[:, "strata"])

                    # Process each stratum
                    for stratum in strata
                        # Filter the DataFrame for the current stratum
                        stratum_df = period_data[period_data[:, "strata"].==stratum, :]

                        # Get unique clusters within the stratum
                        clusters = unique(stratum_df[:, "cluster"])

                        # Sample clusters with replacement
                        sampled_clusters = sample(clusters, length(clusters), replace=true)

                        # Combine the sampled clusters into the bootstrap sample
                        for cluster in sampled_clusters
                            append!(b_sample, stratum_df[stratum_df[:, "cluster"].==cluster, :])
                        end
                    end

                    # elseif occursin("SCF", source)
                    # if yr < 1983
                    # b_sample = period_data[boot_indices[:, s], :]

                    # else
                    #     # Generalize to be more about the data consisting of several imputations
                    #     b_sample = deepcopy(period_data)
                    #     b_sample[!, "newid"] = [parse(Int, string(row[:id])[1:end-1]) for row in eachrow(b_sample)]
                    #     rep_draw = select(filter(row -> !ismissing(row["wgtI95W95_imp1_$s"]), rep_w), ["id", "wgtI95W95_imp1_$s"])
                    #     b_sample = leftjoin(b_sample, rep_draw, on=:newid => :id, indicator=:source) # fully non-parametric bootstrap

                    #     # filter!(row -> row.impnum == 1, b_sample)
                    #     select!(b_sample, Not([:weight, :source]))
                    #     rename!(b_sample, Symbol("wgtI95W95_imp1_$s") => :weight)
                    # end
                else
                    b_sample = period_data[boot_indices[:, s], :]
                end

                # Some quick edits 
                b_sample = coalesce.(b_sample, NaN)
                filter!("weight" => !isnan, b_sample)

                obs_measures = []

                for (m, meas) in enumerate(measures)
                    # First check that the measure is actually in the dataset                     
                    if meas ∉ names(b_sample)
                        nothing
                        # sub_boot_dict[meas]["quantiles"][:, s, count] .= NaN                        
                        # coef_boot[pcf_rows[m], count, s]               .= NaN

                    else
                        # Check that the measure exists for this period 
                        non_missing = filter(meas => !isnan, b_sample)

                        if nrow(non_missing) == 0
                            nothing
                            # println("The measure $rv is full of NaNs")
                            # sub_boot_dict[meas]["quantiles"][:, s, count] .= NaN                        
                            # coef_boot[pcf_rows[m], count, s]               .= NaN

                        else
                            push!(obs_measures, meas)

                            # Find multiplier s.t. average in FRED is equal to average in the sample
                            avg_aggr = avg_aggr_cache[(yr, actual_period, meas)]

                            if occursin("SCF", source)
                                # Over imputations
                                if yr < 1983
                                    imp_boot = zeros(grid_pcf, 5)

                                    for ip in 1:5
                                        # Get the data
                                        non_missingᵢ = filter(x -> x.impnum == ip, non_missing)

                                        # See if the data is severely off from the aggregate
                                        avg_data = mean(non_missingᵢ[:, meas], weights(non_missingᵢ[:, :weight]))
                                        multiplier = abs.(avg_aggr / avg_data)

                                        if multiplier[1] > 20
                                            nothing
                                        end

                                        non_missingᵢ[!, meas] = non_missingᵢ[!, meas] .* multiplier
                                        tot_scale = avg_aggr .- mean(non_missingᵢ[:, meas], weights(non_missingᵢ[:, :weight])) # should be zero if data is all positive 
                                        non_missingᵢ[!, meas] .= non_missingᵢ[!, meas] .+ tot_scale # the average is corrected and ranks don't change.

                                        # sort 
                                        sort!(non_missingᵢ, meas)

                                        # First, transform series 
                                        t_rv = inverse_hyperbolic_sine(non_missingᵢ[:, meas] ./ avg_aggr)

                                        # Estimate weights 
                                        imp_boot[:, ip] .= series_estimator(t_rv, non_missingᵢ[:, :weight], grid_pcf - 1)
                                    end

                                    # Average over coefficients -> generate quantiles
                                    coef_boot[pcf_rows[m], count, s] .= mean(imp_boot, dims=2)

                                    coefs_i = @view coef_boot[pcf_rows[m], count, s]
                                    for i in 1:integral_pcf_grid
                                        integral_val = 0.0
                                        @inbounds for k in 1:n_gl
                                            qf_val = 0.0
                                            for j in 1:grid_pcf
                                                qf_val += gl_basis[i][k, j] * coefs_i[j]
                                            end
                                            integral_val += gl_weights[i][k] * reverse_inverse_hyperbolic_sine(qf_val)
                                        end
                                        sub_boot_dict[meas]["quantiles"][i, s, count] = integral_val * avg_aggr / (intervals[i+1] - intervals[i])
                                    end
                                else
                                    avg_data = mean(non_missing[:, meas], weights(non_missing[:, :weight]))
                                    multiplier = abs.(avg_aggr / avg_data)

                                    if multiplier[1] > 20
                                        nothing
                                    end

                                    non_missing[!, meas] = non_missing[!, meas] .* multiplier
                                    tot_scale = avg_aggr .- mean(non_missing[:, meas], weights(non_missing[:, :weight])) # should be zero if data is all positive 
                                    non_missing[!, meas] .= non_missing[!, meas] .+ tot_scale # the average is corrected and ranks don't change.     

                                    # sort 
                                    sort!(non_missing, meas)

                                    # First, transform series 
                                    t_rv = inverse_hyperbolic_sine(non_missing[:, meas] ./ avg_aggr)

                                    # Estimate weights 
                                    coef_boot[pcf_rows[m], count, s] .= series_estimator(t_rv, non_missing[:, :weight], grid_pcf - 1)


                                    # Generates the integral over the percentile function
                                    coefs_i = @view coef_boot[pcf_rows[m], count, s]
                                    for i in 1:integral_pcf_grid
                                        integral_val = 0.0
                                        @inbounds for k in 1:n_gl
                                            qf_val = 0.0
                                            for j in 1:grid_pcf
                                                qf_val += gl_basis[i][k, j] * coefs_i[j]
                                            end
                                            integral_val += gl_weights[i][k] * reverse_inverse_hyperbolic_sine(qf_val)
                                        end
                                        sub_boot_dict[meas]["quantiles"][i, s, count] = integral_val * avg_aggr / (intervals[i+1] - intervals[i])
                                    end
                                end
                            else
                                avg_data = mean(non_missing[:, meas], weights(non_missing[:, :weight]))

                                multiplier = abs.(avg_aggr / avg_data)

                                if multiplier[1] > 20
                                    nothing
                                end

                                non_missing[!, meas] = non_missing[!, meas] .* multiplier
                                tot_scale = avg_aggr .- mean(non_missing[:, meas], weights(non_missing[:, :weight])) # should be zero if data is all positive 
                                non_missing[!, meas] .= non_missing[!, meas] .+ tot_scale # the average is corrected and ranks don't change.     

                                # sort 
                                sort!(non_missing, meas)

                                # First, transform series 
                                t_rv = inverse_hyperbolic_sine(non_missing[:, meas] ./ avg_aggr)

                                # Estimate weights 
                                coef_boot[pcf_rows[m], count, s] .= series_estimator(t_rv, non_missing[:, :weight], grid_pcf - 1)


                                # Generates the integral over the percentile function
                                coefs_i = @view coef_boot[pcf_rows[m], count, s]
                                for i in 1:integral_pcf_grid
                                    try
                                        integral_val = 0.0
                                        @inbounds for k in 1:n_gl
                                            qf_val = 0.0
                                            for j in 1:grid_pcf
                                                qf_val += gl_basis[i][k, j] * coefs_i[j]
                                            end
                                            integral_val += gl_weights[i][k] * reverse_inverse_hyperbolic_sine(qf_val)
                                        end
                                        sub_boot_dict[meas]["quantiles"][i, s, count] = integral_val * avg_aggr / (intervals[i+1] - intervals[i])
                                    catch e
                                        sub_boot_dict[meas]["quantiles"][i, s, count] = NaN
                                    end
                                end
                            end
                        end
                    end
                end

                # Generate copulas
                obs_dims = length(obs_measures)
                # XX no longer needed — copula_cdf_estimator now uses array indices directly
                # XX = obs_dims == 2 ? [[int_points[i], int_points[j]] for i in id_x, j in id_x] : obs_dims == 3 ? [[int_points[i], int_points[j], int_points[k]] for i in id_x, j in id_x, k in id_x] : [1]
                XX = nothing  # kept for generate_copula_density signature compatibility

                if dim == obs_dims
                    if occursin("SCF", source)
                        if yr < 1983
                            temp_cop = zeros(grid_cop^dim, 5)

                            for i in 1:5
                                period_dataᵢ = filter(row -> row.impnum == i, b_sample)
                                temp_cop[:, i] = get_copulas(period_dataᵢ, measures, obs_measures, estimator; with_immutable=true)
                            end

                            # Compute average over temp_cop 
                            cop_weights = mean(temp_cop, dims=2)
                            cop_weights = reshape(cop_weights, cop_size)

                            sub_boot_dict["copula"][:, s, count] .= generate_copula_density(cop_weights, XX, id_x, integrals_pre, measures, obs_measures, grid_size_data_cop)[:]
                            coef_boot[cop_rows, count, s] .= remove_immutable(cop_weights)

                        else
                            try
                                cop_weights = get_copulas(b_sample, measures, obs_measures, estimator; with_immutable=true)
                                cop_weights = reshape(cop_weights, cop_size)

                                sub_boot_dict["copula"][:, s, count] .= generate_copula_density(cop_weights, XX, id_x, integrals_pre, measures, obs_measures, grid_size_data_cop)[:]
                                coef_boot[cop_rows, count, s] .= remove_immutable(cop_weights)

                            catch e
                                sub_boot_dict["copula"][:, s, count] .= NaN
                                coef_boot[cop_rows, count, s] .= NaN
                            end

                        end
                    else
                        try
                            cop_weights = get_copulas(b_sample, measures, obs_measures, estimator; with_immutable=true)
                            cop_weights = reshape(cop_weights, cop_size)
                            sub_boot_dict["copula"][:, s, count] .= generate_copula_density(cop_weights, XX, id_x, integrals_pre, measures, obs_measures, grid_size_data_cop)[:]
                            coef_boot[cop_rows, count, s] .= remove_immutable(cop_weights)
                        catch e
                            sub_boot_dict["copula"][:, s, count] .= NaN
                            coef_boot[cop_rows, count, s] .= NaN
                        end
                    end
                elseif obs_dims < dim && obs_dims >= 2
                    if typeof(estimator) <: SeriesEstimator # TODO: use multiple dispatch 
                        if occursin("SCF", source)
                            if yr < 1983
                                temp_cop = zeros(grid_cop^dim, 5)

                                for i in 1:5
                                    period_dataᵢ = filter(row -> row.impnum == i, b_sample)
                                    temp_cop[:, i] = get_copulas(period_dataᵢ, measures, obs_measures, estimator; with_immutable=true)
                                end

                                # Compute average over temp_cop 
                                cop_weights = mean(temp_cop, dims=2)
                                cop_weights = reshape(cop_weights, cop_size)

                                sub_boot_dict["copula"][:, s, count] .= generate_copula_density(cop_weights, XX, id_x, integrals_pre, measures, obs_measures, grid_size_data_cop)[:]
                                coef_boot[cop_rows, count, s] .= remove_immutable(cop_weights)

                            else
                                cop_weights = get_copulas(b_sample, measures, obs_measures, estimator; with_immutable=true)
                                cop_weights = reshape(cop_weights, cop_size)

                                sub_boot_dict["copula"][:, s, count] .= generate_copula_density(cop_weights, XX, id_x, integrals_pre, measures, obs_measures, grid_size_data_cop)[:]

                                coef_boot[cop_rows, count, s] .= remove_immutable(cop_weights)
                            end
                        else
                            try
                                cop_weights = get_copulas(b_sample, measures, obs_measures, estimator; with_immutable=true)
                                cop_weights = reshape(cop_weights, cop_size)
                                sub_boot_dict["copula"][:, s, count] .= generate_copula_density(cop_weights, XX, id_x, integrals_pre, measures, obs_measures, grid_size_data_cop)[:]
                                coef_boot[cop_rows, count, s] .= remove_immutable(cop_weights)
                            catch e
                                sub_boot_dict["copula"][:, s, count] .= NaN
                                coef_boot[cop_rows, count, s] .= NaN
                            end
                        end

                    else
                        correction = sqrt(grid_cop)^(dim - obs_dims)
                        cop_with_imm = get_copulas(b_sample, measures, obs_measures, estimator; with_immutable=true) .* correction
                        cop_mat = reshape(cop_with_imm, cop_size)

                        coef_boot[cop_rows, count, s] .= remove_immutable(cop_mat ./ correction)

                        # Get indices to undo treatment 
                        cop_ind = copula_case(obs_measures, measures)
                        vec_of_indices = cop_ci[cop_ind...][:]
                        id_of_cop = findall(x -> x in vec_of_indices, cop_ci[:])

                        # In this case, the idct will return NaNs, so, we need to perform the idct like so:
                        sub_boot_dict["copula"][id_of_cop, s, count] = undo_copula_treatment(cop_mat[cop_ind...], estimator)[:]  #TODO: this only works if the copula has NaNs on the other slices -> otherwise, it would be incorrect 
                    end

                elseif obs_dims <= 1
                    # coef_boot[cop_rows, count, s] .= NaN
                    nothing
                end
            end
            count += 1
        end
    end

    # Remove surplus columns (i.e., draws) which have only zeros 
    clean_sub_boot_dict!(sub_boot_dict)

    return sub_boot_dict, coef_boot
end


function eval_quantile_function(coefficients, order, u)
    val = 0.0
    @inbounds for j in 0:order
        val += coefficients[j+1] * Q_m(j, u)
    end
    return val
end


function clean_sub_boot_dict!(sub_boot_dict) # TODO: do i even need this
    measures = setdiff(collect(keys(sub_boot_dict)), ["copula"])

    for meas in measures
        for o in collect(keys(sub_boot_dict[meas]))
            # for c in axes(sub_boot_dict[meas][o], 3)
            # println(c)
            # println(sub_boot_dict[meas][o])
            sub_boot_dict[meas][o] = sub_boot_dict[meas][o][:, mapslices(col -> all(col .!== 0), sub_boot_dict[meas][o][:, :, 1], dims=1)[:], :] # c = date, and if it's missing in the first draw ,its missing in all of them 
            # sub_boot_dict[meas][o] = sub_boot_dict[meas][o][:, vec(mapslices(col -> all(col .!==NaN), sub_boot_dict[meas][o][:, :, 1], dims = 1)),:] # c = date, and if it's missing in the first draw ,its missing in all of them 
            # end
        end
    end
end

function any_vectors_nan(vectors)
    return any(all(isnan, vector) for vector in vectors)
end


function cdf_to_pdf(cdf_matrix::Array{Float64,2})
    # Initialize the PDF matrix with zeros
    pdf_matrix = zeros(size(cdf_matrix))

    # Compute differences between adjacent cells
    for i in 1:size(cdf_matrix, 1)
        for j in 1:size(cdf_matrix, 2)
            left = j > 1 ? cdf_matrix[i, j-1] : 0
            above = i > 1 ? cdf_matrix[i-1, j] : 0
            left_above = (i > 1 && j > 1) ? cdf_matrix[i-1, j-1] : 0
            pdf_matrix[i, j] = cdf_matrix[i, j] - left - above + left_above
        end
    end

    # correct pdf_matrix 
    pdf_matrix[pdf_matrix.<0] .= 0
    pdf_matrix .= pdf_matrix ./ sum(pdf_matrix)

    return pdf_matrix
end

function cdf_to_pdf(cdf_matrix::Array{Float64,3})
    # Initialize the PDF matrix with zeros
    pdf_matrix = zeros(size(cdf_matrix))

    # Compute differences between adjacent cells in 3D
    for i in 1:size(cdf_matrix, 1)
        for j in 1:size(cdf_matrix, 2)
            for k in 1:size(cdf_matrix, 3)
                left = k > 1 ? cdf_matrix[i, j, k-1] : 0
                above = j > 1 ? cdf_matrix[i, j-1, k] : 0
                back = i > 1 ? cdf_matrix[i-1, j, k] : 0
                left_above = (j > 1 && k > 1) ? cdf_matrix[i, j-1, k-1] : 0
                left_back = (i > 1 && k > 1) ? cdf_matrix[i-1, j, k-1] : 0
                above_back = (i > 1 && j > 1) ? cdf_matrix[i-1, j-1, k] : 0
                left_above_back = (i > 1 && j > 1 && k > 1) ? cdf_matrix[i-1, j-1, k-1] : 0

                pdf_matrix[i, j, k] = cdf_matrix[i, j, k] - left - above - back + left_above + left_back + above_back - left_above_back
            end
        end
    end

    # correct pdf_matrix 
    pdf_matrix[pdf_matrix.<0] .= 0
    pdf_matrix .= pdf_matrix ./ sum(pdf_matrix)

    return pdf_matrix
end


function undo_logit_transformation(x)
    return exp.(x) ./ (1 .+ exp.(x))
end


function perform_proof_of_concept_Γ_comparison(MV, df_vec, gdp_series, time_p, freq_type, time_dict, model_options, means, stds, trend)

    @unpack lags, measures, freq, plot_proof, pca_perspective, estimator, tag = model_options
    @unpack grid_pcf, grid_cop = estimator
    @unpack tmin, tmax, tot_periods, year_vec = time_p

    dimension = length(measures)

    # Step 1: Run PCA with PSID data
    PSID_id = findall(x -> x == "PSID", df_vec.df_names)[1]
    PSID = MV[PSID_id]

    time_dict_k = time_dict[PSID_id]

    # Perform PCA
    proj, pcs, M, _ = perform_pca(PSID, measures, :functional_data, tag)

    # Step 2: Run PCA with overlayed data
    pcs2, proj2 = run_TW_algorithm(MV, 0.8, "second tallest")

    # Step 3: Reconstruct data 
    X̃ = proj * pcs
    X̃2 = proj2 * pcs2

    # Keep the PSID data for X̃2
    X̃2 = X̃2[:, 1:size(X̃, 2)]

    # Step 4: Start reconstructing 
    # Add mean and std back 
    ΓF_σ = add_variance(estimator, X̃, stds, measures)
    ΓF_σ2 = add_variance(estimator, X̃2, stds, measures)
    X1 = ΓF_σ .+ means[PSID_id]
    X2 = ΓF_σ2 .+ means[PSID_id]

    # Find IDs related to PSID 
    # estimation_id = generate_correct_indices(year_vec[PSID_id], freq_type[PSID_id], freq, tmin, time_dict_k)  # for adding back trend to the correct index
    estimation_id_full = []
    for i in axes(MV[PSID_id], 2)
        if sum(isnan.(MV[PSID_id][:, i])) == 0
            push!(estimation_id_full, i)
        end
    end


    for (i, t) in enumerate(estimation_id_full)
        X1[:, i] .= X1[:, i] .+ trend[PSID_id][:, t]
        X2[:, i] .= X2[:, i] .+ trend[PSID_id][:, t]
    end

    d_container = add_multidimensional_immutable(estimator, X1, grid_cop, measures)
    r_container_all = add_multidimensional_immutable(estimator, X2, grid_cop, measures)

    filtered_dict = Dict(k => v for (k, v) in time_dict[PSID_id] if k >= 1999)
    d_data_vector = undo_functional_treatment_validation(estimator, d_container, grid_pcf, grid_cop, measures, gdp_series, filtered_dict, "quarter")
    r_data_vector = undo_functional_treatment_validation(estimator, r_container_all, grid_pcf, grid_cop, measures, gdp_series, filtered_dict, "quarter")

    # subsetting gdp series 
    gdp_series[!, "date"] = QuarterlyDate.(gdp_series[!, "time"])
    filter!(row -> row.date >= QuarterlyDate(tmin["year"], tmin["quarter"]), gdp_series)
    filter!(row -> row.date <= QuarterlyDate(tmax["year"], tmax["quarter"]), gdp_series)


    # From weights, generate data 
    if typeof(estimator) <: SeriesEstimator
        for dv in [d_data_vector, r_data_vector]
            @unpack integral_pcf_grid, integral_cop_grid = estimator

            grid_size_data_pcf = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf
            grid_size_data_cop = typeof(estimator) <: SeriesEstimator ? integral_cop_grid : grid_cop

            # Generate container to store the data of choice 
            T = size(dv[2], 2)
            new_data_pcf = [zeros(grid_size_data_pcf, T) for _ in 1:dimension]

            # Fill 'new_data_pcf' with NaN
            for m in eachindex(new_data_pcf)
                new_data_pcf[m] .= NaN
            end

            grid_points_pcf = select_grid_points(grid_size_data_pcf)
            intervals = vcat([0.0] .+ 1e-6, grid_points_pcf)

            # split the pcfs by measure
            split_pcfs = [dv[2][I, :] for I in Iterators.partition(axes(dv[2], 1), grid_pcf)]

            # Agg data 
            correction_vec = zeros(T, dimension)

            # Storing real gdp pc data
            col_names = sort([meas * "_per_hh" for meas in measures])

            if freq_type[PSID_id] == "year"
                # Collects real gdp pc for all years, not just observed     
                for (i, y) in enumerate(sort(collect(keys(filtered_dict))))
                    p_data = filter(row -> row.date == QuarterlyDate(y, 4), gdp_series)
                    correction_vec[i, :] = Matrix(p_data[:, col_names])[:] #float.()[1]
                end

            elseif freq_type[PSID_id] == "quarter"
                count = 1
                for y in sort(collect(keys(filtered_dict))) # 
                    for p in sort(filtered_dict[y])
                        p_data = filter(row -> row.date == QuarterlyDate(y, p), gdp_series)
                        correction_vec[count, :] = Matrix(p_data[:, col_names])[:]
                        count += 1
                    end
                end
            end

            # Integration 
            integrate_quantile_functions!(new_data_pcf, split_pcfs, grid_pcf, intervals, correction_vec)
            dv[2] = vcat([new_data_pcf[m] for m in eachindex(new_data_pcf)]...)

            # Integrate densities 
            dv[1] = generate_copula_densities(deepcopy(dv[1]), measures, grid_size_data_cop) #FIXME: 
        end
    end

    # Generate shares and levels
    d_levels, d_shares = generate_shares_levels(d_data_vector[2], model_options, gdp_series)
    r_levels, r_shares = generate_shares_levels(r_data_vector[2], model_options, gdp_series)

    # Create dictionary of all data 
    d_data_dict = create_time_series_dictionary([d_data_vector..., d_levels, d_shares], estimator, measures)
    r_data_dict = create_time_series_dictionary([r_data_vector..., r_levels, r_shares], estimator, measures)

    # Generate proof of concept figures for percentile functions 
    gen_proof_of_concept_figure_Γ_comparison(d_data_dict, r_data_dict, model_options, time_p, estimation_id_full)
end


function gen_proof_of_concept_figure_Γ_comparison(d_data_dict, r_data_dict, model_options, time_p, estimation_id)
    @unpack measures, equivalized, bottom_coded, blind_to, estimator = model_options
    @unpack tmin, tmax = time_p

    if typeof(estimator) <: SeriesEstimator
        @unpack integral_pcf_grid = estimator
    else
        @unpack grid_pcf, grid_cop = estimator
    end
    grid_choice = typeof(estimator) <: SeriesEstimator ? integral_pcf_grid : grid_pcf

    # Plots reflect same scale (ex-post)
    @unpack time_dict, tmin, tmax = time_p

    # Get some observations 
    random_measure = first(measures)
    all_periods = axes(d_data_dict[random_measure]["quantiles"]["data"], 2)

    xaxis = collect(1:1:length(all_periods))
    objects = sort(["quantiles"])
    dts = QuarterlyDate(tmin["year"], tmin["quarter"]):Quarter(1):QuarterlyDate(tmax["year"], tmax["quarter"])
    dts_to_keep = dts[estimation_id]
    fully_obs_id = estimation_id #.- estimation_id[1] .+ 1

    # Define path to save figures
    init_path = BASE_PATH
    m_label = measures_folder(measures)
    path = init_path * "/7_Results/proof_of_concept/$m_label/"

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
    interval_dts = measures == ["consum", "income"] ? 3 : 2
    for (n, m) in enumerate(measures)
        M = uppercasefirst(m)
        for (c, o) in enumerate(objects)
            sequences = define_sequences(grid_choice) # sequences = [1, 2:5, 6, 7:9, 10]
            color_choices = [:blue, :green, :sienna, :black, :red]
            file_tags = ["bottom", "middle", "top"]
            circlemaker(x, y, r) = Plots.Shape(r * sind.(0:10:360) .+ x, r * cosd.(0:10:360) .+ y)

            for (s, sequence) in enumerate(sequences)
                Plots.plot()
                dataₒ = scale_transformation(o, d_data_dict[m][o]["data"][sequence, :])
                # cond   = length(sequence) == 1 ? vec(.!any(isnan.(dataₒ), dims=2)) : vec(.!any(isnan.(dataₒ), dims=1))
                cond = length(sequence) == 1 ? vec(.!any(isnan.(r_data_dict[m][o]["data"][sequence..., :]), dims=2)) : vec(.!any(isnan.(r_data_dict[m][o]["data"][sequence, :]), dims=1))

                println(r_data_dict[m][o]["data"][sequence, :])
                println(d_data_dict[m][o]["data"][sequence, :])
                println(cond)
                println(length(cond))
                println(sequence)
                println(dataₒ)

                sxaxis = xaxis[cond]
                sdata = length(sequence) == 1 ? dataₒ[cond] : dataₒ[:, cond]
                l_seq = length(sequence)

                for (i, j) in enumerate(sequence)
                    # sci_l  = scale_transformation(o, ci_l[m][o][j, :][cond])
                    # sci_u  = scale_transformation(o, ci_u[m][o][j, :][cond])

                    indices = l_seq == 1 ? (:, :) : (i, :)
                    Plots.plot!(
                        sxaxis,
                        sdata[indices...],
                        ylabel=M == "Consum" ? L"\textrm{Consumption\,\,in\,\, 10k\, USD}" : L"\textrm{%$(M)\,\,in\,\, 10k\, USD}",
                        xformatter=:latex,
                        yformatter=:latex,
                        xticks=(sxaxis[1:interval_dts:end], [L"%$(date)"[1:5] * "\$" for (_, date) in enumerate(dts_to_keep[cond][1:interval_dts:end])]),
                        label="",
                        linewidth=3,
                        ls=:dot,
                        xtickfontsize=10,
                        ytickfontsize=10,
                        legendfontsize=10,
                        guidefontsize=14,
                        lc=:gray, #color_choices[i],
                        dpi=500,
                    )
                    fs = 12
                    if M == "Consum"
                        if j <= 5
                            if j == 1
                                Plots.annotate!(sxaxis[end-6], sdata[indices...][end-6] .* 1.10, Plots.text(label_q(j), fs, :bottom))
                            else
                                Plots.annotate!(sxaxis[end-6], sdata[indices...][end-6] .* 1.06, Plots.text(label_q(j), fs, :bottom))
                            end
                        elseif j >= 6 && j <= 8
                            Plots.annotate!(sxaxis[end-6], sdata[indices...][end] .* 0.93, Plots.text(label_q(j), fs, :bottom))
                        elseif j == 9
                            Plots.annotate!(sxaxis[end-6], sdata[indices...][end] .* 0.94, Plots.text(label_q(j), fs, :bottom))
                        else
                            Plots.annotate!(sxaxis[end-5], sdata[indices...][end-3] .* 1.015, Plots.text(label_q(j), fs, :bottom))
                        end
                    elseif M == "Income"
                        if j <= 5
                            Plots.annotate!(sxaxis[end-4], sdata[indices...][end-3] * 1.06, Plots.text(label_q(j), fs, :bottom))
                        elseif j >= 6 && j <= 9
                            Plots.annotate!(sxaxis[end-4], sdata[indices...][end-3] .* 1.06, Plots.text(label_q(j), fs, :bottom))
                        else
                            Plots.annotate!(sxaxis[end-5], sdata[indices...][end-4] .* 1.06, Plots.text(label_q(j), fs, :bottom))
                        end
                    elseif M == "Wealth"
                        if j <= 5
                            if j ∈ [4, 5]
                                Plots.annotate!(sxaxis[end-6], sdata[indices...][end-5] .* 1.06, Plots.text(label_q(j), fs, :bottom))
                            elseif j == 2
                                Plots.annotate!(sxaxis[end-2], sdata[indices...][end-5] .* 0.01, Plots.text(label_q(j), fs, :bottom))
                            elseif j == 3
                                Plots.annotate!(sxaxis[end-6], sdata[indices...][end-5] .* 0.3, Plots.text(label_q(j), fs, :bottom))
                            else
                                Plots.annotate!(sxaxis[end-2], sdata[indices...][end-2] .* 0.8, Plots.text(label_q(j), fs, :bottom))
                            end
                        elseif j >= 6 && j <= 9
                            if j == 6
                                Plots.annotate!(sxaxis[end-7], sdata[indices...][end-7] .* 0.5, Plots.text(label_q(j), fs, :bottom))
                            else
                                Plots.annotate!(sxaxis[end-3], sdata[indices...][end] .* 0.8, Plots.text(label_q(j), fs, :bottom))
                            end
                        else
                            Plots.annotate!(sxaxis[end-5], sdata[indices...][end-5] .* 1.05, Plots.text(label_q(j), fs, :bottom))
                        end
                    end

                    if i != l_seq
                        Plots.scatter!(sxaxis, sdata[indices...], marker=:square, mc=:white, msc=:black, msw=1, label="", dpi=500)
                    else
                        Plots.scatter!(sxaxis, sdata[indices...], marker=:square, mc=:white, msc=:black, msw=1, label="", dpi=500)
                        Plots.scatter!([], [], marker=:square, mc=:white, msc=:black, msw=1, label=L"\textrm{\Gamma}", dpi=500)
                    end

                end
                Plots.scatter!(
                    sxaxis,
                    length(sequence) == 1 ? scale_transformation(o, r_data_dict[m][o]["data"][sequence, :])[cond] : scale_transformation(o, r_data_dict[m][o]["data"][sequence, :])[:, cond]',
                    label=["" "" "" ""],
                    dpi=500,
                    msc=:black,
                    mc=:white,
                    ms=3,
                )

                Plots.scatter!(
                    sxaxis,
                    scale_transformation(o, r_data_dict[m][o]["data"][sequence[end], :])[cond],
                    label=L"\Gamma_{\textrm{aug}}",
                    # marker=:square,
                    msc=:black,
                    mc=:white,
                    ms=3,
                    dpi=500,
                )
                Plots.pdf(path * m * "_" * o * "_$(file_tags[s])_quantiles" * ".pdf")
            end
        end
    end
end


function irf_wold(L::AbstractMatrix, h::Integer)
    n = size(L, 1)
    Ψ = Array{eltype(L)}(undef, n, n, h + 1)
    Ψ[:, :, 1] .= I(n)              # zeroth-period response
    Ψ[:, :, 2] .= L                 # 1-step
    for ℓ = 2:h
        @views Ψ[:, :, ℓ+1] .= L * Ψ[:, :, ℓ]
    end
    return Ψ
end

function fevd(L, Σv, h; method=:trace)
    n = size(L, 1)
    # orthogonalised shocks  (use your favourite identification here)
    W = cholesky(Σv).L           # Σv = W W'
    Ψ = irf_wold(L, h)           # Wold IRFs

    # structural IRFs  θ[ℓ] = Ψ[ℓ] * W
    θ = similar(Ψ)
    for ℓ = 1:h+1
        @views θ[:, :, ℓ] .= Ψ[:, :, ℓ] * W
    end

    # cumulative sum of squared IRFs (IEC notation)
    # dims = 3 → sum over lags 0…ℓ
    θ2cum = cumsum(abs2.(θ); dims=3)

    if method == :diag
        numer = permutedims(diag.(eachslice(θ2cum, dims=2)), (2, 1, 3)) # n×n×h+1
        denom = sum(numer; dims=2)                                 # n×1×h+1
        return numer ./ denom
    elseif method == :trace
        numer = sum(θ2cum; dims=(1, 2))           # 1×1×h+1  per shock
        denom = sum(numer; dims=2)               # 1×1×h+1  total
        return dropdims(numer ./ denom; dims=(1, 2))
    else
        throw(ArgumentError("method must be :diag or :trace"))
    end
end

