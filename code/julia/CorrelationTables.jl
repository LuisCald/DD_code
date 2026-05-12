# Make correlation matrix to export 

# Import files that correspond to the different external sources 
function export_all_correlations(m_label::AbstractString="income_and_wealth", tag::AbstractString="")
init_path    = BASE_PATH
base_dir     = init_path * "/7_Results/" * m_label * tag * "/from_mcmc/plots/correlations/"
mkpath(base_dir)
est_s_dict   = Dict()
data_sources = sort(["ACS", "CPS", "WID", "DFA"])

main_df = DataFrame()
for (s,ext_s) in enumerate(data_sources)
    est_s_dict[ext_s] = CSV.read(base_dir * "$ext_s" * "_correlations.csv", DataFrame)
    if s == 1
       main_df =  est_s_dict[ext_s]
    else
        append!(main_df, est_s_dict[ext_s], cols = :union)
    end
end
# Sort the rows by cycle and raw data 
sorted_column_names = sort(names(main_df), by = custom_sort)
cycle_df            = select(main_df, Symbol.(sorted_column_names[occursin.("cycle", sorted_column_names)]))
non_cycle_df        = select(main_df, Not(Symbol.(sorted_column_names[occursin.("cycle", sorted_column_names)])))

# For each dataframe, round the values to 2 decimals
for df in [cycle_df, non_cycle_df]
    for c in names(df)
        df[!, c] = round.(df[!, c], digits = 2)
    end
end

# Transpose the dataframes and create 5 columns: objects, ACS, CPS, DFA, WID
new_cycle_df     = DataFrame()
new_non_cycle_df = DataFrame()

for (df, new_df) in zip([cycle_df, non_cycle_df], [new_cycle_df, new_non_cycle_df])
    new_df[:, :objects] = LaTeXString.(join.(split.(names(df), "_"), " "))
    df = Matrix(df)
    new_df[!, :ACS]     = df[1, :]
    new_df[!, :CPS]     = df[2, :]
    new_df[!, :DFA]     = df[3, :]
    new_df[!, :WID]     = df[4, :]
end

CSV.write(base_dir * "cycle_correlations.csv", cycle_df, writeheader = true)
CSV.write(base_dir * "non_cycle_correlations.csv", non_cycle_df, writeheader = true)
end
copy_to_clipboard(true)
latexify(new_non_cycle_df, latex=false, env=:table, alignment=:l, header = ["Objects", "ACS", "CPS", "DFA", "WID"], 
    # formatters = [x -> replace(x, "_", " ")], 
    # escapefunc = x -> replace(x, "_", " "), 
    caption = "Correlations between the different data sources for the cycle variables", 
    label = "tab:cycle_correlations")


# Custom sorting function
function custom_sort(column_name)
    if occursin("cycle", column_name)
        return 0  # Sorts column names containing "cycle" first
    else
        return 1  # Sorts other column names afterwards
    end
end

