init_path                   = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
cex                     = jldopen(init_path * "/2_Data_processing/confidence_intervals/ci_consum_and_income_and_wealth_deciles_cexnpimp" * ".jld2", "r")
cex_confidence_intervals    = cex["ci"]
dist_dict = Dict("bottom" => 1:5, "middle" => 6:9, "top" => 10)

for j in ["income", "consum"]
    a = cex_confidence_intervals["CEX"]["ci_u"][j]["quantiles"]
    b = cex_confidence_intervals["CEX"]["ci_l"][j]["quantiles"]

    for (s,k) in dist_dict
        Plots.plot()
        for i in k
            Plots.plot!(1:152, a[i, :], fillrange = b[i, :], fillalpha = 0.2, label = "")
            Plots.plot!(1:152, b[i, :], fillalpha = 0.2, label = "")
        end

        Plots.savefig("ci_$(j)_$(s).pdf")
    end
end

aci_dist_old = jldopen("/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/confidence_intervals/ci_consum_and_income_and_wealth_deciles_CEX_and_CPS_and_CPS2_and_PSID_and_SCFnpimp_all2.jld2", "r")
aci_dist_old = aci_dist_old["ci"]

# Save the data from cex 
@unpack func_dict = func_struct

# Save func_dict["CEX"] in jld2 format
init_path = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
save_path = init_path * "/2_Data_processing/data_consum_and_income_and_wealth_deciles_cexnpimp" * ".jld2"
JLD2.save(save_path, "CEX", func_dict["CEX"])

# import CEX all confidence intervals 
cex = jldopen("/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/confidence_intervals/ci_consum_and_income_and_wealth_deciles_npimp.jld2", "r")
cex = cex["ci"]


# import liquid confidence intervals 
w_cex = jldopen("/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/confidence_intervals/ci_consum_and_income_and_wealth_ventiles_series_cex_and_CPS_and_CPS2_and_PSID_and_SCF_all.jld2", "r")
w_cex = w_cex["ci"]

ci_dist_old = jldopen("/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/confidence_intervals/ci_consum_and_income_and_wealth_deciles_CEX_and_CPS_and_CPS2_and_PSID_and_SCFnpimp_all.jld2", "r")
ci_dist_old = ci_dist_old["ci"]
ci_dist_old["CEX"] = w_cex["CEX"]

init_path = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
save_path = init_path * "/2_Data_processing/confidence_intervals/ci_consum_and_income_and_wealth_deciles_cex_and_CPS_and_CPS2_and_PSID_and_SCFnpimp_all.jld2"
JLD2.save(save_path, "ci", ci_dist_old)

# import PSID confidence intervals 
PSID = jldopen("/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/confidence_intervals/ci_consum_and_income_and_wealth_deciles_PSIDnpimp.jld2", "r")
PSID = PSID["ci"]

ci_dist_old = jldopen("/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/confidence_intervals/ci_consum_and_income_and_wealth_deciles_CEX_and_CPS_and_CPS2_and_PSID_and_SCFnpimp_all.jld2", "r")
ci_dist_old = ci_dist_old["ci"]
ci_dist_old["PSID"] = PSID["PSID"]

# Save the new noise distribution file
init_path = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
save_path = init_path * "/2_Data_processing/confidence_intervals/ci_consum_and_income_and_wealth_deciles_CEX_and_CPS_and_CPS2_and_PSID_and_SCFnpimp_all.jld2"
JLD2.save(save_path, "ci", ci_dist_old)

# import CEX intervals, change something and save 
CEX = jldopen("/Users/lc/Dropbox/Distributional_Dynamics/5_Code/ci_draws_consum_and_income_and_wealth_deciles_series_cex_all.jld2", "r")
CEX = CEX["ci"]
init_path = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
save_path = init_path * "/5_Code/ci_draws_consum_and_income_and_wealth_deciles_series_cex_all_test.jld2"

for m in ["income", "consum", "wealth"]
    for o in ["quantiles", "shares", "levels"]
        replace!(CEX[m][o], 0.0=>NaN)
    end
end


jldopen(save_path, "w") do file
    file["ci"] = CEX
end

#################################################################################
# import the noise distributions 
noise_dist = jldopen("/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/noise_distributions/noise_draws_consum_and_income_and_wealth_deciles_series_cex_all.jld2", "r")
noise_dist = noise_dist["noise"]

noise_dist = jldopen("/home/luisc/Distributional_Dynamics/noise_distributions/noise_draws_consum_and_income_and_wealth_deciles_series_cex_all.jld2", "r")

# Replace 0.0 with NaN 
replace!(noise_dist, 0.0=>NaN)

# Save the new noise distribution file
init_path = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
save_path = init_path * "/5_Code/noise_draws_consum_and_income_and_wealth_deciles_series_cex_all.jld2"
jldopen(save_path, "w") do file
    file["noise"] = noise_dist
end


noise_dist = jldopen("/home/luisc/Distributional_Dynamics/noise_distributions/noise_draws_consum_and_income_and_wealth_deciles_series_cex_all.jld2", "r")
noise_dist = noise_dist["noise"]
noise_dist[1]["CEX"]


# import the old noise distribution file and replace some of the entries with the new ones
noise_dist_old = jldopen("/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/noise_distributions/noise_consum_and_income_and_wealth_deciles_CEX_and_CPS_and_CPS2_and_PSID_and_SCFnpimp_all.jld2", "r")
noise_dist_old = noise_dist_old["noise"]
noise_dist_old[1]["CEX"] = noise_dist[1]["CEX"]
noise_dist_old[2]["CEX"] = noise_dist[2]["CEX"]


# Save the new noise distribution file
init_path = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
save_path = init_path * "/2_Data_processing/noise_distributions/noise_consum_and_income_and_wealth_deciles_CEX_and_CPS_and_CPS2_and_PSID_and_SCFnpimp_new.jld2"
JLD2.save(save_path, "noise", noise_dist_old)

# Now repeat the same process, but for the confidence intervals 
ci_dist = jldopen("/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/confidence_intervals/ci_consum_and_income_and_wealth_deciles_CEXnpimp.jld2", "r")
ci_dist = ci_dist["ci"]

ci_dist_old = jldopen("/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/confidence_intervals/ci_consum_and_income_and_wealth_deciles_CEX_and_CPS_and_CPS2_and_PSID_and_SCFnpimp_all.jld2", "r")
ci_dist_old = ci_dist_old["ci"]
ci_dist_old["CEX"] = ci_dist["CEX"]

# Save the new noise distribution file
init_path = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
save_path = init_path * "/2_Data_processing/confidence_intervals/ci_consum_and_income_and_wealth_deciles_CEX_and_CPS_and_CPS2_and_PSID_and_SCFnpimp_all2.jld2"
JLD2.save(save_path, "ci", ci_dist_old)

# Now do the same thing for the confidence intervals for the cex data
ci_dist = jldopen("/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/confidence_intervals/ci_consum_and_income_and_wealth_deciles_cexnpimp_new.jld2", "r")
ci_dist = ci_dist["ci"]

ci_dist_old = jldopen("/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/confidence_intervals/ci_consum_and_income_and_wealth_deciles_CEX_and_cex_and_PSIDnpimp.jld2", "r")
ci_dist_old = ci_dist_old["ci"]
ci_dist["CEX"] = ci_dist_old["CEX"]

# Save the new noise distribution file
init_path = dirname(pwd())[end-7:end] == "Dynamics" ? dirname(pwd()) : pwd()
save_path = init_path * "/2_Data_processing/confidence_intervals/ci_consum_and_income_and_wealth_deciles_cexnpimp_new.jld2"
JLD2.save(save_path, "ci", ci_dist)


# Compare treated_dfs to DCT_boot
noise_scf = jldopen("/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/noise_distributions/noise_draws_consum_and_income_and_wealth_ventiles_series_PSID_all.jld2", "r")
noise_scf2 = noise_scf["noise"];
noise_scf2[109:150, :, 1]
treated_dfs[4][109:150, :]

# import the confidence intervals 
psid_int = jldopen("/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/confidence_intervals/ci_draws_consum_and_income_and_wealth_ventiles_series_SCF_all.jld2", "r")
psid_int = psid_int["ci"];

@unpack estimator = model_options
clean_sub_boot_dict!(psid_int)
ci_u, ci_l = construct_confidence_intervals(psid_int, .025, .975, ["levels", "quantiles", "shares"], ["income", "wealth", "consum"], length(year_vec[4]), "PSID", estimator)
ci_u["income"]["quantiles"]

# Observed 15 times? 
e = psid_int["income"]["quantiles"]
cond                   = mapslices(col -> all(col .!==NaN), e, dims = 1)[:]
findall(x -> x == true, cond)
clean_sub_boot_dict!(psid_int)

# Open this file /Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/noise_distributions/noise_draws_consum_and_income_and_wealth_ventiles_series_cex_all.jld2
noise_scf = jldopen("/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/noise_distributions/noise_draws_consum_and_income_and_wealth_ventiles_series_PSID_all.jld2", "r")
noise_scf2 = noise_scf["noise"];

sigma = transform_DCT_boot(noise_scf2, time_p, year_vec[4], freq, freq_type[4], time_dict[4], estimator, 3, measures)
sigma == Σ̂⁻¹²[1063:1416, 1063:1416]

for i in axes(noise_scf2, 3)
    println(sum((isnan).(noise_scf2[73:108, :, i])))
end


s = size(noise_scf2)
DCT_boot_reshaped = reshape(noise_scf2, s[1], s[2] * s[3])

Σ                 = nancov(DCT_boot_reshaped, dims=2)
replace!(Σ, -0.0=>0.0)

Σ̂⁻¹²[1:36, 1:36]
Σ̂⁻¹²[73:108, 73:108]
Σ̂⁻¹²[325:354, 325:354]

treated_dfs[4][73:108, :]
noise_scf2[73:108, :, 1]
treated_dfs[5][73:108, :]
noise_scf = jldopen("/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/noise_distributions/noise_draws_consum_and_income_and_wealth_ventiles_series_SCF_all.jld2", "r")
noise_scf3 = noise_scf["noise"];
noise_scf3[73:108, :, 1]

pool[73:108, 93:152]

@unpack    boot_noise_processes, y, MV = model_elements
MV[1]

block_size = 354
block_matrices = [
    boot_noise_processes[(i-1)*block_size+1:i*block_size, (i-1)*block_size+1:i*block_size]
    for i in 1:5
]

block_matrices[4][325:354, 325:354]
block_matrices[5][325:354, 325:354]

ys               = [y[I, :] for I in Iterators.partition(axes(y, 1), 354)]

i = 4
cond                   = mapslices(col -> any(col .!==NaN), ys[i][73:108, :], dims = 1)[:]
ys[i][73:108, cond]



cond                   = mapslices(col -> any(col .!==NaN), MV[4][335:end, :], dims = 1)[:]
MV[4][335:end, cond]

cond                   = mapslices(col -> any(col .!==NaN), MV[5][335:end, :], dims = 1)[:]
MV[5][335:end, cond]


confidence_intervals["PSID"]["ci_l"]["wealth"]["quantiles"][:, 24]
confidence_intervals["PSID"]["ci_u"]["income"]["quantiles"][:, 24]


@unpack confidence_intervals, func_dict = func_data
c = confidence_intervals["PSID"]["ci_u"]["wealth"]["quantiles"]
d = func_dict["PSID"]["wealth"]["quantiles"]["data"]

cond                   = mapslices(col -> all(col .!==NaN), c, dims = 1)[:]
cond                   = mapslices(col -> all(col .!==NaN), d, dims = 1)[:]
findall(x -> x == true, cond)

if df_name == "SCF"
    temp_cop  = zeros(size(raw_copulas[:, q + p]))
    temp_cop2 = zeros(size(copulas[:, q + p]))
    
    for i in 1:5
        period_dataᵢ           = filter(row -> row.impnum == i, period_data)
        some_cop               = get_copulas(period_dataᵢ, measures, obs_meas, estimator) 
        temp_cop             .+= some_cop[:]
        
        # Treat the copula
        se                     = typeof(estimator) <: SeriesEstimator ? period_dataᵢ : false
        temp_cop2            .+= treat_copula(estimator, some_cop, obs_meas, measures, grid_choice_cop, se) # Since the series estimators has to construct all the sub-copulas, we need all the data 
    end
    raw_copulas[:, q + p] .= temp_cop ./ 5
    copulas[:, q + p] .= temp_cop2 ./ 5
else
    some_cop               = get_copulas(period_data, measures, obs_meas, estimator) 
    raw_copulas[:, q + p] .= some_cop[:]
    
    # Treat the copula
    se                     = typeof(estimator) <: SeriesEstimator ? period_data : false
    copulas[:, q + p]     .= treat_copula(estimator, some_cop, obs_meas, measures, grid_choice_cop, se) # Since the series estimators has to construct all the sub-copulas, we need all the data 
end
end 

@unpack boot_noise_processes = model_elements 


scf_a = boot_noise_processes[end-353:end, end-353:end]

scf_a[73:108, 73:108]


b = boot_noise_processes .== 3.74553e-5 
sum(b)
c = pool[:, 121:152]
nan_rows = all(isnan.(c), dims=2)[:]

# Subset to non-NaN rows
non_nan_matrix = c[.!nan_rows, :]


Σ̂⁻¹²[1:324, 1:324]

Σ̂⁻¹²[325:356, 325:356]

a = rand(4,4)
a[1, :] .=NaN
nancov(a, dims=2)

# Trying to see whether or not estimates are within the bounds 

# SCF copulas 
@unpack confidence_intervals = func_data
cop_estimates = dv["SCF"][1]
confidence_intervals["SCF"]["ci_l"]["copula"]
a= (!isnan).(confidence_intervals["SCF"]["ci_u"]["copula"][:, 1])
confidence_intervals["SCF"]["ci_u"]["copula"][a, 1]

dct(reshape(confidence_intervals["PSID"]["ci_u"]["copula"][a, 236], (10,10))) ./ sqrt(10)

# Actual intervals 
psid_int = jldopen("/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/confidence_intervals/ci_draws_consum_and_income_and_wealth_ventiles_series_SCF_all.jld2", "r")
psid_int = psid_int["ci"];

a = (!isnan).(psid_int["copula"][:, 1, 1])
dct(reshape(psid_int["copula"][a, 1, 1], (10, 10))) ./ sqrt(10)

# the plan now is to re-do the confidence intervals s.t. they are reflecting the actual object. This means I have DCT them and do the scale correction for the copulas based on what is observed or not
datasets = ["PSID", "SCF", "CEX", "CPS2", "CPS"]
datasets = ["SCF"]
measures = ["consum", "income", "wealth"]
grid_size_data_cop = 10

for j in datasets
    # import the intervals 
    cis = jldopen("/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/confidence_intervals/ci_draws_consum_and_income_and_wealth_deciles_series_$(j)_all.jld2", "r")["ci"]
    
    new_sub_boot_dict = Dict()
    cop_n                   = grid_size_data_cop^3
    T                       = size(cis["copula"], 3)
    draws                   = size(cis["copula"], 2)
    cop_dens                = fill(NaN, (cop_n, draws, T))

    objects = ["copula", "consum", "income", "wealth"]
    for obj in objects
        new_sub_boot_dict[obj] = Dict()
        
        if obj != "copula"
            for series in ["levels", "quantiles", "shares"]
                new_sub_boot_dict[obj][series] = deepcopy(cis[obj][series])
            end
        else
            # Create density and then store it 
            for i in 1:draws
                X        = cis[obj][:, i, :]

                # Reshape 'X' into a 3D array
                
                if any(isfinite, X)
                    X        = reshape(X, (11, 11, 11, T))
                    cop_dens[:, i, :] = generate_copula_densities(X, measures, grid_size_data_cop)[:]
                else
                    break
                end
            end
            finish!(p)
            new_sub_boot_dict[obj] = cop_dens
        end
    end
    save_path = "/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/confidence_intervals/ci_draws_consum_and_income_and_wealth_deciles_series_$(j)_all_test.jld2"
    jldopen(save_path, "w") do file
        file["ci"] = new_sub_boot_dict
    end
    # JLD2.save(save_path, "ci", new_ci)
end


    # # for each point in time and for each draw 
    # for a in axes(cis["copula"], 3)
    #     for b in axes(cis["copula"], 2)
    #         # find the non-NaN values 
    #         c = (!isnan).(cis["copula"][:, b, a])
    #         if sum(c) == 100
    #             # DCT the values 
    #             cis["copula"][c, b, a] .= (dct(reshape(cis["copula"][c, b, a], (10, 10))) ./ sqrt(10))[:]
    #         elseif sum(c) == 1000
    #             cis["copula"][c, b, a] .= dct(reshape(cis["copula"][c, b, a], (10, 10, 10)))[:] 
    #         end
    #     end
    # end

    # Save the intervals 
    save_path = "/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/confidence_intervals/ci_draws_consum_and_income_and_wealth_ventiles_series_$(j)_all_test.jld2"
    JLD2.save(save_path, "ci", new_ci)
# end

cis["copula"] 
cis["wealth"]["quantiles"]

cis = jldopen("/Users/lc/Dropbox/Distributional_Dynamics/2_Data_processing/confidence_intervals/ci_draws_consum_and_income_and_wealth_ventiles_series_cex_all_test.jld2", "r")["ci"]
d = (!isnan).(cis["copula"][:, 1, 1])
reshape(cis["copula"][a, 2, 1], (10, 10))

a = (!isnan).(dv["SCF"][1][:,:,:,238])

b = reshape(dv["SCF"][1][a,238], (10, 10))
c = b ./ maximum(b)

@unpack confidence_intervals = func_data
a = (!isnan).(confidence_intervals["SCF"]["ci_u"]["copula"][:, 1])
b = (!isnan).(confidence_intervals["SCF"]["ci_l"]["copula"][:, 1])
reshape(confidence_intervals["SCF"]["ci_u"]["copula"][a, 229], (10, 10))
reshape(confidence_intervals["SCF"]["ci_l"]["copula"][b, 229], (10, 10))


# Why is SCF income and wealth noisy? 
@unpack MV, y = model_elements
scf_y  = y[end-353:end, :]
psid_y = y[end-(353+353):end-353, :] 

cond1                   = mapslices(col -> all(col .!==NaN), MV[4], dims = 1)[:]
MV[4][73:108, cond1]

cond2                   = mapslices(col -> any(col .!==NaN), MV[5], dims = 1)[:]
MV[5][73:108, cond2]

psid_y[73:108, cond1]
scf_y[73:108, cond2]

a = rand(Normal(.542, 10), 10000)

# plot 
p = barhist(a)
Plots.savefig(p, "test.pdf")