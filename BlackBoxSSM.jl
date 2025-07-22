function run_black_box_opt(SSM, param_vector, param_sizes, priors, measures)
    step_ranges = define_step_ranges(param_sizes, priors, param_vector)
    opttime = length(param_vector) * 7
    res = bboptimize(SSM, param_vector; SearchRange=step_ranges, Method=:adaptive_de_rand_1_bin_radiuslimited, MaxTime=opttime, TraceMode=:compact, TraceInterval=60)
    return best_candidate(res)
end

function run_black_box_hyperopt(objective_function, hyperpriors)
    step_ranges = [(hyperpriors[i].μ - 20 * hyperpriors[i].σ, hyperpriors[i].μ + 20 * hyperpriors[i].σ) for i in eachindex(hyperpriors)]
    func_evals = 10000
    res = bboptimize(objective_function; SearchRange=step_ranges, Method=:adaptive_de_rand_1_bin, MaxFuncEvals=func_evals, TraceMode=:compact, TraceInterval=60)

    return best_candidate(res)
end


function define_step_ranges(param_sizes, priors, param_vector)
    step_ranges = Vector{Tuple{Float64,Float64}}(undef, length(param_vector))

    # # For A, B
    # for i in eachindex(priors[1].μ)
    #     step_ranges[i] = (priors[1].μ[i] - priors[1].Σ[i,i], priors[1].μ[i] + priors[1].Σ[i,i]) 
    # end
    # For A 
    l_A = prod(param_sizes[1])
    l_B = prod(param_sizes[2])
    l_C = prod(param_sizes[3])
    l_D = param_sizes[4][1]
    l_Ωf = param_sizes[1][1]
    l_Ωy = param_sizes[2][2]
    l_Σ = param_sizes[6][1]

    for i in 1:l_A
        # step_ranges[i] = (-0.99, 0.99) 
        σ_a = priors[1].Σ[i, i]
        step_ranges[i] = (priors[1].μ[i] - 3 * σ_a, priors[1].μ[i] + 3 * σ_a)
    end

    # For B
    for i in (l_A+1):(l_A+l_B)
        σ_b = priors[1].Σ[i, i]
        step_ranges[i] = (priors[1].μ[i] - 3 * σ_b, priors[1].μ[i] + 3 * σ_b)
    end

    # For C
    for i in (l_A+l_B+1):(l_A+l_B+l_C)
        σ_c = priors[1].Σ[i, i]
        step_ranges[i] = (priors[1].μ[i] - 3 * σ_c, priors[1].μ[i] + 3 * σ_c)
    end

    # For D
    for i in (l_A+l_B+l_C+1):(l_A+l_B+l_C+l_D)
        σ_d = priors[1].Σ[i, i]
        step_ranges[i] = (priors[1].μ[i] - 3 * σ_d, priors[1].μ[i] + 3 * σ_d)
    end

    # For Ωf
    for i in (l_A+l_B+l_C+l_D+1):(l_A+l_B+l_C+l_D+l_Ωf)
        step_ranges[i] = (-3 * priors[2].σ + priors[2].μ, 3 * priors[2].σ + priors[2].μ)
    end

    # For Ωy
    for i in (l_A+l_B+l_C+l_D+l_Ωf+1):(l_A+l_B+l_C+l_D+l_Ωf+l_Ωy)
        step_ranges[i] = (-3 * priors[3].σ + priors[3].μ, 3 * priors[3].σ + priors[3].μ)
    end

    # For Σ
    for (q, i) in enumerate((l_A+l_B+l_C+l_D+l_Ωf+l_Ωy+1):(l_A+l_B+l_C+l_D+l_Ωf+l_Ωy+l_Σ))
        step_ranges[i] = (-3 * priors[3+q].σ + priors[3+q].μ, 3 * priors[3+q].σ + priors[3+q].μ)
    end

    return step_ranges
end
