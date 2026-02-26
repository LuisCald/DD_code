# General helper that works even if start/stop aren't on Q1/Q4 boundaries
# This is for the CEX, every 4 years
function muted_quarters_between(start::QuarterlyDate, stop::QuarterlyDate)
    y0, y1 = year(start), year(stop)

    # keep every 4th year counting from the start year (y0)
    is_kept_year(y) = (y - y0) % 4 == 0

    # for partial edge years (if any)
    qstart(y) = (y == y0) ? quarter(start) : 1
    qstop(y) = (y == y1) ? quarter(stop) : 4

    muted_years = filter(y -> !is_kept_year(y), y0:y1)

    [QuarterlyDate(y, q)
     for y in muted_years
     for q in qstart(y):qstop(y)]
end

GLOBAL_ECON = "1"


# File of Structures 
function retrieve_data_files()
    local files::Dict
    base_path = DATA_PROCESSING * "/"

    files = Dict(
        "HANK a $(GLOBAL_ECON)" => joinpath(base_path, "HANK_PSID_$(GLOBAL_ECON).csv"),
        "HANK b $(GLOBAL_ECON)" => joinpath(base_path, "HANK_CPS_$(GLOBAL_ECON).csv"),
        "HANK c $(GLOBAL_ECON)" => joinpath(base_path, "HANK_CEX_$(GLOBAL_ECON).csv"),
        "HANK d $(GLOBAL_ECON)" => joinpath(base_path, "HANK_SCF_$(GLOBAL_ECON).csv"),
    )

    # "PSID" => joinpath(base_path, "PSID.csv"),
    # # "CEX" => joinpath(base_path, "CEX_all.csv"), # batching all to 4th quarter
    # # "CEX_all_q" => joinpath(base_path, "CEX_all_q.csv"),
    # # "CEX" => joinpath(base_path, "CEX.csv"),
    # # "CPS" => joinpath(base_path, "CPS.csv"),
    # # "CPS2" => joinpath(base_path, "CPS2.csv"),
    # # "SIPP1" => joinpath(base_path, "SIPP1.csv"),
    # # "SIPP2" => joinpath(base_path, "SIPP2.csv"),

    return files
end


function retrieve_aggregate_data(sheet)
    local aggregate_data::DataFrame
    base_path = DATA_PROCESSING * "/"
    
    # Read as .csv if .xlsx fails
    aggregate_data = CSV.read(joinpath(base_path, "HANK_shocks_economy_$(GLOBAL_ECON).csv"), DataFrame)
    # aggregate_data = DataFrame(XLSX.readtable(joinpath(base_path, "aggregates_HHs_NPs.XLSX"), sheet, header=true,))
    
    return aggregate_data

end


function retrieve_rgdp()
    local aggregate_data
    base_path = DATA_PROCESSING * "/"

    # aggregate_data = DataFrame(XLSX.readtable(joinpath(base_path, "inflation_corrected_correction_series.XLSX"), "data", header=true,))
    aggregate_data = CSV.read(joinpath(base_path, "HANK_truth_$(GLOBAL_ECON).csv"), DataFrame)

    return aggregate_data
end


function retrieve_data(files::Dict)
    data = Vector{DataFrame}(undef, length(keys(files)))
    df_names = String[]

    for (i, k) in enumerate(sort(collect(keys(files))))
        # data[i] = XLSX.readxlsx(files[k])
        push!(df_names, k)
        data[i] = CSV.read(files[k], DataFrame)
    end
    return (data=data, df_names=df_names)
end


# Create an abstract type for all estimators 
abstract type AbstractEstimator end
abstract type AbstractPrior end
abstract type FakeEstimator <: AbstractEstimator end

# Data 
@with_kw struct ObservedData{S<:String,D<:DataFrame,NT<:NamedTuple} # X<:XLSX.XLSXFile}
    files::Dict{S,S} = retrieve_data_files()
    # agg_data::D = retrieve_aggregate_data("stationary_series")
    agg_data::D = retrieve_aggregate_data("HANK_shocks_economy_$(GLOBAL_ECON)")
    gdp_series::D = retrieve_rgdp()  # not just gdp, but many measure specific correction series
    df_vec::NT = retrieve_data(files)
end


@with_kw mutable struct HistogramEstimator{I<:Integer,S<:String} <: AbstractEstimator
    grid_pcf::I
    grid_cop::I
    grid_type_pcf::S
    grid_type_cop::S
end


@with_kw mutable struct KernelEstimator{I<:Integer,S<:String} <: AbstractEstimator
    grid_pcf::I
    grid_cop::I
    grid_type_pcf::S
    grid_type_cop::S
end


@with_kw mutable struct SeriesEstimator{I<:Integer} <: AbstractEstimator
    grid_pcf::I # equivalent to order here!! if 'data_driven_order' is true, set to zero. order will be maximum order of the polynomial # TODO: a fixed order for now
    grid_cop::I # equivalent to order here!! 
    integral_pcf_grid::I
    integral_cop_grid::I
end


struct FakeHistEstimator <: FakeEstimator end # for confidence_intervals 
struct FakeKernelEstimator <: FakeEstimator end # for confidence_intervals 


# order for copula is 'grid_cop' + 1
@with_kw mutable struct ModelOptions{AE<:AbstractEstimator,AP<:AbstractPrior,B<:Bool,S<:String,I<:Integer,VS<:Vector{String},D<:Dict,DS<:Dict{String,String}} #,
    estimator::AE = SeriesEstimator(grid_pcf=11 + 1, grid_cop=11 + 1, integral_pcf_grid=5, integral_cop_grid=5) # 11 + 1 = order 11
    # estimator::AE = SeriesEstimator(grid_pcf=11 + 1, grid_cop=11 + 1, integral_pcf_grid=10, integral_cop_grid=10) # or "KDE" or "series" ,     
    prior::AP = Minnesota(hyperparameters=[0.05, 0.1, 0.5, 0.1, 2.0, 0.9, 0.9]) # We use: 1,2,4,6
    measures::VS = sort(["consum", "wealth", "income"])
    number_of_dfs::I = 4
    plot_proof::B = false # to plot the proof of concept
    case::S = "A non-diag" # ["diag", "A non-diag", "A, Σ non-diag"]  
    blind_to::D = Dict() #Dict("PSID" => ["wealth"])
    estimation_object::S = "copulas and percentile functions" #"copulas and percentile functions"
    information::S = "partial"  # OR "partial"
    lags::I = 1
    agg_lags::I = 4
    freq::I = 4  # the frequency you would like e.g., quarterly
    agg_freq::I = 4
    constant::B = false  # this is a time trend and a quadratic trend. 
    measurement_error::S = "one per object, per dataset" #"one per income quantile, per dataset" #"one per measurement" # "one per object" # "one per object, per dataset
    collapse::B = true  # # TODO: see where relevant  
    errors_process::S = "one per object, per dataset" #"average" # for the measurement error stipulation 
    pre_multiply::B = true
    pca_perspective::S = "frequentist" # or "bayesian"
    best_aggs::B = false


    # Data Treatment 
    equivalized::B = false
    bottom_coded::Vector{Any} = []
    rm_seasonality::B = true # TODO: basically if we use the entire CEX

    data_cutoffs::DS = Dict("begin" => "", "end" => "") # Define boundaries (extensive margin)
    data_to_mute = Dict("begin" => "", "end" => "")
    # data_to_mute = Dict("CEX" => muted_quarters_between(QuarterlyDate(1984, 1), QuarterlyDate(2021, 4))) # Dict("begin" => "", "end" => "")
    logit_transform::B = false # TODO: doesnt work anyway
    compare_to_other_est::B = false # compare with other models
    tag::S = " HANK $(GLOBAL_ECON)"  # descriptive tag for the model
    # tag::S = " additional factors"  # descriptive tag for the model

    # tag::S = " HANK" # in process
    # tag::S = " excluding housing cycle" # in process
    # tag::S = " every 4 years" # TODO: change data_to_mute
    # tag::S = " excluding recent 20 quarters"

    # Models to do: changing distributional factors
    ## 3 dist and 10 aggs = " less DF and AF" (in process)
    ## 6 factors, 25 aggs = " 6 factors" (in process)
    ## 7 factors, 25 aggs = " 7 factors" (in process)

    # Models to do: changing aggs
    ## 15 aggregate factors = " less AF" (in process)
    ## 20 aggregate factors = " more AF" (in process)

    # Other Models to do:
    ## " excluding housing cycle short" (in process) # Dict("begin" => QuarterlyDate(2007, 4), "end" => QuarterlyDate(2011, 4))
    ## " Γ all" = 10 factors (in process)
    ## " Γ all 85" = 12 factors (in process)

    # Once done, do:
    ## set compare_to_other_est to true and then run export_table_to_tex_with_strings(measures, :from_mcmc) from the baseline
    ## do historical decomposition
    ## do variance decomposition
    ## do " best aggs" case
end


# List of models: PP CEX is for adding a projection factor from CEX
# " PP CEX excluding housing cycle"
# " PP CEX excluding housing cycle short"
# " PP CEX excluding recent 20 quarters"
# " PP CEX every 4 years"

# " PP SCF"
# " Γ estimated"
# " all AF"
# " less AF" (less aggregate factors) #TODO: in process

# " Γ all 85"
# " Γ all"
# " more AF"
# " less DF and AF" (less distributional and aggregate factors) #TODO: in process

# " 6 factors"
# " 7 factors"
# " additional factors" 
# " excluding housing cycle wealth" # Dict("begin" => QuarterlyDate(2004, 1), "end" => QuarterlyDate(2009, 4))
# " every 4 years" # data_to_mute = Dict("CEX" => muted_quarters_between(QuarterlyDate(1984, 1), QuarterlyDate(2021, 4)))
# " excluding recent 20 quarters" # #Dict("begin" => QuarterlyDate(2020, 1), "end" => QuarterlyDate(2024, 1)) 

# " less AF"
# # " excluding housing cycle" # Dict("begin" => QuarterlyDate(2004, 4), "end" => QuarterlyDate(2009, 4))
# # " excluding housing cycle short" # Dict("begin" => QuarterlyDate(2007, 4), "end" => QuarterlyDate(2011, 4))
# # " higher order15" 

# check compare_to_other_est, tag, data_to_mute, data above


@with_kw struct MCMCOptions{B<:Bool,S<:String,I<:Int64,F<:Float64}
    nsave::I = 100000  # MCMC
    nburn::I = 200000   # MCMC 
    chains::I = 4       # MCMC
    adaptive_rwmh::B = true
    compute_hessian::B = true   # not necessary with Gibbs 
    mhscale::F = 0.35
    sampler::S = "MH"
    mcmc_jsd_draws::I = 100    # for state trajectory in kalman  
    n_jsd_draws::I = 10    # for state trajectory in kalman-gibbs
    thinning_step::I = 1  # 1 = no thinning -- controversial ... 
end

@with_kw mutable struct Minnesota{VF<:Vector{Float64}} <: AbstractPrior
    hyperparameters::VF = [0.2, 0.3, 0.01, 5, 2.0, 0.90, 0.90] # trying [0.2, 0.3, .01, 5, 2.0, 0.95]
end


# Structure of Parameters 
# @with_kw mutable struct MethodOptions{B<:Bool, S<:String, I<:Integer, F<:Float64, D<:Dict, DS<:Dict{String, String}}
#     case::S                     = "A non-diag" # ["diag", "A non-diag", "A, Σ non-diag"]  
#     equivalized::B              = false
#     bottom_coded::Vector{Any}   = []
#     grid::I                     = 10   # only a single number since axes have to be equal since uniform distributions  
#     grid_type::S                = "uniform"  # or "chebyshev"
#     functional_data::B          = false
#     measures::Vector{S}         = sort(["income", "consum", "wealth"]) # sort(["income", "consum", "liquid"])
#     blind_to::D                 = Dict() #Dict("PSID" => ["wealth"])
#     estimation_object::S        = "copulas and percentile functions" #"copulas and percentile functions"
#     information::S              = "partial"  # OR "partial"
#     percentile_estimation::S    = "naive" #TODO: delete 
#     minnesota_params::Vector{F} = [0.2, 0.3, .01, 5, 2.0, 0.90] # trying [0.2, 0.3, .01, 5, 2.0, 0.95]
#     lags::I                     = 1
#     freq::I                     = 4  # the frequency you would like e.g., quarterly
#     agg_freq::I                 = 4
#     constant::B                 = false  # this is a time trend and a quadratic trend. 
#     number_of_dfs::I            = 5
#     sampler::S                  = "MH"
#     n_jsd_draws::I              = 10    # for state trajectory in kalman-gibbs
#     # var_method::S               = "Bayesian"  # or "LS" for least squares 
#     pca_perspective::S          = "Frequentist"
#     pca_method::S               = "linear"
#     std_method::S               = "z-score"  # or "norm"
#     measurement_error::S        = "one per object, per dataset" #"one per income quantile, per dataset" #"one per measurement" # "one per object" # "one per object, per dataset
#     collapse::B                 = true  
#     plot_proof::B               = false
#     reconstruction_to_show      = "SCF" #TODO: remove
#     rm_seasonality::B           = true # TODO: basically if we use the entire CEX
#     errors_process::S           = "one per object, per dataset" #"average" # 
#     scf_ci::S                   = "npimp"
#     pre_multiply::B             = true
#     data_cutoffs::DS            = Dict("begin" => "", "end" => "") # data cutoff. Estimation will always run until 2021Q4 or whenever the end of the aggregates is 
#     data_to_mute::DS             = Dict("begin" => "", "end" => "") # Dict("begin" => QuarterlyDate(2017, 1), "end" => QuarterlyDate(2021, 4))
#     logit_transform::B          = false # TODO: doesnt work anyway 
# end

mutable struct TimeParams{D<:Dict,VV<:Vector{Vector{Int64}},I<:Int64,VS<:Vector{String}}
    year_vec::VV
    tmin::D
    tmax::D
    tot_years::I
    tot_periods::I
    time_dict::D
    freq_type::VS
end

# TODO: Norm Explanation
# Normalization is a good technique to use when you do not know the distribution of your 
# data or when you know the distribution is not Gaussian (a bell curve). 
# Normalization is useful when your data has varying scales and the algorithm 
# you are using does not make assumptions about the distribution of your data, 
# such as k-nearest neighbors and artificial neural networks. Since PCA is non-parametric, this is ok!

# Standardization assumes that your data has a Gaussian (bell curve) distribution. 
# This does not strictly have to be true, but the technique is more effective if 
# your attribute distribution is Gaussian. Standardization is useful when your data 
# has varying scales and the algorithm you are using does make assumptions about your 
# data having a Gaussian distribution, such as linear regression, logistic regression,
# and linear discriminant analysis.



# The theory behinds Markov Chain Monte Carlo (the family of algorithms that includes the Metropolis
#  algorithm as a special case) guarantees that (under certain conditions) Metropolis will give you the right answer. 
#  The acceptance rate is related only to the number of steps you will need to find the correct answer
# (with a certain level of precision). If the acceptance rate is too low or too high it will take longer 
# time, but eventually you will get the correct answer. The intuition is as follows; you will get high 
# acceptance when you make small movements not too far away of your current position, but if you move 
# slowly it will take a lot of time to sample all the distribution. On the contrary, if you try to 
# make larger movements chances are high you propose movements to regions of low probability and 
# then you will get a lot of rejections. So in practice even when you will "eventually" get the 
# right answer is a good idea to tune the acceptance rate to provide better performance in this 
# paper by Gelman they provide evidence that a acceptance rate of ~0.23 is optimal. PyMC3 has an 
# option to tune a Metropolis algorithm based on monitoring the acceptance rate.

@with_kw struct DiagnosticOptions{F<:Float64,B<:Bool}
    alpha::F = 0.05   # to estimate the upper limits of scale reduction factors, level of significance 
    mpsrf::B = false  # multivariate potential scale reduction factor (joint convergence). This factor will not be calculable if any one of the parameters in the output is a linear combination of others.
    transform::B = false  # produce output that is more normally distributed, an assumption of the PSRF formulations.
    first_opt::F = 0.1    # proportion of iterations to include in the first window.
    last_opt::F = 0.5    # proportion of iterations to include in the last window.
    etype = :imse  # method for computing Monte Carlo standard errors. Default: initial monotone sequence estimator
    # q::
    # r::
    # s::
    # eps:: 
end

# Structure of objects necessary for estimation 
mutable struct StateSpaceModel{A<:Array{Float64},V<:Vector{Float64},VV<:Vector{Vector{Float64}},M<:Matrix{Float64},VM<:Vector{Matrix{Float64}},I<:Int64}
    #    years::I
    MV::VM            # Vector of measurements, split by dataset 
    y::M              # Matrix of all measurements 
    G::VM             # Vector of selection matrices, one per period
    Gⱼ::VM             # Selection matrix for the j-th measurement
    u::M
    pcs::M
    means::VV
    stds::V
    proj::M
    u_proj::M
    factor_count::I
    agg_count::I
    n_less_than_one::I
    βs::A
    trend::VM
    boot_noise_processes
end

mutable struct FunctionalData{D<:Dict,VS<:Vector{String}}
    func_dict::D
    confidence_intervals::D
    year_vec
    data_sources::VS
end

# Structure for defining the prior 
#mutable struct Prior{VS<:Vector{Distribution{Multivariate, Continuous}}, V<:Vector{AbstractMatrix{Float64}}}
mutable struct Prior
    priors
    param_vector
end


# function set_prior(model_elements::StateSpaceModel, method_options::MethodOptions)
#     @unpack exogenous, constant, lags, number_of_dfs, grid, measures = method_options
#     @unpack  factor_count, agg_count = model_elements

#     # State Equation
#     state_count = factor_count + agg_count
#     exo_count  = size(exogenous, 1)   
#     var_count  = state_count * lags + exo_count + sum(constant) 

#     # Measurement Equation
#     meas_count = agg_count + number_of_dfs * (grid^length(measures) + grid * length(measures))  # number of possible measurements in the measurement eq. + length of Diagonal(Σ)
#     noisy_meas = meas_count - agg_count

#     priors   = Vector{Sampleable}(undef, 2)
#     param_vector = zeros(length(prior_A) + length(S_prior) - )

#     # # For the Normal Distribution 
#     # A_prior  = zeros(state_count, var_count)
#     # V_prior  = zeros(length(prior_A), length(prior_A))

#     # # For the inverse wishart of Ω and Σ
#     # Ω_prior  = zeros(state_count, state_count)
#     # Σ_prior  = zeros(noisy_meas, noisy_meas)

#     return Prior(priors, A_prior, var_Ω, V_prior, param_vector)
# end

# function set_up_mcmc(method_options::MethodOptions)
#     """Pre-allocates parameter matrices for MCMC."""

#     @unpack grid, measures, aggregates, var_count, state_count, nsave, chains, number_of_dfs = method_options
#     m_size       = (grid^length(measures) + grid * length(measures)) * number_of_dfs + length(aggregates) - 1 # TODO: number_of_dfs??
#     ALPHA        = zeros(Float64, var_count * state_count)
#     alpha_draws  = zeros(Float64, nsave, var_count * state_count, chains)
#     OMEGA_draws  = zeros(Float64, nsave, state_count * state_count, chains)
#     R_draws      = zeros(Float64, nsave, m_size * m_size, chains)

#     return PosteriorDistribution(ALPHA, alpha_draws, OMEGA_draws, R_draws)
# end

# Store results from the model estimation 
struct SmootherResults{MF<:Matrix{Float64},VMF<:Vector{Matrix{Float64}}}
    x_filtered::MF
    sigma_filtered::VMF
    x_updated::MF
    sigma_updated::VMF
    x_smoothed::MF
    sigma_smoothed::VMF
    cross_covariances::VMF
    dε_smoothed::MF
end


# Posterior Distribution 
mutable struct PosteriorDistribution{F<:Float64}
    alpha::Vector{F}
    alpha_draws::Array{F}
    OMEGA_draws::Array{F}
    SIGMA_draws::Array{F}
end


# Posterior Results
mutable struct GibbsPosteriorResults{A<:Array{Float64,3},V<:Vector{Float64}}
    A_chains::A
    Ω_chains::A
    Σ_chains::A
    par_vec::V

    function GibbsPosteriorResults(A_chains, Ω_chains, Σ_chains)
        A_summary = summarize(A_chains)[:, :mean]
        Ω_summary = summarize(Ω_chains)[:, :mean]
        Σ_summary = summarize(Σ_chains)[:, :mean]
        par_vec = vcat(A_summary, Ω_summary, Σ_summary)  # they are already vectors 

        new{A,V}(A_chains, Ω_chains, Σ_chains, par_vec)
    end
end



mutable struct RWMHPosteriorResults{A<:Array{Float64,3},VV<:Vector{Vector{Float64}}}
    parameter_chain::A
    par_mcmc::VV
    bounds::VV

    function RWMHPosteriorResults(parameter_chain)
        par_mcmc = []

        for ch in axes(parameter_chain, 3)
            push!(par_mcmc, mean(Chains(parameter_chain[:, :, ch]))[:, 2])
        end
        c = Chains(parameter_chain)
        # par_mcmc    = mean(c)[:, 2]
        quants = quantile(c, q=[0.025, 0.975])
        quant_lb = quants[:, 2]
        quant_ub = quants[:, 3]
        new{Array{Float64,3},Vector{Vector{Float64}}}(parameter_chain, par_mcmc, [quant_lb, quant_ub])
    end
end

mutable struct DiagnosticResults{NM<:NamedTuple,V<:Vector{Float64},AV<:AbstractVector{ChainDataFrame{NamedTuple{(:parameters, :zscore, :pvalue),Tuple{Vector{Symbol},Vector{Float64},Vector{Float64}}}}}}
    acceptance_rate::V  # the acceptance rate for each chain 
    ess::V              # the effective sample size for each parameter 
    gelman::NM          # R̂
    geweke::AV           # To see if the burn-in is too short, checks stationarity by comparing two samples. 
end


# Figures 
@with_kw struct PlotParams{I<:Int}
    years::Vector{I} = [1971, 1992, 2010]
    grid::I = 10
    draws::I = 2000
end


# Calling Structs 
const plot_params = PlotParams()
# const method_options         = MethodOptions()
const model_options = ModelOptions()
const diagnostics_options = DiagnosticOptions()
const mcmc_options = MCMCOptions()
const obs_data = ObservedData()


# # Once everything is finished 
# mutable struct Results{F<:Float64}
#     x_updated::Matrix{F}
#     sigma_updated::Array{F}
#     x_smoothed::Matrix{F}
#     sigma_smoothed::Array{F}
#     posterior::F
#     accept_rate::F
#     par_final::Vector{F}
#     draws::Matrix{F}
#     inv_hessian::Matrix{F} 
# end

