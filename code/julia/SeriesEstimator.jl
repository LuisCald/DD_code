include(joinpath(@__DIR__, "DistributionalDynamics.jl"))

using Integrals

# Output directory for Figure 2 plots (Overleaf)
const PLOT_DIR = "/Users/lc/Dropbox/Apps/Overleaf/Distributional Dynamics/Plots/proof_of_concept"
const RESULTS_DIR = "/Users/lc/Dropbox/Distributional_Dynamics/7_Results/order_analysis"

# Plan of attack: create 3D object <-> 2D object mapping
# Step 1: estimate the 3D copula
# Step 2: find coefficients for that 
# Step 3: estimate the 2D copula
# Step 4: find coefficients for that
# Step 5: compare coefficients -> figure out mapping 

grid = 10
dom = [1.0, 0.0] # upper bound first -> package says so
p = nodes(grid, :chebyshev_nodes, dom)

# Get PSID data 
@unpack df_vec = obs_data
df = df_vec[1][1]
df = df[df[:, :year].==2019, :]

# Cleaning out NaNs
non_missing = coalesce.(df, NaN)
non_missing = filter("income" => !isnan, select(df, ["income", "wealth", "consum", "weight", "id"]))
filter!("weight" => !isnan, non_missing)
filter!("income" => !isnan, non_missing)
filter!("wealth" => !isnan, non_missing)
filter!("consum" => !isnan, non_missing)

df_3D = select(non_missing, ["consum", "income", "wealth", "weight"])
df_2D = select(non_missing, ["income", "wealth", "weight"])



# Define the orthonormal shifted Legendre polynomials as in the paper 
function legendre_polynomial(m, x)
    if m == 0
        return 1.0
    elseif m == 1
        return x
    else
        P_prev_prev = 1.0
        P_prev = x
        P_current = 0.0

        for n in 2:m
            P_current = ((2n - 1) * x * P_prev - (n - 1) * P_prev_prev) / n
            P_prev_prev, P_prev = P_prev, P_current
        end

        return P_current
    end
end



# Define the Legendre polynomial of degree m
function Q_m(m, x)
    L_m = legendre_polynomial(m, 2x - 1)

    return sqrt(2m + 1) * L_m
end


function integrate_legendre_polynomial(m, u)
    if m == 0
        return u
    else
        integral_cop, _ = quadgk(u -> Q_m(m, u), 0, u, rtol=1e-3)

        return integral_cop
    end
end

# Function to integrate the shifted Legendre polynomial
function I_m(m, u)
    return integrate_legendre_polynomial(m, u)
end




# Estimate the copula coefficients ρ_m
function estimate_rho(R, weights, m)
    n, d = size(R)
    rho_m = 0.0

    if all(iszero, m)
        return 1.0
    end

    # If d-1 elements of m are zero, then the product is zero
    if sum(m .== 0) >= d - 1
        return 0.0
    end

    Threads.@threads for i in 1:n
        product = 1.0

        for j in 1:d
            # product *= weights[i] * Q_m(m[j], R[i, j])
            product *= Q_m(m[j], R[i, j])
        end

        rho_m += weights[i] * product
        # rho_m += product
    end

    return rho_m / sum(weights)
    # return rho_m / n
end

# Helper function to compute ranks
function rankdata(a)
    order = sortperm(a)
    ranks = similar(order)
    ranks[order] .= 1:length(a)
    return ranks
end

# Construct the N-th order estimator for copula density
function copula_density_estimator(X, weights, N, u)
    d = size(X, 2)

    ranges = [(0:N[j]) for j in 1:d]
    cl = length(collect(Iterators.product(ranges...)))
    c_N = 0.0


    # Multi-thread this 
    Threads.@threads for (xx, m) in collect(enumerate(Iterators.product(ranges...)))
        m = Tuple(m)
        rho_m = estimate_rho(X, weights, m)
        product = 1.0

        for j in 1:d
            product *= Q_m(m[j], u[j])
        end

        c_N += rho_m * product
    end

    # Apply positivity correction
    c_N = maximum([0, c_N])

    return c_N
end

function copula_cdf_estimator(X, weights, N, u)
    d = size(X, 2)

    ranges = [(0:N[j]) for j in 1:d]
    C_N = 0.0

    Threads.@threads for m in collect(Iterators.product(ranges...))
        m = Tuple(m)
        rho_m = estimate_rho(X, weights, m)
        product = 1.0
        for j in 1:d
            product *= I_m(m[j], u[j])
        end
        C_N += rho_m * product
    end

    return C_N
end


function copula_function(cop_coefs, N, u)
    d = length(N)

    ranges = [(0:N[j]) for j in 1:d]
    cl = length(collect(Iterators.product(ranges...)))
    c_N = 0.0


    # Multi-thread this 
    Threads.@threads for (xx, m) in collect(enumerate(Iterators.product(ranges...)))
        m = Tuple(m)
        product = 1.0

        for j in 1:d
            product *= Q_m(m[j], u[j])
        end

        c_N += cop_coefs[xx] * product
    end

    # Apply positivity correction
    c_N = maximum([0, c_N])

    return c_N
end


function get_copula_coefficients(X, W, N)
    d = size(X, 2)

    ranges = [(0:N) for j in 1:d] # TODO: can be made more flexible
    cl = length(collect(Iterators.product(ranges...)))
    c_N = 0.0
    rho_m = zeros(cl)

    # Multi-thread this 
    Threads.@threads for (xx, m) in collect(enumerate(Iterators.product(ranges...)))
        tup_m = Tuple(m)
        rho_m[xx] = estimate_phi(X, W, tup_m)
    end

    return rho_m
end


# Example usage with a sample data X
X = Matrix(select(df_3D, ["consum", "income", "wealth"]))  # Example bivariate data
X2 = Matrix(select(df_2D, ["income", "wealth"]))  # Example bivariate data

# rank the data 
for i in 1:size(X, 2)
    X[:, i] = rankdata(X[:, i]) / (size(X, 1) + 1)
end

for i in 1:size(X2, 2)
    X2[:, i] = rankdata(X2[:, i]) / (size(X2, 1) + 1)
end


N = (10, 10, 10)  # Example truncation order
N2 = (10, 10)  # Example truncation order

function copula_cdf_estimator(X, u)
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
    Threads.@threads for ci in CartesianIndices(m_combos)
        m = Tuple(m_combos[ci])
        rho_m = X[ci]
        product = 1.0

        for j in 1:d
            product *= I_m(m[j], u[j])
        end
        C_N += rho_m * product
    end

    return C_N
end


d = length(size(X)) - 1

x = select_grid_points(10)
x[end] = x[end] - 1e-6
XX = [[x[i], x[j]] for i in eachindex(x), j in eachindex(x)]



X = Matrix(select(df_3D, ["consum", "income", "wealth", "weight"]))  # Example bivariate data
m_combos = generate_unique_combinations(["income", "consum", "wealth"])
filter!(x -> length(x) == 2, m_combos)

for combo in m_combos
    # Example usage with a sample data X
    X = Matrix(select(df_3D, combo))  # Example bivariate data

    # rank the data 
    for i in 1:size(X, 2)
        X[:, i] = rankdata(X[:, i]) / (size(X, 1))
    end

    rho = get_copula_coefficients(X, df_3D[:, :weight], 20)
    rho_mat = reshape(rho, Int(sqrt(size(rho, 1))), Int(sqrt(size(rho, 1))))
    rho_mat = rho_mat[2:end, 2:end]

    tag = join(combo, "_")

    Plots.surface(1:20, 1:20, rho_mat, xlabel=L"\textrm{Order}", ylabel=L"\textrm{Order}", xformatter=:latex, yformatter=:latex, zformatter=:latex, zlabel="",
        camera=(30, 10),
        size=(500, 500),
        color=:winter,
        legend=false,
        guidefontsize=18, tickfontsize=16, legendfontsize=16,
        display_option=Plots.GR.OPTION_SHADED_MESH)
    outpath = joinpath(PLOT_DIR, "copula_weight_$tag.pdf")
    Plots.savefig(outpath)
    run(`pdfcrop $outpath $outpath`)
end



# For the percentile functions 

# Function to compute Legendre polynomials up to a given order
function legendre_polynomials(x, order)
    P = [1, x]
    for n in 2:order
        Pn = ((2n - 1) * x * P[end] - (n - 1) * P[end-1]) / n
        push!(P, Pn)
    end
    return P
end


function series_estimator(data, weights, order, p_values)
    # Estimate Phi using the legendre polynomials
    Phi = zeros(length(data), order + 1)

    s_weights = cumsum(weights) / sum(weights)

    for i in eachindex(data)
        for j in 0:order
            # Phi[i, j+1]          = Q_m(j, normalize_to_0_1(data[i], min_data, max_data)) # Q_m places [0,1] data internally to [-1, 1]
            Phi[i, j+1] = Q_m(j, s_weights[i]) # Q_m places [0,1] data internally to [-1, 1]
        end
    end
    # Estimate the coefficients, (Phi' * Phi)^(-1) * Phi' * y
    # coefficients = inv(transpose(Phi) * Phi) * transpose(Phi) * data
    # Incorporate weights into the Phi matrix
    W = Diagonal(sqrt.(weights))  # Use square root of weights for correct weighting
    Phi_weighted = W * Phi

    # More stable and efficient way to estimate the coefficients
    # Weighted least squares solution
    # coefficients = (transpose(Phi_weighted) * Phi_weighted) \ (transpose(Phi_weighted) * (W * data))
    # coefficients = (transpose(Phi_weighted) * (W * data))

    coefficients = zeros(order + 1)
    for j in 1:order+1
        for i in 1:length(data)
            coefficients[j] += weights[i] * Phi[i, j] * data[i]
        end
        coefficients[j] /= sum(weights)
    end
    # cdf_est = Phi * coefficients

    basis = zeros(length(p_values), order + 1)
    for (i, p) in enumerate(p_values)
        for j in 0:order
            basis[i, j+1] = Q_m(j, p)
        end
    end

    quants = basis * coefficients

    # return cdf_est
    return quants
end

function get_quantile_weights(data, weights, order)
    # Estimate Phi using the legendre polynomials
    Phi = zeros(length(data), order + 1)
    min_data = nanminimum(data)
    max_data = nanmaximum(data)
    s_weights = cumsum(weights) / sum(weights)

    for i in eachindex(data)
        for j in 0:order
            Phi[i, j+1] = Q_m(j, s_weights[i]) # Q_m places [0,1] data internally to [-1, 1]
        end
    end

    # Incorporate weights into the Phi matrix
    W = Diagonal(sqrt.(weights))  # Use square root of weights for correct weighting
    Phi_weighted = W * Phi

    coefficients = zeros(order + 1)
    for j in 1:order+1
        for i in 1:length(data)
            coefficients[j] += weights[i] * Phi[i, j] * data[i]
        end
        coefficients[j] /= sum(weights)
    end

    return coefficients
end

function eval_quantile_function(coefficients, order, u)

    basis = zeros(1, order + 1)

    for j in 0:order
        basis[j+1] = Q_m(j, u)
    end

    quants = basis * coefficients

    return quants
end

# Example usage:
measure_of_choice = :wealth
sort!(df_3D, measure_of_choice)
data = df_3D[!, measure_of_choice][:]
weights = df_3D[!, :weight][:]

@unpack gdp_series = obs_data
gdp_series[!, "date"] = QuarterlyDate.(gdp_series[!, "time"])
a = filter(x -> x.date == QuarterlyDate(2019, 4), gdp_series)
data = inverse_hyperbolic_sine(data ./ a[1, "$(String(measure_of_choice))_per_hh"])

# Estimate the CDF
p_values = collect(0.02:0.01:1)
p_values[end] = 0.9999

# Generate the quantile function
emp_est = []
s_weights = cumsum(weights) / sum(weights)

for p in p_values
    est = [meas for (meas, cdf) in zip(df_3D[!, measure_of_choice], s_weights) if cdf >= p][1]
    push!(emp_est, est)
end

emp_est = inverse_hyperbolic_sine(emp_est)
orders_to_estimate = collect(3:20)
mse = zeros(length(orders_to_estimate))
mse_bot = zeros(length(orders_to_estimate))
mse_mid = zeros(length(orders_to_estimate))
mse_top = zeros(length(orders_to_estimate))

# Plot the results
for (mse_i, i) in enumerate(orders_to_estimate)
    Plots.plot()
    q_est = series_estimator(data, weights, i, p_values)
    q_est = inverse_hyperbolic_sine(reverse_inverse_hyperbolic_sine(q_est) .* a[1, "$(String(measure_of_choice))_per_hh"])
    Plots.plot!(axes(q_est), q_est, color=:black, ls=:dash, xlabel="x", ylabel="CDF", title="CDF Estimation using Polynomial Basis")
    # Plots.plot!(axes(data), data, color=:red, ls=:dot, label="Empirical CDF")
    Plots.plot!(axes(emp_est), emp_est, color=:red, ls=:dot, label="Empirical CDF")
    Plots.savefig("cdf_est$i.pdf")
    # Compute the MSE between empirical and estimated quantiles
    mse[mse_i] = mean((emp_est .- q_est) .^ 2)

    # Compute the MSE for the bottom, middle, and top quantiles
    mse_bot[mse_i] = mean((emp_est[1:50] .- q_est[1:50]) .^ 2)
    mse_mid[mse_i] = mean((emp_est[51:90] .- q_est[51:90]) .^ 2)
    mse_top[mse_i] = mean((emp_est[91:end] .- q_est[91:end]) .^ 2)
end

# Concatenate the MSE results
mse_results = hcat(mse, mse_bot, mse_mid, mse_top)

# Plot the MSE
Plots.plot(orders_to_estimate, mse_results, xlabel=L"\textrm{Order}", ylabel=L"\textrm{MSE}",
    xformatter=:latex, yformatter=:latex,
    color=:auto, marker=:circle, markersize=5, label=[L"\textrm{All}" L"\textrm{Bottom}" L"\textrm{Middle}" L"\textrm{Top}"],
    legend=:topright,
    grid=true,
    guidefontsize=22, tickfontsize=20, legendfontsize=20,
    bottom_margin=5Plots.mm,
    size=(800, 400))
Plots.savefig(joinpath(RESULTS_DIR, "mse_qf_$(String(measure_of_choice)).pdf"))



# Integration 
n_coefs = 21
q_weights = get_quantile_weights(data, weights, n_coefs)
integral, err = quadgk(u -> reverse_inverse_hyperbolic_sine(eval_quantile_function(q_weights, n_coefs, u)) .* a[1, :income_per_hh], 0.1, 0.2, rtol=1e-8)

to_store = rand(n_coefs)
for (i, j) in enumerate(0.0:0.01:0.099)
    to_store[i] = eval_quantile_function(q_weights, n_coefs, j)[1]
end


# Plotting the coefficients 

for (i, m) in enumerate(sort(["income", "consum", "wealth"]))
    sort!(df_3D, m)
    data = df_3D[!, m][:]
    weights = df_3D[!, :weight][:]

    @unpack gdp_series = obs_data
    gdp_series[!, "date"] = QuarterlyDate.(gdp_series[!, "time"])
    a = filter(x -> x.date == QuarterlyDate(2019, 4), gdp_series)

    data = inverse_hyperbolic_sine(data ./ a[1, m*"_per_hh"])
    w = get_quantile_weights(data, weights, n_coefs)

    Plots.plot(axes(w), w, color=:blue, xformatter=:latex, yformatter=:latex, guidefontsize=22, xtickfontsize=20, ytickfontsize=20, legendfontsize=20, bottom_margin=5Plots.mm, ls=:dot, label="", xlabel=L"\textrm{Order}", ylabel="", size=(700, 400))
    Plots.scatter!(axes(w), w, color=:blue, xformatter=:latex, yformatter=:latex, ls=:dot, markersize=4, label=i == 1 ? L"\textrm{Legendre \,\, Coefficients}" : "", xlabel=L"\textrm{Order}", ylabel="")
    Plots.savefig(joinpath(PLOT_DIR, "Legendre_coefs_" * "$m" * "_pcfs.pdf"))
end
