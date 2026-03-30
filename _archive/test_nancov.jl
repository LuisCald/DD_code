using Statistics, LinearAlgebra, Random

function nancov_simple(X; dims=2)
    if dims == 2
        n = size(X, 1)
        m = size(X, 2)
        Σ = zeros(n, n)
        for i in 1:n
            for j in i:n
                s_xy = 0.0; s_x = 0.0; s_y = 0.0; k = 0
                for t in 1:m
                    xi = X[i, t]; xj = X[j, t]
                    if !isnan(xi) && !isnan(xj)
                        s_x += xi; s_y += xj; s_xy += xi * xj; k += 1
                    end
                end
                if k > 1
                    Σ[i, j] = (s_xy - s_x * s_y / k) / (k - 1)
                    Σ[j, i] = Σ[i, j]
                end
            end
        end
        return Σ
    end
end

n = 2000
T = 40 * 200

println("Matrix size: $n × $T  (scaled down for timing)")
println("="^60)

for miss_pct in [0.0, 0.5, 0.8, 0.9]
    X = randn(n, T)
    
    if miss_pct > 0
        n_nan_rows = round(Int, n * miss_pct)
        nan_rows = randperm(n)[1:n_nan_rows]
        X[nan_rows, :] .= NaN
    end
    
    # Warmup
    if miss_pct == 0.0
        _ = nancov_simple(randn(10, 100), dims=2)
        _ = cov(randn(10, 100), dims=2)
    end
    # t0 = @elapsed Σ0 = nancov(X, dims=2)  # This will error if there are NaNs, but we just want the timing for the non-NaN case
    t1 = @elapsed Σ1 = nancov_simple(X, dims=2)
    
    X2 = copy(X)
    replace!(X2, NaN => 0.0)
    t2 = @elapsed Σ2 = cov(X2, dims=2)
    
    non_zero = diag(Σ1) .> 1e-12
    if any(non_zero)
        # sub0 = Σ0[non_zero, non_zero]
        sub1 = Σ1[non_zero, non_zero]
        sub2 = Σ2[non_zero, non_zero]
        # max_diff0 = maximum(abs.(sub0 .- sub2))
        max_diff1 = maximum(abs.(sub1 .- sub2))

    else
        max_diff0 = NaN
        max_diff1 = NaN
    end
    
    println("\n$(round(Int, miss_pct*100))% rows NaN:")
    # println(" real nancov:  $(round(t0, digits=3))s")
    println("  nancov:      $(round(t1, digits=3))s")
    println("  replace+cov: $(round(t2, digits=3))s")
    println("  speedup:     $(round(t1/t2, digits=1))x")
    println("  max |diff|:  $max_diff1")
    # println("  max |diff| (real nancov): $max_diff0")
end

# Extrapolate to full size
println("\n" * "="^60)
println("Extrapolation to full size (n=2746, T=40*999):")
println("  nancov scales as O(n²·T): ~$(round(2.5 * (2746/500)^2 * (40*999)/(40*200), digits=0))s")
println("  cov scales as O(n²·T) via BLAS: ~$(round(1.5 * (2746/500)^2 * (40*999)/(40*200), digits=1))s")
