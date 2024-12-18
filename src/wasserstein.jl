# Function to repeat values based on corresponding repeat counts
rep(values, repeats) = reduce(vcat, map((v, r) -> repeat([v], r), values, repeats))

# Function to calculate the Wasserstein distance for 1D distributions
function wasserstein1d(a::Vector, b::Vector; p=1, wa=nothing, wb=nothing)
    m = length(a)
    n = length(b)
    @assert m > 0 && n > 0 "Both input vectors must have positive length."

    if m == n && wa === nothing && wb === nothing
        return mean(abs.(sort(b) - sort(a)).^p)^(1/p)
    end

    @assert wa === nothing || length(wa) == m "Weights wa must be of length m."
    @assert wb === nothing || length(wb) == n "Weights wb must be of length n."

    if wa === nothing
        wa = ones(m)
    else
        valid_indices = wa .> 0
        wa = wa[valid_indices]
        a = a[valid_indices]
        m = length(a)
    end

    if wb === nothing
        wb = ones(n)
    else
        valid_indices = wb .> 0
        wb = wb[valid_indices]
        b = b[valid_indices]
        n = length(b)
    end

    orda = sortperm(a)
    ordb = sortperm(b)
    a = a[orda]
    b = b[ordb]
    wa = wa[orda]
    wb = wb[ordb]

    ua = wa / sum(wa)
    ub = wb / sum(wb)
    ua = ua[1:end-1]
    ub = ub[1:end-1]
    cua = cumsum(ua)
    cub = cumsum(ub)

    arep = fit(Histogram, cub, vcat(-Inf, cua, Inf), closed=:left).weights .+ 1
    brep = fit(Histogram, cua, vcat(-Inf, cub, Inf), closed=:left).weights .+ 1

    aa = rep(a, arep)
    bb = rep(b, brep)

    uu = sort(vcat(cua, cub))
    uu0 = vcat(0, uu)
    uu1 = vcat(uu, 1)

    areap = sum((uu1 - uu0) .* abs.(bb - aa).^p)^(1/p)
    return areap
end

# Function to get a random projection vector
function get_random_projection(n::Int)
    v = randn(n)
    return v / norm(v)
end

# Function to project a sample onto a given axis
function project_sample(x::Vector{V}, axis::Vector{W}) where {V<:Real, W<:Real}
    return dot(x, axis)
end

# Function to project samples onto a given axis
function project_samples(X, axis::Vector{W}) where {W<:Real}
    return project_sample.(X, Ref(axis))
end

# Function to calculate the sliced Wasserstein distance
function get_sliced_wasserstein_distance(sample1::DensitySampleVector, sample2::DensitySampleVector; L=1000, p=1, N=0, parallel=true)
    sample1, sample2 = prepare_twosample_dsv(sample1, sample2, N=N)

    dim = length(unshaped.(sample1).v[1])
    ds = []
    if !parallel || Threads.nthreads() < 2
        for i in 1:L
            z = get_random_projection(dim)
            W = project_samples(Array{Array{Float64}}(unshaped.(sample1).v), z)
            V = project_samples(Array{Array{Float64}}(unshaped.(sample2).v), z)
            if is_weighted(sample1) || is_weighted(sample2)
                push!(ds, wasserstein1d(V, W, wa=sample1.weight, wb=sample2.weight))
            else
                push!(ds, wasserstein1d(V, W))
            end
        end
        return norm(ds, p) / (length(ds)^(1/p))
    else
        ds = Vector{Float64}(undef, L)
        Threads.@threads for i in 1:L
            z = get_random_projection(dim)
            W = project_samples(Array{Array{Float64}}(unshaped.(sample1).v), z)
            V = project_samples(Array{Array{Float64}}(unshaped.(sample2).v), z)
            if is_weighted(sample1) || is_weighted(sample2)
                ds[i] = wasserstein1d(V, W, wa=sample1.weight, wb=sample2.weight)
            else
                ds[i] = wasserstein1d(V, W)
            end
        end
    end
    return norm(ds, p) / (length(ds)^(1/p))
end

# Helper function to calculate the sliced Wasserstein distance in parallel
function i_get_sliced_wasserstein_distance(dim, sample1, sample2)
    z = get_random_projection(dim)
    W = project_samples(Array{Array{Float64}}(unshaped.(sample1).v), z)
    V = project_samples(Array{Array{Float64}}(unshaped.(sample2).v), z)
    if is_weighted(sample1) || is_weighted(sample2)
        return wasserstein1d(V, W, wa=sample1.weight, wb=sample2.weight)
    else
        return wasserstein1d(V, W)
    end
end

# Function to calculate the sliced Wasserstein distance between two distributions
"""
    get_sliced_wasserstein_distance(dist1, dist2; d=50, L=1000, p=1, N=100_000)

Compute the sliced Wasserstein distance between two probability distributions.

## Arguments
- `dist1`: The first probability distribution.
- `dist2`: The second probability distribution.
- `d`: The dimension of the samples. Default is 50.
- `L`: The number of random projections. Default is 1000.
- `p`: The norm order. Default is 1.
- `N`: The number of samples. Default is 100_000.

## Returns
The sliced Wasserstein distance between `dist1` and `dist2`.
"""

function get_sliced_wasserstein_distance(dist1, dist2; d=50, L=1000, p=1, N=100_000)
    samples_A = [rand(dist1, d) for i in 1:N]
    samples_B = [rand(dist2, d) for i in 1:N]

    ds = []

    for i in 1:L
        v = get_random_projection(d)
        Y_A = project_samples(samples_A, v)
        Y_B = project_samples(samples_B, v)
        distance = wasserstein1d(Y_A, Y_B)
        push!(ds, distance)
    end

    return norm(ds, p) / (length(ds)^(1/p))
end
