function get_mmd(s1::DensitySampleVector, s2::DensitySampleVector; g=0, N=0)
    s1, s2 = prepare_twosample_dsv(s1, s2,N=N)
    x = Matrix{Float64}(hcat(unshaped.(s1).v...))
    y = Matrix{Float64}(hcat(unshaped.(s2).v...))
    get_mmd(x, y, g=g)
end

function get_mmd(x::AbstractArray, y::AbstractArray; g=0)
    if g != 0
        return mmd(GaussianKernel(g), x, y)
    end
    return mmd(GaussianKernel(compute_bandwidth(x, y)), x, y)
end

function compute_bandwidth(x, y)
    combined_data = [x y]'

    # Compute pairwise squared distances
    D = pairwise(SqEuclidean(), combined_data; dims=1)

    # Extract upper triangle without diagonal
    upper_triangle = D[triu(ones(Bool, size(D, 1), size(D, 2)), 1)]

    γ = median(upper_triangle)
    return γ
end

function compute_bandwidth(s1::DensitySampleVector, s2::DensitySampleVector)
    s1, s2 = prepare_twosample_dsv(s1, s2)
    x = Matrix{Float64}(hcat(unshaped.(s1).v...))
    y = Matrix{Float64}(hcat(unshaped.(s2).v...))
    compute_bandwidth(x, y)
end