abstract type AbstractKernel end
const MetricOrFun = Union{PreMetric,Function}

pairwisel2(x::Matrix, y::Matrix) = pairwise(SqEuclidean(), x, y, dims=2)
pairwisel2(x::AbstractMatrix) = pairwisel2(x,x)


function kernelsum(k::AbstractKernel, x::AbstractMatrix, y::AbstractMatrix, dist::MetricOrFun)
    sum(k(dist(x,y))) / (size(x,2) * size(y,2))
end

function kernelsum(k::AbstractKernel, x::AbstractMatrix{T}, dist::MetricOrFun) where T
    l = size(x,2)
    (sum(k(dist(x,x))) - l*k(T(0)))/(l^2 - l)
end

kernelsum(k::AbstractKernel, x::AbstractVector, dist::MetricOrFun) = zero(eltype(x))

struct GaussianKernel{T} <: AbstractKernel
	γ::T
end

(m::GaussianKernel)(x::Number) = exp(-m.γ * x)
(m::GaussianKernel)(x::AbstractArray) = exp.(-m.γ .* x)

struct MMD{K<:AbstractKernel,D<:MetricOrFun} <: PreMetric
    kernel::K
    dist::D
end

function (m::MMD)(x::AbstractArray, y::AbstractArray)
    xx = kernelsum(m.kernel, x, m.dist)
    yy = kernelsum(m.kernel, y, m.dist)
    xy = kernelsum(m.kernel, x, y, m.dist)
    xx + yy - 2xy
end

mmd(k::AbstractKernel, x::AbstractArray, y::AbstractArray, dist=pairwisel2) = MMD(k, dist)(x, y)
mmd(k::AbstractKernel, x::AbstractArray, y::AbstractArray, n::Int, dist=pairwisel2) = mmd(k, samplecolumns(x,n), samplecolumns(y,n), dist)

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