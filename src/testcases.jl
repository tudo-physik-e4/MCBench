# Abstract type for test cases
abstract type AbstractTestcase end
abstract type Target end

# Struct for test cases with specified types for distributions, bounds, dimensions, and additional information
struct Testcases{
    D<:Union{Distribution,Target},
    B<:NamedTupleDist,
    A<:Any,
    N<:Int,
} <: AbstractTestcase
    f::D          # Distribution
    bounds::B     # Bounds as NamedTupleDist
    dim::N        # Dimension
    info::A       # Additional info
end

# Constructor for test cases without bounds
function Testcases(f::D, dim::N, info::A) where {D <: Union{Distribution,Target}, N <: Int, A <: Any}
    bounds = NamedTupleDist(x = fill(-10..10, dim))
    Testcases(f, bounds, dim, info)
end

# Indicate that the given type is a density
@inline DensityInterface.DensityKind(::Testcases) = IsDensity()

# IID sampling from the distribution of the test case
# Wrapping the samples into DensitySampleVector
function sample(t::Testcases; n_steps=10^5)
    s = rand(t.f, n_steps)  # Generate random samples from the distribution
    lgd = logpdf(t.f, s)  # Compute the log densities of the samples
    make_dsv(s, lgd)
end

function sample(t::Testcases, n::Int)
    sample(t, n_steps=n)
end

function sample(t::Testcases, s::IIDSamplingAlgorithm; n_steps=10^5)
    n = n_steps
    s = rand(t.f, n)  # Generate random samples from the distribution
    lgd = logpdf(t.f, s)  # Compute the log densities of the samples
    if typeof(s) == Vector{Float64}
        s = reshape(s, (1, length(s)))  # Reshape the samples if needed
    end
    DensitySampleVector([x = s[:, i] for i in 1:size(s, 2)], lgd)
end

#Sampling using a file-based sampler
function sample(t::AT, s::FBA; n_steps=10^4) where {AT <: AbstractTestcase, FBA <: AbstractFileBasedSampler}
    sample(s,t=t,n_steps=n_steps)
end

function sample(s::FileBasedSampler; t=0, n_steps=10^4)
    samples = [parse.(Float64, split(read_sample!(s), ",")) for i in 1:n_steps]
    samples = hcat(samples...)
    lgd = isa(t,Testcases) ? logpdf(t.f, samples) : ones(size(samples,2))
    if typeof(samples) == Vector{Float64} 
        samples = reshape(samples, (1, length(samples)))  # Reshape the samples if needed
    end
    DensitySampleVector([x = samples[:, i] for i in 1:size(samples, 2)], lgd)
end

function sample(fbs::CsvBasedSampler, n::Int)
    samples = Vector{Float64}[]
    for i in 1:n
        push!(samples, [parse(Float64,i) for i in split(read_sample!(fbs), ",")[fbs.mask]])
    end
    return make_dsv(samples)
end

function sample(s::CsvBasedSampler; t=0, n_steps=10^4)
    samples = [parse.(Float64, split(read_sample!(s), ",")[s.mask]) for i in 1:n_steps]
    samples = hcat(samples...)
    lgd = isa(t,Testcases) ? logpdf(t.f, samples) : ones(size(samples,2))
    if typeof(samples) == Vector{Float64} 
        samples = reshape(samples, (1, length(samples)))  # Reshape the samples if needed
    end
    DensitySampleVector([x = samples[:, i] for i in 1:size(samples, 2)], lgd)
end

function sample(t::AT, s::CsvBasedSampler; n_steps=10^4) where {AT <: AbstractTestcase}
    sample(s,t=t,n_steps=n_steps)
end

function sample(s::DsvSampler; t=0, n_steps=10^4)
    if(n_steps > s.neff[s.current_dsv_index])
        println("WARNING: Number of steps is greater than the number of effective samples. Resampling to the number of effective samples.")
        n_steps = Int(floor(s.neff[s.current_dsv_index]))
    end
    resample_dsv(s.dsvs[s.current_dsv_index],n_steps)
end


# Struct for test cases using DSV to use as IID for twosampleteststatics from precalculated samples
struct DsvTestcase{
    DS<:DsvSampler,
    A<:Any,
    N<:Int,
} <: AbstractTestcase
    sampler::DS   # Sampler
    dim::N        # Dimension
    info::A       # Additional info
end

function DsvTestcase(s::DS; n=0, info="DsvTestcase") where {DS <: DsvSampler}
    n = length(s.dsvs[1].v[1])
    DsvTestcase(s, n, info)
end

function sample(t::DsvTestcase, n::Int)
    sample(t.sampler, t=t, n_steps=n)
end
function sample(t::DsvTestcase; n_steps=10^5)
    sample(t.sampler, t=t, n_steps=n_steps)
end
function sample(t::DsvTestcase, s::SA; n_steps=10^5) where {SA <: SamplingAlgorithm}
    sample(t.sampler, t=t, n_steps=n_steps)
end
function sample(t::DsvTestcase, s::FBA; n_steps=10^5) where {FBA <: FileBasedSampler}
    sample(s, t=t, n_steps=n_steps)
end

