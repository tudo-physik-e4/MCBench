"""
    abstract type Target

Targets are used to represent the target distributions of test cases.
Make sure to define the `Base.rand`, `Base.length` and `Distributions.logpdf` methods for the target distribution.
See the posterior database example for more details.
"""

"""
    abstract type AbstractTestcase
An abstract type for any test cases.
"""
abstract type Target end
abstract type AbstractTestcase end

"""
    struct Testcases <: AbstractTestcase

A struct representing a test case to be used in the framework of MCBench.
Testcases must consisit of a distribution or Target that is sampleable and a set of bounds.

# Fields
- `f::D`: The distribution or Target of the test case.
- `bounds::B`: The bounds of the test case.
- `dim::N`: The dimension of the test case.
- `info::A`: Additional information about the test case.

# Constructors
- `Testcases(; fields...)`: Creates a test case with the given fields.
- `Testcases(f::D, bounds::B, dim::N, info::A)`: Creates a test case with the given distribution or target, bounds, dimension and additional information.
- `Testcases(f::D, dim::N, info::A)`: Creates a test case with the given distribution or target, dimension and additional information. Bounds are set to `[-10..10]` for each dimension.

"""
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


"""
    sample(t::Testcases; n_steps=10^5)::DensitySampleVector
    sample(t::Testcases, n::Int)::DensitySampleVector
    sample(t::Testcases, s::IIDSamplingAlgorithm; n_steps=10^5)::DensitySampleVector

IID sampling from the distribution of the test case `t` with `n_steps` or `n` samples.
When integrating a custom sampling algorithm, the `sample` method should be overloaded for the new sampling algorithm type for `s`.
Returns a density sample vector with the samples and the log densities.
"""
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

"""

    sample(t<:AbsatractTestcase, s<:AbstractFileBasedSampler; n_steps=10^5)::DensitySampleVector

Sampling from the distribution of the test case `t` using a file-based sampler `s` with `n_steps` samples.
The test case is not used in the sampling process, but it is used to calculate the log densities of the samples. 
Returns a density sample vector with the samples and the log densities.

"""
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


"""
    struct DsvTestcase <: AbstractTestcase

A struct representing a test case based on DensitySampleVectors.
This is meant to be used as IIDs when the samples are precalculated and stored in a file.
In that case samples should be read from the file converted to a DensitySampleVector and used as the test case.

# Fields
- `sampler::DS`: The sampler that generates the samples.
- `dim::N`: The dimension of the test case.
- `info::A`: Additional information about the test case.

# Constructors
- `DsvTestcase(s::DS, n::Int, info::A)`: Creates a test case with the given sampler, dimension and additional information.
- `DsvTestcase(s::DS, info::A)`: Creates a test case with the given sampler and additional information. The dimension is set to the dimension of the samples.

"""
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


"""

    sample(t::DsvTestcase; n_steps=10^5)::DensitySampleVector
    sample(t::DsvTestcase, s::SamplingAlgorithm; n_steps=10^5)::DensitySampleVector

The `sample` methods using `DsvTestcase` using the same logic as for `Testcases` but using the precalculated samples no matter the sampler used. 

"""
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