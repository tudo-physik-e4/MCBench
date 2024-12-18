# Useful functions to manipulate samples

# Function to create a DensitySampleVector (DSV) from a matrix of samples
function make_dsv(samples::Matrix{<:Real}, logdensities::AbstractVector{<:Real}; weights::AbstractVector{<:Real} = ones(length(logdensities)))
    return DensitySampleVector([x = samples[:, i] for i in 1:size(samples, 2)], logdensities, weight=weights)
end
function make_dsv(samples::Matrix{<:Real}, logdensities::AbstractVector{<:Real}, infos; weights::AbstractVector{<:Real} = ones(length(logdensities)))
    return DensitySampleVector([x = samples[:, i] for i in 1:size(samples, 2)], logdensities, weight=weights,info=infos)
end

# Function to create a DSV from a vector of samples
function make_dsv(samples::Vector{<:Real}, logdensities::AbstractVector{<:Real}; weights::AbstractVector{<:Real} = ones(length(logdensities)))
    samples = collect(hcat(samples...))
    return DensitySampleVector([x = samples[:, i] for i in 1:size(samples, 2)], logdensities, weight=weights)
end

# Function to create a DSV from a vector of vector samples
function make_dsv(samples::Vector{Vector{T}}, logdensities::AbstractVector{<:Real}; weights::AbstractVector{<:Real} = ones(length(logdensities))) where {T <: Real}
    samples = collect(hcat(samples...))
    return DensitySampleVector([x = samples[:, i] for i in 1:size(samples, 2)], logdensities, weight=weights)
end

#
function make_dsv(samples::Vector{Vector{T}}; logdensities::AbstractVector{<:Real} = ones(length(samples)), weights::AbstractVector{<:Real} = ones(length(logdensities))) where {T <: Real}
    make_dsv(samples, logdensities, weights=weights)
end
function make_dsv(samples::Vector{<:Real}; logdensities::AbstractVector{<:Real} = ones(length(samples)), weights::AbstractVector{<:Real} = ones(length(logdensities)))
    make_dsv(samples, logdensities, weights=weights)
end

# Function to get the effective sample size of a DSV
function get_effective_sample_size(dsv::DensitySampleVector)
    return BAT.bat_eff_sample_size(dsv, BAT.KishESS()).result
end

# Overloaded function to get the effective sample size of a DSV with an integer parameter
function get_effective_sample_size(dsv::DensitySampleVector, s::Int)
    return BAT.bat_eff_sample_size(dsv, BAT.KishESS()).result
end

# Overloaded function to get the effective sample size of a DSV with a sampling algorithm
function get_effective_sample_size(dsv::DensitySampleVector, s::SamplingAlgorithm)
    if isa(s, MCMCSamplingAlgorithm)
        ess = BAT.bat_eff_sample_size(dsv, BAT.KishESS())
        #ess = BAT.bat_eff_sample_size(dsv, BAT.EffSampleSizeFromAC())
    else
        ess = BAT.bat_eff_sample_size(dsv, BAT.KishESS())
    end
    if isa(ess.result, Vector{<:Real}) 
        return minimum(ess.result)
    else 
        return minimum(collect(ess.result...))
    end
    return ess
end

# Function to condense a DSV by removing repeated samples
function condense_dsv(dsv::DensitySampleVector)
    @assert dsv.weight == ones(length(dsv.weight)) "Already condensed"
    v = unshaped.(dsv.v)
    logval = dsv.logd
    info = dsv.info
    aux = dsv.aux
    idxs, weight = BAT.repetition_to_weights(v)
    #return DensitySampleVector(ArrayOfSimilarArrays(v[idxs]), logval[idxs], weight=weight, info=info[idxs], aux=aux[idxs])
    return DensitySampleVector((v[idxs]), logval[idxs], weight=weight, info=info[idxs], aux=aux[idxs])
end

# Function to resample a DSV to its effective sample size using a sampling algorithm
function resample_dsv_to_ess(dsv::DensitySampleVector, s::SamplingAlgorithm)
    N_eff = Int(floor(get_effective_sample_size(dsv, s)))
    return bat_sample(dsv, RandResampling(nsamples=N_eff)).result
end

# Overloaded function to resample a DSV to its effective sample size
function resample_dsv_to_ess(dsv::DensitySampleVector)
    N_eff = Int(floor(get_effective_sample_size(dsv)))
    return bat_sample(dsv, RandResampling(nsamples=N_eff)).result
end

# Function to resample a DSV to a specified number of samples
function resample_dsv(dsv::DensitySampleVector, n::Int)
    return bat_sample(dsv, RandResampling(nsamples=n)).result
end

# Function to check if a DSV is weighted
function is_weighted(dsv::DensitySampleVector)
    return any(dsv.weight .!= 1)
end

# Function to prepare two sample DSVs for comparison, ensuring they have the same effective sample size
function prepare_twosample_dsv(dsv1::DensitySampleVector, dsv2::DensitySampleVector; N=0)
    if length(dsv1) != length(dsv2)
        N = Int(floor(min(get_effective_sample_size(dsv1), get_effective_sample_size(dsv2))))
        println("WARNING: Samples must have the same length, setting n = Neff(min(sample1, sample2)) = $N")
    elseif N != 0 && (N > get_effective_sample_size(dsv1) || N > get_effective_sample_size(dsv2))
        N = Int(floor(min(get_effective_sample_size(dsv1), get_effective_sample_size(dsv2))))
        println("WARNING: N set too high, setting n = Neff(min(sample1, sample2)) = $N")
    end

    if N != 0
        @assert N <= get_effective_sample_size(dsv1) && N <= get_effective_sample_size(dsv2)
        dsv1 = resample_dsv(dsv1, N)
        dsv2 = resample_dsv(dsv2, N)
    end
    return dsv1, dsv2
end

# include("/ceph/groups/e4/users/slacagnina/MCBench/scripts/BAT_utils.jl")

# # Function to calculate the sliced Wasserstein distance for chains of samples
# function sliced_wasserstein_distance_chains(samples; N=0, L=1000)
#     cids = _unique_chainids(samples)
#     fss = [filter_by_chainID(samples, i) for i in cids]
#     els = vcat([(i, j) for i in 1:length(cids), j in 1:length(cids)]...)
#     perms = [x for x in els if x[1] >= x[2]]

#     return Folds.collect(get_sliced_wasserstein_distance(fss[x[1]], fss[x[2]], N=N, L=L) for x in perms)
# end
