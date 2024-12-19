# Abstract types for metrics
abstract type TestMetric end
abstract type TwoSampleMetric <: TestMetric end

# Kolmogorov-Smirnov (KS) metric struct
struct ks{
    V<:Real,
} <: TwoSampleMetric
    val::V
end

# Calculation of metric for a test case and two density sample vectors
function calc_metric(t::AT, s1::DensitySampleVector, s2::DensitySampleVector, m::TM) where {TM <: TestMetric, AT <: AbstractTestcase}
    ms1, ms2 = marginalize_with_binning(s1, s2)
    return [calc_metric(ms1[i], ms2[i], m) for i in 1:length(ms1)]
end

# Marginal mean metric struct and constructors
struct marginal_mean{
    V<:Real,
    A<:Any,
} <: TestMetric
    val::V
    info::A
end

function marginal_mean() marginal_mean(0.0, "Mean") end
function marginal_mean(val::V) where {V<:Real} marginal_mean(val, "Mean") end

# Calculation of marginal mean metric
function calc_metric(t::AT, s::DensitySampleVector, m::marginal_mean) where {AT <: AbstractTestcase}
    vals = isa(BAT.mean(s), NamedTuple) ? Vector{Float64}(collect(values(BAT.mean(s)))[1]) : BAT.mean(s)
    return [marginal_mean(i) for i in vals]
end

# Marginal variance metric struct and constructors
struct marginal_variance{
    V<:Real,
    A<:Any,
} <: TestMetric
    val::V
    info::A
end

function marginal_variance() marginal_variance(0.0, "Variance") end
function marginal_variance(val::V) where {V<:Real} marginal_variance(val, "Variance") end

# Calculation of marginal variance metric
function calc_metric(t::AT, s::DensitySampleVector, m::marginal_variance) where {AT <: AbstractTestcase}
    vals = isa(BAT.var(s), NamedTuple) ? Vector{Float64}(collect(values(BAT.var(s)))[1]) : BAT.var(s)
    return [marginal_variance(i) for i in vals]
end

# Global mode metric struct and constructors
struct global_mode{
    V<:Real,
    A<:Any,
} <: TestMetric
    val::V
    info::A
end

function global_mode() global_mode(0.0, "Globalmode") end
function global_mode(val::V) where {V<:Real} global_mode(val, "Globalmode") end

# Calculation of global mode metric
function calc_metric(t::AT, s::DensitySampleVector, m::global_mode) where {AT <: AbstractTestcase}
    vals = isa(BAT.mode(s), NamedTuple) ? Vector{Float64}(collect(values(BAT.mode(s)))[1]) : BAT.mode(s)
    return [global_mode(i) for i in vals]
end

# Marginal mode metric struct and constructors
struct marginal_mode{
    V<:Real,
    A<:Any,
} <: TestMetric
    val::V
    info::A
end

function marginal_mode() marginal_mode(0.0, "Marginalmode") end
function marginal_mode(val::V) where {V<:Real} marginal_mode(val, "Marginalmode") end

# Calculation of marginal mode metric
function calc_metric(t::AT, s::DensitySampleVector, m::marginal_mode) where {AT <: AbstractTestcase}
    vals = isa(BAT.bat_marginalmode(s).result, NamedTuple) ? Vector{Float64}(collect(values(BAT.bat_marginalmode(s).result))[1]) : Vector{Float64}(BAT.bat_marginalmode(s).result)
    return [marginal_mode(i) for i in vals]
end

# Marginal skewness metric struct and constructors
struct marginal_skewness{
    V<:Real,
    A<:Any,
} <: TestMetric
    val::V
    info::A
end

function marginal_skewness() marginal_skewness(0.0, "Skewness") end
function marginal_skewness(val::V) where {V<:Real} marginal_skewness(val, "Skewness") end

# Calculation of marginal skewness metric
function calc_metric(t::AT, s::DensitySampleVector, m::marginal_skewness) where {AT <: AbstractTestcase}
    return [marginal_skewness(skewness([i[idim] for i in unshaped.(s).v], FrequencyWeights(s.weight))) for idim in 1:t.dim]
end

# Marginal kurtosis metric struct and constructors
struct marginal_kurtosis{
    V<:Real,
    A<:Any,
} <: TestMetric
    val::V
    info::A
end

function marginal_kurtosis() marginal_kurtosis(0.0, "Kurtosis") end
function marginal_kurtosis(val::V) where {V<:Real} marginal_kurtosis(val, "Kurtosis") end

# Calculation of marginal kurtosis metric
function calc_metric(t::AT, s::DensitySampleVector, m::marginal_kurtosis) where {AT <: AbstractTestcase}
    return [marginal_kurtosis(kurtosis([i[idim] for i in unshaped.(s).v], FrequencyWeights(s.weight))) for idim in 1:t.dim]
end


# Wasserstein 1D
struct wasserstein_1d{
    V<:Real,
    A<:Any,
} <: TwoSampleMetric
    val::V
    info::A
end

function wasserstein_1d() wasserstein_1d(0.0, "Wasserstein") end
function wasserstein_1d(val::V) where {V<:Real} wasserstein_1d(val, "Wasserstein") end

# Calculation of marginal kurtosis metric
function calc_metric(t::AT, s1::DensitySampleVector, s2::DensitySampleVector, m::wasserstein_1d) where {AT <: AbstractTestcase}
    us1 = hcat(unshaped.(s1).v...)
    us2 = hcat(unshaped.(s2).v...)
    return [wasserstein_1d(wasserstein1d(Vector{Float64}(us1[idim,:]), Vector{Float64}(us2[idim,:]), wa=s1.weight, wb=s2.weight)) for idim in 1:t.dim]
end

# Sliced Wasserstein distance metric struct and constructors
struct sliced_wasserstein_distance{
    V<:Real,
    I<:Int,
    A<:Any,
    P<:Any
} <: TwoSampleMetric
    val::V
    N::I
    info::A
    proc::P
end

# Helper function for initializing Sliced Wasserstein distance metric
function wd_pint()
    if Threads.nthreads() < 5
        return ""
    else
        return ""
    end
end

function sliced_wasserstein_distance() sliced_wasserstein_distance(0.0, 10^5, "SlicedWasserstein", wd_pint()) end
function sliced_wasserstein_distance(val::V) where {V<:Real} sliced_wasserstein_distance(val, 10^5, "SlicedWasserstein", "") end
function sliced_wasserstein_distance(val::V, N::Int) where {V<:Real} sliced_wasserstein_distance(val, N, "SlicedWasserstein", wd_pint()) end

# Calculation of Sliced Wasserstein distance metric
function calc_metric(t::AT, s1::DensitySampleVector, s2::DensitySampleVector, m::sliced_wasserstein_distance) where {AT<:AbstractTestcase}
    return [(sliced_wasserstein_distance(get_sliced_wasserstein_distance(s1, s2, N=m.N)))]
end


# Maximum Mean Discrepancy (MMD) metric struct and constructors
struct maximum_mean_discrepancy{
    V<:Real,
    I<:Int,
    A<:Any,
    P<:Any
} <: TwoSampleMetric
    val::V
    N::I
    info::A
    proc::P
end


function maximum_mean_discrepancy() maximum_mean_discrepancy(0.0, 10^4, "MaximumMeanDiscrepancy", wd_pint()) end
function maximum_mean_discrepancy(val::V) where {V<:Real} maximum_mean_discrepancy(val, 10^4, "MaximumMeanDiscrepancy", "") end
function maximum_mean_discrepancy(val::V, N::Int) where {V<:Real} maximum_mean_discrepancy(val, N, "MaximumMeanDiscrepancy", wd_pint()) end

# Calculation of Sliced Wasserstein distance metric
function calc_metric(t::AT, s1::DensitySampleVector, s2::DensitySampleVector, m::maximum_mean_discrepancy) where {AT<:AbstractTestcase}
    return [(maximum_mean_discrepancy(get_mmd(s1, s2, N=m.N)))]
end

#Chi^2 test metric struct and constructors
struct chi_squared{
    V<:Real,
    A<:Any,
} <: TwoSampleMetric
    val::V
    info::A
end


function chi_squared() chi_squared(0.0, "Chi-Squared") end
function chi_squared(val::V) where {V<:Real} chi_squared(val,"Chi-Squared") end

# Calculation of Sliced Wasserstein distance metric
function calc_metric(t::AT, s1::DensitySampleVector, s2::DensitySampleVector, m::chi_squared) where {AT<:AbstractTestcase}
    return chisq_test(s1, s2)
end

function chisq_test(s1::DensitySampleVector, s2::DensitySampleVector)
    #s1 = iid
    s1, s2 = prepare_twosample_dsv(s1, s2)
    x = Matrix{Float64}(hcat(unshaped.(s1).v...))
    y = Matrix{Float64}(hcat(unshaped.(s2).v...))
    mx = [mean(i) for i in eachrow(x)]
    sx = [var(i) for i in eachrow(x)]
    [chi_squared(sum(y[i,:] .- mx[i])/sx[i]) for i in 1:size(x,1)]
end
