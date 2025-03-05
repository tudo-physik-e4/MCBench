"""
    abstract type TestMetric

An abstract type for any test metrics. Each metric is its own struct and must implement the `calc_metric` function.
Metrics are designed to be both structs that hold the metric value and determine which version of `calc_metric` to dispatch.

    abstract type TwoSampleMetric

An abstract type for any two-sample test metrics. Such as Wasserstein, MMD, etc.
"""
abstract type TestMetric end
abstract type TwoSampleMetric <: TestMetric end

"""
    calc_metric(t::AbstractTestcase, 
    s::DensitySampleVector, 
    m::TestMetric)

    Calculate the metric `m` for the test case `t` and the samples `s`.
    This function is used to calculate the metric value for a given test case and samples for any metric type.
    `TwoSampleMetric` metrics can also be used, however the `calc_metric` then samples a second set of IID samples to compare against. 
"""
"""
    struct marginal_mean{V<:Real,A<:Any} <: TestMetric

    # Fields
    - `val::V`: The value of the metric.
    - `info::A`: Information about the metric.

    # Constructors
    - `marginal_mean(; fields...)`
    - `marginal_mean(val::Real)` : Creates a `marginal_mean` metric with the given metric value.

    A struct for the mean value for each dimension of the samples as a metric.
    
"""
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

"""
    struct marginal_variance{V<:Real,A<:Any} <: TestMetric

    # Fields
    - `val::V`: The value of the metric.
    - `info::A`: Information about the metric.

    # Constructors
    - `marginal_variance(; fields...)`
    - `marginal_variance(val::Real)` : Creates a `marginal_variance` metric with the given metric value.

    A struct for the variance for each dimension of the samples as a metric.
    
"""
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

"""
    struct global_mode{V<:Real,A<:Any} <: TestMetric

    # Fields
    - `val::V`: The value of the metric.
    - `info::A`: Information about the metric.

    # Constructors
    - `global_mode(; fields...)`
    - `global_mode(val::Real)` : Creates a `global_mode` metric with the given metric value.

    A struct for the global mode of the samples as a metric. This will return the mode as a vector of the mode for each dimension.
    
"""
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

"""
    struct marginal_mode{V<:Real,A<:Any} <: TestMetric

    # Fields
    - `val::V`: The value of the metric.
    - `info::A`: Information about the metric.

    # Constructors
    - `marginal_mode(; fields...)`
    - `marginal_mode(val::Real)` : Creates a `marginal_mode` metric with the given metric value.

    A struct for the marginal_mode of the samples as a metric. Marginal modes are calculated for each dimension of the samples and returned as a vector of `marginal_mode` values.
    
"""
struct marginal_mode{
    V<:Real,
    A<:Any,
} <: TestMetric
    val::V
    info::A
end

function marginal_mode() marginal_mode(0.0, "Marginalmode") end
function marginal_mode(val::V) where {V<:Real} marginal_mode(val, "Marginalmode") end

function calc_metric(t::AT, s::DensitySampleVector, m::marginal_mode) where {AT <: AbstractTestcase}
    vals = isa(BAT.bat_marginalmode(s).result, NamedTuple) ? Vector{Float64}(collect(values(BAT.bat_marginalmode(s).result))[1]) : Vector{Float64}(BAT.bat_marginalmode(s).result)
    return [marginal_mode(i) for i in vals]
end

"""
    struct marginal_skewness{V<:Real,A<:Any} <: TestMetric

    # Fields
    - `val::V`: The value of the metric.
    - `info::A`: Information about the metric.

    # Constructors
    - `marginal_skewness(; fields...)`
    - `marginal_skewness(val::Real)` : Creates a `marginal_skewness` metric with the given metric value.

    A struct for the marginal_skewness of the samples as a metric. Marginal skewness are calculated for each dimension of the samples and returned as a vector of `marginal_skewness` values.
    
"""
struct marginal_skewness{
    V<:Real,
    A<:Any,
} <: TestMetric
    val::V
    info::A
end

function marginal_skewness() marginal_skewness(0.0, "Skewness") end
function marginal_skewness(val::V) where {V<:Real} marginal_skewness(val, "Skewness") end

function calc_metric(t::AT, s::DensitySampleVector, m::marginal_skewness) where {AT <: AbstractTestcase}
    return [marginal_skewness(skewness([i[idim] for i in unshaped.(s).v], FrequencyWeights(s.weight))) for idim in 1:t.dim]
end

"""
    struct marginal_kurtosis{V<:Real,A<:Any} <: TestMetric

    # Fields
    - `val::V`: The value of the metric.
    - `info::A`: Information about the metric.

    # Constructors
    - `marginal_kurtosis(; fields...)`
    - `marginal_kurtosis(val::Real)` : Creates a `marginal_kurtosis` metric with the given metric value.

    A struct for the marginal_kurtosis of the samples as a metric. Marginal kurtosis are calculated for each dimension of the samples and returned as a vector of `marginal_kurtosis` values.
    
"""
struct marginal_kurtosis{
    V<:Real,
    A<:Any,
} <: TestMetric
    val::V
    info::A
end

function marginal_kurtosis() marginal_kurtosis(0.0, "Kurtosis") end
function marginal_kurtosis(val::V) where {V<:Real} marginal_kurtosis(val, "Kurtosis") end

function calc_metric(t::AT, s::DensitySampleVector, m::marginal_kurtosis) where {AT <: AbstractTestcase}
    return [marginal_kurtosis(kurtosis([i[idim] for i in unshaped.(s).v], FrequencyWeights(s.weight))) for idim in 1:t.dim]
end


"""
    struct wasserstein_1d{V<:Real,A<:Any} <: TwoSampleMetric

    # Fields
    - `val::V`: The value of the metric.
    - `info::A`: Information about the metric.

    # Constructors
    - `wasserstein_1d(; fields...)`
    - `wasserstein_1d(val::Real)` : Creates a `wasserstein_1d` metric with the given metric value.

    A struct for the wasserstein_1d of the samples as a metric. The Wasserstein distance is calculated for each dimension separately and returned as a vector of `wasserstein_1d` values.
    
"""
struct wasserstein_1d{
    V<:Real,
    A<:Any,
} <: TwoSampleMetric
    val::V
    info::A
end

function wasserstein_1d() wasserstein_1d(0.0, "Wasserstein") end
function wasserstein_1d(val::V) where {V<:Real} wasserstein_1d(val, "Wasserstein") end

function calc_metric(t::AT, s1::DensitySampleVector, s2::DensitySampleVector, m::wasserstein_1d) where {AT <: AbstractTestcase}
    us1 = hcat(unshaped.(s1).v...)
    us2 = hcat(unshaped.(s2).v...)
    return [wasserstein_1d(wasserstein1d(Vector{Float64}(us1[idim,:]), Vector{Float64}(us2[idim,:]), wa=s1.weight, wb=s2.weight)) for idim in 1:t.dim]
end

"""
    struct sliced_wasserstein_distance{V<:Real,A<:Any} <: TwoSampleMetric

    # Fields
    - `val::V`: The value of the metric.
    - `N::I`: The number of points used to calculate the sliced Wasserstein distance.
    - `info::A`: Information about the metric.
    - `proc::P`: The processor information for the metric.

    # Constructors
    - `sliced_wasserstein_distance(; fields...)`
    - `sliced_wasserstein_distance(val::Real)` : Creates a `sliced_wasserstein_distance` metric with the given metric value.

    A struct for the sliced_wasserstein_distance of the samples as a metric. 
    The sliced Wasserstein distance is calculated by projecting the samples onto a random direction and calculating the Wasserstein distance in that direction.
    
"""
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

function calc_metric(t::AT, s1::DensitySampleVector, s2::DensitySampleVector, m::sliced_wasserstein_distance) where {AT<:AbstractTestcase}
    return [(sliced_wasserstein_distance(get_sliced_wasserstein_distance(s1, s2, N=m.N)))]
end

"""
    struct maximum_mean_discrepancy{V<:Real,A<:Any} <: TwoSampleMetric

    # Fields
    - `val::V`: The value of the metric.
    - `N::I`: The number of points used to calculate the maximum mean discrepancy.
    - `info::A`: Information about the metric.
    - `proc::P`: The processor information for the metric.

    # Constructors
    - `maximum_mean_discrepancy(; fields...)`
    - `maximum_mean_discrepancy(val::Real)` : Creates a `maximum_mean_discrepancy` metric with the given metric value.
    - `maximum_mean_discrepancy(val::Real, N::Int)` : Creates a `maximum_mean_discrepancy` metric with the given metric value and number of points.

    A struct for the maximum_mean_discrepancy of the samples as a metric. 
    Per default the number of points used to calculate the maximum mean discrepancy is 10^4.
    
"""
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

function calc_metric(t::AT, s1::DensitySampleVector, s2::DensitySampleVector, m::maximum_mean_discrepancy) where {AT<:AbstractTestcase}
    return [(maximum_mean_discrepancy(get_mmd(s1, s2, N=m.N)))]
end

"""
    struct chi_squared{V<:Real,A<:Any} <: TwoSampleMetric

    # Fields
    - `val::V`: The value of the metric.
    - `info::A`: Information about the metric.

    # Constructors
    - `chi_squared(; fields...)`
    - `chi_squared(val::Real)` : Creates a `chi_squared` metric with the given metric value.

    A struct for the chi_squared of the samples as a metric. 
    The chi-squared is calculated for each dimension separately using an unbinned chi-squared. 
    
"""
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
