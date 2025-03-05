"""
    run_teststatistic(
    t::AbstractTestcase, 
    samples1::DensitySampleVector,
    samples2::DensitySampleVector,
    m::TwoSampleMetric,
    s::SamplingAlgorithm;)

    run_teststatistic(
    t::AbstractTestcase, 
    samples::DensitySampleVector,
    m::TwoSampleMetric, s::AnySampler;)

    run_teststatistic(
    t::AbstractTestcase, 
    samples::DensitySampleVector,
    m::TwoSampleMetric, s::Int;)

    run_teststatistic(
    t::AbstractTestcase, m::TwoSampleMetric,
    s::AnySampler; n_steps::Int=10^5)

    run_teststatistic(
    t::AbstractTestcase, m::TwoSampleMetric;
    n_steps::Int=10^5)

Functions to run a two-sample test statistic on a given testcase and metric.
When the samples are provided, the function calculates the metric value using the samples.
When the sampling algorithm is provided, the function samples using the algorithm and the IID samples and calculates the metric value.
When no samples and sampling algorithm are provided, the function samples using the testcase and calculates the metric value.

# Arguments
- `t::AbstractTestcase`: The test case, is a subtype of `AbstractTestcase`.
- `samples1::DensitySampleVector`: The first density sample vector.
- `samples2::DensitySampleVector`: The second density sample vector.
- `m::TwoSampleMetric`: The metric, is a subtype of `TwoSampleMetric`.
- `s::SamplingAlgorithm`: The sampling algorithm used for calculations.
- `n_steps::Int`: The number of steps to be generated.

# Returns
- `TwoSampleMetric`: The calculated metric value.

"""
function run_teststatistic(
    t::TSM, 
    samples1::DensitySampleVector,
    samples2::DensitySampleVector,
    m::TM,
    s::SamplingAlgorithm;
    ) where {TSM <: AbstractTestcase, TM <: TwoSampleMetric}
    mval = calc_metric(t,samples1,samples2,m)
end

function run_teststatistic(
    t::TSM, 
    samples::DensitySampleVector,
    m::TM,
    s::AS;
    ) where {TSM <: AbstractTestcase, TM <: TwoSampleMetric, AS <: AnySampler}
    iids1 = sample(t, n_steps=length(samples))
    mval = calc_metric(t,iids1,samples,m)
end

function run_teststatistic(
    t::TSM, 
    samples::DensitySampleVector,
    m::TM,
    s::Int;
    ) where {TSM <: AbstractTestcase, TM <: TwoSampleMetric}
    iids1 = sample(t, n_steps=length(samples))
    mval = calc_metric(t,iids1,samples,m)
end

function run_teststatistic(
    t::TSM, 
    m::TM,
    s::AS;
    n_steps::Int=10^5,
    ) where {TSM <: AbstractTestcase, TM <: TwoSampleMetric, AS <: AnySampler}
    iids1 = sample(t, n_steps=n_steps)
    samples = sample(t,s, n_steps=n_steps)
    mval = calc_metric(t,iids1,samples,m)
end

function run_teststatistic(
    t::TSM, 
    m::TM;
    n_steps::Int=10^5
    ) where {TSM <: AbstractTestcase, TM <: TwoSampleMetric}
    iids1 = sample(t, n_steps=n_steps)
    iids2 = sample(t, n_steps=n_steps)
    mval = calc_metric(t,iids1,iids2,m)
end