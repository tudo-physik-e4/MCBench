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