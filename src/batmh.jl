abstract type MCMCSamplingAlgorithm <: SamplingAlgorithm end
# Struct for BAT Metropolis-Hastings (MH) sampler
struct BATMH{
    SA <: BAT.AbstractSamplingAlgorithm,
    A <: Any,
} <: MCMCSamplingAlgorithm
    sampler::SA
    info::A
end

# Constructor for posterior measure
function buildBATPosterior(testcase::Testcases)
    # Build the BAT posterior using the distribution and bounds from the test case
    BAT.PosteriorMeasure(testcase, testcase.bounds)
end

# Log density of the test case for BAT.jl usage
function DensityInterface.logdensityof(d::Testcases, x)
    logpdf(d.f, x.x)[1]
end

# Constructor for BAT Metropolis-Hastings (MH) sampler
function BATMH(; n_steps=10^5, nchains=10)
    BATMH(MCMCSampling(mcalg = MetropolisHastings(), nsteps=n_steps, nchains=nchains), "BAT-MH")
end

# Sample function for test cases using BAT Metropolis-Hastings (MH) sampler
function sample(t::Testcases, s::BATMH; n_steps::Int=10^5, nchains::Int=10)
    sampler = (n_steps != 10^5 || nchains != 10) ? BATMH(n_steps=n_steps, nchains=nchains).sampler : s.sampler
    bat_sample(buildBATPosterior(t), sampler).result
end