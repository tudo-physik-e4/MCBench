function sample_eight_schools_accept_reject(y,cov_matrix,maxval;trafo=false)
    mu = rand(Normal(0,5))
    tau = abs.(rand(Cauchy(0,5)))
    theta =  [rand(Normal(mu,tau)) for i in 1:8]
    lgd = logpdf(MvNormal(theta, cov_matrix), y)
    if lgd > log(rand(Uniform(0,maxval)))
        tau = trafo ? log.(tau) : tau
        return (v=(theta=theta,mu=mu,tau=tau),lgd=lgd)
    else
        return []
    end
end

function sample_eight_schools_accept_reject(n,data;trafo=false)
    y = Vector{Float64}(data["y"])
    sigma = data["sigma"]
    cov_matrix = Diagonal(sigma .^ 2)
    maxval = pdf(MvNormal(zeros(8), cov_matrix),zeros(8))
    result = Vector{NamedTuple}()
    logvals = Vector{Float64}()
    while length(result) < n
        res = sample_eight_schools_accept_reject(y,cov_matrix,maxval,trafo=trafo)
        if length(res) > 0
            push!(result,res.v)
            push!(logvals,res.lgd)
        end
    end
    return (result=result,logvals=logvals)
end


struct EightSchoolsAcceptReject{
    D<:Dict{String,Any}
} <: Target
    data::D
    trafo::Bool
end
function EightSchoolsAcceptReject(data)
    return EightSchoolsAcceptReject(data,false)
end
function Base.rand(d::EightSchoolsAcceptReject, n::Int)
    res = sample_eight_schools_accept_reject(n,d.data,trafo=d.trafo)
    x = zeros(10,n)
    for i in axes(x,2)
        x[:,i] .= collect(Tuple(Iterators.flatten(res.result[i])))
    end 
    x
end

Base.length(d::EightSchoolsAcceptReject) = 10
function Distributions.logpdf(d::EightSchoolsAcceptReject, x::AbstractArray{<:Real})
    #global a = x
    trafo = d.trafo
    y = Vector{Float64}(d.data["y"])
    sigma = d.data["sigma"]
    cov_matrix = Diagonal(sigma .^ 2)
    mu = x[9,:]
    tau = x[10,:]
    theta =  x[1:8,:]
    lgmu = [logpdf(Normal(0, 5), mu_i) for mu_i in mu]      # Prior for mu
    lgtau = 0
    lgtheta = 0
    if trafo
        tau = exp.(tau)
        lgtau = [logpdf(Cauchy(0, 5), tau_i) + log(tau_i) for tau_i in tau]  # Prior for tau
    else
        lgtau = [logpdf(Cauchy(0, 5), tau_i) for tau_i in tau]  # Prior for tau
    end
    lgtheta = [sum(logpdf(Normal(mu_i, tau_i), theta_i)) for (mu_i,tau_i,theta_i) in zip(mu,tau,eachcol(theta))]
    lgdata = [logpdf(MvNormal(theta_i, cov_matrix), y) for theta_i in eachcol(theta)]
    lgd = lgdata .+ lgmu .+ lgtau .+ lgtheta
    return lgd
end


bounds = NamedTupleDist(x = [i==10 ? 10^-100..100 : -100..100 for i in 1:10])
stan_data = JSON.parsefile("/ceph/groups/e4/users/slacagnina/MCBench/posteriordb/tmp.json")
eight_schools_testcase = Testcases(EightSchoolsAcceptReject(stan_data,false),bounds,10,"EightSchoolsAcceptReject")
eight_schools_testcase_trafo = Testcases(EightSchoolsAcceptReject(stan_data,true),bounds,10,"EightSchoolsAcceptReject")

# d = ESTest.f
# x = a
# #ESTest = Testcases(EightSchoolsAcceptReject(stan_data),10,"EightSchoolsAcceptReject")
# ESTest = Testcases(EightSchoolsAcceptReject(stan_data),bounds,10,"EightSchoolsAcceptReject")
# pwd()

# m = [marginal_mean(),marginal_variance(),global_mode(),marginal_skewness(),marginal_kurtosis(),sliced_wasserstein_distance(0.0,10^5)]
# t = ESTest
# s = BATMH(MCMCSampling(mcalg = MetropolisHastings(), nsteps=10^6, nchains=10,init = MCMCChainPoolInit(), burnin = MCMCMultiCycleBurnin(max_ncycles=10000)),"BAT-MH")

# build_teststatistic(t,m,n=1,n_steps=10^6,n_samples=10^5,par=true,s=s,clean=true)


# get_effective_sample_size(sbat,s)
# sbat = sample(ESTest,s)
# sbat = sample(ESTest,BATMH(n_steps=10^5,nchains=20))
# siid = sample(ESTest,10^6)




# plot(sbat,globalmode=true)
# plot(siid,globalmode=true)
# mode(sbat)
# mode(siid)


# plot(sample(ESTest,10^5))

# posterior = BAT.PosteriorMeasure(ESTest, bounds)
# posterior.function
# sam = bat_sample(posterior, BATMH().sampler).result
# plot(sam)


# x = rand(bounds,100)
# sh = varshape(bounds)
# x.x

# using InverseFunctions
# noshape = inverse(sh)
# xx = collect(hcat(x.x...))


# logpdf(ESTest.f,xx)
# logpdf(posterior,xx)


# rand(EightSchoolsAcceptReject(stan_data),12)

# s = sample(ts,10^5)
# plot(s)

# using BenchmarkTools
# @btime sample(ts,10^5)
# @btime sample(normal_100d_weakly_correlated,10^5)


# @btime sample_eight_schools_accept_reject(10^5,stan_data)

# @profview sample_eight_schools_accept_reject(10^5,stan_data)
# @profview sample(ts,10^5)

# d = EightSchoolsAcceptReject(stan_data)

# isless(Normal,Distribution)

# typeof(EightSchoolsAcceptReject) <: Distribution

# Distribution
# sample_eight_schools_accept_reject(1,stan_data).result
# collect(Tuple(Iterators.flatten(sample_eight_schools_accept_reject(1,stan_data).result[1])))



# # function rand(d::EightSchoolsAcceptReject, x::AbstractArray{<:Real})
# #     res = sample_eight_schools_accept_reject(size(x,2),d.data)
# #     indx = eachcol(x)
# #     for i in axes(x,2)
# #         x[:,i] .= collect(Tuple(Iterators.flatten(res.result[i])))[1]
# #     end 
# #     x
# # end