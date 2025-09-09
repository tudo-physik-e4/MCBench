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
stan_data = JSON.parsefile(joinpath(@__DIR__, "..", "examples/eight_schools.json"))
eight_schools_testcase = Testcases(EightSchoolsAcceptReject(stan_data,false),bounds,10,"EightSchoolsAcceptReject")
eight_schools_testcase_trafo = Testcases(EightSchoolsAcceptReject(stan_data,true),bounds,10,"EightSchoolsAcceptReject")
