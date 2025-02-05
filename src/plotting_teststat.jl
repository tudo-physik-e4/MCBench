"""
    plot_metrics(t::AbstractTestcase, normvals::Vector{NamedTuple}; plotalpha=0.2, infos=[], s=[]) 
    plot_teststatistic(t::AbstractTestcase, m::TM) where {TM <: TestMetric}
    plot_teststatistic(t::AbstractTestcase, m::TM, s::AnySampler; nbins=32, same_bins=true,sampler_bins=false)

Functions to plot test statistics and overview plots for a given test case and metrics.
These functions are subject to change and changed for recepies and are not guaranteed to be stable

# Arguments
- `t::AbstractTestcase`: The open file stream from which the line is read.
- `m::TM`: The metric for which the test statistic is plotted.
- `s::AnySampler`: The sampling algorithm for which the test statistic is plotted.
- `nbins::Int`: The number of bins for the histogram.
- `same_bins::Bool`: Whether to use the same bins for the test statistic and the sampling algorithm.
- `sampler_bins::Bool`: Whether to use the bins of the sampling algorithm for the test statistic or from the IID samples. This should be set to `true` as the distribution of the test statistics for the sampler i svery likely have a wider range than the IID samples.


# Example
```julia
t = MyTestcase()
m = [MyMetric_1(), MyMetric_2()]
s = MySampler()

plot_teststatistic(t, m, s; nbins=32, same_bins=true, sampler_bins=true)
plot_metrics(t, normvals; plotalpha=0.2, infos=[], s=[])

```
"""

# Function to plot the test statistic for a given test case and metric
function plot_teststatistic(t::AbstractTestcase, m::TM) where {TM <: TestMetric}
    mvals = read_teststatistic(t, m)
    dims = t.dim
    for idim in 1:dims
        histogram(mvals[idim, :], st=:stephist, title=string(t.info, "-", m.info, "-x", idim),
                  xlabel=m.info, label=string("n = ", length(mvals[idim, :])), ylabel="Entries",
                  size=(500, 2*500/3))
        savefig(string("./teststatistics/", t.info, "-", m.info, "-x", idim, ".pdf"))
    end
end

# Overloaded function to plot test statistic with additional sampling algorithm
function plot_teststatistic(t::AbstractTestcase, m::TM, s::AnySampler; nbins=32, same_bins=true,sampler_bins=false) where {TM <: TestMetric}
    mkpath(string("./", t.info))
    mvals = read_teststatistic(t, m)
    svals = read_teststatistic(t, m, s)
    dims = t.dim
    if isa(m,TwoSampleMetric)
        dims=1  
    end
    for idim in 1:dims
        hbins = fit(Histogram, Vector{Float64}(svals[idim, :]), nbins=nbins).edges[1]
        h = sampler_bins ? fit(Histogram, Vector{Float64}(mvals[idim, :]), hbins) : fit(Histogram, Vector{Float64}(mvals[idim, :]), nbins=nbins)
        sh = same_bins ? fit(Histogram, Vector{Float64}(svals[idim, :]), h.edges[1]) : fit(Histogram, Vector{Float64}(svals[idim, :]), nbins=nbins)
        h = normalize(h)
        sh = normalize(sh)
        plot(h, st=:step, label=string("IID, n = ", length(mvals[idim, :])), title=string(t.info, "-", m.info, "-x", idim), xlabel=m.info)
        plot!(ylabel="Entries", size=(500, 2*500/3), title=string(t.info, "-", m.info, "-x", idim), xlabel=m.info)
        plot!(sh, st=:step, xlabel=m.info, label=string(s.info, ", n = ", length(svals[idim, :])))
        savefig(string("./", t.info, "/", t.info, "-", m.info, "", s.info, "-x", idim, ".pdf"))
    end
end


# Core function to plot normalized metric values
function plot_metrics(t::AbstractTestcase, normvals::Vector{NamedTuple}; plotalpha=0.2, infos=[], s=[]) where {TM <: TestMetric}
    normvals=reverse(normvals)
    mkpath(string("./", t.info))
    ylen = length(normvals) * 30
    plot((0, 0), size=(550, ylen+100), legend=:topleft, xrange=(-3, 3), label="", bottom_margin=3Plots.mm, left_margin=3Plots.mm, framestyle=:box, dpi=400)
    xlabel!("metric - mean(metric) / std(metric)")
    
    yvals = [10 + 30*(i-1) for i in 1:length(normvals)]
    xmax = 3
    for i in 1:length(normvals)
        if findall(x -> x == :std, keys(normvals[i])) != []
            scatter!((normvals[i].val, yvals[i]), xerr=normvals[i].std, label="", color=:black)
            #plot!((normvals[i].val, yvals[i]), xerr=normvals[i].std, label="", color=:black)
        else
            scatter!((normvals[i].val, yvals[i]), label="", color=:black)
        end
        xmax = max(xmax, abs(normvals[i].val))
    end
    
    vline!([0], label="", color=:black)
    plot!([-1, 1], [ylen*2, ylen*2], fillrange=[-ylen, -ylen], label="", color=:green, fillalpha=plotalpha, lw=0)
    plot!([-2, -1], [ylen*2, ylen*2], fillrange=[-ylen, -ylen], label="", color=:yellow, fillalpha=plotalpha, lw=0)
    plot!([1, 2], [ylen*2, ylen*2], fillrange=[-ylen, -ylen], label="", color=:yellow, fillalpha=plotalpha, lw=0)
    plot!([-3, -2], [ylen*2, ylen*2], fillrange=[-ylen, -ylen], label="", color=:red, fillalpha=plotalpha, lw=0)
    plot!([2, 3], [ylen*2, ylen*2], fillrange=[-ylen, -ylen], label="", color=:red, fillalpha=plotalpha, lw=0)
    
    yticks!(yvals, [normvals[i].name for i in 1:length(yvals)], size=(550, ylen), xrange=((-xmax-0.1, xmax+0.1)), yrange=((-10, ylen+10)))
    if infos != []
        title!(string(t.info, "-"), size=(550, ylen+100))
    else
        title!(string(t.info), size=(550, ylen+100))
    end
    if s != []
        title!(string(t.info, "-", s.info), size=(550, ylen+100))
        savefig(string("./", t.info, "/", t.info, "-", s.info, "-metrics", ".pdf"))
        savefig(string("./", t.info, "/", t.info, "-", s.info, "-metrics", ".png"))
    else
        savefig(string("./", t.info, "/", t.info, "-metrics", ".pdf"))
        savefig(string("./", t.info, "/", t.info, "-metrics", ".png"))
    end
end

# Function to plot metrics for a given test case, metrics, and sampling algorithm
function plot_metrics(t::AbstractTestcase, ms::Vector{TM}, s::AnySampler; names=[]) where {TM <: TestMetric}
    normvals = Vector{NamedTuple}(undef, 0)
    for m in ms
        mvals = read_teststatistic(t, m)
        svals = read_teststatistic(t, m, s)
        dims = size(mvals)[1]
        for idim in 1:dims
            midval = mean(mvals[idim, :])
            stdval = std(mvals[idim, :])
            smidval = mean(svals[idim, :])
            sstdval = std(svals[idim, :])
            mvalnorm = (smidval - midval) / stdval
            mvalstdnorm = sstdval / stdval
            if names == []
                push!(normvals, (name=string(m.info, "-x", idim), val=mvalnorm, std=mvalstdnorm))
            else
                push!(normvals, (name=name=string(m.info, "-", names[idim]), val=mvalnorm, std=mvalstdnorm))
            end
            #push!(normvals, (name=string(m.info, "-x", idim), val=mvalnorm, std=mvalstdnorm))
        end
    end
    plot_metrics(t, normvals, s=s)
end