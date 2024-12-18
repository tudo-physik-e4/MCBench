# Overloaded function to run a test statistic with a test case, density sample vector, and sampling algorithm
function run_teststatistic(t::AT, samples::DensitySampleVector, m::TM, s::AS) where {TM <: TestMetric, AT <: AbstractTestcase, AS <: AnySampler}
    mval = calc_metric(t, samples, m)
end

# Overloaded function to run a test statistic with a test case, density sample vector, and integer
function run_teststatistic(t::Testcases, samples::DensitySampleVector, m::TM, s::Int) where {TM <: TestMetric}
    mval = calc_metric(t, samples, m)
end

# Function to build a test statistic
function build_teststatistic(t::T, m::Vector{TM}; s=IIDSampler(), n::Int=10^2, n_steps::Int=10^5, n_samples::Int=10^5, par::Bool=true, clean::Bool=false, unweight::Bool=true, use_sampler::Bool=true, iid=false) where {T<:AbstractTestcase, TM <: TestMetric}
    fnm = [string("./teststatistics/", string(t.info, "-", im.info, ".txt")) for im in m]
    if !isa(s, IIDSamplingAlgorithm) && !iid
        fnm = [string("./teststatistics_sampler/", string(t.info, "-", im.info, "-", s.info, ".txt")) for im in m]
    end
    if clean
        fs = [open(f, "w+") for f in fnm]
        close.(fs)
    end

    fs = [open(f, "a+") for f in fnm]
    try
        build_teststat_reshuffle(t, n, m, fs, s=s, n_steps=n_steps, n_samples=n_samples, unweight=unweight, par=par, use_sampler=use_sampler)
    catch
        close.(fs)
    end
end

# Function to build reshuffled test statistics
function build_teststat_reshuffle(t::T, n::Int, m::Vector{TM}, fnm::Vector{IOStream}; 
    s=IIDSampler(), n_steps=10^5, n_samples=0, unweight=true, par=false, use_sampler=true) where {TM <: TestMetric, T<:AbstractTestcase}
    v = 0
    j = 0
    while j < n
        v = use_sampler ? sample(t, s) : sample(t, s, n_steps=n_steps)
        if unweight && v.weight != ones(length(v))
            dsv_unw = v
            dsv_unw = resample_dsv_to_ess(dsv_unw, s)
            v = dsv_unw
        end
        if n_samples <= 0
            for i in 1:length(m)
                write(fnm[i], string([j.val for j in run_teststatistic(t, v, m[i], s)], "\n"))
            end
            a = j % (Int(floor((n_steps/100)))) == 0 ? println(j) : nothing
            j += 1
            continue
        end
        if n_samples > 0
            nv = v
            rv = 0
            while length(nv) >= 0 && !isempty(rv) && j < n
                while length(nv) < n_samples
                    rv = use_sampler ? sample(t, s) : sample(t, s, n_steps=n_steps)
                    if unweight && rv.weight != ones(length(v))
                        dsv_unw = rv
                        dsv_unw = resample_dsv_to_ess(dsv_unw, s)
                        rv = dsv_unw
                    end
                    nv = vcat(nv, rv)
                end
                while length(nv) >= n_samples && j < n
                    nvcalc = nv[1:n_samples]
                    if par
                        Folds.collect(write(fnm[i], string([j.val for j in run_teststatistic(t, nvcalc, m[i], s)], "\n")) for i in 1:length(m))
                    else
                        for i in 1:length(m)
                            write(fnm[i], string([j.val for j in run_teststatistic(t, nvcalc, m[i], s)], "\n"))
                        end
                    end
                    nv = length(nv) > n_samples ? nv[n_samples+1:end] : sample(t, s, n_steps=n_steps)
                    a = j % (Int(floor((n_steps/100)))) == 0 ? println(j) : nothing
                    j += 1
                end
            end
        end
    end
    close.(fnm)
end


# Function to read test statistic
function read_teststatistic(t::AbstractTestcase, m::TM) where {TM <: TestMetric}
    filename = string(t.info, "-", m.info, ".txt")
    parse_teststatistic(string("./teststatistics/", filename))
end

# Overloaded function to read test statistic with a sampling algorithm
function read_teststatistic(t::AbstractTestcase, m::TM, s::AnySampler) where {TM <: TestMetric}
    filename = string(t.info, "-", m.info, "-", s.info, ".txt")
    parse_teststatistic(string("./teststatistics_sampler/", filename))
end

# Function to parse test statistic from a file
function parse_teststatistic(filename::String)
    f = open(string(filename), "r")
    mvals = JSON.parse.(readlines(f))
    close(f)
    reshape(vcat(mvals...), length(mvals[1]), length(mvals))
end

# Function to parse a line from a file
function parseline(f::IOStream)
    iol = readline(f)
    if iol == ""
        []
    else
        return Vector{Float64}(JSON.parse(iol))
    end
end
