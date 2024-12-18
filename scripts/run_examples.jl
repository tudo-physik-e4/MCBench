include("../MCBench.jl") 
mkpath("teststatistics_sampler")
mkpath("teststatistics")


BSam = BATMH(MCMCSampling(mcalg = MetropolisHastings(),nsteps = 10^5,nchains = 10,init = MCMCChainPoolInit(),strict=false,burnin = MCMCMultiCycleBurnin(max_ncycles = 30)),"BAT-MH")
#custom settings for BAT sampler as example using the multimodal results in non-convergence 
function run_BATMH(testcase)
    @sync begin
        @async build_teststatistic(testcase,m,n=100,n_steps=10^5,n_samples=10^5,par=true,clean=true)
        @async build_teststatistic(testcase,m,n=100,n_steps=10^5,n_samples=10^5,par=true,clean=false,s=BSam)
    end
    plot_metrics(testcase,m,BSam)
    plot_teststatistic(testcase,m[1],BSam)
end

m = [marginal_mean(),marginal_variance(),global_mode(),marginal_skewness(),marginal_kurtosis(),sliced_wasserstein_distance(0.0,10^5),maximum_mean_discrepancy(0.0,10^4)]

run_BATMH(normal_3d_uncorrelated)
run_BATMH(normal_3d_multimodal_10std)

m = [marginal_mean(),marginal_variance(),sliced_wasserstein_distance(0.0,10^5),maximum_mean_discrepancy(0.0,10^4)]

plot_metrics(normal_3d_uncorrelated,m,BSam)
plot_metrics(normal_3d_multimodal_10std,m,BSam)

plot_teststatistic(normal_3d_uncorrelated,m[1],BATMH(),sampler_bins=true,nbins=16)
plot_teststatistic(normal_3d_uncorrelated,m[3],BATMH(),sampler_bins=true,nbins=16)

plot_teststatistic(normal_3d_multimodal_10std,m[1],BATMH(),sampler_bins=true,nbins=16)
plot_teststatistic(normal_3d_multimodal_10std,m[3],BATMH(),sampler_bins=true,nbins=16)