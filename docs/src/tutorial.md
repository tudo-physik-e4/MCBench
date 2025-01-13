# Using MCBench
This is a simple example of how to use the MCBench package. 

## Specify the test case
```
f = MvNormal(zeros(3), I(3))
bounds = NamedTupleDist(x = [-10..10 for i in 1:3])
Standard_Normal_3D_Uncorrelated = Testcases(f,bounds,3,"Normal-3D-Uncorrelated")
```
 
## Selecting the metrics to be applied
```
metrics = [marginal_mean(), marginal_variance(), sliced_wasserstein_distance(), maximum_mean_discrepancy()]
```
 
## Load the external MC samples to be tested
```
sampler = FileBasedSampler("samples_from_my_algorithm.csv")
```

## Generating the teststatistics 
Evaluate the metrics both

- for the IID samples (IID samples are generated automatically in the background):
```
teststatistics_IID = build_teststatistic(Standard_Normal_3D_Uncorrelated, metrics,
n=100, n_steps=10^5, n_samples=10^5)
```
- and for the MC samples to be tested:
```
teststatistics_my_samples = build_teststatistic(Standard_Normal_3D_Uncorrelated, metrics,
n=100, n_steps=10^5, n_samples=10^5, s=sampler)
```
## Generating comparison plots
- Overview plot of all selected metrics
```
plot_metrics(Standard_Normal_3D_Uncorrelated, metrics, sampler)
```
<img src="../images/Normal-3D-Uncorrelated-metrics.svg" width="480"/>

- Individual metrics
```
plot_teststatistic(Standard_Normal_3D_Uncorrelated, marginal_mean(), sampler, nbins=20)
```
<img src="../images/Normal-3D-Uncorrelated-SlicedWasserstein.svg" width="480"/>
