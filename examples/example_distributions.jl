#Normal Distributions

f = Normal()
bounds = NamedTupleDist(x = [-10..10])
normal_1d_uncorrelated = Testcases(f,bounds,1,"Normal-1D-Uncorrelated")

f = MvNormal(zeros(2), I(2))
bounds = NamedTupleDist(x = [-10..10 for i in 1:2])
normal_2d_uncorrelated = Testcases(f,bounds,2,"Normal-2D-Uncorrelated")

f = MvNormal(zeros(3), I(3))
bounds = NamedTupleDist(x = [-10..10 for i in 1:3])
normal_3d_uncorrelated = Testcases(f,bounds,3,"Normal-3D-Uncorrelated")

f = MvNormal(zeros(10), I(10))
bounds = NamedTupleDist(x = [-10..10 for i in 1:10])
normal_10d_uncorrelated = Testcases(f,bounds,10,"Normal-10D-Uncorrelated")

f = MvNormal(zeros(100), I(100))
bounds = NamedTupleDist(x = [-10..10 for i in 1:100])
normal_100d_uncorrelated = Testcases(f,bounds,100,"Normal-100D-Uncorrelated")

f = MvNormal(zeros(2), [1.0 0.3; 0.3 1.0])
bounds = NamedTupleDist(x = [-10..10 for i in 1:2])
normal_2d_weakly_correlated = Testcases(f,bounds,2,"Normal-2D-Weakly-Correlated")

f = MvNormal(zeros(2), [1.0 0.9; 0.9 1.0])
bounds = NamedTupleDist(x = [-10..10 for i in 1:2])
normal_2d_strongly_correlated = Testcases(f,bounds,2,"Normal-2D-Strongly-Correlated")

f = MvNormal(zeros(10), ones(10,10)*0.2 + I(10)*0.8)
bounds = NamedTupleDist(x = [-10..10 for i in 1:10])
normal_10d_weakly_correlated = Testcases(f,bounds,10,"Normal-10D-Weakly-Correlated")

f = MvNormal(zeros(10), ones(10,10)*0.9 + I(10)*0.1)
bounds = NamedTupleDist(x = [-10..10 for i in 1:10])
normal_10d_strongly_correlated = Testcases(f,bounds,10,"Normal-10D-Strongly-Correlated")

f = MvNormal(zeros(100), ones(100,100)*0.2 + I(100)*0.8)
bounds = NamedTupleDist(x = [-10..10 for i in 1:100])
normal_100d_weakly_correlated = Testcases(f,bounds,100,"Normal-100D-Weakly-Correlated")

f = MvNormal(zeros(100), ones(100,100)*0.9 + I(100)*0.1)
bounds = NamedTupleDist(x = [-10..10 for i in 1:100])
normal_100d_strongly_correlated = Testcases(f,bounds,100,"Normal-100D-Strongly-Correlated")

f = MixtureModel([Normal(2,1), Normal(-2,1)], [0.5, 0.5])
bounds = NamedTupleDist(x = [-10..10])
normal_1d_multimodal_4std = Testcases(f,bounds,1,"Normal-1D-Multimodal-4std")

f = MixtureModel([Normal(10,1), Normal(-10,1)], [0.5, 0.5])
bounds = NamedTupleDist(x = [-20..20])
normal_1d_multimodal_20std = Testcases(f,bounds,1,"Normal-1D-Multimodal-20std") 

f = MixtureModel([Normal(2,1), Normal(-2,1)], [0.25, 0.75])
bounds = NamedTupleDist(x = [-10..10])
normal_1d_multimodal_4std_1to3 = Testcases(f,bounds,1,"Normal-1D-Multimodal-4std-1to3")

f = MixtureModel([Normal(10,1), Normal(-10,1)], [0.25, 0.75])
bounds = NamedTupleDist(x = [-20..20])
normal_1d_multimodal_20std_1to3 = Testcases(f,bounds,1,"Normal-1D-Multimodal-20std-1to3")


##Cauchy Distributions
f = Cauchy()
bounds = NamedTupleDist(x = [-10..10])
cauchy_1d = Testcases(f,bounds,1,"Cauchy-1D")


##Multimodal Mixture 
r = 5
f1 = MvNormal(r*ones(10), ones(10,10)*0.9 + I(10)*0.1)
f2 = MvNormal(-r*ones(10), ones(10,10)*0.9 + I(10)*0.1)
f = MixtureModel([f1,f2], [0.25, 0.75])
bounds = NamedTupleDist(x = [-100..100 for i in 1:10])
normal_10d_multimodal_10std = Testcases(f,bounds,10,"Normal-10D-Multimodal-10std")

r = 5
f1 = MvNormal(r*ones(3), ones(3,3)*0.9 + I(3)*0.1)
f2 = MvNormal(-r*ones(3), ones(3,3)*0.9 + I(3)*0.1)
f = MixtureModel([f1,f2], [0.25, 0.75])
bounds = NamedTupleDist(x = [-100..100 for i in 1:3])
normal_3d_multimodal_10std = Testcases(f,bounds,3,"Normal-3D-Multimodal-10std")

