using Distributions 
using Statistics 
using ValueShapes
using LinearAlgebra
using Folds
using IntervalSets
using BAT
using StatsBase
using Plots
using JSON  #for pasting teststatistic quickly 
using StructArrays
using HypothesisTests
using DensityInterface
import DensityInterface: logdensityof
using IPMeasures
using Distances 
using DataFrames
using CSV

include("src/samplers.jl")
include("src/testcases.jl")
include("src/mmd.jl")
include("src/testmetrics.jl")
include("src/wasserstein.jl")
include("src/twosampleteststatics.jl")
include("src/marginal_utils.jl")
include("src/teststatistic.jl")
include("src/sample_utils.jl")
include("src/plotting_teststat.jl")
include("src/batmh.jl")
include("examples/example_distributions.jl")
include("examples/example_posteriordb.jl")