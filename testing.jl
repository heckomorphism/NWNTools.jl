using NWNTools, StaticArrays, BenchmarkTools, LinearAlgebra, SparseArrays, LightGraphs, SimpleWeightedGraphs
dims = SA[20.0,20.0]
ρₙ = 0.4
l = 7.0
ρₑ = 22.2
D = 0.02
V₊ = 5.0
Rⱼ = 11.0
Rₑ = 7.0

nwn = rand_network_JDA(dims, round(Int,prod(dims)*ρₙ), l, ρₑ, D, V₊, Rⱼ, Rₑ)

e = ensemble_JDA(100, dims, ρₙ, l, ρₑ, D, V₊, Rⱼ, Rₑ)