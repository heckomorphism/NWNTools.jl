using NWNTools, StaticArrays, BenchmarkTools, LinearAlgebra, SparseArrays, LightGraphs, SimpleWeightedGraphs
dims = SA[1.0,1.0]
lines = [Line(SA[0.0,0.5],SA[1.0,0.5]),Line(SA[0.5,0.0],SA[0.5,0.5])]
wprops = repeat([WireProp(22.0,0.02)],length(lines))
elecs = [Line(SA[0.0,0.0],SA[0.0,1.0]),Line(SA[0.0,0.0],SA[1.0,0.0])]
volts = [5.0,3.3]
grnds = [Line(SA[1.0,0.0],SA[1.0,1.0])]
nwn = NWN_JDA(lines, wprops, elecs, volts, grnds, dims, 7.0,11.0)
