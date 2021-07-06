using NWNTools, StaticArrays, BenchmarkTools, LinearAlgebra, SparseArrays, LightGraphs, SimpleWeightedGraphs
dims = SA[50.0,50.0]
wprop = NWNTools.WireProp(22.0,0.02)
lines = rand_lines(dims, 1250, 7.0)
elecs = [Line(SA[0.0,0.0],SA[0.0,50.0])]
grnds = [Line(SA[50.0,0.0],SA[50.0,50.0])]
nwn = NWN_JDA(lines, repeat([wprop],1250), elecs, [5.0], grnds, dims, 7.0, 11.0);
A,z,inds=JDA_eqs(nwn)
