using NWNTools, StaticArrays, BenchmarkTools, LinearAlgebra, SparseArrays, LightGraphs, SimpleWeightedGraphs
dims = SA[20.0,20.0];
ρₙ = 0.4;
l = 7.0;
ρₑ = 22.2/1000; #put in micro Ohm metres
D = 0.02;
V₊ = 5.0;
Rⱼ = 11.0;
Rₑ = 7.0;
lines = rand_lines(dims, round(Int, prod(dims)*ρₙ), l);
props = repeat([WireProp(ρₑ,D)],length(lines));
grnds = [Line(SA[dims[1],0.0],dims)];
elecs = [Line(SA[0.0,0.0],SA[0.0,dims[2]])];

nwn = NWN_circuit{:MNR}(lines, props, elecs, [V₊], grnds, dims, Rⱼ, Rₑ);
