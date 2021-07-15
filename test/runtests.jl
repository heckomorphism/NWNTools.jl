using NWNTools
using Test
using DelimitedFiles, StaticArrays

@testset "NWNTools.jl" begin
    @testset "Benchmark1" begin
        data = readdlm("benchmark1.txt",',')
        dims = SA[3.0,5.0]
        elecs = [Line(SVector{2}(data[1,:]),SVector{2}(data[2,:]))]
        grnds = [Line(SVector{2}(data[3,:]),SVector{2}(data[4,:]))]

        lines = Line{Float64,2}[]
        for i ∈ 5:2:size(data)[1]
            push!(lines, Line(SVector{2}(data[i,:]),SVector{2}(data[i+1,:])))
        end
        wprops = repeat([WireProp(22.63676,0.06)],length(lines))

        nwn = NWN_JDA(lines, wprops, elecs, [5.0], grnds, dims, 20.0,20.0)
        R = sheet_resistance_JDA(nwn)
        @test R≈(160/3)
    end
end
