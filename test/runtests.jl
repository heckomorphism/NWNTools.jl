using NWNTools
using Test
using DelimitedFiles, StaticArrays

@testset "NWNTools.jl" begin
    @testset "Benchmark1" begin
        data = readdlm("benchmark1.txt",',', BigFloat)
        dims = SA[BigFloat("3.0"),BigFloat("5.0")]
        elecs = [Line(SVector{2}(data[1,:]),SVector{2}(data[2,:]))]
        grnds = [Line(SVector{2}(data[3,:]),SVector{2}(data[4,:]))]
        Rⱼ = BigFloat("20.0")
        V₊ = BigFloat("5.0")
        volts = [V₊]
        ρₑ = BigFloat("22.63676")/BigFloat("1000")
        D = BigFloat("0.06")

        lines = Line{BigFloat,2}[]
        for i ∈ 5:2:size(data)[1]
            push!(lines, Line(SVector{2}(data[i,:]),SVector{2}(data[i+1,:])))
        end
        props = repeat([WireProp(ρₑ,D)],length(lines))

        nwn_JDA = NWN_circuit{:JDA}(lines, props, elecs, volts, grnds, dims, Rⱼ,Rⱼ)
        R_JDA = sheet_resistance(nwn_JDA)
        nwn_MNR = NWN_circuit{:MNR}(lines, props, elecs, volts, grnds, dims, Rⱼ,Rⱼ)
        R_MNR = sheet_resistance(nwn_MNR)
        @test R_JDA≈(8*Rⱼ/3)
        numerator = 72*(3+5*sqrt(big(2)))*ρₑ^2+2*D^2*big(π)*Rⱼ*((81+69*sqrt(big(2)))*ρₑ+20*D^2*big(π)*Rⱼ)
        denominator = 15*D^2*big(π)*(2*(1+sqrt(big(2)))*ρₑ+D^2*big(π)*Rⱼ)
        R_MNR_analytic = numerator/denominator
        @test R_MNR≈R_MNR_analytic
    end
end
