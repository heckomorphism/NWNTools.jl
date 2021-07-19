export Line, Wire, NWN, WireProp, NWN_JDA, len

using StaticArrays, LinearAlgebra, LightGraphs, SimpleWeightedGraphs


"""
A structure representing a line with end points 
p₁ and p₂. Parameterized by the datatype of the 
coordinates and the length (i.e. dimension) N. 
"""
struct Line{T<:Real,N}
    p₁::SVector{N,T}
    p₂::SVector{N,T}
end

struct WireProp{T<:Number}
    ρ::T # Resistivity
    D::T # Wire diameter
end

mutable struct NWN{S<:Real,T<:Number,N,M,P}
    lines::Vector{Line{S,N}}
    props::Vector{WireProp{T}}
    elecs::Vector{Line{S,N}}
    volts::Vector{T}
    grnds::Vector{Line{S,N}}
    dims::SVector{N,S}
    Rⱼ::T
    Rₑ::T
    graph::SimpleWeightedGraph{Int, T}
    function NWN{S,T,N,M,P}(
        lines::Vector{Line{S,N}},props::Vector{WireProp{T}},
        elecs::Vector{Line{S,N}},volts::Vector{T},grnds::Vector{Line{S,N}},
        dims::SVector{N,S},Rⱼ::T,Rₑ::T
        ) where {S,T,N,M,P}

        nwn = new{S,T,N,M,P}(lines, props, elecs, volts, grnds, dims, Rⱼ, Rₑ)
        nwn.graph = graph(nwn)
        nwn
    end
end

function graph(nwn::NWN{S,T,N,:JDA,P}) where {S,T,N,P}
    ww_src, ww_dst = connections(nwn.lines, nwn.dims)
    we_src, we_dst = connections(nwn.lines, nwn.elecs, nwn.dims)
    wg_src, wg_dst = connections(nwn.lines, nwn.grnds, nwn.dims)
    nₗ, nₑ = length(nwn.lines), length(nwn.elecs)
    mₗ, mₑ, m₀ = length(ww_src), length(we_src), length(wg_src)
    src = vcat(ww_src, we_src, wg_src)
    dst = vcat(ww_dst, we_dst.+nₗ, wg_dst.+(nₗ+nₑ))
    wgt = vcat(ones(mₗ)./nwn.Rⱼ,ones(mₑ+m₀)./nwn.Rₑ)
    g = SimpleWeightedGraph(src, dst, wgt)
    g
end
