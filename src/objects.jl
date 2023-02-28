export Line, Wire, StickNetwork, WireProp, NWN_circuit, len

using StaticArrays, LinearAlgebra, Graphs, SimpleWeightedGraphs


"""
```julia
Line(p₁,p₂)
```
Represents a line segment with end points `p₁` and `p₂` given 
by `SVector{N,T}`. Parameterized by the datatype of the 
coordinates, `T`, and the length (i.e. dimension) `N`. 
"""
struct Line{T<:Real,N}
    p₁::SVector{N,T}
    p₂::SVector{N,T}
end

"""
```julia
WireProp(ρ,D)
```
Contains the resistivity, `ρ`, and diameter, `D`, data for a nanowire.
"""
struct WireProp{T<:Number}
    ρ::T # Resistivity
    D::T # Wire diameter
end

"""
```julia
StickNetwork(lines,dims)
```
Stores a network of sticks, `lines`, bounded in the region `[0,0]` to `dims`.
"""
struct StickNetwork{T<:Real,N}
    lines::Vector{Line{T,N}}
    dims::SVector{N,T}
end

"""
```julia
NWN_circuit{M}(lines::Vector{Line{T,N}},props::Vector{WireProp{T}},elecs::Vector{Line{T,N}},volts::Vector{T},grnds::Vector{Line{T,N}},dims::SVector{N,T},Rⱼ::T,Rₑ::T)
```
Contains the data to perform circuit analysis of a nanowire 
network including all geometric, electrical, and structural 
data.
"""
mutable struct NWN_circuit{M,T<:Number,N} # T is numeric type, N is number of spatial dimensions, M is Model parameter
    lines::Vector{Line{T,N}}
    props::Vector{WireProp{T}}
    elecs::Vector{Line{T,N}}
    volts::Vector{T}
    grnds::Vector{Line{T,N}}
    dims::SVector{N,T}
    Rⱼ::T
    Rₑ::T
    graph::SimpleWeightedGraph{Int, T}
    ls::Vector{Int}
    ns::Vector{Int}
    function NWN_circuit{M}(lines::Vector{Line{T,N}},props::Vector{WireProp{T}},elecs::Vector{Line{T,N}},volts::Vector{T},grnds::Vector{Line{T,N}},dims::SVector{N,T},Rⱼ::T,Rₑ::T) where {M,T,N}
        update_graph(new{M,T,N}(lines,props,elecs,volts,grnds,dims,Rⱼ,Rₑ))
    end
end

NWN_circuit{M2}(o::NWN_circuit{M1}) where {M1,M2} = NWN_circuit{M2}(o.lines, o.props, o.elecs, o.volts, o.grnds, o.dims, o.Rⱼ, o.Rₑ)



# Conversion methods
import Base.convert
Line(p₁::Vector{T},p₂::Vector{S}) where {T,S} = Line(promote(SVector{length(p₁)}(p₁), SVector{length(p₂)}(p₂))...)
convert(::Type{Line{T,N}},x::Line{S,N}) where {T,S,N} = Line(convert(SVector{N,T},x.p₁), convert(SVector{N,T},x.p₂))

"""
```julia
len(ℓ::Line{T,N})
```
Determines the Euclidean length of the line segment ℓ.
"""
len(ℓ::Line{T,N}) where {T,N} = sqrt(sum((ℓ.p₁[i]-ℓ.p₂[i])^2 for i ∈ 1:N))