export Line, Wire, NWN, WireProp, NWN_circuit, len

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

"""
A structure to hold all of the information relevant to a wire.
"""
struct Wire{T<:Real,N}
    line::Line{T,N} # Wire geomoetry
    # length::Float64 # Wire length
    ρ::Float64 # Resistivity
    D::Float64 # Wire diameter
    # A::Float64 # Cross sectional area
    # Rℓ⁻¹::Float64 # Resistance per length
end

struct WireProp{T<:Number}
    ρ::T # Resistivity
    D::T # Wire diameter
end

"""
A structure to store a collection of wires in a certain size. 
"""
struct NWN{T<:Real,N}
    wires::Array{Wire{T,N},1}
    dims::SVector{N,T}
end

mutable struct NWN_circuit{M,T<:Number,N} # T is numeric type, N is number of spatial dimensions, M is Model parameter
    lines::Array{Line{T,N},1}
    props::Array{WireProp{T},1}
    elecs::Array{Line{T,N},1}
    volts::Array{T,1}
    grnds::Array{Line{T,N},1}
    dims::SVector{N,T}
    Rⱼ::T
    Rₑ::T
    graph::SimpleWeightedGraph{Int, T}
    ls::Vector{Int}
    ns::Vector{Int}
    function NWN_circuit{M}(lines::Array{Line{T,N},1},props::Array{WireProp{T},1},elecs::Array{Line{T,N},1},volts::Array{T,1},grnds::Array{Line{T,N},1},dims::SVector{N,T},Rⱼ::T,Rₑ::T) where {M,T,N}
        update_graph(new{M,T,N}(lines,props,elecs,volts,grnds,dims,Rⱼ,Rₑ))
    end
end



# Conversion methods
import Base.convert
Line(p₁::Array{T,1},p₂::Array{S,1}) where {T,S} = Line(promote(SVector{length(p₁)}(p₁), SVector{length(p₂)}(p₂))...)
convert(::Type{Line{T,N}},x::Line{S,N}) where {T,S,N} = Line(convert(SVector{N,T},x.p₁), convert(SVector{N,T},x.p₂))

"""
Determines the Euclidean length of the line segment ℓ.
"""
len(ℓ::Line{T,N}) where {T,N} = sqrt(sum((ℓ.p₁[i]-ℓ.p₂[i])^2 for i ∈ 1:N))
len(w::Wire{T,N}) where {T,N} = len(w.line)