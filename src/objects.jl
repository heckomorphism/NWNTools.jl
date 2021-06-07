export Line, Wire, NWN, len

using StaticArrays


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

"""
A structure to store a collection of wires in a certain size. 
"""
struct NWN{T<:Real,N}
    wires::Array{Wire{T,N},1}
    dims::SVector{N,T}
end



import Base.convert
Line(p₁::Array{T,1},p₂::Array{S,1}) where {T,S} = Line(promote(SVector{length(p₁)}(p₁), SVector{length(p₂)}(p₂))...)
convert(::Type{Line{T,N}},x::Line{S,N}) where {T,S,N} = Line(convert(SVector{N,T},x.p₁), convert(SVector{N,T},x.p₂))

"""
Determines the Euclidean length of the line segment ℓ.
"""
len(ℓ::Line{T,N}) where {T,N} = sqrt(sum((ℓ.p₁[i]-ℓ.p₂[i])^2 for i ∈ 1:N))
len(w::Wire{T,N}) where {T,N} = len(w.line)