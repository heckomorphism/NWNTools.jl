export rand_lines, rand_wires, rand_network, rand_network_JDA, ensemble

using Random
using Distributions: Normal

"""
Generates a random collection of lines with initial point 
laying within the rectangle between the origin and `dims`. 
Determines second point of line to be distance 'l' and random 
direction (as sampled unifromly from the N-1 sphere) from 
the initial point. 

Parameters:
    dims:   Dimensions of rectangle containing lines.
    n:      Number of lines.
    l:      Length of lines.
"""
function rand_lines(dims::SVector{N,T},n::Int,l::Float64) where {N,T}
    x₀s = map(x->x.*dims,rand(SVector{N,T},n))
    A = rand(Normal(),N,n)
    # normalize columns of A to leave columns 
    # as uniformly sampled on the N-1 sphere
    A .= A./sqrt.(sum(A.^2,dims=1)).*l
    x₁s = [SVector{N}(c) for c ∈ eachcol(A)]   
    Line.(x₀s,x₀s+x₁s)
end

"""
Generates a random collection of wires using `rand_lines` 
and adds the properties of resistivity and wire diameter.

Parameters:
    dims:   Dimensions of rectangle containing lines.
    n:      Number of lines.
    l:      Length of lines.
    ρ:      Wire resistivity.
    D:      Wire diameter.
"""
function rand_network(dims::SVector{N,T},n,l,ρ,D) where {N,T}
    lines = rand_lines(dims::SVector{N,T},n::Int,l::Float64)
    wires = Wire.(lines, ρ, D)
    NWN(wires, dims)
end

function rand_network_JDA(dims::SVector{2,T},n,l,ρ,D,V₊,Rⱼ,Rₑ) where {T}
    lines = rand_lines(dims::SVector{2,T},n::Int,l::Float64)
    wprops = repeat([WireProp(ρ,D)],n)
    elec = Line(SA[0.0,0.0],SA[0.0,dims[2]])
    grnd = Line(SA[dims[1],0.0],dims)
    NWN_JDA(lines, wprops, [elec], [V₊], [grnd], dims, Rⱼ, Rₑ)
end

"""
Creates an ensemble of nanowire networks with the given 
paramaters. 
"""
function ensemble(n, dims::SVector{N,T},ρₙ,l,ρ,D) where {N,T}
    [rand_network(dims, round(Int, ρₙ*prod(dims)), l, ρ, D) for i∈1:n]
end

ensemble(n, params) = ensemble(n, params...)

