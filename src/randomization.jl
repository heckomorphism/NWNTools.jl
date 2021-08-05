export rand_lines, rand_network, rand_network_circuit, ensemble, ensemble_circuit

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
function rand_network(dims::SVector{N,T},n,l) where {N,T}
    lines = rand_lines(dims::SVector{N,T},n::Int,l::Float64)
    StickNetwork(lines, dims)
end

function rand_network_circuit(dims::SVector{2,T},n,l,ρ,D,V₊,Rⱼ,Rₑ,M) where {T}
    lines = rand_lines(dims::SVector{2,T},n::Int,l::Float64)
    wprops = repeat([WireProp(ρ,D)],n)
    elec = Line(SA[eps(T),zero(T)],SA[eps(T),dims[2]])
    grnd = Line(SA[dims[1]-eps(T),zero(T)],SA[dims[1]-eps(T),dims[2]])
    NWN_circuit{M}(lines, wprops, [elec], [V₊], [grnd], dims, Rⱼ, Rₑ)
end

function rand_network_circuit(dims::SVector{2,T},n,l,ρ,D,V₊,Rⱼ,Rₑ,M,::Val{:conducting}) where {T}
    lines = rand_lines(dims::SVector{2,T},n::Int,l::Float64)
    wprops = repeat([WireProp(ρ,D)],n)
    elec = Line(SA[eps(T),zero(T)],SA[eps(T),dims[2]])
    grnd = Line(SA[dims[1]-eps(T),zero(T)],SA[dims[1]-eps(T),dims[2]])
    nwn = NWN_circuit{M}(lines, wprops, [elec], [V₊], [grnd], dims, Rⱼ, Rₑ)
    while !is_conducting(nwn)
        lines = rand_lines(dims::SVector{2,T},n::Int,l::Float64)
        nwn = NWN_circuit{M}(lines, wprops, [elec], [V₊], [grnd], dims, Rⱼ, Rₑ)
    end
    nwn
end

"""
Creates an ensemble of nanowire networks with the given 
paramaters. 
"""
function ensemble(n,dims::SVector{N,T},ρₙ,l) where {N,T}
    [rand_network(dims, round(Int, ρₙ*prod(dims)), l) for i∈1:n]
end

function ensemble_circuit(n,dims::SVector{2,T},ρₙ,l,ρ,D,V₊,Rⱼ,Rₑ,M,args...) where {T}
    [rand_network_circuit(dims,round(Int,ρₙ*prod(dims)),l,ρ,D,V₊,Rⱼ,Rₑ,M,args...) for i∈1:n]
end

