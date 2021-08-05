export rand_lines, rand_network, rand_network_circuit, ensemble, ensemble_circuit

using Random
using Distributions: Normal

"""
```julia
rand_lines(dims,n,l)
```
Generates a random collection of lines with initial point 
laying within the rectangle between the origin and `dims`. 
Determines second point of line to be distance 'l' and random 
direction (as sampled unifromly from the N-1 sphere) from 
the initial point.
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
```julia
rand_network(dims,n,l)
```
Generates a random `StickNetwork` using `rand_lines(dims,n,l)` 
and adds the properties of resistivity and wire diameter.
"""
function rand_network(dims::SVector{N,T},n,l) where {N,T}
    lines = rand_lines(dims::SVector{N,T},n::Int,l::Float64)
    StickNetwork(lines, dims)
end

"""
```julia
rand_network_circuit(dims,n,l,ρ,D,V₊,Rⱼ,Rₑ,M)
```
Generates a random nanowire network using model `M` with dimensions `dims`, 
number of wires `n`, wire resistivity `ρ`, wire diameter `D`, voltage `V₊`, 
junction resistance `Rⱼ`, and electrode resistance `Rₑ`. Electrode is along 
the leftmost extent of the region and ground along the rightmost extent. 
Electrode and ground are positioned `eps(T)` in from the edges to minimize 
rejected intersections.
"""
function rand_network_circuit(dims::SVector{2,T},n,l,ρ,D,V₊,Rⱼ,Rₑ,M) where {T}
    lines = rand_lines(dims::SVector{2,T},n::Int,l::Float64)
    wprops = repeat([WireProp(ρ,D)],n)
    elec = Line(SA[eps(T),zero(T)],SA[eps(T),dims[2]])
    grnd = Line(SA[dims[1]-eps(T),zero(T)],SA[dims[1]-eps(T),dims[2]])
    NWN_circuit{M}(lines, wprops, [elec], [V₊], [grnd], dims, Rⱼ, Rₑ)
end

"""
```julia
rand_network_circuit(dims,n,l,ρ,D,V₊,Rⱼ,Rₑ,M,Val(:conducting))
```
Generates a random nanowire network using model `M` which is gaurenteed to conduct. Otherwise the same as `rand_network_circuit`.
"""
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
```julia
ensemble(n,dims,ρₙ,l)
```
Creates an ensemble of `StickNetwork`s with dimensions `dims`, stick densites `ρₙ`, and stick lengths `l`.
"""
function ensemble(n,dims::SVector{N,T},ρₙ,l) where {N,T}
    [rand_network(dims, round(Int, ρₙ*prod(dims)), l) for i∈1:n]
end

"""
```julia
ensemble_circuit(n,dims,ρₙ,l,ρ,D,V₊,Rⱼ,Rₑ,M,args...)
```
Generates an ensemble of nanowire networks by calling `rand_network_circuit(dims,round(Int,ρₙ*prod(dims)),l,ρ,D,V₊,Rⱼ,Rₑ,M,args...)` `n` times.
"""
function ensemble_circuit(n,dims::SVector{2,T},ρₙ,l,ρ,D,V₊,Rⱼ,Rₑ,M,args...) where {T}
    [rand_network_circuit(dims,round(Int,ρₙ*prod(dims)),l,ρ,D,V₊,Rⱼ,Rₑ,M,args...) for i∈1:n]
end

