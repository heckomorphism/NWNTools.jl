export graph_JDA, NWN_JDA, eqs_JDA, sheet_resistance_JDA

using SparseArrays

swap_index(a) = collect(broadcast(x->getindex(x,i),a) for i ∈ eachindex(a[1]))

"""
Generates the graph of the nanowire network.
"""
function graph_JDA(lines, elecs, grnds, dims, Rⱼ, Rₑ)
    ww_src, ww_dst = swap_index(connections(lines, dims))
    we_src, we_dst = swap_index(connections(lines, elecs, dims))
    wg_src, wg_dst = swap_index(connections(lines, grnds, dims))
    nₗ, nₑ = length(lines), length(elecs)
    mₗ, mₑ, m₀ = length(ww_src), length(we_src), length(wg_src)
    src = vcat(ww_src, we_src, wg_src)
    dst = vcat(ww_dst, we_dst.+nₗ, wg_dst.+(nₗ+nₑ))
    wgt = vcat(ones(mₗ)./Rⱼ,ones(mₑ+m₀)./Rₑ)
    g = SimpleWeightedGraph(src, dst, wgt)
    g
end

# a constuctor for NWN_JDA that makes the graph
function NWN_JDA(lines::Array{Line{S,N},1}, props::Array{WireProp{T},1}, 
    elecs::Array{Line{S,N},1}, volts::Array{T,1}, grnds::Array{Line{S,N},1}, 
    dims::SVector{N,S}, Rⱼ::T, Rₑ::T) where {S<:Real,T<:Number,N}
    graph = graph_JDA(lines, elecs, grnds, dims, Rⱼ, Rₑ)
    NWN_JDA(lines, props, elecs, volts, grnds, dims, Rⱼ, Rₑ, graph)
end

"""
Creates the system of equations which solves the 
electrical problem of the nanowire network under 
the JDA. Removes ill conditioned nodes (those not 
connected through the circuit to a voltage source 
or ground).
"""
function eqs_JDA(nwn::NWN_JDA{S,T,N}) where {S,T,N}
    nₗ, nₑ, n₀ = length(nwn.lines), length(nwn.elecs), length(nwn.grnds)

    # include only nodes part of the actual circuit
    con_comps = connected_components(nwn.graph)
    connected = vcat(con_comps[.!(con_comps .⊆ [1:nₗ])]...)

    # remove grounds from equations as MNA requires
    inds = setdiff(connected, nₗ+nₑ+1:nₗ+nₑ+n₀)

    # view of Laplacian of graph with indices needed
    L = -laplacian_matrix(nwn.graph)[inds,inds]

    # generates B matrix
    B = spzeros(length(inds),nₑ)
    for i ∈ 1:nₑ
        B[findfirst(inds.==(nₗ+i)),i] = 1
    end

    # generates system of equations in the form Ax=z
    A = [L B;
        transpose(B) spzeros(nₑ,nₑ)]
    z = [spzeros(length(inds));
        sparse(nwn.volts)]
    A,z,inds
end

function sheet_resistance_JDA(nwn)
    @assert length(nwn.elecs)==1 "Sheet resistance is only defined for one electrode and one ground."
    A, z, inds = JDA_eqs(nwn)
    sol = A\collect(z)
    nwn.volts[1]/sol[end]
end