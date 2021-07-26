export graph_JDA, NWN_JDA, eqs_JDA, sheet_resistance_JDA, is_conducting, swap_index

using SparseArrays

"""
Takes an iterable collection of iterable collections of the 
form (a,b,c,...) and returns an object with "reversed" nesting 
of the form ((a[1],b[1],c[1],...),(a[2],b[2],c[2],...),...).
"""
function swap_index(a)
    # if a is empty then set n=0 otherwise n=length(a[1])
    n = length(a) 
    n==0 || (n=length(a[1])) 
    collect(broadcast(x->getindex(x,i),a) for i ∈ 1:n)
end

"""
Generates the graph of the nanowire network. Requires 
at least one of each of wire-wire, wire-electrode, 
and wire-ground intersections.
"""
function graph_JDA(lines, elecs, grnds, dims, Rⱼ, Rₑ)
    ww_src, ww_dst = connections(lines, dims)
    we_src, we_dst = connections(lines, elecs, dims)
    wg_src, wg_dst = connections(lines, grnds, dims)
    nₗ, nₑ, n₀ = length(lines), length(elecs), length(grnds)
    mₗ, mₑ, m₀ = length(ww_src), length(we_src), length(wg_src)
    src = vcat(ww_src, we_src, wg_src)
    dst = vcat(ww_dst, we_dst.+nₗ, wg_dst.+(nₗ+nₑ))
    wgt = vcat(ones(mₗ)./Rⱼ,ones(mₑ+m₀)./Rₑ)
    sp = sparse(vcat(src,dst), vcat(dst,src), vcat(wgt,wgt),nₗ+nₑ+n₀,nₗ+nₑ+n₀,+)
    SimpleWeightedGraph(sp)
end

# a constuctor for NWN_JDA that makes the graph
function NWN_JDA(lines::Array{Line{S,N},1}, props::Array{WireProp{T},1}, 
    elecs::Array{Line{S,N},1}, volts::Array{T,1}, grnds::Array{Line{S,N},1}, 
    dims::SVector{N,S}, Rⱼ::T, Rₑ::T) where {S<:Real,T<:Number,N}
    graph = graph_JDA(lines, elecs, grnds, dims, Rⱼ, Rₑ)
    NWN_JDA(lines, props, elecs, volts, grnds, dims, Rⱼ, Rₑ, graph)
end

function graph_MNR(lines, wprops, elecs, grnds, dims, Rⱼ, Rₑ)
    nₗ, nₑ, n₀ = length(lines), length(elecs), length(grnds)
    n = nₗ + nₑ + n₀
    line_ip = SVector{n}([Float64[] for i ∈ 1:n])
    juncts = NamedTuple{(:R,:c₁,:c₂),Tuple{Float64,Tuple{Int,Int},Tuple{Int,Int}}}[]
    each_intersect(lines, dims, []) do i₁,i₂,arr,dims,ps,p,arg
        push!(line_ip[i₁],ps[1])
        push!(line_ip[i₂],ps[2])
        push!(juncts,(R=Rⱼ,c₁=(i₁,length(line_ip[i₁])),c₂=(i₂,length(line_ip[i₂]))))
    end
    each_intersect(lines, vcat(elecs,grnds), dims, []) do i₁,i₂,arr₁,arr₂,dims,ps,p,arg
        i₂ += nₗ
        push!(line_ip[i₁],ps[1])
        push!(line_ip[i₂],ps[2])
        push!(juncts,(R=Rₑ,c₁=(i₁,length(line_ip[i₁])),c₂=(i₂,length(line_ip[i₂]))))
    end
    ls = length.(line_ip)
    ns = cumsum(ls).-ls
    nᵢ = ns[nₗ]+ls[nₗ]-sum(ls[1:nₗ].>1)
    nⱼ = length(juncts)
    perms = sortperm.(line_ip)
    invperms = invperm.(perms)
    permute!.(line_ip,perms)

    src, dst, wgt = zeros(Int,nᵢ+nⱼ), zeros(Int,nᵢ+nⱼ), zeros(Float64,nᵢ+nⱼ)
    for i ∈ 1:nₗ
        js = 1:(ls[i]-1)
        slice = js.+(ns[i]-i+1)
        Rₗ = len(lines[i])*wprops[i].ρ/(π*wprops[i].D*wprops[i].D)
        src[slice] = js.+ns[i]
        dst[slice] = js.+(1+ns[i]) 
        wgt[slice] = 1 ./(Rₗ.*(diff(line_ip[i])))
    end
    for (i,ju) ∈ enumerate(juncts)
        src[i+nᵢ] = ns[ju.c₁[1]]+invperms[ju.c₁[1]][ju.c₁[2]]
        dst[i+nᵢ] = ns[ju.c₂[1]]+invperms[ju.c₂[1]][ju.c₂[2]]
        wgt[i+nᵢ] = 1/ju.R
    end
    sp = sparse(vcat(src,dst), vcat(dst,src), vcat(wgt,wgt),ns[end]+ls[end],ns[end]+ls[end],+)
    SimpleWeightedGraph(sp), ls
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

"""
Calculates the sheet resistance of a NWN object. 
Solves the network using MNA and returns the ratio 
of the voltage to the current. Requires the network 
to have only one electrode.
"""
function sheet_resistance_JDA(nwn::NWN_JDA{S,T,N}) where {S,T,N}
    @assert length(nwn.elecs)==1 "Sheet resistance is only defined for one electrode and one ground."
    A, z, inds = eqs_JDA(nwn)
    sol = A\collect(z)
    nwn.volts[1]/sol[end]
end

"""
Returns a boolean value for if the given network 
has a path from at least one electrode and ground.
"""
function is_conducting(nwn::NWN_JDA{S,T,N}) where {S,T,N}
    nₗ, nₑ = length(nwn.lines), length(nwn.elecs)
    for i ∈ 1:length(nwn.elecs)
        for j ∈ 1:length(nwn.grnds)
            if has_path(nwn.graph,i+nₗ,j+nₗ+nₑ)
                return true
            end
        end
    end
    return false
end
