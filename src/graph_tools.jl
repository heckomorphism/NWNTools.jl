export update_graph, circuit_eqs, sheet_resistance, is_conducting, swap_index

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
Generates the graph of the nanowire network.
"""
function update_graph(nwn::NWN_circuit{:JDA,T,N}) where {T,N}
    ww_src, ww_dst = connections(nwn.lines, nwn.dims)
    we_src, we_dst = connections(nwn.lines, nwn.elecs, nwn.dims)
    wg_src, wg_dst = connections(nwn.lines, nwn.grnds, nwn.dims)
    nₗ, nₑ, n₀ = length(nwn.lines), length(nwn.elecs), length(nwn.grnds)
    n = nₗ + nₑ + n₀
    mₗ, mₑ, m₀ = length(ww_src), length(we_src), length(wg_src)
    src = vcat(ww_src, we_src, wg_src)
    dst = vcat(ww_dst, we_dst.+nₗ, wg_dst.+(nₗ+nₑ))
    wgt = vcat(ones(mₗ)./nwn.Rⱼ,ones(mₑ+m₀)./nwn.Rₑ)
    sp = sparse(vcat(src,dst), vcat(dst,src), vcat(wgt,wgt),n,n,+)
    nwn.graph=SimpleWeightedGraph(sp)
    nwn.ls = ones(Int,n)
    nwn.ns = cumsum(nwn.ls)
    nwn
end

function update_graph(nwn::NWN_circuit{:MNR,T,N}) where {T,N}
    nₗ, nₑ, n₀ = length(nwn.lines), length(nwn.elecs), length(nwn.grnds)
    n = nₗ + nₑ + n₀
    line_ip = SVector{n}([T[] for i ∈ 1:n])
    juncts = NamedTuple{(:R,:c₁,:c₂),Tuple{T,Tuple{Int,Int},Tuple{Int,Int}}}[]
    each_intersect(nwn.lines, nwn.dims, []) do i₁,i₂,arr,dims,ps,p,arg
        push!(line_ip[i₁],ps[1])
        push!(line_ip[i₂],ps[2])
        push!(juncts,(R=nwn.Rⱼ,c₁=(i₁,length(line_ip[i₁])),c₂=(i₂,length(line_ip[i₂]))))
    end
    each_intersect(nwn.lines, vcat(nwn.elecs,nwn.grnds), nwn.dims, []) do i₁,i₂,arr₁,arr₂,dims,ps,p,arg
        i₂ += nₗ
        push!(line_ip[i₁],ps[1])
        push!(line_ip[i₂],ps[2])
        push!(juncts,(R=nwn.Rₑ,c₁=(i₁,length(line_ip[i₁])),c₂=(i₂,length(line_ip[i₂]))))
    end
    nwn.ls = length.(line_ip) # number of intersections on each wire
    nwn.ns = cumsum(nwn.ls) # number of nodes up to wire i
    # to count number of nodes involved in inner-wire resistances:
    # if ls[i]==0 then line i does not add to the node count or add a ressitance
    # if ls[1]==1 then line i has no inner wire resistances; -1 to nodes in inner-wires
    # if ls[1]>1 then line i has i-1 inner wire resistances; -1 to nodes in inner-wires
    # ⟹ if ls[i]>0 then -1 to nodes in inner-wires
    ms = vcat([0],cumsum(nwn.ls[1:nₗ].>0)) 
    nᵢ = nwn.ns[nₗ]-ms[end]
    nⱼ = length(juncts)
    perms = sortperm.(line_ip)
    invperms = invperm.(perms)
    permute!.(line_ip,perms)

    src, dst, wgt = zeros(Int,nᵢ+nⱼ), zeros(Int,nᵢ+nⱼ), zeros(T,nᵢ+nⱼ)
    for i ∈ 1:nₗ
        js = 1:(nwn.ls[i]-1)
        slice = js.+(nwn.ns[i]-nwn.ls[i]-ms[i])
        Rₗ = 4*len(nwn.lines[i])*nwn.props[i].ρ/(π*nwn.props[i].D*nwn.props[i].D)
        src[slice] = js.+(nwn.ns[i]-nwn.ls[i])
        dst[slice] = js.+(1+nwn.ns[i]-nwn.ls[i]) 
        wgt[slice] = 1 ./(Rₗ.*(diff(line_ip[i])))
    end
    for (i,ju) ∈ enumerate(juncts)
        src[i+nᵢ] = nwn.ns[ju.c₁[1]]-nwn.ls[ju.c₁[1]]+invperms[ju.c₁[1]][ju.c₁[2]]
        dst[i+nᵢ] = nwn.ns[ju.c₂[1]]-nwn.ls[ju.c₂[1]]+invperms[ju.c₂[1]][ju.c₂[2]]
        wgt[i+nᵢ] = 1/ju.R
    end
    sp = sparse(vcat(src,dst), vcat(dst,src), vcat(wgt,wgt),nwn.ns[end],nwn.ns[end],+)
    nwn.graph = SimpleWeightedGraph(sp)
    nwn
end




"""
Creates the system of equations which solves the 
electrical problem of the nanowire network under 
the JDA. Removes ill conditioned nodes (those not 
connected through the circuit to a voltage source 
or ground).
"""
function circuit_eqs(nwn::NWN_circuit{M,T,N}) where {M,T,N}
    nₗ, nₑ, n₀ = length(nwn.lines), length(nwn.elecs), length(nwn.grnds)
    nₗₙ = nwn.ns[nₗ] # number of line related nodes
    nₑₙ = nwn.ns[nₗ + nₑ] - nₗₙ
    # n₀ₙ = nwn.ns[nₗ + nₑ + n₀] - nₗₙ - nₑₙ

    # include only nodes part of the actual circuit
    con_comps = connected_components(nwn.graph)
    connected = vcat(con_comps[.!(con_comps .⊆ [1:nₗₙ])]...)

    # remove grounds from equations as MNA requires
    inds = setdiff(connected, (nwn.ns[nₗ+nₑ]+1):nwn.ns[nₗ+nₑ+n₀])

    # view of Laplacian of graph with indices needed
    L = -laplacian_matrix(nwn.graph)[inds,inds]

    # generates B matrix
    B = spzeros(T,length(inds),nₑₙ)
    vs = spzeros(T,nₑₙ)
    for i ∈ 1:nₑ
        js = (nwn.ns[nₗ+i]-nwn.ls[nₗ+i]+1):nwn.ns[nₗ+i]
        B[CartesianIndex.(findall(inds.∈[js]),js.-nₗₙ)] .= one(T)
        vs[js.-nₗₙ] .= nwn.volts[i]
    end

    # generates system of equations in the form Ax=z
    A = [L B;
        transpose(B) spzeros(T,nₑₙ,nₑₙ)]
    z = [spzeros(T,length(inds));
        vs]
    A,z,inds
end

"""
Calculates the sheet resistance of a NWN object. 
Solves the network using MNA and returns the ratio 
of the voltage to the current. Requires the network 
to have only one electrode.
"""
function sheet_resistance(nwn::NWN_circuit)
    @assert length(nwn.elecs)==1 "Sheet resistance is only defined for one electrode and one ground."
    A, z, inds = circuit_eqs(nwn)
    sol = A\collect(z)
    nwn.volts[1]/sum(sol[length(inds)+1:end])
end

function sheet_resistance(nwn::NWN_circuit{M,BigFloat,N}) where {M,N}
    @assert length(nwn.elecs)==1 "Sheet resistance is only defined for one electrode and one ground."
    A, z, inds = circuit_eqs(nwn)
    sol = collect(A)\collect(z)
    nwn.volts[1]/sol[end]
end

"""
Returns a boolean value for if the given network 
has a path from at least one electrode and ground.
"""
function is_conducting(nwn::NWN_circuit)
    nₗ, nₑ = length(nwn.lines), length(nwn.elecs)
    for i ∈ (nwn.ns[nₗ]+1):nwn.ns[nₗ+nₑ]
        for j ∈ (nwn.ns[nₗ+nₑ]+1):nwn.ns[end]
            if has_path(nwn.graph,i,j)
                return true
            end
        end
    end
    return false
end
