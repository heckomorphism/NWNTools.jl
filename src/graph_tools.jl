export JDA_graph1, JDA_graph2, JDA_graph3

using LightGraphs, SimpleWeightedGraphs, SparseArrays

swap_index(a) = collect(broadcast(x->getindex(x,i),a) for i ∈ eachindex(a[1]))

# Test JDA Laplacian with no terminals

function JDA_graph1(nwn::NWN{T,N}, Rⱼ) where {T,N}
    n = length(nwn.wires)
    g = SimpleWeightedGraph(n)
    for i ∈ 2:(n)
        for j ∈ 1:(i-1)
            ps = intersection_params(nwn.wires[i].line,nwn.wires[j].line)
            if (0.0 ≤ ps[1] ≤ 1.0) && (0.0 ≤ ps[2] ≤ 1.0)
                p = line_segment_point(ps[1],nwn.wires[i].line)
                if sum( zeros(SVector{N,T}) .≤ p .≤ nwn.dims ) == N
                    add_edge!(g,i,j,1/Rⱼ)
                end
            end
        end
    end
    return laplacian_matrix(g)
end

function JDA_graph2(nwn::NWN{T,N}, Rⱼ) where {T,N}
    connects = find_connect(nwn)
    g = SimpleWeightedGraph(swap_index(connects)..., repeat([1/Rⱼ],length(connects)))
    return laplacian_matrix(g)
end

function JDA_graph3(nwn::NWN{T,N}, Rⱼ) where {T,N}
    n = length(nwn.wires)
    L = spzeros(n,n)
    for i ∈ 2:(n)
        for j ∈ 1:(i-1)
            ps = intersection_params(nwn.wires[i].line,nwn.wires[j].line)
            if (0.0 ≤ ps[1] ≤ 1.0) && (0.0 ≤ ps[2] ≤ 1.0)
                p = line_segment_point(ps[1],nwn.wires[i].line)
                if sum( zeros(SVector{N,T}) .≤ p .≤ nwn.dims ) == N
                    L[i,j] = L[j,i] = 1/Rⱼ
                end
            end
        end
    end

    for i ∈ 1:n
        L[i,i] = -sum(L[:,i])
    end
    return L
end

# Try implementing terminals

function JDA_graph_elecs1(nwn::NWN_electrodes)
    ww_src, ww_dst = swap_index(connections(nwn.lines, nwn.dims))
    we_src, we_dst = swap_index(connections(nwn.lines, nwn.elecs, nwn.dims))
    wg_src, wg_dst = swap_index(connections(nwn.lines, nwn.grnds, nwn.dims))
    nₗ, nₑ = length(nwn.lines), length(nwn.elecs)
    mₗ, mₑ, m₀ = length(ww_src), length(we_src), length(wg_src)
    src = vcat(ww_src, we_src, wg_src)
    dst = vcat(ww_dst, we_dst.+nₗ, wg_dst.+nₗ.+nₑ)
    wgt = 1 ./ vcat(repeat([nwn.Rⱼ],mₗ),repeat([nwn.Rₑ],mₑ+m₀))
    g = SimpleWeightedGraph(src, dst, wgt)
    laplacian_matrix(g)
end

    