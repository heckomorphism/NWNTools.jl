export JDA_graph1, JDA_graph2, JDA_graph3

using LightGraphs, SimpleWeightedGraphs, SparseArrays

swap_index(a) = collect(broadcast(x->getindex(x,i),a) for i ∈ eachindex(a[1]))

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