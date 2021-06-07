export intersection_params, line_segment_point, find_juncts, find_connect

"""
Given two lines computes the parameters `[t,s]` so that 
the intersection of the lines occurs at 
`(w₁.p₂-w₁.p₁)t+w₁.p₁==(w₂.p₂-w₂.p₁)s+w₂.p₁`. If the 
lines are parallel, one is a point and not a line, or 
the lines are colinear the system of equations will 
not be consitent and so will return some solution 
like `[NaN,NaN]` or `[Inf,Inf]`. The lines intersect on 
the segments given whenever both `0≤t≤1` and `0≤s≤1`.
"""
function intersection_params(w₁::Line, w₂::Line)
    hcat(w₁.p₂-w₁.p₁,w₂.p₁-w₂.p₂)\(w₂.p₁-w₁.p₁)
end

"""
Computes the point on the line w at parameter `t` where 
`t` is the parameter of the line so that `p₂*t+(1-t)p₁` 
is the parametric equation of the line. 
"""
function line_segment_point(t::T, w::Line) where T<:Real 
    w.p₂*t+(1-t)*w.p₁
end

"""
Given an NWN object finds all points of 
intersection between wires. 
"""
function find_juncts(nwn::NWN{T,N}) where {T,N}
    juncts = Array{Tuple{SVector{N,T},Tuple{Int,Int}},1}()
    n = length(nwn.wires)
    for i ∈ 2:(n)
        for j ∈ 1:(i-1)
            ps = intersection_params(nwn.wires[i].line,nwn.wires[j].line)
            if (0.0 ≤ ps[1] ≤ 1.0) && (0.0 ≤ ps[2] ≤ 1.0)
                p = line_segment_point(ps[1],nwn.wires[i].line)
                if sum( SA[0.0,0.0] .≤ p .≤ nwn.dims ) == N
                    push!(juncts,(p,(i,j)))
                end
            end
        end
    end
    juncts
end

"""
Given an NWN object finds all connections between wires. 
"""
function find_connect(nwn::NWN{T,N}) where {T,N}
    connects = Array{Tuple{Int,Int},1}()
    n = length(nwn.wires)
    for i ∈ 2:(n)
        for j ∈ 1:(i-1)
            ps = intersection_params(nwn.wires[i].line,nwn.wires[j].line)
            if (0.0 ≤ ps[1] ≤ 1.0) && (0.0 ≤ ps[2] ≤ 1.0)
                p = line_segment_point(ps[1],nwn.wires[i].line)
                if sum( SA[0.0,0.0] .≤ p .≤ nwn.dims ) == N
                    push!(connects,(i,j))
                end
            end
        end
    end
    connects
end