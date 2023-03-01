export intersection_params, line_segment_point, each_intersect, find_juncts, find_connect, connections

"""
```julia
intersection_params(w₁::Line, w₂::Line)
```
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
```julia
line_segment_point(t::T, w::Line)
```
Computes the point on the line `w` at parameter `t` where 
`t` is the parameter of the line so that `p₂*t+(1-t)p₁` 
is the parametric equation of the line. 
"""
function line_segment_point(t::T, w::Line) where T<:Real 
    w.p₂*t+(1-t)*w.p₁
end

"""
```julia
each_intersect(f!,arr₁,arr₂,dims,arg)
```
    Executes `f!(i₁,i₂,arr₁,arr₂,dims,ps,p,arg)` for any intersection between pairs of 
    lines in `arr₁` and `arr₂`. `i₁` and `i₂` are the indices of the lines in `arr₁` 
    and `arr₂` respectively, `ps` is the vector of line intersection parameters as 
    computed by `intersection_params(arr₁[i₁],arr₂[i₂])`, and `p` is the point at which 
    the intersection occurs.
"""
function each_intersect(f!, arr₁::Array{Line{S,N}}, arr₂::Array{Line{S,N}}, dims::SVector{N,S}, arg) where {S,N}
    for i₁ ∈ eachindex(arr₁)
        for i₂ ∈ eachindex(arr₂)
            ps = intersection_params(arr₁[i₁],arr₂[i₂])
            if (0.0 ≤ ps[1] ≤ 1.0) && (0.0 ≤ ps[2] ≤ 1.0)
                p = line_segment_point(ps[1],arr₁[i₁])
                if all( zero(SVector{N, S}) .≤ p .≤ dims )
                    f!(i₁,i₂,arr₁,arr₂,dims,ps,p,arg)
                end
            end
        end
    end
    arg
end

"""
```julia
each_intersect(f!,arr,dims,arg)
```
    Executes `f!(i₁,i₂,arr,dims,ps,p,arg)` for any intersection between pairs of 
    lines in `arr`. `i₁` and `i₂` are the indices of the lines in `arr`, `ps` is 
    the vector of line intersection parameters as computed by 
    `intersection_params(arr[i₁],arr[i₂])`, and `p` is the point at which the 
    intersection occurs.
"""
function each_intersect(f!, arr::Array{Line{S,N}}, dims::SVector{N,S}, arg) where {S,N}
    n = length(arr)
    for i₁ ∈ 2:(n)
        for i₂ ∈ 1:(i₁-1)
            ps = intersection_params(arr[i₁],arr[i₂])
            if (0.0 ≤ ps[1] ≤ 1.0) && (0.0 ≤ ps[2] ≤ 1.0)
                p = line_segment_point(ps[1],arr[i₁])
                if all( zero(SVector{N, S}) .≤ p .≤ dims )
                    f!(i₁,i₂,arr,dims,ps,p,arg)
                end
            end
        end
    end
    arg
end

"""
```julia
find_juncts(net)
```
Given an `StickNetwork` object finds all points of intersection between lines. 
"""
function find_juncts(net::StickNetwork{T,N}) where {T,N}
    each_intersect(net.lines, net.dims, Array{Tuple{SVector{N,T},Tuple{Int,Int}},1}()) do i₁,i₂,arr,dims,ps,p,juncts
        push!(juncts,(p,(i₁,i₂)))
    end
end

"""
```julia
find_connect(net)
```
Given an `StickNetwork` object finds all pairs `(i₁,i₂)` with an intersection between lines `i₁` and `i₂` in `net`.
"""
function find_connect(net::StickNetwork{T,N}) where {T,N}
    each_intersect(net.lines, net.dims, Array{Tuple{Int,Int},1}()) do i₁,i₂,arr,dims,ps,p,connects
        push!(connects,(i₁,i₂))
    end
end

"""
```julia
connections(arr₁,arr₂,dims)
```
Finds all connections between arrays of `Line`s `arr₁` and `arr₂`. 
Returns `Tuple{Vector{Int},Vector{Int}}` with each vector representing the source and destination line indices for each connection.
"""
function connections(arr₁::Array{Line{S,N}}, arr₂::Array{Line{S,N}}, dims::SVector{N,S}) where {S,N}
    each_intersect(arr₁, arr₂, dims, (Int[],Int[])) do i₁,i₂,arr₁,arr₂,dims,ps,p,connects
        push!(connects[1],i₁)
        push!(connects[2],i₂)
    end
end

"""
```julia
connections(arr,dims)
```
Finds all connections between `Line`s in `arr`. 
Returns `Tuple{Vector{Int},Vector{Int}}` with each vector representing the source and destination line indices for each connection.
"""
function connections(arr::Array{Line{S,N}}, dims::SVector{N,S}) where {S,N}
    each_intersect(arr, dims, (Int[],Int[])) do i₁,i₂,arr,dims,ps,p,connects
        push!(connects[1],i₁)
        push!(connects[2],i₂)
    end
end