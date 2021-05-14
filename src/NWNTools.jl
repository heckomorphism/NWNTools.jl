module NWNTools

    export Line, Wire, NWN
    export len, rand_wires, rand_network, ensemble
    export intersection_params, line_segment_point, find_juncts
    export ρ_junct, ρ_junct_avg, ρ_junct_mean_std


    using Random
    using StaticArrays
    using SparseArrays
    using Distributions: Normal
    using Statistics: mean, var, std

    """
    A structure representing a line with end points 
    p₁ and p₂. Parameterized by the datatype of the 
    coordinates and the length (i.e. dimension) N. 
    """
    struct Line{T<:Real,N}
        p₁::SVector{N,T}
        p₂::SVector{N,T}
    end

    import Base.convert
    Line(p₁::Array{T,1},p₂::Array{S,1}) where {T,S} = Line(promote(SVector{length(p₁)}(p₁), SVector{length(p₂)}(p₂))...)
    convert(::Type{Line{T,N}},x::Line{S,N}) where {T,S,N} = Line(convert(SVector{N,T},x.p₁), convert(SVector{N,T},x.p₂))

    """
    Determines the Euclidean length of the line segment ℓ.
    """
    len(ℓ::Line{T,N}) where {T,N} = sqrt(sum((ℓ.p₁[i]-ℓ.p₂[i])^2 for i ∈ 1:N))

    """
    A structure to hold all of the information relevant to a wire.
    """
    struct Wire{T<:Real,N}
        line::Line{T,N}
        length::Float64
        # d::Float64
        # ρ::Float64
    end

    """
    Constructor the Wire struct which only requires a Line object.
        """
    Wire(ℓ) = Wire(ℓ,len(ℓ))

    """
    A structure to store a collection of wires in a certain size. 
    """
    struct NWN{T<:Real,N}
        wires::Array{Wire{T,N},1}
        dims::SVector{N,T}
    end


    """
    Generates a random collection of wires with initial point 
    laying within the rectangle between the origin and `dims`. 
    Determines second point of line to be distance l and random 
    direction (as sampled unifromly from the N-1 sphere) from 
    the initial point. 
    """
    function rand_wires(dims::SVector{N,T},n::Int,l) where {N,T}
        x₀s = map(x->x.*dims,rand(SVector{N,T},n))
        A = rand(Normal(),N,n)
        # normalize columns of A to leave columns 
        # as uniformly sampled on the N-1 sphere
        A .= A./sqrt.(sum(A.^2,dims=1)).*l
        x₁s = [SVector{N}(c) for c ∈ eachcol(A)]   
        ls = Line.(x₀s,x₀s+x₁s)
        Wire.(ls)
    end

    function rand_network(dims::SVector{N,T},ρ,l) where {N,T}
        NWN(rand_wires(dims, round(Int,ρ*prod(dims)),l),dims)
    end

    """
    Creates an ensemble of nanowire networks with the given 
    paramaters. 
    """
    function ensemble(n, dims::SVector{N,T},ρ,l) where {N,T}
        [rand_network(dims, ρ, l) for i∈1:n]
    end
    # function ensemble(n, dims::SVector{N,T},ρ,l) where {N,T}
    #     # number of wires in each network
    #     m = round(Int,ρ*prod(dims)) 
    #     # generate m*n wires
    #     wires = rand_wires(dims, round(Int,ρ*prod(dims))*n, l)
    #     # make n networks with m wires each
    #     [NWN(collect(ws),dims) for ws in Iterators.partition(wires,m)]

    # end

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
    Given an array of line objects finds all points of 
    intersection between lines. 
    """
    # function find_juncts(lines::Array{Line{T,N},1}) where {T,N}
    #     juncts = Array{Tuple{SVector{N,T},Tuple{Int,Int}},1}()
    #     inds = eachindex(lines)
    #     for i ∈ inds, j ∈ inds
    #         ps = intersection_params(lines[i],lines[j])
    #         if (0.0 ≤ ps[1] ≤ 1.0) && (0.0 ≤ ps[2] ≤ 1.0)
    #             push!(juncts,(line_segment_point(ps[1],lines[i]),(i,j)))
    #         end
    #     end
    #     juncts
    # end
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

    function ρ_junct(nwn::NWN{T,N}) where {T,N}
        n = length(nwn.wires)
        m = 0
        for i ∈ 2:(n)
            for j ∈ 1:(i-1)
                ps = intersection_params(nwn.wires[i].line,nwn.wires[j].line)
                if (0.0 ≤ ps[1] ≤ 1.0) && (0.0 ≤ ps[2] ≤ 1.0)
                    p = line_segment_point(ps[1],nwn.wires[i].line)
                    if sum( SA[0.0,0.0] .≤ p .≤ nwn.dims ) == N
                        m += 1
                    end
                end
            end
        end
        m/prod(nwn.dims)
    end

    ρ_junct_avg(dims, ρ, l; N=32) = mean(ρ_junct.(ensemble(N,dims, ρ, l)))
    function ρ_junct_mean_std(dims, ρ, l; N=32)
        e = ensemble(N,dims, ρ, l)
        ρ_js = ρ_junct.(e)
        mean(ρ_js), std(ρ_js)
    end

end