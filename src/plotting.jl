using RecipesBase

@recipe function f(lines::Vector{Line{T,2}}) where {T}
    xs = hcat(([line.p₁[1],line.p₂[1]] for line ∈ lines)...)
    ys = hcat(([line.p₁[2],line.p₂[2]] for line ∈ lines)...)
    linecolor --> :black
    legend := false
    xs,ys
end

@recipe function f(nwn::NWN_JDA{S,T,N}) where {S,T,N}
    legend := false
    grid := false
    seriestype := :line
    @series begin
        linecolor := :black
        nwn.lines
    end
    @series begin
        linecolor := :red
        linewidth --> 2
        nwn.elecs
    end
    @series begin
        linecolor := :blue
        linewidth --> 2
        nwn.grnds
    end
end