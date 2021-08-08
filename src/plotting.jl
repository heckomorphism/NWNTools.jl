using RecipesBase

@recipe function f(lines::Vector{Line{T,2}}) where {T}
    xs = hcat(([line.p₁[1],line.p₂[1]] for line ∈ lines)...)
    ys = hcat(([line.p₁[2],line.p₂[2]] for line ∈ lines)...)
    linecolor --> :black
    legend := false
    xs,ys
end

@recipe function f(nwn::NWN_circuit)
    legend -> false
    grid -> false
    seriestype := :line
    @series begin
        linecolor := :black
        label := "Wires"
        nwn.lines
    end
    @series begin
        linecolor := :red
        label := "Electrode"
        linewidth --> 2
        nwn.elecs
    end
    @series begin
        linecolor := :blue
        label := "Ground"
        linewidth --> 2
        nwn.grnds
    end
end

@recipe function f(net::StickNetwork)
    legend -> false
    grid -> false
    seriestype := :line
    @series begin
        linecolor := :black
        net.lines
    end
end