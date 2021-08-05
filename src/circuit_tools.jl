export circuit_eqs, sheet_resistance, currents

"""
```julia
circuit_eqs(nwn::NWN_circuit{M,T,N})
```
Creates the system of equations which solves the 
electrical problem of the nanowire network under 
the given model `M`. Removes ill conditioned nodes (those not 
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
```julia
sheet_resistance(nwn::NWN_circuit)
```
Calculates the sheet resistance of a nanowire network object. 
Solves the network using MNA and returns the ratio of the voltage 
to the current. Requires the network to have only one electrode.
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
    # For BigFloat the matrix must be dense
    sol = collect(A)\collect(z) 
    nwn.volts[1]/sol[end]
end

"""
```julia
currents(nwn::NWN_circuit{M,T,N})
```
Returns a list of currents as iterated by `edges(nwn.graph)` and using thhe node labeling as defined by the model `M`.
"""
function currents(nwn::NWN_circuit{M,T,N}) where {M,T,N}
	A,z,inds = circuit_eqs(nwn)
	sol = A\collect(z)
	Vs = zeros(T,nv(nwn.graph))
	Vs[inds] = sol[1:length(inds)]
	es = edges(nwn.graph)
	# (Vₒᵤₜ-Vᵢₙ)/R
	(Vs[dst.(es)].-Vs[src.(es)]).*weight.(es)
end