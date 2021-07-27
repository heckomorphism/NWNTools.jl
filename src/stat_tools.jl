export emsd_f, emsd_f_JDA, ρ_junct, ρ_junct_avg, ρ_junct_mean_std

using Statistics: mean, var, std

"""
Computes the mean and standard deviation of an observable 
f(nw) for nw a NWN.

Parameters:
    f:      Function, takes a NWN as the first argument 
            and returns an object which the mean and 
            standard deviation make sense for.
    params: Dimensions, wire density, and wire length.
    N:      Number of networks in the ensemble.
"""
function emsd_f(f, params; N=32)
    e = ensemble(N, params)
    res = f.(e)
    mean(res), std(res)
end

function emsd_f_circuit(f, params, M, args...; N=32)
    e = ensemble_circuit(N, params..., M, args...)
    res = f.(e)
    mean(res), std(res)
end

"""
Computes junction density in a simiar algorithm to `find_juncts`.
"""
function ρ_junct(nwn::NWN{T,N}) where {T,N}
    n = length(nwn.wires)
    m = 0
    for i ∈ 2:(n)
        for j ∈ 1:(i-1)
            ps = intersection_params(nwn.wires[i].line,nwn.wires[j].line)
            if (0.0 ≤ ps[1] ≤ 1.0) && (0.0 ≤ ps[2] ≤ 1.0)
                p = line_segment_point(ps[1],nwn.wires[i].line)
                if sum( zeros(SVector{N,T}) .≤ p .≤ nwn.dims ) == N
                    m += 1
                end
            end
        end
    end
    m/prod(nwn.dims)
end

"""
Computes the mean and standard deviation of the junction 
density for a certain configuration of NWN parameters.
"""
ρ_junct_mean_std(params; N=32) = emsd_f(ρ_junct, params; N=N)
ρ_junct_mean_std(params...; N=32) = ρ_junct_mean_std(params; N=N)

"""
Deprecated, use `ρ_junct_mean_std` instead.
"""
ρ_junct_avg(dims, ρₙ, l, ρ, D; N=32) = emsd_f(ρ_junct, [dims, ρₙ, l, ρ, D], N=N)[1]