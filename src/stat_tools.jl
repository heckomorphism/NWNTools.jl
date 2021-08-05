export emsd_f, emsd_f_circuit, ρ_junct, ρ_junct_avg, ρ_junct_mean_std

using Statistics: mean, var, std

"""
```julia
emsd_f(f, params; N=32)
```
Computes the mean and standard deviation of an ensemble of `StickNetwork`s determined by `params` with ensemble size `N`.
"""
function emsd_f(f, params; N=32)
    e = ensemble(N, params...)
    res = f.(e)
    mean(res), std(res)
end

"""
```julia
emsd_f_circuit(f, params, M, args...; N=32)
```
Computes the mean and standard deviation of an ensemble of `N` nanowire networks 
with model `M` and parameters 'params'. `args` may specify additonal parameters to 
`ensemble_circuit` such as the `Val(:conducting)` property.
"""
function emsd_f_circuit(f, params, M, args...; N=32)
    e = ensemble_circuit(N, params..., M, args...)
    res = f.(e)
    mean(res), std(res)
end

"""
```julia
ρ_junct(net)
```
Computes junction density of `StickNetwork` `net`.
"""
function ρ_junct(net::StickNetwork{T,N}) where {T,N}
    n = length(net.lines)
    m = 0
    for i ∈ 2:(n)
        for j ∈ 1:(i-1)
            ps = intersection_params(net.lines[i],net.lines[j])
            if (0.0 ≤ ps[1] ≤ 1.0) && (0.0 ≤ ps[2] ≤ 1.0)
                p = line_segment_point(ps[1],net.lines[i])
                if sum( zeros(SVector{N,T}) .≤ p .≤ net.dims ) == N
                    m += 1
                end
            end
        end
    end
    m/prod(net.dims)
end

"""
```julia
ρ_junct_mean_std(params; N=32)
```
Computes the mean and standard deviation of the junction 
density of an ensemble of `StickNetwork`s given by `params`.
"""
ρ_junct_mean_std(params; N=32) = emsd_f(ρ_junct, params; N=N)
ρ_junct_mean_std(params...; N=32) = ρ_junct_mean_std(params; N=N)

"""
Deprecated, use `ρ_junct_mean_std` instead.
"""
ρ_junct_avg(dims, ρₙ, l, ρ, D; N=32) = emsd_f(ρ_junct, [dims, ρₙ, l, ρ, D], N=N)[1]