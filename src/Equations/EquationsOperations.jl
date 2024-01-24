using FillArrays
using LinearAlgebra
using Parameters

val(x) = x
val(x::Gridap.Fields.ForwardDiff.Dual) = x.value
conv(u, ∇u) = (∇u') ⋅ u

"""
    τsu(u, h, ν::Real, dt::Real)

Stabilization parameter for SUPG formulation.
Janssens, B. (2014). Numerical modeling and experimental investigation of ﬁne particle coagulation and dispersion in dilute ﬂows.
"""
function τsu(u, h, ν::Real, dt::Real)

    function τsu(u, h)
        r = 1
        τ₂ = h^2 / (4 * ν)
        τ₃ = dt / 2

        u = val(norm(u))
        if iszero(u)
            return (1 / τ₂^r + 1 / τ₃^r)^(-1 / r)
        end
        τ₁ = h / (2 * u)
        return (1 / τ₁^r + 1 / τ₂^r + 1 / τ₃^r)^(-1 / r)

    end

    return τsu ∘ (u, h)
end

"""
    τb(u, h, ν::Real, dt::Real)

Stabilization parameter for SUPG formulation.
Janssens, B. (2014). Numerical modeling and experimental investigation of ﬁne particle coagulation and dispersion in dilute ﬂows.
"""
function τb(u, h, ν::Real, dt::Real)
    return (u ⋅ u) * τsu(u, h, ν, dt)
end


"""
    τm(uu, params::Dict{Symbol,Any})

Stabilization parameter for VMS formulation.
Bazilevs, Y., Calo, V. M., Cottrell, J. A., Hughes, T. J. R., Reali, A., & Scovazzi, G. (2007). Variational multiscale residual-based turbulence modeling for large eddy simulation of incompressible flows. Computer Methods in Applied Mechanics and Engineering, 197(1–4), 173–201. https://doi.org/10.1016/j.cma.2007.07.016
"""
function τm(uu, params::Dict{Symbol,Any})
    
    @unpack  G,GG,gg,Cᵢ,ν,dt = params

    function τm(uu, G, GG)
        τ₁ = Cᵢ[1] * (2 / dt)^2 #Here, you can increse the 2 if CFL high
        τ₃ = Cᵢ[2] * (ν^2 * GG)




        uu_new = VectorValue(val.(uu)...)

        if iszero(norm(uu_new))
            return (τ₁ .+ τ₃) .^ (-1 / 2)
        end

        τ₂ = uu_new ⋅ G ⋅ uu_new
        return (τ₁ .+ τ₂ .+ τ₃) .^ (-1 / 2)
    end

    return τm ∘ (uu, G, GG)
end


"""
    τc(uu, params::Dict{Symbol,Any})

Stabilization parameter for VMS formulation.
Bazilevs, Y., Calo, V. M., Cottrell, J. A., Hughes, T. J. R., Reali, A., & Scovazzi, G. (2007). Variational multiscale residual-based turbulence modeling for large eddy simulation of incompressible flows. Computer Methods in Applied Mechanics and Engineering, 197(1–4), 173–201. https://doi.org/10.1016/j.cma.2007.07.016
"""
function τc(uu, params::Dict{Symbol,Any})
 @unpack   gg = params
    return 1 / (τm(uu,params) ⋅ gg)
end