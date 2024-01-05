val(x) = x
val(x::Gridap.Fields.ForwardDiff.Dual) = x.value
conv(u, ∇u) = (∇u') ⋅ u

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

function τb(u, h, ν::Real, dt::Real)
    return (u ⋅ u) * τsu(u, h, ν, dt)
end


function τm(uu, G, GG, ν, dt)

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



function τc(uu, gg, G, GG, ν, dt)
    return 1 / (τm(uu, G, GG, ν, dt) ⋅ gg)
end