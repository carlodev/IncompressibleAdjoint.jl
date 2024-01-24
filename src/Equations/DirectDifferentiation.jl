function eq_direct_differentiation_steady(params::Dict{Symbol,Any})
    @unpack method, D,dΩ= params
    if method == :VMS
        #
    elseif method == :SUPG
        res_prim = direct_differentiation_SUPG(params)
    end

    rhs((v, q)) = ∫(VectorValue(zeros(D)...) ⋅ v)dΩ

    return res_prim, rhs
end

function direct_differentiation_SUPG(params::Dict{Symbol,Any})

    @unpack ν, dt, dΩ, D, Ω, θ,uh = params
    h = h_param(Ω, D)
    updatekey(params, :h,h)
    
    Rm(ub, pb) = transpose(∇(ub))⋅uh+∇(uh)⋅ub+∇(pb)

    a((ub, pb), (v, q)) = ∫(v⊙Rm(ub, pb))dΩ + ∫(q * (∇ ⋅ ub) )dΩ + ∫(ν * ∇(v) ⊙ ∇(ub) )dΩ

    astab((u, p), (v, q)) = ∫( τsu(uh, h,ν,dt) ⋅ (uh ⋅ ∇(v) + ∇(q)) ⊙ Rm(u, p))dΩ + ∫( τb(uh, h,ν,dt) ⋅ (∇ ⋅ v) ⊙ (∇ ⋅ u))dΩ

    res_prim((u, p), (v, q)) = a((u, p), (v, q)) + astab((u, p), (v, q))

    return res_prim

end