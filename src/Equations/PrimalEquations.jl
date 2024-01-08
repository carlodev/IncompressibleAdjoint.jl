####################################################################
#STEADY PRIMAL
####################################################################
function eq_primal_steady(params::Dict{Symbol,Any})
    @unpack method, D,dΩ= params
    if method == :VMS
        res_prim = primal_steady_VMS(params)
    elseif method == :SUPG
        res_prim = primal_steady_SUPG(params)
    end

    rhs((v, q)) = ∫(VectorValue(zeros(D)...) ⋅ v)dΩ

    return res_prim, rhs
end


function primal_steady_SUPG(params::Dict{Symbol,Any})

    @unpack ν, dt, dΩ, D, Ω, θ,uh = params
    h = h_param(Ω, D)
    updatekey(params, :h,h)



    a((u, p), (v, q)) = ∫(ν * ∇(v) ⊙ ∇(u) + ∇(p) ⊙ v + q * (∇ ⋅ u))dΩ + ∫(v ⊙ (conv ∘ (uh, ∇(u))))dΩ

    Rm(u, p) = conv ∘ (uh, ∇(u)) + ∇(p)
    astab((u, p), (v, q)) = ∫( τsu(uh, h,ν,dt) ⋅ (uh ⋅ ∇(v) + ∇(q)) ⊙ Rm(u, p))dΩ + ∫( τb(uh, h,ν,dt) ⋅ (∇ ⋅ v) ⊙ (∇ ⋅ u))dΩ

    res_prim((u, p), (v, q)) = a((u, p), (v, q)) + astab((u, p), (v, q))

    return res_prim

end

function primal_steady_VMS(params::Dict{Symbol,Any})
    @error("Not implemented yet")

end




####################################################################
#UNSTEADY PRIMAL
####################################################################
function eq_primal_unsteady(params::Dict{Symbol,Any})
    @unpack method, D,dΩ = params
    if method == :VMS
        m, res_prim = primal_unsteady_VMS(params)
    elseif method == :SUPG
        m, res_prim = primal_unsteady_SUPG(params)
    end

    rhs(t, (v, q)) = ∫(VectorValue(zeros(D)...) ⋅ v)dΩ

    return m, res_prim, rhs
end


function primal_unsteady_SUPG(params::Dict{Symbol,Any})

    @unpack ν, dt, dΩ, D, Ω, θ,uh = params
    h = h_param(Ω, D)
    updatekey(params, :h,h)
    θvp = get_θvp(θ)
    
    a(t, (u, p), (v, q)) = ∫(ν * ∇(v) ⊙ ∇(u) + ∇(p) ⊙ v + q * (∇ ⋅ u))dΩ + ∫(v ⊙ (conv ∘ (uh, ∇(u))))dΩ

    Rm(u, p) = conv ∘ (uh, ∇(u)) + ∇(p)
    astab(t, (u, p), (v, q)) = ∫( τsu(uh, h,ν,dt) ⋅ (uh ⋅ ∇(v) + ∇(q)) ⊙ Rm(u, p))dΩ + ∫( τb(uh, h,ν,dt) ⋅ (∇ ⋅ v) ⊙ (∇ ⋅ u))dΩ

    res_prim(t, (u, p), (v, q)) = a(t, (u, p), (v, q)) + astab(t, (u, p), (v, q))

    m(t, (u, p), (v, q)) = ∫(u ⋅ v)dΩ + ∫( τsu(uh, h,ν,dt) ⋅ (uh ⋅ ∇(v) + θvp * ∇(q)) ⊙ u)dΩ


    return m, res_prim

end

function primal_unsteady_VMS(params::Dict{Symbol,Any})
    @error("Not implemented yet")
end

function get_θvp(θ)
    θvp = 1.0

    if length(θ) == 2
        θv, θp = θ
        θvp = θv / θp
    end
    return θvp

end

