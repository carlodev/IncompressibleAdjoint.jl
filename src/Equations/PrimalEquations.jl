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


"""
    primal_steady_SUPG(params::Dict{Symbol,Any})

Navier-Stokes SUPG stabilized set of equations
"""
function primal_steady_SUPG(params::Dict{Symbol,Any})

    @unpack ν, dt, dΩ, D, Ω, θ,uh = params
    h = h_param(Ω, D)
    updatekey(params, :h,h)
    

    a((u, p), (v, q)) = ∫(ν * ∇(v) ⊙ ∇(u) + ∇(p) ⊙ v + q * (∇ ⋅ u))dΩ + ∫(v ⊙ (conv ∘ (uh, ∇(u))))dΩ

    Rm(u, p) = conv ∘ (uh, ∇(u)) + ∇(p)
    astab((u, p), (v, q)) = ∫( τsu(uh, h,ν,dt) ⋅ ((conv ∘ (uh, ∇(v))) + ∇(q)) ⊙ Rm(u, p))dΩ + ∫( τb(uh, h,ν,dt) ⋅ (∇ ⋅ v) ⊙ (∇ ⋅ u))dΩ

    res_prim((u, p), (v, q)) = a((u, p), (v, q)) + astab((u, p), (v, q))

    return res_prim

end

"""
    primal_steady_VMS(params::Dict{Symbol,Any})

Navier-Stokes VMS stabilized set of equations
"""
function primal_steady_VMS(params::Dict{Symbol,Any})
    @unpack ν, dt, dΩ, D, Ω, θ,uh,Cᵢ = params
    G = compute_G(Ω,params)
    GG = compute_GG(Ω,params)
    gg = compute_gg(Ω,params)

    updatekey(params, :G,G)
    updatekey(params, :GG,GG)
    updatekey(params, :gg,gg)

    Rm(u, p) = conv ∘ (uh, ∇(u)) + ∇(p)

    TRm(u, p) = τm(uh,params) * Rm(u, p)

    Bᴳ((u, p), (v, q)) = ∫((u ⋅ ∇(u)) ⋅ v)dΩ - ∫((∇ ⋅ v) * p)dΩ + ∫((q * (∇ ⋅ u)))dΩ + ν * ∫(∇(v) ⊙ ∇(u))dΩ
    B_SUPG((u, p), (v, q)) = ∫((uh ⋅ ∇(v) + ∇(q)) ⊙ TRm(u, p))dΩ + ∫((∇ ⋅ v) ⊙ (τc(uh, params) * (∇ ⋅ u) ))dΩ
    B_VMS1((u, p), (v, q)) = ∫((uh ⋅ (∇(v))') ⊙ TRm(u, p))dΩ


    res_prim((u, p), (v, q)) = Bᴳ((u, p), (v, q)) +B_SUPG((u, p), (v, q)) + B_VMS1((u, p), (v, q))

    return res_prim

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

"""
    primal_steady_SUPG(params::Dict{Symbol,Any})

Navier-Stokes unsteady SUPG stabilized set of equations
"""
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

"""
    primal_steady_VMS(params::Dict{Symbol,Any})

Navier-Stokes unsteady VMS stabilized set of equations
"""
function primal_unsteady_VMS(params::Dict{Symbol,Any})
    @unpack ν, dt, dΩ, D, Ω, θ,uh,Cᵢ = params
    G = compute_G(Ω,params)
    GG = compute_GG(Ω,params)
    gg = compute_gg(Ω,params)

    updatekey(params, :G,G)
    updatekey(params, :GG,GG)
    updatekey(params, :gg,gg)

    Rm(u, p) = conv ∘ (uh, ∇(u)) + ∇(p)

    TRm(u, p) = τm(uh,params) * Rm(u, p)

    Bᴳ(t,(u, p), (v, q)) = ∫(v ⊙ (conv ∘ (uh, ∇(u))))dΩ - ∫((∇ ⋅ v) * p)dΩ + ∫((q * (∇ ⋅ u)))dΩ + ν * ∫(∇(v) ⊙ ∇(u))dΩ
    B_SUPG(t,(u, p), (v, q)) = ∫((uh ⋅ ∇(v) + ∇(q)) ⊙ TRm(u, p))dΩ + ∫((∇ ⋅ v) ⊙ (τc(uh, params) * (∇ ⋅ u) ))dΩ
    B_VMS1(t,(u, p), (v, q)) = ∫((uh ⋅ (∇(v))') ⊙ TRm(u, p))dΩ

    m(t, (u, p), (v, q)) =   ∫(u ⋅ v)dΩ +  ∫( τm(uh,params) ⋅ (uh ⋅ ∇(v) + (uh ⋅ (∇(v))') + ∇(q)) ⋅ u)dΩ


    res_prim(t,(u, p), (v, q)) = Bᴳ(t,(u, p), (v, q)) +B_SUPG(t,(u, p), (v, q)) + B_VMS1(t,(u, p), (v, q))

    return m,res_prim
end

function get_θvp(θ)
    θvp = 1.0

    if length(θ) == 2
        θv, θp = θ
        θvp = θv / θp
    end
    return θvp

end

