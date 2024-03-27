######################################################
#Steady Direct Differentiation
#######################################################
function eq_direct_differentiation_steady(params::Dict{Symbol,Any})
    @unpack method, D,dΩ= params
    if method == :VMS
        res_dd = direct_differentiation_steady_VMS(params)
    elseif method == :SUPG
        res_dd = direct_differentiation_steady_SUPG(params)
    end

    rhs((v, q)) = ∫(VectorValue(zeros(D)...) ⋅ v)dΩ

    return res_dd, rhs
end

"""
    direct_differentiation_SUPG(params::Dict{Symbol,Any})

Formulation from
Janssens -P Vandenschrick -K Stevens -G Alessi, B. (n.d.). THE CONTINUOUS ADJOINT APPROACH APPLIED TO THE STABILIZED FINITE-ELEMENT FORMULATION OF THE INCOMPRESSIBLE NAVIER-STOKES EQUATIONS. www.euroturbo.eu
"""
function direct_differentiation_steady_SUPG(params::Dict{Symbol,Any})

    @unpack ν, dt, dΩ, D, Ω, θ,uh = params
    h = h_param(Ω, D)
    updatekey(params, :h,h)

    ### ub ∂u/∂bi
    ### pb ∂p/∂bi
    Rm(ub, pb) = transpose(∇(ub))⋅uh+∇(uh)⋅ub+∇(pb)

    a((ub, pb), (v, q)) = ∫(v⊙Rm(ub, pb))dΩ + ∫(q * (∇ ⋅ ub) )dΩ + ∫(ν * ∇(v) ⊙ ∇(ub) )dΩ

    astab((u, p), (v, q)) = ∫( τsu(uh, h,ν,dt) ⋅ (uh ⋅ ∇(v) + ∇(q)) ⊙ Rm(u, p))dΩ + ∫( τb(uh, h,ν,dt) ⋅ (∇ ⋅ v) ⊙ (∇ ⋅ u))dΩ

    res_dd((u, p), (v, q)) = a((u, p), (v, q)) + astab((u, p), (v, q))

    return res_dd

end



function direct_differentiation_steady_VMS(params::Dict{Symbol,Any})

    @unpack ν, dt, dΩ, D, Ω, θ,uh = params
    G = compute_G(Ω,params)
    GG = compute_GG(Ω,params)
    gg = compute_gg(Ω,params)

    updatekey(params, :G,G)
    updatekey(params, :GG,GG)
    updatekey(params, :gg,gg)

    ### ub ∂u/∂bi
    ### pb ∂p/∂bi
    Rm(ub, pb) = transpose(∇(ub))⋅uh+∇(uh)⋅ub+∇(pb)
    TRm(u, p) = τm(uh,params) * Rm(u, p)

    Bᴳ((ub, pb), (v, q)) = ∫(v⊙Rm(ub, pb))dΩ + ∫(q * (∇ ⋅ ub) )dΩ + ∫(ν * ∇(v) ⊙ ∇(ub) )dΩ

    B_SUPG((u, p), (v, q)) = ∫((uh ⋅ ∇(v) + ∇(q)) ⊙ TRm(u, p))dΩ + ∫((∇ ⋅ v) ⊙ (τc(uh, params) * (∇ ⋅ u) ))dΩ
    B_VMS1((u, p), (v, q)) = ∫((uh ⋅ (∇(v))') ⊙ TRm(u, p))dΩ

    # m((u, p), (v, q)) =   ∫(u ⋅ v)dΩ +  ∫( τm(uh,params) ⋅ (uh ⋅ ∇(v) + (uh ⋅ (∇(v))') + ∇(q)) ⋅ u)dΩ

    res_dd((u, p), (v, q)) = Bᴳ((u, p), (v, q)) +B_SUPG((u, p), (v, q)) + B_VMS1((u, p), (v, q))

    return res_dd

end



######################################################
#UnSteady Direct Differentiation
#######################################################
function eq_direct_differentiation_unsteady(params::Dict{Symbol,Any})
    @unpack method, D,dΩ= params
    if method == :VMS
        m,res_dd = direct_differentiation_unsteady_VMS(params)
    elseif method == :SUPG
        m,res_dd = direct_differentiation_unsteady_SUPG(params)
    end

    rhs(t, (v, q)) = ∫(VectorValue(zeros(D)...) ⋅ v)dΩ

    return m, res_prim, rhs
end




"""
    direct_differentiation_unsteady_SUPG(params::Dict{Symbol,Any})

"""
function direct_differentiation_unsteady_SUPG(params::Dict{Symbol,Any})
    

    @unpack ν, dt, dΩ, D, Ω, θ,uh = params
    h = h_param(Ω, D)
    updatekey(params, :h,h)
    θvp = get_θvp(θ)

    ### ub ∂u/∂bi
    ### pb ∂p/∂bi
    Rm(ub, pb) = transpose(∇(ub))⋅uh+∇(uh)⋅ub+∇(pb)

    a((ub, pb), (v, q)) = ∫(v⊙Rm(ub, pb))dΩ + ∫(q * (∇ ⋅ ub) )dΩ + ∫(ν * ∇(v) ⊙ ∇(ub) )dΩ

    astab((u, p), (v, q)) = ∫( τsu(uh, h,ν,dt) ⋅ (uh ⋅ ∇(v) + ∇(q)) ⊙ Rm(u, p))dΩ + ∫( τb(uh, h,ν,dt) ⋅ (∇ ⋅ v) ⊙ (∇ ⋅ u))dΩ

    res_dd(t,(u, p), (v, q)) = a((u, p), (v, q)) + astab((u, p), (v, q))

    m(t, (u, p), (v, q)) = ∫(u ⋅ v)dΩ + ∫( τsu(uh, h,ν,dt) ⋅ (uh ⋅ ∇(v) + θvp * ∇(q)) ⊙ u)dΩ


    return res_dd

end

"""
    direct_differentiation_unsteady_VMS(params::Dict{Symbol,Any})

"""
function direct_differentiation_unsteady_VMS(params::Dict{Symbol,Any})
    @unpack ν, dt, dΩ, D, Ω, θ,uh = params
    G = compute_G(Ω,params)
    GG = compute_GG(Ω,params)
    gg = compute_gg(Ω,params)

    updatekey(params, :G,G)
    updatekey(params, :GG,GG)
    updatekey(params, :gg,gg)

    ### ub ∂u/∂bi
    ### pb ∂p/∂bi
    Rm(ub, pb) = transpose(∇(ub))⋅uh+∇(uh)⋅ub+∇(pb)
    TRm(u, p) = τm(uh,params) * Rm(u, p)

    Bᴳ((ub, pb), (v, q)) = ∫(v⊙Rm(ub, pb))dΩ + ∫(q * (∇ ⋅ ub) )dΩ + ∫(ν * ∇(v) ⊙ ∇(ub) )dΩ

    B_SUPG((u, p), (v, q)) = ∫((uh ⋅ ∇(v) + ∇(q)) ⊙ TRm(u, p))dΩ + ∫((∇ ⋅ v) ⊙ (τc(uh, params) * (∇ ⋅ u) ))dΩ
    B_VMS1((u, p), (v, q)) = ∫((uh ⋅ (∇(v))') ⊙ TRm(u, p))dΩ

    res_dd(t,(u, p), (v, q)) = Bᴳ((u, p), (v, q)) +B_SUPG((u, p), (v, q)) + B_VMS1((u, p), (v, q))
    
    m(t, (u, p), (v, q)) =   ∫(u ⋅ v)dΩ +  ∫( τm(uh,params) ⋅ (uh ⋅ ∇(v) + (uh ⋅ (∇(v))') + ∇(q)) ⋅ u)dΩ

    return m,res_dd

end

