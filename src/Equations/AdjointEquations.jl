####################################################################
#STEADY ADJOINT
####################################################################

"""
    adjoint_conservation(params)

Steady conservation equations of the adjoint problem 
"""
function adjoint_conservation(params)
    @unpack uh = params
    Rmadj(u, p) =  - ∇(p) - transpose(∇(u)) ⋅ uh - ∇(u)⋅uh
    Rcadj(u) = -1 .* (∇ ⋅ (u))
    return Rmadj, Rcadj
end

function eq_adjoint_steady(params::Dict{Symbol,Any})
    @unpack method,D,dΩ = params
    if method==:VMS
        res_adj = adjoint_steady_VMS(params)
    elseif method==:SUPG
        res_adj = adjoint_steady_SUPG(params)
    end

    rhs((v, q))= ∫(VectorValue(zeros(D)...) ⋅ v)dΩ

    return res_adj,rhs
end


"""
    adjoint_steady_SUPG(params::Dict{Symbol,Any})

It provides the set of Adjoint Equations
Symbol convention: O. Soto, & R. Lohner. (2004). On the Boundary Computation of Flow Sensitivities. https://doi.org/10.2514/6.2004-112
SUPG Stabilization: Srinath, D. N., & Mittal, S. (2010). An adjoint method for shape optimization in unsteady viscous flows. Journal of Computational Physics, 229(6), 1994–2008. https://doi.org/10.1016/j.jcp.2009.11.019
"""
function adjoint_steady_SUPG(params::Dict{Symbol,Any})
    
    @unpack ν, dt, dΩ, D, Ω, θ,uh = params
    h = h_param(Ω, D)
    updatekey(params,:h,h)   


    Rmadj, Rcadj = adjoint_conservation(params)


 

    a_adj((u, p), (v, q)) = ∫(ν * ∇(v) ⊙ ∇(u) - q * (∇ ⋅ u))dΩ + ∫(v ⊙ Rmadj(u, p))dΩ


    a_adj_stab((u, p), (v, q)) = ∫( τsu(uh, h,ν,dt) ⋅ (- transpose(∇(v)) ⋅ uh - ∇(v)⋅uh- ∇(q)) ⊙ Rmadj(u, p))dΩ +
                                 -1 * ∫( τb(uh, h,ν,dt) ⋅ (∇ ⋅ v) ⊙ Rcadj(u))dΩ

    res_adj((u, p), (v, q)) = a_adj((u, p), (v, q)) + a_adj_stab((u, p), (v, q)) #+ bc((u, p), (v, q))

    return res_adj

end

"""
    adjoint_steady_VMS(params::Dict{Symbol,Any})

It provides the set of Adjoint Equations
Symbol convention: O. Soto, & R. Lohner. (2004). On the Boundary Computation of Flow Sensitivities. https://doi.org/10.2514/6.2004-112
"""
function adjoint_steady_VMS(params::Dict{Symbol,Any})
    @unpack ν, dt, dΩ, D, Ω, θ,uh = params
    
    G, GG, gg = G_params(Ω,params)
    updatekey(params,:G,G)   
    updatekey(params,:GG,GG)   
    updatekey(params,:gg,gg)   

    Rmadj, Rcadj = adjoint_conservation(params)

    TRm(u, p) = τm(uh,params) * Rmadj(u, p)
    
    ADJBᴳ((u, p), (v, q)) = ∫(ν * ∇(v) ⊙ ∇(u) - q * (∇ ⋅ u))dΩ + ∫(v ⊙ Rmadj(u, p))dΩ
    
    ADJB_SUPG((u, p), (v, q)) = ∫( (- transpose(∇(v)) ⋅ uh - ∇(v)⋅uh- ∇(q)) ⊙ TRm(u, p))dΩ +
            -1 * ∫( τc(uh, params) ⋅ (∇ ⋅ v) ⊙ Rcadj(u))dΩ

    ADJB_VMS1((u, p), (v, q)) = ∫((-uh ⋅ (∇(v))') ⊙ TRm(u, p))dΩ


    res_adj( (u, p), (v, q)) = ADJBᴳ((u, p), (v, q)) + ADJB_SUPG((u, p), (v, q)) + ADJB_VMS1((u, p), (v, q))
    return res_adj
end




####################################################################
#UNSTEADY ADJOINT
####################################################################
function eq_adjoint_unsteady(params::Dict{Symbol,Any})
    @unpack method,D,dΩ = params
    if method==:VMS
        m,res_adj = adjoint_unsteady_VMS(params)
    elseif method==:SUPG
        m,res_adj = adjoint_unsteady_SUPG(params)
    end

    rhs(t, (v, q)) = ∫(VectorValue(zeros(D)...) ⋅ v)dΩ

    return m,res_adj,rhs
end


"""
    adjoint_unsteady_SUPG(params::Dict{Symbol,Any})

It provides the set of Adjoint Equations
Symbol convention: O. Soto, & R. Lohner. (2004). On the Boundary Computation of Flow Sensitivities. https://doi.org/10.2514/6.2004-112
SUPG Stabilization: Srinath, D. N., & Mittal, S. (2010). An adjoint method for shape optimization in unsteady viscous flows. Journal of Computational Physics, 229(6), 1994–2008. https://doi.org/10.1016/j.jcp.2009.11.019
"""
function adjoint_unsteady_SUPG(params::Dict{Symbol,Any})
    
    @unpack ν, dt, dΩ, D, Ω, θ,uh = params
    h = h_param(Ω, D)
    updatekey(params,:h,h)   

    Rmadj, Rcadj = adjoint_conservation(params)


    a_adj((u, p), (v, q)) = ∫(ν * ∇(v) ⊙ ∇(u) - q * (∇ ⋅ u))dΩ + ∫(v ⊙ Rmadj(u, p))dΩ
    a_adj_stab((u, p), (v, q)) = ∫( τsu(uh, h,ν,dt) ⋅ (-uh ⋅ ∇(v) - ∇(v) ⋅ uh - ∇(q)) ⊙ Rmadj(u, p))dΩ +
                                 -1 * ∫( τb(uh, h,ν,dt) ⋅ (∇ ⋅ v) ⊙ Rcadj(u))dΩ

    time_sign = -1
    m(t, (u, p), (v, q)) =   time_sign *∫(u ⋅ v)dΩ +  time_sign *∫( τsu(uh, h,ν,dt) ⋅ (-uh ⋅ ∇(v) - ∇(v) ⋅ uh - ∇(q)) ⋅ u)dΩ


    res_adj(t, (u, p), (v, q)) = a_adj((u, p), (v, q)) + a_adj_stab((u, p), (v, q))

    return m,res_adj

end

"""
    adjoint_unsteady_VMS(params::Dict{Symbol,Any})

It provides the set of Adjoint Equations
Symbol convention: O. Soto, & R. Lohner. (2004). On the Boundary Computation of Flow Sensitivities. https://doi.org/10.2514/6.2004-112
SUPG Stabilization: Srinath, D. N., & Mittal, S. (2010). An adjoint method for shape optimization in unsteady viscous flows. Journal of Computational Physics, 229(6), 1994–2008. https://doi.org/10.1016/j.jcp.2009.11.019
"""
function adjoint_unsteady_VMS(params::Dict{Symbol,Any})
    @unpack ν, dt, dΩ, D, Ω, θ,uh = params
    
    G, GG, gg = G_params(Ω,params)
    updatekey(params,:G,G)   
    updatekey(params,:GG,GG)   
    updatekey(params,:gg,gg)   

    Rmadj, Rcadj = adjoint_conservation(params)

    TRm(u, p) = τm(uh,params) * Rmadj(u, p)
    
    ADJBᴳ((u, p), (v, q)) = ∫(ν * ∇(v) ⊙ ∇(u) - q * (∇ ⋅ u))dΩ + ∫(v ⊙ Rmadj(u, p))dΩ
    ADJB_SUPG((u, p), (v, q)) = ∫( (- transpose(∇(v)) ⋅ uh - ∇(v)⋅uh- ∇(q)) ⊙ TRm(u, p))dΩ +
            -1 * ∫( τc(uh, params) ⋅ (∇ ⋅ v) ⊙ Rcadj(u))dΩ
            
    ADJB_VMS1((u, p), (v, q)) = ∫((-uh ⋅ (∇(v))') ⊙ TRm(u, p))dΩ

    time_sign = -1
    ADJm(t, (u, p), (v, q)) =   time_sign *∫(u ⋅ v)dΩ +  time_sign *∫( τm(uh,params) ⋅ (-uh ⋅ ∇(v) - ∇(v) ⋅ uh - ∇(q)) ⋅ u)dΩ


    res_adj(t, (u, p), (v, q)) = ADJBᴳ((u, p), (v, q)) + ADJB_SUPG((u, p), (v, q)) + ADJB_VMS1((u, p), (v, q))

    return ADJm,res_adj

end

