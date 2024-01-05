####################################################################
#STEADY ADJOINT
####################################################################
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


function adjoint_steady_SUPG(params::Dict{Symbol,Any})
    
    @unpack ν, dt, dΩ, D, Ω, θ,uh = params
    h = h_param(Ω, D)
    updatekey(params,:h,h)   
    Rmadj(u, p) = -uh ⋅ ∇(u) - ∇(u) ⋅ uh - ∇(p)
    Rcadj(u) = -1 .* (∇ ⋅ (u))

    a_adj((u, p), (v, q)) = ∫(ν * ∇(v) ⊙ ∇(u) - q * (∇ ⋅ u))dΩ + ∫(v ⊙ Rmadj(u, p))dΩ
    a_adj_stab((u, p), (v, q)) = ∫( τsu(uh, h,ν,dt) ⋅ (-uh ⋅ ∇(v) - ∇(v) ⋅ uh - ∇(q)) ⊙ Rmadj(u, p))dΩ +
                                 -1 * ∫( τb(uh, h,ν,dt) ⋅ (∇ ⋅ v) ⊙ Rcadj(u))dΩ



    res_adj((u, p), (v, q)) = a_adj((u, p), (v, q)) + a_adj_stab((u, p), (v, q))

    return res_adj

end

function adjoint_steady_VMS(params::Dict{Symbol,Any})

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


function adjoint_unsteady_SUPG(params::Dict{Symbol,Any})
    
    @unpack ν, dt, dΩ, D, Ω, θ,uh = params
    h = h_param(Ω, D)
    updatekey(params,:h,h)   
    Rmadj(u, p) = -uh ⋅ ∇(u) - ∇(u) ⋅ uh - ∇(p)
    Rcadj(u) = -1 .* (∇ ⋅ (u))

    a_adj((u, p), (v, q)) = ∫(ν * ∇(v) ⊙ ∇(u) - q * (∇ ⋅ u))dΩ + ∫(v ⊙ Rmadj(u, p))dΩ
    a_adj_stab((u, p), (v, q)) = ∫( τsu(uh, h,ν,dt) ⋅ (-uh ⋅ ∇(v) - ∇(v) ⋅ uh - ∇(q)) ⊙ Rmadj(u, p))dΩ +
                                 -1 * ∫( τb(uh, h,ν,dt) ⋅ (∇ ⋅ v) ⊙ Rcadj(u))dΩ

    time_sign = -1
    m(t, (u, p), (v, q)) =   time_sign *∫(u ⋅ v)dΩ +  time_sign *∫( τsu(uh, h,ν,dt) ⋅ (-uh ⋅ ∇(v) - ∇(v) ⋅ uh - ∇(q)) ⋅ u)dΩ


    res_adj(t, (u, p), (v, q)) = a_adj((u, p), (v, q)) + a_adj_stab((u, p), (v, q))

    return m,res_adj

end

function adjoint_unsteady_VMS(params::Dict{Symbol,Any})

end

