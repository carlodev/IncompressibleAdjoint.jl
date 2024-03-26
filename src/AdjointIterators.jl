###################################################################################
#OBJECTIVE FUNCTIONS
##################################################################################
"""
    compute_drag(uh,ph,nΓ,dΓ,params)

Objective function minimizing the Drag
"""
function compute_drag(uh,ph,nΓ,dΓ,params)
    CD, _ = compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)
    return CD, CD
end


"""
    compute_lift(uh,ph,nΓ,dΓ,params)

Objective function, it wants the CL to reach `CLtarget`. It return CL so it can be monitored.
"""
function compute_lift(uh,ph,nΓ,dΓ,params; CLtarget=0.75)
    _, CL = compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)
    return 0.5*(CL-CLtarget)^2, CL
end

###################################################################################
#Iterators
##################################################################################
"""
    compute_airfoil_forces(uh,ph,nΓ,dΓ,params::Dict{Symbol,Any})

It computes Drag and Lift over params[:tagname] boundary. It takes into account pressure and velocity gradient.
It needs the normals pointings outward respect to the body.
"""
function compute_airfoil_forces(uh,ph,nΓ,dΓ,params::Dict{Symbol,Any})
    @unpack tagname,ν = params
    IForce = ∫( -ph ⋅ nΓ+ ν* transpose(∇(uh)) ⋅ nΓ)dΓ #+ 
    D,L = sum(IForce)
    return D,L
end

"""
    compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params::Dict{Symbol,Any})

Compute the normaization of airfoil forces, obtaining CD and CL
"""
function compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params::Dict{Symbol,Any})
    @unpack chord, u_in = params
    q = 0.5 * chord*u_in
    D,L = compute_airfoil_forces(uh,ph,nΓ,dΓ,params)
    CD = D/q
    CL = L/q
    return CD,CL
end

"""
    obj_fun(model::DiscreteModel, params::Dict{Symbol,Any}, uh,ph, fun)

Wrapper that evaluates `fun` the objective function. `fun` is a user-defined function.
See eg. `compute_drag`, `compute_lift` functions. It return fitnessval,E. 
`fitnessval` is the function value that need to be minimized. `E` is a value that we want to monitor.
"""
function obj_fun(model::DiscreteModel, params::Dict{Symbol,Any}, uh,ph, fun)
    @unpack tagname= params
    Γ = BoundaryTriangulation(model; tags=tagname)
    dΓ = Measure(Γ, 8)
    nΓ = -1 .* get_normal_vector(Γ)

    fitnessval,E = fun(uh,ph,nΓ,dΓ,params)

    return fitnessval, E 
end

function rotation(n::VectorValue{2,Float64})
    n1, n2 = [n...]
    VectorValue(-n2, n1)
end

"""
    compute_sensitivity(model::DiscreteModel, params::Dict{Symbol,Any}, uh0,ph0,ϕu0, ϕp0; objective_function=compute_drag)

From the solution of the primal flow `uh0` `ph0`, and the adjoint flow `ϕu0` `ϕp0` it computes the senstivities according to the objective function.
J2 contributions are splitted and the gradients computed individually to avoid numerical cancellation
"""
function  compute_sensitivity(model::DiscreteModel, params::Dict{Symbol,Any}, uh0,ph0,ϕu0, ϕp0; objective_function=compute_drag)
    @unpack ν,tagname,i =params
    @unpack D,order,t_endramp,t0,tf,θ,dt,u_in, d_boundary = params

    V,Q = create_primal_spaces(model,params)

    u0 = D == 2 ? VectorValue(u_in, 0.0) :  VectorValue(u_in, 0.0, 0.0)
    u_walls = D == 2 ? VectorValue(0.0, 0.0) :  VectorValue(0.0, 0.0, 0.0)
    p0 = 0.0

    U = TrialFESpace(V, [u0, u0, u_walls])
    P = TrialFESpace(Q, p0)


    V_adj,Q_adj = create_adjoint_spaces(model, params)
    U_adj = TrialFESpace(V_adj, [d_boundary, VectorValue(0, 0),VectorValue(0, 0)])
    P_adj = TrialFESpace(Q_adj)
    
    
    uh = FEFunction(U,uh0.free_values)
    ph = FEFunction(P,ph0.free_values)
    ϕu = FEFunction(U_adj,ϕu0.free_values)
    ϕp = FEFunction(P_adj,ϕp0.free_values)

    Ω = Triangulation(model)
    dΩ = Measure(Ω, 8)
    J1, _ = obj_fun(model,params, uh,ph, objective_function)

    Γ = BoundaryTriangulation(model; tags=tagname)
    dΓ = Measure(Γ, 8)
    nΓ =  -1 .*get_normal_vector(Γ) #-1, pointing outward

    J2_1 = sum(∫(ϕu⋅(transpose(∇(uh))⋅uh ))dΩ) 
    J2_2 = sum(∫(ϕu⋅(∇(ph)))dΩ)
    J2_3 = sum(∫(ϕp⋅(∇⋅(uh)))dΩ)
    J2_4 = ν*sum(∫((∇(ϕu)⊙∇(uh)))dΩ)
    J2_5 = -ν* sum(∫(ν*ϕu⋅transpose(∇(uh))⋅ nΓ)dΓ)
   
    writevtk(Ω, "MeshPerturb/model-$(100+i)", cellfields=["uh"=>uh,"ph"=>ph, "uadj"=>ϕu,"padj"=>ϕp])

    return J1,(J2_1,J2_2,J2_3,J2_4,J2_5)
end



###################################################################################
#ITERATORS
##################################################################################
function iterate_optimization(des_params, params; iter=0, objective_function=compute_drag)
    
    modelname = create_msh(des_params; iter = iter,  mesh_ref=4)
    model = GmshDiscreteModel(modelname)
    uh,ph = solve_inc_primal(model, tagname; filename=joinpath("Results_primal", "res-$(iter)"))
    _, S = obj_fun(model,params,uh, ph,objective_function)

return model, S, (uh,ph)
end



function iterate_optimization(uh,ph,model, des_params, params; iter=0, detail=false, δ=0.02, α=0.01, objective_function=compute_drag)

ϕu, ϕp = solve_inc_adj_s(model, uh,ph, params; filename="adj-$iter")#joinpath("Results","inc-adj-res-iter-$iter")

J1ref,J2ref= compute_sensitivity(model,params, uh,ph,ϕu, ϕp; objective_function=objective_function)
modelgrid0 = copy(model.grid.node_coordinates)
Nc = get_designparameters_number(des_params)
tag = get_designparameters_tags(des_params)

vv = (tag .== "top") .- (tag .== "bottom")
J1=zeros(Nc)
J2=zeros(Nc)

shift = vv.*δ
for (i,ss) in enumerate(shift)
        des_tmp = perturb_DesignParameters(des_params, i, ss)
        modelname = create_msh(des_tmp; iter = i+100, mesh_ref=4,AoA=4.0)
        model_tmp = GmshDiscreteModel(modelname)
        model.grid.node_coordinates .= model_tmp.grid.node_coordinates 
       J1tmp,J2tmp=compute_sensitivity(model,params, uh,ph,ϕu, ϕp; objective_function=objective_function)
       J1[i] = J1tmp
       J2[i] = J2tmp
       model.grid.node_coordinates .= modelgrid0
end



J1s = (J1 .- J1ref)./ (δ) #Geometric gradient
J2s =  (J2 .- J2ref)./(δ)
Jtot = J1s +J2s #Total Gradient

#Step, -1 because opposite direction of the gradient
ΔD = ShiftUpdate(-vv.*α.*Jtot, tag)

contr_new = perturb_DesignParameters(des_params, ΔD)

modelname = create_msh(contr_new; iter = iter+1, mesh_ref=4, AoA=4.0)
model = GmshDiscreteModel(modelname)
model.grid.node_coordinates .= model.grid.node_coordinates 

uh,ph = solve_inc_primal_s(model, params; filename=joinpath("Results_primal", "res-$(iter+1)"))

fitnessval, S = obj_fun(model,params,uh, ph,objective_function)

if detail
    return model, (J1s,J2s,Jtot), (fitnessval,S), contr_new, (uh,ph)
else
    return model, Jtot, S, contr_new, (uh,ph)
end


end


