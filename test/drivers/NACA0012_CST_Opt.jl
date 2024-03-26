module NACA0012_CST_Opt

using IncompressibleAdjoint
using IncompressibleAdjoint.Geometry
using IncompressibleAdjoint.IncompressibleSolvers
using IncompressibleAdjoint.Utils

using Gridap,GridapGmsh
using Parameters
using NLopt
using Test

"""
NACA0012 - Reynolds 1000, AoA 2.5, CL adjoint optimization, aiming a targhet CL=0.35.
The primal flow is solved unsteady, the adjoint is solved steadily from the average of the primal one
The parametrization is CST as β shape parameters
"""


Reynolds = 1000

params=Dict(
    :chord => 1.0,
    :u_in=>1.0,
    :D=>2,
    :Re => Reynolds,
    :dt => 0.05,
    :ν =>1 / Reynolds,
    :D =>2,
    :order => 1,
    :method=>:VMS,
    :tagname=>"airfoil",
    :t0=>0.0,
    :tf=>10.0,
    :t_endramp=>0.0,
    :t_primal_start_avg=>5.0,
    :θ=>1.0,
    :d_boundary=>VectorValue(-1.0,0.0), # It is the boundary condition for the adjoint problem at tag `:tagname`
    :AoA=>2.5,
    :Cᵢ=>[4,36],
    :mesh_ref=>1.0,
)
mesh_folder = joinpath("test", "TestMesh")


@unpack AoA,ν,tagname,mesh_ref = params

airfoil_cst = AirfoilCST(CST_NACA0012(;N=15,t= 0.0), 0.0)
xx = collect(0:0.0001:1)
des_points = AirfoilCSTDesign(airfoil_cst,xx)

modelname =create_msh(des_points; AoA=AoA, mesh_ref=mesh_ref,folder=mesh_folder)
model = GmshDiscreteModel(modelname)

@test typeof(model) <: Gridap.Geometry.UnstructuredDiscreteModel

#Store the node_coordinates of the original model
modelgrid0 = copy(model.grid.node_coordinates)

Ndes = get_designparameters_number(des_points)
des_tag = get_designparameters_tags(des_points)

@test Ndes == 30
@test typeof(des_tag) == Vector{String}


### Deformation used to compute the gradients with the adjoint
δ=0.01
des_bool = (des_tag .== "top") .- (des_tag .== "bottom")
shift = des_bool.*δ


"""
Optimization function, it aims for a CL=0.35.
Also the CL is obtained, to check its evolution
This is a model on how to write an optimization function.
It return fitnessval, val
where 
fitnessval: value that the algorithm is minimizing
val: some output that is interesting to monitor
"""
function compute_CL(uh,ph,nΓ,dΓ,params)
    _, CL = compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)
    return 0.5*(CL-0.35)^2, CL
end


#Get Airfoil characteristics
airfoil_nodes_top,airfoil_nodes_bottom,_,_ = get_airfoil_characteristics(model,params)

airfoil_nodes_top

@test typeof(airfoil_nodes_top) <: Vector
@test typeof(airfoil_nodes_bottom) <: Vector

### Initialize output Vectors. 20 maximum iters

#Store the shape parameters
βv= [Vector{Float64}(undef, Ndes) for _ = 1:20]

#Store the airfoil nodes and pressure distribution
nodes_top = [Vector{Float64}(undef, length(airfoil_nodes_top)) for _ = 1:20] 
nodes_bottom = [Vector{Float64}(undef, length(airfoil_nodes_bottom)) for _ = 1:20] 
cp_top = [Vector{Float64}(undef, length(airfoil_nodes_top)) for _ = 1:20] 
cp_bottom = [Vector{Float64}(undef, length(airfoil_nodes_bottom)) for _ = 1:20] 




merge!(params,Dict(:iter=>0,:i=>0, :fitnessval=>zeros(100), :CL=>zeros(100),:βv=>copy(βv),:fgrad=>copy(βv),
:des_points=>des_points,:Ndes=>Ndes,:shift=>shift,:obj_fun=>compute_CL,
:nodes_top=>nodes_top,:nodes_bottom=>nodes_bottom,:cp_top=>cp_top,:cp_bottom=>cp_bottom))

### Solve the primal unsteady flow
(uh,duhdt, UH, DUHDT), (ph, PH) = solve_inc_primal_u(model, params; filename="res_unsteady")

typeof(uh)

@test typeof(uh) <: Gridap.FESpaces.SingleFieldFEFunction
@test typeof(ph) <: Gridap.FESpaces.SingleFieldFEFunction

@test typeof(UH) <: Vector
@test typeof(PH) <: Vector



CL_vec = zeros(length(PH))

for (i,(uhval,phval)) in enumerate(zip(UH,PH))
    copyto!(uh.free_values, uhval)
    copyto!(ph.free_values, phval)
    fitnessval, CL  = IncompressibleAdjoint.obj_fun(model, params, uh,ph, compute_CL)
    CL_vec[i] = copy(CL)
end

CL_history = [Vector{Float64}(undef, length(PH)) for _ = 1:20] 
#It stores the CL at each time step for every iteration optimization cycle
CL_history[1] = copy(CL_vec)




### using uh, ph as initial solution for the new unsteady simulations
merge!(params,Dict(:uh00=>uh,:ph00=>ph, :UH=>UH,:PH=>PH, :CL_history=>CL_history))



"""
    optimization_loop(w::Vector,grad::Vector)

It receives:
-w: CST weights 
-grad: gradient dI/dw_i

It is using the solution at the previous iteration on the new geometry in order to reduce the initial transient
"""
function optimization_loop(w::Vector,grad::Vector)
 @unpack des_points,Ndes,shift,AoA, obj_fun, iter, 
            mesh_ref,uh00,ph00,UH,PH=params
    des_points = update_CST_weights(w,des_points)


    modelname =create_msh(des_points; AoA=AoA, mesh_ref=mesh_ref, folder=mesh_folder)
    model = GmshDiscreteModel(modelname)
    
    uh = uh00
    ph=ph00
    
    (uh,duhdt, UH, DUHDT), (ph, PH) = solve_inc_primal_u(model, params; filename="inc-results-$iter", uh00=uh00,ph00=ph00)
    
    
    
    #### Compute Unsteady CL
    CL_vec = zeros(length(PH))

    for (i,(uhval,phval)) in enumerate(zip(UH,PH))
    copyto!(uh.free_values, uhval)
    copyto!(ph.free_values, phval)
        fitnessval, CL  = IncompressibleAdjoint.obj_fun(model, params, uh,ph, compute_CL)
        CL_vec[i] = copy(CL)
    end

    params[:CL_history][iter+1] = copy(CL_vec)

    CD,CL = average_CD_CL(model,params,(uh,UH),(ph,PH))


#### Adjoint Boundary Conditions
    params[:d_boundary] = VectorValue(0.0, (0.35-CL))
    

    ### run a steady adjoint on the average field
    average_field!(uh,UH[end-50:end])
    average_field!(ph,PH[end-50:end])

    #### average Cp
    airfoil_nodes_top, airfoil_nodes_bottom,_,_ = get_airfoil_characteristics(model, params; tag="airfoil")
    cp_top,cp_bottom = get_aerodynamic_features(params,model, uh,ph)
    params[:nodes_top][iter+1] = copy(getindex.(airfoil_nodes_top,1))
    params[:nodes_bottom][iter+1] =  copy(getindex.(airfoil_nodes_bottom,1))
    params[:cp_top][iter+1] =  copy(cp_top)
    params[:cp_bottom][iter+1] =  copy(cp_bottom)




    uhadj,phadj = solve_inc_adj_s(model, uh, ph,params; filename="res-adj-steady")

    J1ref,(J2_1ref,J2_2ref,J2_3ref,J2_4ref,J2_5ref)= IncompressibleAdjoint.compute_sensitivity(model,params, uh,ph,uhadj,phadj; objective_function=obj_fun)

    fitnessval, CL  = IncompressibleAdjoint.obj_fun(model, params, uh,ph, obj_fun)


    params[:iter]= iter +1 
    params[:βv][iter+1] = copy(w)

    params[:CL][iter+1] = copy(CL)
    params[:fitnessval][iter+1]= copy(fitnessval)
    J1=zeros(Ndes)
    J2_1=zeros(Ndes)
    J2_2=zeros(Ndes)
    J2_3=zeros(Ndes)
    J2_4=zeros(Ndes)
    J2_5=zeros(Ndes)

    copyto!(params[:uh00].free_values, uh.free_values)
    copyto!(params[:ph00].free_values, ph.free_values)


    for (i,ss) in enumerate(shift)

        des_temp = perturb_DesignParameters(des_points, i, ss)
        modelname_tmp = create_msh(des_temp; AoA=AoA, mesh_ref=mesh_ref,folder=mesh_folder)
        model_tmp = GmshDiscreteModel(modelname_tmp)

        J1tmp,(J2_1tmp,J2_2tmp,J2_3tmp,J2_4tmp,J2_5tmp)= IncompressibleAdjoint.compute_sensitivity(model_tmp,params, uh,ph,uhadj,phadj; objective_function=compute_CL)


        J1[i] = J1tmp
        J2_1[i] = J2_1tmp
        J2_2[i] = J2_2tmp
        J2_3[i] = J2_3tmp
        J2_4[i] = J2_4tmp
        J2_5[i] = J2_5tmp
        println(i)
    end
    
    J1s = (J1 .- J1ref)./δ
    J2s1 = (J2_1 .- J2_1ref)./δ
    J2s2 = (J2_2 .- J2_2ref)./δ
    J2s3 = (J2_3 .- J2_3ref)./δ
    J2s4 = (J2_4 .- J2_4ref)./δ
    J2s5 = (J2_5 .- J2_5ref)./δ
    
    
    J2s = J2s1 + J2s2 + J2s3 + J2s4 + J2s5
    Jtot = J1s+ J2s
    

    grad .= J2s.*des_bool

    params[:fgrad][iter+1] = copy(grad)
    println("CL = $CL")
    # jldsave("results/NACA0012_1000_2p5/parameters_2.jld2"; params)
    return J1ref

end



## Set the optimization algorithm
opt = Opt(:LD_LBFGS, Ndes)

## Set the bounds
opt.lower_bounds = [-0.1 .* ones(Int(Ndes/2));-1 .* ones(Int(Ndes/2))]
opt.upper_bounds = [ones(Int(Ndes/2)); 0.1 .* ones(Int(Ndes/2))]

## Set the tolerance
opt.xtol_rel = 2.5e-2

opt.min_objective = optimization_loop

#Initialization
w_init,_ = get_CST_values(des_points) 

@test typeof(opt) == Opt

###### Resolution
(minf,minx,ret) = optimize(opt,w_init)


@unpack iter,CL = params

@test isapprox(CL[iter],0.35,atol=0.2 )


#Delete Mesh created
rm(mesh_folder; force=true, recursive=true)

end