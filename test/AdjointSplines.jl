using IncompressibleAdjoint
using Gridap,GridapGmsh
using Parameters
using Statistics

Reynolds = 1_000

params=Dict(
    :chord => 1.0,
    :u_in=>1.0,
    :D=>2,
    :Re => Reynolds,
    :dt => 0.05,
    :ν =>1 / Reynolds,
    :D =>2,
    :order => 1,
    :method=>:SUPG,
    :tagname=>"airfoil",
    :t0=>0.0,
    :tf=>10.0,
    :t_endramp=>0.0,
    :t_primal_start_avg=>5.0,
    :θ=>1.0,
    :d_boundary=>VectorValue(0.0,-1.0),

)



@unpack order,tagname,ν,dt = params

des_points = IncompressibleAdjoint.initialize_control_points(xx; chord = 1.0, initfun=NACA00)


modelname =create_msh(des_points; AoA=4.0, mesh_ref=4.0)
model = GmshDiscreteModel(modelname)
writevtk(model,"model-spline")

iter=0
CLvec= Float64[]

# for iter = 0:1:20

writevtk(model, "model-$iter")

(uh0,UH),(ph0,PH)=solve_inc_primal_u(model, params; filename="res-unsteady")

CDavg,CLavg = average_CD_CL(model, params, (uh0,UH),(ph0,PH))
params[:d_boundary] = VectorValue(0.0,CLavg)


uhadj,phadj = solve_inc_adj_s(model, uh0, ph0,params; filename="res-adj-steady")

Uhadj,Phadj = solve_inc_adj_u(model, (uh0,UH),(ph0,PH), (uhadj,phadj),params; filename="res-adj-unsteady")

average_field!(uh0,UH[50:end])
average_field!(ph0,PH[50:end])

average_field!(uhadj,Uhadj[50:150])
average_field!(phadj,Phadj[50:150])
fitnessval, S = IncompressibleAdjoint.obj_fun(model,params,uh0, ph0,compute_lift)
push!(CLvec,S)

writevtk(params[:Ω], "Avg", cellfields=["uh"=>uh0,"ph"=>ph0,"uhadj"=>uhadj,"phadj"=>phadj])




J1ref,J2ref= IncompressibleAdjoint.compute_sensitivity(model,params, uh0,ph0,uhadj, phadj; objective_function=compute_lift)
modelgrid0 = copy(model.grid.node_coordinates)
Nc = IncompressibleAdjoint.get_designparameters_number(des_points)
tag = IncompressibleAdjoint.get_designparameters_tags(des_points)

vv = (tag .== "top") .- (tag .== "bottom")
J1=zeros(Nc)
J2=zeros(Nc)
δ=0.01
shift = vv.*δ
for (i,ss) in enumerate(shift)
        des_tmp = IncompressibleAdjoint.perturb_DesignParameters(des_points, i, ss)
        modelname = create_msh(des_tmp; iter = i+100, mesh_ref=4,AoA=4.0)
        model_tmp = GmshDiscreteModel(modelname)
        model.grid.node_coordinates .= model_tmp.grid.node_coordinates 
       J1tmp,J2tmp= IncompressibleAdjoint.compute_sensitivity(model,params, uh0,ph0,uhadj,phadj; objective_function=compute_lift)
       J1[i] = J1tmp
       J2[i] = J2tmp
       model.grid.node_coordinates .= modelgrid0
end


J1s = (J1 .- J1ref)./ (δ) #Geometric gradient
J2s =  (J2 .- J2ref)./(δ)
Jtot = J1s +J2s #Total Gradient


α=0.5/(iter+1)
#Step, -1 because opposite direction of the gradient
ΔD = IncompressibleAdjoint.ShiftUpdate(-vv.*α.*Jtot, tag)
iter = iter+1
contr_new = IncompressibleAdjoint.perturb_DesignParameters(des_points, ΔD)
modelname = create_msh(contr_new; iter = iter+1, mesh_ref=4, AoA=4.0)
model = GmshDiscreteModel(modelname)
des_points = contr_new

# end
CLvec



# using IncompressibleAdjoint

function compute_lift(uh,ph,nΓ,dΓ,params)
    _, CL = compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)
    return 0.5*(CL-0.75)^2, CL
end

function compute_drag(uh,ph,nΓ,dΓ,params)
    CD, _ = compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)
    return CD, CD
end

function compute_efficiency(uh,ph,nΓ,dΓ,params)
    CD, CL = compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)
    return Drag - 0.1*Lift, CL/CD
end


using IncompressibleAdjoint

CDv,CLv =istantaneus_CD_CL(model, params, (uh0,UH),(ph0,PH))
