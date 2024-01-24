using IncompressibleAdjoint
using Gridap,GridapGmsh
using Parameters

Reynolds = 10e3

params=Dict(
    :chord => 1.0,
    :D=>2,
    :Re => Reynolds,
    :dt => 0.05,
    :ν =>1 / Reynolds,
    :D =>2,
    :order => 1,
    :method=>:VMS,
    :tagname=>"airfoil",
    :t0=>0.0,
    :tf=>20.0,
    :t_endramp=>0.0,
    :θ=>1.0,
    :u_in=>1.0,
    :d_boundary=>VectorValue(-1.0,0.0),
    :Cᵢ=>[4,36],
)

xx = collect(0:0.001:1)
airfoil_cst = AirfoilCST(CST_NACA0012(;N=15,t= 0.00126), 0.00126)
des_points = AirfoilCSTDesign(airfoil_cst,xx)


modelname =create_msh(des_points; AoA=0.0, mesh_ref=1.0)
model = GmshDiscreteModel(modelname)


function compute_drag(uh,ph,nΓ,dΓ,params)
    CD, _ = compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)
    return CD, CD
end



params[:method]=:VMS
params[:dt]=0.01
uh,ph = solve_inc_primal_s(model, params; filename="inc-steady-VMS")


params[:method]=:SUPG
params[:dt]=:0.01
params[:tf]=:5.0

(uh,duhdt, UH, DUHDT), (ph, PH)=solve_inc_primal_u(model, params; filename="res-unsteady")

CDavg,CLavg = average_CD_CL(model, params, (uh,UH),(ph,PH))
params[:d_boundary] = VectorValue(-CDavg,0.0)

UH
IDXAVG = collect(200:500)
average_field!(uh,UH[IDXAVG])
average_field!(ph,PH[IDXAVG])



uhadj,phadj = solve_inc_adj_s(model, uh, ph,params; filename="res-adj-steady")

Uhadj,Phadj = solve_inc_adj_u(model, (uh,UH),(ph,PH), (uhadj,phadj),params; filename="res-adj-unsteady")

average_field!(uhadj,reverse(Uhadj)[10:300])
average_field!(phadj,reverse(Phadj)[10:300])

writevtk(params[:Ω], "Average", cellfields=["uh"=>uh, "ph"=>ph,"uhadj"=>uhadj, "phadj"=>phadj ])


params[:method]=:SUPG
uhadj,phadj = solve_inc_adj_s(model, uh, ph,params; filename="res-adj-steady")


modelgrid0 = copy(model.grid.node_coordinates)
Nc = IncompressibleAdjoint.get_designparameters_number(des_points)
tag = IncompressibleAdjoint.get_designparameters_tags(des_points)


vv = ones(Nc) #(tag .== "top") .- (tag .== "bottom")

J1=zeros(Nc)
J2=zeros(Nc)

using IncompressibleAdjoint

δ=0.000001
shift = vv.*δ
J1ref,J2ref= IncompressibleAdjoint.compute_sensitivity(model,params, uh,ph,uhadj,phadj; objective_function=compute_drag)

for (i,ss) in enumerate(shift)
    des_tmp = IncompressibleAdjoint.perturb_DesignParameters(des_points, i, ss)
    modelname = create_msh(des_tmp; iter = i+100, mesh_ref=1.0,AoA=0.0)
    model_tmp = GmshDiscreteModel(modelname)
    model.grid.node_coordinates .= model_tmp.grid.node_coordinates 
    J1tmp,J2tmp= IncompressibleAdjoint.compute_sensitivity(model,params, uh,ph,uhadj,phadj; objective_function=compute_drag)
    writevtk(model, "MeshPerturb/model-$(100+i)")

    J1[i] = J1tmp
    J2[i] = J2tmp

    model.grid.node_coordinates .= modelgrid0
end
J1s = (J1 .- J1ref)./shift
J2s = (J2 .- J2ref)./shift

Jtot = J1s + J2s 


using Plots
plot(J1s)
plot!(J2s)

plot(Jtot)


CDFD=zeros(Nc)
uh0s,ph0s = solve_inc_primal_s(model, params; filename=nothing)

_,CDref = IncompressibleAdjoint.obj_fun(model, params, uh0s,ph0s, compute_drag; degree=4)

for (i,ss) in enumerate(shift)
    des_tmp = IncompressibleAdjoint.perturb_DesignParameters(des_points, i, ss)
    modelname = create_msh(des_tmp; iter = i+100, mesh_ref=1.0,AoA=0.0)
    model_tmp = GmshDiscreteModel(modelname)
    model.grid.node_coordinates .= model_tmp.grid.node_coordinates 
    uh0tmp,ph0tmp = solve_inc_primal_s(model, params; filename=nothing)
    _,CDtmp = IncompressibleAdjoint.obj_fun(model, params, uh0tmp,ph0tmp, compute_drag; degree=4)

    CDFD[i] = CDtmp

    model.grid.node_coordinates .= modelgrid0
end
CDFD

GradFD = (CDFD .- CDref)./shift 

Jtot

plot(GradFD,seriestype=:scatter)
plot!(J2s*20,seriestype=:scatter)
plot!(J1s,seriestype=:scatter)


