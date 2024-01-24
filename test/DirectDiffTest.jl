using IncompressibleAdjoint
using Gridap,GridapGmsh
using Parameters

Reynolds = 20_000

params=Dict(
    :chord => 1.0,
    :D=>2,
    :Re => Reynolds,
    :dt => 0.05,
    :ν =>1 / Reynolds,
    :D =>2,
    :order => 1,
    :method=>:SUPG,
    :tagname=>"airfoil",
    :t0=>0.0,
    :tf=>20.0,
    :t_endramp=>0.0,
    :θ=>1.0,
    :u_in=>1.0,
    :d_boundary=>VectorValue(-1.0,0.0),
    :Cᵢ=>[4,36],
)

xx = collect(0:0.0002:1.0)

airfoil_cst = AirfoilCST(CST_NACA0012(;N=15,t= 0.00126), 0.00126)
des_points = AirfoilCSTDesign(airfoil_cst,xx)


modelname =create_msh(des_points; AoA=0.0, mesh_ref=8.0)
model = GmshDiscreteModel(modelname)


function compute_drag(uh,ph,nΓ,dΓ,params)
    CD, _ = compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)
    return CD, CD
end



params[:dt]=0.05
params[:tf]=5.0

uh,ph = solve_inc_primal_s(model, params; filename="inc-steady-SUPG")

(uh,duhdt, UH, DUHDT), (ph, PH)=solve_inc_primal_u(model, params; filename="res-unsteady")

writevtk(params[:Ω], "Finaluhph", cellfields=["uh"=>uh,"ph"=>ph])

CD0,_=IncompressibleAdjoint.obj_fun(model,params, uh,ph,compute_drag)
using IncompressibleAdjoint



@unpack ν=params
Γtagname = BoundaryTriangulation(model; tags="airfoil")
nΓtagname = -1 .* get_normal_vector(Γtagname)
dΓtagname = Measure(Γtagname, 4)

cf = Dict("airfoil_nΓ" => nΓtagname)
ttrian = nΓtagname.trian
f = (reffe) -> Gridap.Geometry.UnstructuredGrid(reffe)
ref_grids = map(f, Gridap.Geometry.get_reffes(ttrian))
visgrid = Gridap.Visualization.VisualizationGrid(ttrian, ref_grids)
pdata = Gridap.Visualization._prepare_pdata(ttrian, cf, visgrid.cell_to_refpoints)

ttrian = Γtagname.trian
ref_grids = map(f, Gridap.Geometry.get_reffes(ttrian))
visgrid = Gridap.Visualization.VisualizationGrid(ttrian, ref_grids)
airfoil_points = visgrid.sub_grid.node_coordinates

nv_airfoil = pdata["airfoil_nΓ"]

NV = nv_airfoil[1:50:800]
dJdb = zeros(length(NV))
dir = VectorValue(1,0)


for (i,nv) in enumerate(NV)
    println(i)
    ub,pb=solve_inc_direct_differentiation_s(model,uh,params, nv)
    -sum(∫(pb⋅nΓtagname⋅ dir)dΓtagname)
    dJdb[i] = -sum(∫(pb⋅nΓtagname⋅ dir)dΓtagname)
end
writevtk(params[:Ω], "DirDiff", cellfields=["ub"=>ub,"pb"=>pb])

params[:d_boundary]=VectorValue(-1.0,0.0)

uhadj,phadj = solve_inc_adj_s(model, uh, ph,params; filename="res-adj-steady")

params[:tf] = 1.0
Uhadj,Phadj = solve_inc_adj_u(model, (uh,UH),(ph,PH), (uhadj,phadj),params; filename="res-adj-unsteady")

copyto!(uhadj.free_values, Uhadj[end])
copyto!(phadj.free_values, Phadj[end])

modelgrid0 = copy(model.grid.node_coordinates)
Nc = IncompressibleAdjoint.get_designparameters_number(des_points)
tag = IncompressibleAdjoint.get_designparameters_tags(des_points)

vv = ones(Nc) #(tag .== "top") .- (tag .== "bottom")

J1=zeros(Nc)
J2=zeros(Nc)
J2extra=zeros(Nc)


using IncompressibleAdjoint

δ=0.0001
shift = vv.*δ
J1ref,(J2ref,J2refExtra)= IncompressibleAdjoint.compute_sensitivity(model,params, uh,ph,uhadj,phadj; objective_function=compute_drag)

for (i,ss) in enumerate(shift[1:15])
    des_tmp = IncompressibleAdjoint.perturb_DesignParameters(des_points, i, ss)
    modelname = create_msh(des_tmp; iter = i+100, mesh_ref=8.0,AoA=0.0)
    model_tmp = GmshDiscreteModel(modelname)
    model.grid.node_coordinates .= model_tmp.grid.node_coordinates 
    J1tmp,(J2tmp,J2tmpextra)= IncompressibleAdjoint.compute_sensitivity(model,params, uh,ph,uhadj,phadj; objective_function=compute_drag)
    writevtk(model, "MeshPerturb/model-$(100+i)")

    J1[i] = J1tmp
    J2[i] = J2tmp
    J2extra[i] = J2tmpextra

    model.grid.node_coordinates .= modelgrid0
end

J1s = (J1 .- J1ref)./shift
J2s = (J2 .- J2ref)./shift
J2sextra = (J2extra .- J2refExtra)./shift
Jtot = J1s + J2s +J2sextra


using Plots
plot(J1s)
plot!(J2s)
plot!(J2sextra)
plot!(Jtot)

CDFD=zeros(15)
uh0s,ph0s = solve_inc_primal_s(model, params; filename=nothing)
_,CDref = IncompressibleAdjoint.obj_fun(model, params, uh0s,ph0s, compute_drag; degree=4)

for (i,ss) in enumerate(shift[1:15])
    des_tmp = IncompressibleAdjoint.perturb_DesignParameters(des_points, i, ss)
    modelname = create_msh(des_tmp; iter = i+100, mesh_ref=8.0,AoA=0.0)
    model_tmp = GmshDiscreteModel(modelname)
    model.grid.node_coordinates .= model_tmp.grid.node_coordinates 
    uh0tmp,ph0tmp = solve_inc_primal_s(model, params; filename=nothing)
    _,CDtmp = IncompressibleAdjoint.obj_fun(model, params, uh0tmp,ph0tmp, compute_drag; degree=4)

    CDFD[i] = CDtmp

    model.grid.node_coordinates .= modelgrid0
end
CDFD

GradFD = (CDFD .- CDref)./shift[1:15]

Jtot

plot(getindex.(airfoil_points[1:50:800],1), dJdb,seriestype=:scatter, label="direct differentiation\n normal vector")
plot!(LinRange(0,1.0,15),Jtot[1:15]*50,seriestype=:scatter,label="class shape transformation weights")
plot!(LinRange(0,1.0,15),GradFD*50,seriestype=:scatter,label="class shape transformation weights FD")

plot!(xlabel="")


plot(J2s[1:15],seriestype=:scatter,label="class shape transformation weights")
plot!(J2sextra[1:15],seriestype=:scatter,label="class shape transformation weights")

plot(J1s[1:15],seriestype=:scatter,label="class shape transformation weights")
plot!(GradFD,seriestype=:scatter,label="class shape transformation weights FD")
