include("../src/Morph.jl")
using IncompressibleAdjoint
using Plots
using GridapGmsh
using Parameters
using Revise

Reynolds = 1_000

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
    :AoA=>2.5,
)

xx = vcat(collect(0.0:0.0001:0.1),collect(0.101:0.001:1.0))

airfoil_cst = AirfoilCST(CST_NACA0012(;N=15,t= 0.0), 0.0)
des_points = AirfoilCSTDesign(airfoil_cst,xx)


modelname =create_msh(des_points; AoA=params[:AoA], mesh_ref=4.0)
model = GmshDiscreteModel(modelname)
modelgrid0 = copy(model.grid.node_coordinates)
writevtk(model,"model_original")




function compute_drag(uh,ph,nΓ,dΓ,params)
    CD, _ = compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)
    return CD, CD
end

function compute_lift(uh,ph,nΓ,dΓ,params)
    _, CL = compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)
    return CL,CL
end

params[:d_boundary]=VectorValue(-1.0,0.0)
params[:dt] = 1.0

uh,ph = solve_inc_primal_s(model, params; filename="inc-steady-SUPG")

using IncompressibleAdjoint
uhadj,phadj = solve_inc_adj_s(model, uh, ph,params; filename="res-adj-steady")


include("../src/Morph.jl")

control_points_position = [collect(0.001:0.002:0.01);collect(0.02:0.01:0.1);collect(0.1:0.1:0.9) ]

bnodes_top,bnodes_bottom,bnormals_top,bnormals_bottom = get_control_boundary(control_points_position,model,params)

bnodes=bnodes_top
bnormals=bnormals_top


using Plots


Nc = length(bnodes)
J1=zeros(Nc)
J2_1=zeros(Nc)
J2_2=zeros(Nc)
J2_3=zeros(Nc)
J2_4=zeros(Nc)
J2_5=zeros(Nc)

using IncompressibleAdjoint
updatekey(params,:i,0)

J1ref,(J2_1ref,J2_2ref,J2_3ref,J2_4ref,J2_5ref)= IncompressibleAdjoint.compute_sensitivity(model,params, uh,ph,uhadj,phadj; objective_function=compute_drag)


δ=0.01
for (i,(node,normal)) in enumerate(zip(bnodes,bnormals))
i = 1
node = bnodes[i]
normal = bnormals[i]
    sv = get_shift_vec(model, node, normal.*δ)
    updatekey(params,:i,i)

    model.grid.node_coordinates .= model.grid.node_coordinates .+ sv
    J1tmp,(J2_1tmp,J2_2tmp,J2_3tmp,J2_4tmp,J2_5tmp)= IncompressibleAdjoint.compute_sensitivity(model,params, uh,ph,uhadj,phadj; objective_function=compute_drag)

    J1[i] = J1tmp
    J2_1[i] = J2_1tmp
    J2_2[i] = J2_2tmp
    J2_3[i] = J2_3tmp
    J2_4[i] = J2_4tmp
    J2_5[i] = J2_5tmp

    println(i)
    model.grid.node_coordinates .=modelgrid0
end


diff = model.grid.node_coordinates .- modelgrid0
sum(norm.(diff) .>0)


J1s = (J1 .- J1ref)./δ
J2s1 = (J2_1 .- J2_1ref)./δ
J2s2 = (J2_2 .- J2_2ref)./δ
J2s3 = (J2_3 .- J2_3ref)./δ
J2s4 = (J2_4 .- J2_4ref)./δ
J2s5 = (J2_5 .- J2_5ref)./δ


J2s = J2s1 + J2s2 + J2s3 + J2s4 + J2s5
Jtot = J1s+ J2s
plot(Jtot)


plot(J1s)
plot(J2s1)
plot!(J2s2)
plot!(J2s3)
plot!(J2s4)

plot(J2s)
include("../src/Morph.jl")

### Finite Differences
_,CDref = IncompressibleAdjoint.obj_fun(model, params, uh,ph, compute_drag)
CDFD_F=zeros(Nc)
get_shift_vec(model, bnodes[1], bnormals[1].*δ)

for (i,(node,normal)) in enumerate(zip(bnodes,bnormals))
    sv = get_shift_vec(model, node, normal.*δ)
    
    model.grid.node_coordinates .= model.grid.node_coordinates .+ sv
    uh0tmp,ph0tmp = solve_inc_primal_s(model, params; filename=nothing)
    writevtk(params[:Ω], "MeshPerturb/model-$(100+i)", cellfields=["uh"=>uh0tmp,"ph"=>ph0tmp])

    _,CDtmp = IncompressibleAdjoint.obj_fun(model, params, uh0tmp,ph0tmp, compute_drag)
    CDFD_F[i] = CDtmp
    println(i)
    model.grid.node_coordinates .=modelgrid0
end

CDFD_B=zeros(Nc)

for (i,(node,normal)) in enumerate(zip(bnodes,bnormals))
    sv = get_shift_vec(model, node, -normal.*δ)
    
    model.grid.node_coordinates .= model.grid.node_coordinates .+ sv
    uh0tmp,ph0tmp = solve_inc_primal_s(model, params; filename=nothing)
    writevtk(params[:Ω], "MeshPerturb/model-$(100+i)", cellfields=["uh"=>uh0tmp,"ph"=>ph0tmp])

    _,CDtmp = IncompressibleAdjoint.obj_fun(model, params, uh0tmp,ph0tmp, compute_drag)
    CDFD_B[i] = CDtmp
    println(i)
    model.grid.node_coordinates .=modelgrid0
end


GradFD = (CDFD_F .- CDFD_B)./(2*δ)

GradFD
Jtot



plot(GradFD,seriestype=:scatter,markercolor=:black, label="Finite Difference")
plot!(GradFD,linecolor=:black,label=false)
plot!(Jtot,seriestype=:scatter,markercolor=:red,label="Adjoint")
plot!(Jtot,linecolor=:red,label=false)

