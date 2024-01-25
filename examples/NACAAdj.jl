using IncompressibleAdjoint
using Gridap,GridapGmsh
using Parameters
using Statistics
using Revise
using Plots

"""
Trying to re-create the result obtained:
Castro, C., Lozano, C., Palacios, F., & Zuazua, E. (2007). Systematic continuous adjoint approach to viscous aerodynamic design on unstructured grids. AIAA Journal, 45(9), 2125–2139. https://doi.org/10.2514/1.24859

Sorgiovanni, G., Quadrio, M., & Ponzini, R. (2016). A robust open-source adjoint optimization method for external aerodynamics.

The Gradients are obtained with the classic adjoint FORMULATION

"""

Reynolds = 1000

params=Dict(
    :chord => 1.0,
    :u_in=>1.0,
    :D=>2,
    :Re => Reynolds,
    :dt => 0.25,
    :ν =>1 / Reynolds,
    :D =>2,
    :order => 1,
    :method=>:SUPG,
    :tagname=>"airfoil",
    :t0=>0.0,
    :tf=>15.0,
    :t_endramp=>0.0,
    :t_primal_start_avg=>5.0,
    :θ=>1.0,
    :d_boundary=>VectorValue(-1.0,0.0), # It is the boundary condition for the adjoint problem at tag `:tagname`
    :AoA=>2.5,
)


@unpack AoA,ν,tagname = params


airfoil_cst = AirfoilCST(CST_NACA0012(;N=15,t= 0.0), 0.0)
xx = collect(0:0.0001:1)
des_points = AirfoilCSTDesign(airfoil_cst,xx)

modelname =create_msh(des_points; AoA=AoA, mesh_ref=4.0)
model = GmshDiscreteModel(modelname)

#Store the node_coordinates of the original model
modelgrid0 = copy(model.grid.node_coordinates)

uh,ph = solve_inc_primal_s(model, params; filename="inc-steady-SUPG")

uhadj,phadj = solve_inc_adj_s(model, uh, ph,params; filename="res-adj-steady")


#Get Control points from the model
control_points_position = [collect(0.001:0.002:0.01);collect(0.02:0.01:0.1);collect(0.1:0.1:0.9) ]

bnodes_top,bnodes_bottom,bnormals_top,bnormals_bottom = get_control_boundary(control_points_position,model,params)

bnodes=bnodes_top
bnormals=bnormals_top



updatekey(params,:i,0)

J1ref,(J2_1ref,J2_2ref,J2_3ref,J2_4ref,J2_5ref)= IncompressibleAdjoint.compute_sensitivity(model,params, uh,ph,uhadj,phadj; objective_function=compute_drag)

Nc = length(bnodes_top)
J1=zeros(Nc)
J2_1=zeros(Nc)
J2_2=zeros(Nc)
J2_3=zeros(Nc)
J2_4=zeros(Nc)
J2_5=zeros(Nc)

δ=0.001
for (i,(node,normal)) in enumerate(zip(bnodes_top,bnormals_top))
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
J1s = (J1 .- J1ref)./δ
J2s1 = (J2_1 .- J2_1ref)./δ
J2s2 = (J2_2 .- J2_2ref)./δ
J2s3 = (J2_3 .- J2_3ref)./δ
J2s4 = (J2_4 .- J2_4ref)./δ
J2s5 = (J2_5 .- J2_5ref)./δ


J2s = J2s1 + J2s2 + J2s3 + J2s4 + J2s5
Jtot = J1s+ J2s
plot(Jtot, seriestype=:scatter)