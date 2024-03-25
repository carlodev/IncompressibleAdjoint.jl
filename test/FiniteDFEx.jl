using IncompressibleAdjoint
using Gridap,GridapGmsh
using Parameters
using Statistics
using Revise
using Plots

"""
Finite Difference Validation
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
    :tf=>15.0,
    :t_endramp=>0.0,
    :t_primal_start_avg=>5.0,
    :θ=>1.0,
    :d_boundary=>VectorValue(-1.0,0.0), # It is the boundary condition for the adjoint problem at tag `:tagname`
    :AoA=>2.5,
    :Cᵢ=>[4,36],
)


@unpack AoA,ν,tagname = params


airfoil_cst = AirfoilCST(CST_NACA0012(;N=15,t= 0.0), 0.0)
xx = collect(0:0.0001:1)
des_points = AirfoilCSTDesign(airfoil_cst,xx)

modelname =create_msh(des_points; AoA=AoA, mesh_ref=4.0)
modelFD = GmshDiscreteModel(modelname)
writevtk(modelFD_tmp, "model_original")

#Store the node_coordinates of the original model
modelgrid0FD = copy(modelFD.grid.node_coordinates)




#Get Control points from the model
control_points_position = [collect(0.001:0.002:0.01);collect(0.02:0.01:0.1);collect(0.1:0.1:0.9) ]

bnodes_top,bnodes_bottom,bnormals_top,bnormals_bottom = get_control_boundary(control_points_position,modelFD,params)

bnodes=bnodes_top
bnormals=bnormals_top

using JLD2
load("NACA0012_AoA2p5_1000.jld2")["UH"][end]
load("NACA0012_AoA2p5_1000.jld2")["PH"][end]

params[:dt] = 0.05
params[:tF] = 5.0
uhFD,phFD = solve_inc_primal_s(modelFD, params; filename="res-steady")


uhFD.free_values .= load("NACA0012_AoA2p5_1000.jld2")["UH"][end]
phFD.free_values .= load("NACA0012_AoA2p5_1000.jld2")["PH"][end]


des_tmp = IncompressibleAdjoint.perturb_DesignParameters(des_points, 1, 0.01)

modelname_tmp =create_msh(des_tmp; AoA=AoA, mesh_ref=4.0)
modelFD_tmp = GmshDiscreteModel(modelname_tmp)
writevtk(modelFD_tmp, "model_cst")


_,CLref = IncompressibleAdjoint.obj_fun(modelFD, params, uhFD,phFD, compute_CL)

CLFD_F=zeros(length(bnodes_top))
δFD = 0.0001
for (i,(node,normal)) in enumerate(zip(bnodes_top,bnormals_top))
    sv = get_shift_vec(modelFD, node, normal.*δFD)
    
    modelFD.grid.node_coordinates .= modelFD.grid.node_coordinates .+ sv
    (uh0tmp,duhdt, UH, DUHDT), (ph0tmp, PH) = solve_inc_primal_u(modelFD, params; filename="res-unsteady")

    _,CDtmp = IncompressibleAdjoint.obj_fun(modelFD, params, uh0tmp,ph0tmp, compute_CL)
    CLFD_F[i] = CDtmp
    println(i)
    modelFD.grid.node_coordinates .=modelgrid0FD
end

