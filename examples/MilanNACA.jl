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

The Adjoint is computed as explained in the last paper.

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
model = GmshDiscreteModel(modelname)

#Store the node_coordinates of the original model
modelgrid0 = copy(model.grid.node_coordinates)



uh,ph = solve_inc_primal_s(model, params; filename="inc-steady-VMS")
using JLD2

uh.free_values .= load("NACA0012_AoA2p5_1000.jld2")["UH"][end]
ph.free_values .= load("NACA0012_AoA2p5_1000.jld2")["PH"][end]

params[:d_boundary] = VectorValue(0.0,-1.0)
uhadj,phadj = solve_inc_adj_s(model, uh, ph,params; filename="res-adj-steady")


#Get Control points from the model
control_points_position = [collect(0.001:0.002:0.01);collect(0.02:0.01:0.1);collect(0.1:0.1:0.9) ]

bnodes_top,bnodes_bottom,bnormals_top,bnormals_bottom = get_control_boundary(control_points_position,model,params)




plot(getindex.(bnodes_top,1),getindex.(bnodes_top,2),seriestype=:scatter, label="control points - TOP")
plot!(getindex.(bnodes_bottom,1),getindex.(bnodes_bottom,2),seriestype=:scatter,  label="control points - TOP")

#######################################################
### Adjoint Sensitivity
#######################################################
Γ = BoundaryTriangulation(model; tags="airfoil")
dΓ = Measure(Γ, 8)
nΓ = -1 .* get_normal_vector(Γ)
tΓ = rotation ∘ nΓ #extract tangent

δ=0.001
GvecT = zeros(length(bnodes_top))

for (i,(point,bnorm)) in enumerate(zip(bnodes_top,bnormals_top))
function node_filter(x)
    R = 0.001 #get_radius_shift(bnorm.*δ)
    dist = compute_point_dist(x, point)
   return morph_kernel(dist,R)
end

G = ν *sum(∫(((transpose(∇(uh)) ⋅ nΓ) ⋅ tΓ)* ((transpose(∇(uhadj)) ⋅ nΓ) ⋅ tΓ) * node_filter )dΓ)
GvecT[i] = G
end



GvecB = Float64[]
for (point,bnorm) in zip(bnodes_bottom,bnormals_bottom)

    function node_filter(x)
        R = get_radius_shift(bnorm.*δ)
        dist = compute_point_dist(x, point)
       return δ .*morph_kernel(dist,R)
    end

G = ν *sum(∫(((transpose(∇(uh)) ⋅ nΓ) ⋅ tΓ)* ((transpose(∇(uhadj)) ⋅ nΓ) ⋅ tΓ) * node_filter )dΓ)
push!(GvecB,G)
end


plot(GvecT./ δ, seriestype=:scatter, label = "top points Sensitivities")
plot!(xlabel="desing variable", ylabel="CD gradient")

#######################################################
### Finite Difference Sensitivity
#######################################################modelnameFD =create_msh(des_points; AoA=2.5, mesh_ref=3.0)
#Same model as the Adjoint
modelFD = GmshDiscreteModel(modelname)

modelgrid0FD = copy(modelFD.grid.node_coordinates)
uhFD,phFD = solve_inc_primal_s(modelFD, params; filename="inc-steady-SUPG")


_,CDref = IncompressibleAdjoint.obj_fun(modelFD, params, uhFD,phFD, compute_drag)

CDFD_F=zeros(length(bnodes_top))
δFD = 0.0025
for (i,(node,normal)) in enumerate(zip(bnodes_top,bnormals_top))
    sv = get_shift_vec(modelFD, node, normal.*δFD)
    
    modelFD.grid.node_coordinates .= modelFD.grid.node_coordinates .+ sv
    uh0tmp,ph0tmp = solve_inc_primal_s(modelFD, params; filename=nothing)
    writevtk(params[:Ω], "MeshPerturb/model-$(100+i)", cellfields=["uh"=>uh0tmp,"ph"=>ph0tmp])

    _,CDtmp = IncompressibleAdjoint.obj_fun(modelFD, params, uh0tmp,ph0tmp, compute_drag)
    CDFD_F[i] = CDtmp
    println(i)
    modelFD.grid.node_coordinates .=modelgrid0FD
end



CDFD_B=zeros(length(bnodes_top))

for (i,(node,normal)) in enumerate(zip(bnodes_top,bnormals_top))
    sv = get_shift_vec(modelFD, node, -normal.*δFD)
    
    modelFD.grid.node_coordinates .= modelFD.grid.node_coordinates .+ sv
    uh0tmp,ph0tmp = solve_inc_primal_s(modelFD, params; filename=nothing)
    writevtk(params[:Ω], "MeshPerturb/model-$(100+i)", cellfields=["uh"=>uh0tmp,"ph"=>ph0tmp])

    _,CDtmp = IncompressibleAdjoint.obj_fun(modelFD, params, uh0tmp,ph0tmp, compute_drag)
    CDFD_B[i] = CDtmp
    println(i)
    modelFD.grid.node_coordinates .=modelgrid0FD
end


GradFD = (CDFD_F .- CDFD_B)./(2*δFD)

plot(Gvec./δ, seriestype=:scatter)
plot!(GradFD, seriestype=:scatter)


