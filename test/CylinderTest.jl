using IncompressibleAdjoint
using Gridap,GridapGmsh
using Parameters
using Statistics
using Revise

Reynolds = 50
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
    :d_boundary=>VectorValue(-1.0,0.0),
    :AoA=>0.0,

)



@unpack order,tagname,ν,dt = params
xx = collect(0:0.0001:1)

des_points = IncompressibleAdjoint.initialize_control_points(xx; chord = 1.0, initfun=circle)

modelname =create_msh(des_points; AoA=0.0, mesh_ref=4.0)
model = GmshDiscreteModel(modelname)

modelgrid0 = copy(model.grid.node_coordinates)

function compute_drag(uh,ph,nΓ,dΓ,params)
    CD, _ = compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)
    return CD, CD
end
writevtk(model, "model0")


uh,ph = solve_inc_primal_s(model, params; filename="inc-steady-SUPG")

params[:d_boundary] = VectorValue(-2,0)
uhadj,phadj = solve_inc_adj_s(model, uh, ph,params; filename="res-adj-steady")

Nc = IncompressibleAdjoint.get_designparameters_number(des_points)
tag = IncompressibleAdjoint.get_designparameters_tags(des_points)

vv = ones(Nc) #(tag .== "top") .- (tag .== "bottom")


δ=0.01

updatekey(params,:i,0)
using IncompressibleAdjoint

J1ref,(J2_1ref,J2_2ref,J2_3ref,J2_4ref,J2_5ref)= IncompressibleAdjoint.compute_sensitivity(model,params, uh,ph,uhadj,phadj; objective_function=compute_drag)




xmorph=collect(0.0:0.05:1.0)
ymorph = sqrt.(0.5^2 .-(xmorph .- 0.5).^2)

using LinearAlgebra
bnorm(x,y) = normalize([x-0.5, y])

bnodes = Point[]
bnormals = VectorValue[]

for (x,y) in zip(xmorph,ymorph)
    push!(bnodes, Point(x,y))
    xn,yn = bnorm(x,y)
    push!(bnormals, Point(xn,yn))
end



Ns = length(bnodes)
J1=zeros(Ns)
J2_1=zeros(Ns)
J2_2=zeros(Ns)
J2_3=zeros(Ns)
J2_4=zeros(Ns)
J2_5=zeros(Ns)
include("../src/Morph.jl")
for (i,(node,normal)) in enumerate(zip(bnodes,bnormals))
    updatekey(params,:i,i)
    
    sv = get_shift_vec(model, node, normal.*δ)
    model.grid.node_coordinates .= model.grid.node_coordinates .+ sv
    println(i)

    J1tmp,(J2_1tmp,J2_2tmp,J2_3tmp,J2_4tmp,J2_5tmp)= IncompressibleAdjoint.compute_sensitivity(model,params, uh,ph,uhadj,phadj; objective_function=compute_drag)

    J1[i] = J1tmp
    J2_1[i] = J2_1tmp
    J2_2[i] = J2_2tmp
    J2_3[i] = J2_3tmp
    J2_4[i] = J2_4tmp
    J2_5[i] = J2_5tmp
    model.grid.node_coordinates .=modelgrid0

end

J1s = (J1.- J1ref)./δ
J2s1 = (J2_1 .- J2_1ref)./δ
J2s2 = (J2_2 .- J2_2ref)./δ
J2s3 = (J2_3 .- J2_3ref)./δ
J2s4 = (J2_4 .- J2_4ref)./δ
J2s5 = (J2_5 .- J2_5ref)./δ


J2s = J2s1 + J2s2 + J2s3 + J2s4 + J2s5
Jtot = J1s * 0.5+ J2s




using Plots
plot(xmorph[2:end-1],Jtot[2:end-1])

plot!(xmorph[2:end-1],J2s[2:end-1])


### Using Finite Differences
_,CD0 = IncompressibleAdjoint.obj_fun(model, params, uh,ph, compute_drag)

CDFD = zeros(Ns)
for (i,(node,normal)) in enumerate(zip(bnodes,bnormals))
    updatekey(params,:i,i)
    
    sv = get_shift_vec(model, node, normal.*δ)
    model.grid.node_coordinates .= model.grid.node_coordinates .+ sv
    println(i)

    uh,ph = solve_inc_primal_s(model, params; filename="MeshPerturb/FD$i")
    _,CDtmp = IncompressibleAdjoint.obj_fun(model, params, uh,ph, compute_drag)
    CDFD[i] = CDtmp
    model.grid.node_coordinates .=modelgrid0

end

GradFD = (CDFD .- CD0) ./δ

plot(GradFD[2:end-1],seriestype=:scatter)
plot!(J2s[2:end-1],seriestype=:scatter)
plot!(J1s[2:end-1],seriestype=:scatter)
plot!(Jtot[2:end-1],seriestype=:scatter)

function rotation(n::VectorValue{2,Float64})
    n1, n2 = [n...]
    VectorValue(-n2, n1)
end
Γ = BoundaryTriangulation(model; tags="airfoil")
dΓ = Measure(Γ, 8)
nΓ = -1 .* get_normal_vector(Γ)

tΓ = rotation ∘ nΓ #extract tangent

G = ∫(ν* ((transpose(∇(uh)) ⋅ nΓ) ⋅ tΓ) * ((transpose(∇(uhadj)) ⋅ nΓ) ⋅ tΓ))dΓ

cf = Dict("G" => G)
ttrian = G.trian
f = (reffe) -> Gridap.Geometry.UnstructuredGrid(reffe)
ref_grids = map(f, Gridap.Geometry.get_reffes(ttrian))
visgrid = Gridap.Visualization.VisualizationGrid(ttrian, ref_grids)
pdata = Gridap.Visualization._prepare_pdata(ttrian, cf, visgrid.cell_to_refpoints)




ttrian0 = Γ.trian
ref_grids0 = map(f, Gridap.Geometry.get_reffes(ttrian0))
visgrid0 = Gridap.Visualization.VisualizationGrid(ttrian0, ref_grids0)
airfoil_points = visgrid0.sub_grid.node_coordinates

airfoil_points
pdata["G"]
using UnicodePlots
ps = scatterplot(getindex.(airfoil_points,1), pdata["G"])
scatterplot!(ps,getindex.(bnodes,1), GradFD .*20)

ps = scatterplot(getindex.(bnodes,1),J2s)
scatterplot!(ps,getindex.(bnodes,1),GradFD)
