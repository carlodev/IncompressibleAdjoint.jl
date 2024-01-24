using IncompressibleAdjoint
using Gridap,GridapGmsh
using Parameters
using Statistics
using Revise

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
    :d_boundary=>VectorValue(-1.0,0.0),
    :AoA=>2.5,

)



@unpack order,tagname,ν,dt = params
airfoil_cst = AirfoilCST(CST_NACA0012(;N=15,t= 0.0), 0.0)
xx = collect(0:0.0001:1)
des_points = AirfoilCSTDesign(airfoil_cst,xx)

modelname =create_msh(des_points; AoA=2.5, mesh_ref=4.0)
model = GmshDiscreteModel(modelname)

modelgrid0 = copy(model.grid.node_coordinates)

function compute_drag(uh,ph,nΓ,dΓ,params)
    CD, _ = compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)
    return CD, CD
end

iter = 0.0
writevtk(model, "model$iter")

# (uh,duhdt, UH, DUHDT), (ph, PH)=solve_inc_primal_u(model, params; filename="res-unsteady")

# writevtk(params[:Ω], "uhph", cellfields=["uh"=>uh, "ph"=>ph])

# CDavg,CLavg = average_CD_CL(model, params, (uh,UH),(ph,PH))
# params[:d_boundary] = VectorValue(-1.0,0.0)

params[:d_boundary] = VectorValue(-1.0,0.0)
uh,ph = solve_inc_primal_s(model, params; filename="inc-steady-SUPG")

uhadj,phadj = solve_inc_adj_s(model, uh, ph,params; filename="res-adj-steady")



Nc = IncompressibleAdjoint.get_designparameters_number(des_points)
tag = IncompressibleAdjoint.get_designparameters_tags(des_points)
Nc = length(bnodes)


vv = ones(Nc) #(tag .== "top") .- (tag .== "bottom")

J1=zeros(Nc)
J2_1=zeros(Nc)
J2_2=zeros(Nc)
J2_3=zeros(Nc)
J2_4=zeros(Nc)
J2_5=zeros(Nc)

δ=0.01
shift = vv.*δ
updatekey(params,:i,0)
using IncompressibleAdjoint

J1ref,(J2_1ref,J2_2ref,J2_3ref,J2_4ref,J2_5ref)= IncompressibleAdjoint.compute_sensitivity(model,params, uh,ph,uhadj,phadj; objective_function=compute_drag)


for (i,ss) in enumerate(shift)
    updatekey(params,:i,i)

    des_tmp = IncompressibleAdjoint.perturb_DesignParameters(des_points, i, ss)
    modelname = create_msh(des_tmp; iter = i+100, mesh_ref=4.0,AoA=0.0)
    model_tmp = GmshDiscreteModel(modelname)
    # model.grid.node_coordinates .= model_tmp.grid.node_coordinates 
    # writevtk(model, "model-$(iter+100+i)")

    J1tmp,(J2_1tmp,J2_2tmp,J2_3tmp,J2_4tmp,J2_5tmp)= IncompressibleAdjoint.compute_sensitivity(model_tmp,params, uh,ph,uhadj,phadj; objective_function=compute_drag)

    J1[i] = J1tmp
    J2_1[i] = J2_1tmp
    J2_2[i] = J2_2tmp
    J2_3[i] = J2_3tmp
    J2_4[i] = J2_4tmp
    J2_5[i] = J2_5tmp

end

J1s = (J1 .- J1ref)./shift
J2s1 = (J2_1 .- J2_1ref)./shift
J2s2 = (J2_2 .- J2_2ref)./shift
J2s3 = (J2_3 .- J2_3ref)./shift
J2s4 = (J2_4 .- J2_4ref)./shift
J2s5 = (J2_5 .- J2_5ref)./shift


J2s = J2s1 + J2s2 + J2s3 + J2s4 + J2s5
Jtot = J1s+ J2s


using Plots


J2s
J2s2 + J2s3 + J2s4 +J2s5

plot!(J2s)

function rotation(n::VectorValue{2,Float64})
    n1, n2 = [n...]
    VectorValue(-n2, n1)
end

Γ = BoundaryTriangulation(model; tags="airfoil")
dΓ = Measure(Γ, 8)
nΓ = -1 .* get_normal_vector(Γ)

tΓ = rotation ∘ nΓ #extract tangent



include("../src/Morph.jl")
δ=0.01
GvecT = Float64[]

for point in bnodes
    println(point)

function morph_kernel(x)
    R = get_radius_shift(δ)

    dist = compute_point_dist(x, point)
    δ.*morph_kernel(dist,R)

end

G = ν *sum(∫(((transpose(∇(uh)) ⋅ nΓ) ⋅ tΓ)* ((transpose(∇(uhadj)) ⋅ nΓ) ⋅ tΓ) * morph_kernel )dΓ)
push!(Gvec,G)
end


plot(Gvec./ δ, seriestype=:scatter)


Ω = Triangulation(model)
dΩ = Measure(Ω, 8)





writevtk(Ω,"Kernel", cellfields=["Kernel"=>morph_kernel])
