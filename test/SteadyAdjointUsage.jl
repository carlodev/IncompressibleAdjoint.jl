using IncompressibleAdjoint
using Gridap,GridapGmsh
using Parameters

Reynolds = 1_000

params=Dict(
    :chord => 1.0,
    :D=>2,
    :Re => Reynolds,
    :dt => 0.005,
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
    :d_boundary=>VectorValue(0.0,-1.0),
)

airfoil_cst = AirfoilCST(CST_NACA0012(), 0.0)
xx = collect(0:0.001:1)
des_points = AirfoilCSTDesign(airfoil_cst,xx)


modelname =create_msh(des_points; AoA=4.0, mesh_ref=4.0)
model = GmshDiscreteModel(modelname)


function compute_lift(uh,ph,nΓ,dΓ,params)
    _, CL = compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)

    return 0.5*(CL-0.75)^2, CL
end

uh,ph = solve_inc_primal_s(model, params; filename="inc-steady")

# ϕu, ϕp = solve_inc_adj_s(model, (uh, [uh.free_values]), (ph, [ph.free_values]), params)

@unpack ν=params
Γ = BoundaryTriangulation(model; tags="airfoil")
dΓ = Measure(Γ, 2)
nΓ = -1 .* get_normal_vector(Γ)

IForce = ∫(-ph ⋅ nΓ + ν* transpose(∇(uh)) ⋅ nΓ)dΓ
D,L = sum(IForce)
compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)


IncompressibleAdjoint.obj_fun(model, params, uh,ph, compute_lift; degree=2)


Sv = Float64[]
FVv=Float64[]
Jtotv=Vector[]
CPvec=typeof(des_points)[]
α = 0.1

using IncompressibleAdjoint


for iter=0:1:2
    model, (J1s,J2s,Jtot), (fitnessval,St), des_points, (uh,ph) = iterate_optimization(uh,ph, model, des_points,params; iter=iter,α = α/log(iter+3), δ=5e-4, objective_function=compute_lift,detail=true)
    push!(Sv,St)
    push!(FVv,fitnessval)
    push!(Jtotv,Jtot)
    push!(CPvec,des_points)
end

using Plots
plot(Sv)