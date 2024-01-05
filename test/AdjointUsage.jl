using IncompressibleAdjoint
using Gridap,GridapGmsh
using Parameters

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
)



@unpack order,tagname,ν,dt = params

airfoil_cst = AirfoilCST(CST_NACA0012(), 0.0)
xx = collect(0:0.001:1)
des_points = AirfoilCSTDesign(airfoil_cst,xx)


modelname =create_msh(des_points; AoA=4.0, mesh_ref=4.0)
model = GmshDiscreteModel(modelname)

(uh0,UH),(ph0,PH)=solve_inc_primal(model, params; filename="res-unsteady")


# uhadj,phadj = solve_inc_adj_s(model, (uh0,UH), (ph0,PH),params; filename="res-adj-steady")

# using IncompressibleAdjoint
# Uhadj,Phadj = solve_inc_adj_u(model, (uh0,UH),(ph0,PH), (uhadj,phadj),params; filename="res-adj-unsteady")


function compute_lift(uh,ph,nΓ,dΓ,params)
    _, CL = compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)
    return 0.5*(CL-0.1)^2, CL
end

function compute_drag(uh,ph,nΓ,dΓ,params)
    CD, _ = compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)
    return CD, CD
end

function compute_efficiency(uh,ph,nΓ,dΓ,params)
    CD, CL = compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)
    return Drag - 0.1*Lift, CL/CD
end



istantaneus_CL_CD(model, params, (uh0,UH),(ph0,PH))


