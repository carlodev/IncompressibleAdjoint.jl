using IncompressibleAdjoint
using Gridap,GridapGmsh
using Parameters
using Statistics
using Revise
using Plots

#Number of Process to add
using SharedArrays
using Distributed
# np = 8
# addprocs(np)

"""
Trying to re-create the result obtained:
Castro, C., Lozano, C., Palacios, F., & Zuazua, E. (2007). Systematic continuous adjoint approach to viscous aerodynamic design on unstructured grids. AIAA Journal, 45(9), 2125–2139. https://doi.org/10.2514/1.24859

Sorgiovanni, G., Quadrio, M., & Ponzini, R. (2016). A robust open-source adjoint optimization method for external aerodynamics.

Direct Differentiation as in 
Janssens -P Vandenschrick -K Stevens -G Alessi, B. (n.d.). THE CONTINUOUS ADJOINT APPROACH APPLIED TO THE STABILIZED FINITE-ELEMENT FORMULATION OF THE INCOMPRESSIBLE NAVIER-STOKES EQUATIONS. www.euroturbo.eu
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

modelname =create_msh(des_points; AoA=2.5, mesh_ref=4.0)
model = GmshDiscreteModel(modelname)
modelgrid0 = copy(model.grid.node_coordinates)

uh,ph = solve_inc_primal_s(model, params; filename="inc-steady-SUPG")


control_points_position = [collect(0.001:0.002:0.01);collect(0.02:0.01:0.1);collect(0.1:0.1:0.9) ]

bnodes_top,bnodes_bottom,bnormals_top,bnormals_bottom = get_control_boundary(control_points_position,model,params)

δ = 0.01
updatekey(params,:δ,δ)

Γ = BoundaryTriangulation(model; tags=tagname)
dΓ = Measure(Γ, 8)
nΓ = -1 .* get_normal_vector(Γ)


GradDF = SharedArray{Float64}(length(bnodes_top))


## Parallelized Loop
#@distributed
for i in 1:1:length(bnodes_top)

point = bnodes_top[i]
nb = bnormals_top[i]

uhb,phb = solve_inc_direct_differentiation_s(model,uh, params, point, nb; filename="DirectDifferentiation/inc-direct-diff$i")
#Drag Sensitivity
dJdB = sum(∫( -phb ⋅ nΓ ⋅VectorValue(1.0,0.0))dΓ )
println(dJdB)
GradDF[i] = dJdB
end


GradDF = [-0.2621306345626182
-0.48837494633399264
-0.5483345905356058
-0.5510878895261202
-0.5392293958599353
-0.3455177047900573
-0.2463536263675509
-0.17127723353459423
-0.1294737392751218
-0.09646698639202629
-0.07460575131793601
-0.05814327390057963
-0.04721502479150033
-0.03834860988428915
-0.03834860988428915
-0.00039337665944287457
0.0027710117457030127
0.0019769775600276917
0.0009248128298273695
0.0003330578356437004
9.505630656367683e-5
-1.7568313172498394e-6
-2.260843774552878e-5]


plot(-GradDF, seriestype=:scatter, label="Direct Differentiation")
plot!(Jtot/2, seriestype=:scatter, label="Adjoint")
plot!(xlabel="design variable", ylabel="CD gradient")
savefig("Adjoint_DirDiff.pdf")