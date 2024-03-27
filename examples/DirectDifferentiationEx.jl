"""
Trying to re-create the result obtained:
Castro, C., Lozano, C., Palacios, F., & Zuazua, E. (2007). Systematic continuous adjoint approach to viscous aerodynamic design on unstructured grids. AIAA Journal, 45(9), 2125–2139. https://doi.org/10.2514/1.24859

Sorgiovanni, G., Quadrio, M., & Ponzini, R. (2016). A robust open-source adjoint optimization method for external aerodynamics.

Direct Differentiation as in 
Janssens -P Vandenschrick -K Stevens -G Alessi, B. (n.d.). THE CONTINUOUS ADJOINT APPROACH APPLIED TO THE STABILIZED FINITE-ELEMENT FORMULATION OF THE INCOMPRESSIBLE NAVIER-STOKES EQUATIONS. www.euroturbo.eu
"""

using IncompressibleAdjoint
using IncompressibleAdjoint.Geometry
using IncompressibleAdjoint.IncompressibleSolvers

using Gridap,GridapGmsh
using Parameters
using Statistics
using Revise
using Plots
using JLD2




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
    :tf=>10.0,
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

modelname =create_msh(des_points; AoA=2.5, mesh_ref=1.0)
model = GmshDiscreteModel(modelname)
modelgrid0 = copy(model.grid.node_coordinates)

uh,ph = solve_inc_primal_s(model, params; filename=nothing)

params[:tf]=0.1
(uh,duhdt, UH, DUHDT), (ph, PH) = solve_inc_primal_u(model, params; filename="res-unsteady")

average_field!(uh,UH[end-50:end])
average_field!(ph,PH[end-50:end])
CD,CL = average_CD_CL(model,params,(uh,UH),(ph,PH))


control_points_position = [collect(0.001:0.002:0.01);collect(0.02:0.01:0.1);collect(0.1:0.1:0.9) ]
bnodes_top,bnodes_bottom,bnormals_top,bnormals_bottom = get_control_boundary(control_points_position,model,params)

δ = 0.0025
updatekey(params,:δ,δ)

Γ = BoundaryTriangulation(model; tags=tagname)
dΓ = Measure(Γ, 8)
nΓ = -1 .* get_normal_vector(Γ)


control_nodes = [bnodes_top;bnodes_bottom]
control_normals = [bnormals_top;bnormals_bottom]
GradDDF = zeros(length(control_nodes))
GradLDF = zeros(length(control_normals))


### Direct Differentiation Loop
for i in 1:1:length(control_nodes)

point = control_nodes[i]
nb = control_normals[i]



uhb,phb = solve_inc_direct_differentiation_u(model,(uh,UH),(ph,PH), params, point, nb; filename="DirectDifferentiation/inc-direct-diff$i")
#Drag Sensitivity
dJDdB = sum(∫( -phb ⋅ nΓ ⋅VectorValue(1.0,0.0))dΓ )
dJLdB = sum(∫( -phb ⋅ nΓ ⋅VectorValue(0.0,1.0))dΓ )

println(dJDdB)
GradDDF[i] = dJDdB
GradLDF[i] = dJLdB
jldsave("results/DirectDiff_2.jld2"; GradDDF,GradLDF)

end





####POST PROCS
using Plots
GradDDF = load("results/DirectDiff_2.jld2")["GradDDF"]

plot(GradDDF./2)

plot(GradLDF,seriestype=:scatter)


for d in GradLDF
    println(d)
end

using XLSX, Plots, DataFrames
df = DataFrame(XLSX.readtable("examples/Res.xlsx", "DirectDiff"))

df_top = filter(:tag=> d->d=="top",df)
df_bottom = filter(:tag=> d->d=="bottom",df)

plot(df_top.Lift, seriestype=:scatter, label="Direct Differentiation - top")
plot!(J1s, seriestype=:scatter, label="Adjoint - top")

plot!(xlabel="design variable", ylabel="CD gradient")
savefig("Adjoint_DirDiff.pdf")