using IncompressibleAdjoint
using Gridap,GridapGmsh
using Parameters
using Statistics
using Revise
using Plots
using NLopt

using JLD2,Statistics

"""
NACA0012 - Reynolds 1000, AoA 2.5, gradients validation using Finite Differences.
The parametrization is CST as β shape parameters
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
    :tf=>10.0,
    :t_endramp=>0.0,
    :t_primal_start_avg=>5.0,
    :θ=>1.0,
    :d_boundary=>VectorValue(-1.0,0.0), # It is the boundary condition for the adjoint problem at tag `:tagname`
    :AoA=>2.5,
    :Cᵢ=>[4,36],
    :mesh_ref=>1.0,
)


@unpack AoA,ν,tagname,mesh_ref = params

airfoil_cst = AirfoilCST(CST_NACA0012(;N=15,t= 0.0), 0.0)
xx = collect(0:0.0001:1)
des_points = AirfoilCSTDesign(airfoil_cst,xx)

modelname =create_msh(des_points; AoA=AoA, mesh_ref=mesh_ref)
model = GmshDiscreteModel(modelname)
writevtk(model, "model_original")


#Store the node_coordinates of the original model
modelgrid0 = copy(model.grid.node_coordinates)

Ndes = IncompressibleAdjoint.get_designparameters_number(des_points)
des_tag = IncompressibleAdjoint.get_designparameters_tags(des_points)

### Deformation used to compute the gradients with the adjoint
δ=0.01
des_bool = (des_tag .== "top") .- (des_tag .== "bottom")
shift = des_bool.*δ


"""
Optimization function, it checks the value of CL.
Also the CL is obtained, to check its evolution
"""
function compute_CL(uh,ph,nΓ,dΓ,params)
    _, CL = compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)
    return CL, CL
end



βv= [Vector{Float64}(undef, Ndes) for _ = 1:100]

merge!(params,Dict(:iter=>0,:i=>0, :fitnessval=>zeros(100), :CL=>zeros(100),:βv=>copy(βv),:fgrad=>copy(βv),
:des_points=>des_points,:Ndes=>Ndes,:shift=>shift,:obj_fun=>compute_CL))


uh,ph = solve_inc_primal_s(model, params; filename=nothing)


uh0=load("results/NACA0012_1000_2p5/NACA0012_AoA2p5_1000_ref2.jld2")["uh"]
ph0=load("results/NACA0012_1000_2p5/NACA0012_AoA2p5_1000_ref2.jld2")["ph"]
UH=load("results/NACA0012_1000_2p5/NACA0012_AoA2p5_1000_ref2.jld2")["UH"]
PH=load("results/NACA0012_1000_2p5/NACA0012_AoA2p5_1000_ref2.jld2")["PH"]

copyto!(uh.free_values, uh0.free_values)
copyto!(ph.free_values, ph0.free_values)
### Solve the primal unsteady flow

# (uh,duhdt, UH, DUHDT), (ph, PH) = solve_inc_primal_u(model, params; filename="res-unsteady")
# CD,CL = average_CD_CL(model,params,(uh,UH),(ph,PH))

### Solve the adjoint flow
average_field!(uh,UH[end-50:end])
average_field!(ph,PH[end-50:end])
fitnessval, CL  = IncompressibleAdjoint.obj_fun(model, params, uh,ph, compute_CL)

params[:d_boundary] = VectorValue(0.0, CL)

uhadj,phadj = solve_inc_adj_s(model, uh, ph,params; filename="res-adj-steady")
J1ref,(J2_1ref,J2_2ref,J2_3ref,J2_4ref,J2_5ref)= IncompressibleAdjoint.compute_sensitivity(model,params, uh,ph,uhadj,phadj; objective_function=compute_CL)
J1=zeros(Ndes)
J2_1=zeros(Ndes)
J2_2=zeros(Ndes)
J2_3=zeros(Ndes)
J2_4=zeros(Ndes)
J2_5=zeros(Ndes)


for (i,ss) in enumerate(shift)

    des_temp = IncompressibleAdjoint.perturb_DesignParameters(des_points, i, ss)
    modelname_tmp =create_msh(des_temp; AoA=AoA, mesh_ref=mesh_ref)
    model_tmp = GmshDiscreteModel(modelname_tmp)

    J1tmp,(J2_1tmp,J2_2tmp,J2_3tmp,J2_4tmp,J2_5tmp)= IncompressibleAdjoint.compute_sensitivity(model_tmp,params, uh,ph,uhadj,phadj; objective_function=compute_CL)


    J1[i] = J1tmp
    J2_1[i] = J2_1tmp
    J2_2[i] = J2_2tmp
    J2_3[i] = J2_3tmp
    J2_4[i] = J2_4tmp
    J2_5[i] = J2_5tmp
    println(i)
end

J1s = (J1 .- J1ref)./δ
J2s1 = (J2_1 .- J2_1ref)./δ
J2s2 = (J2_2 .- J2_2ref)./δ
J2s3 = (J2_3 .- J2_3ref)./δ
J2s4 = (J2_4 .- J2_4ref)./δ
J2s5 = (J2_5 .- J2_5ref)./δ


J2s = J2s1 + J2s2 + J2s3 + J2s4 + J2s5
Jtot = J1s+ J2s


adj_grad = J2s.*des_bool


plot(adj_grad)
plot!(Jtot.*des_bool)



"""
Finite Differences Validation
"""


CL_FD = zeros(length(shift))
for (i,ss) in enumerate(shift)

    des_temp = IncompressibleAdjoint.perturb_DesignParameters(des_points, i, ss)
    modelname_tmp =create_msh(des_temp; AoA=AoA, mesh_ref=mesh_ref)
    model_tmp = GmshDiscreteModel(modelname_tmp)

    (uh,duhdt, UH, DUHDT), (ph, PH) = solve_inc_primal_u(model_tmp, params; filename="FD")
    average_field!(uh,UH[end-50:end])
    average_field!(ph,PH[end-50:end])
    fitnessval, CL  = IncompressibleAdjoint.obj_fun(model_tmp, params, uh,ph, compute_CL)
    
    CL_FD[i] = CL
    println("CL $CL, shape params = $i")
    jldsave("results/NACA0012_1000_2p5/CL_FD.jld2"; CL_FD)
end

FD_grad = (CL_FD .- CL) ./ shift

plot(-adj_grad)
plot!(FD_grad)

