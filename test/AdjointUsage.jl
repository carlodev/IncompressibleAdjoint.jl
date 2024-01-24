using IncompressibleAdjoint
using Gridap,GridapGmsh
using Parameters
using Statistics
using Revise

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
    :method=>:VMS,
    :tagname=>"airfoil",
    :t0=>0.0,
    :tf=>30.0,
    :t_endramp=>0.0,
    :t_primal_start_avg=>10.0,
    :θ=>1.0,
    :d_boundary=>VectorValue(0.0,1.0),

)



@unpack order,tagname,ν,dt = params

airfoil_cst = AirfoilCST(CST_NACA0012(;N=6), 0.0)
xx = collect(0:0.001:1)
des_points = AirfoilCSTDesign(airfoil_cst,xx)

spoints = vcat(collect(0.0:0.001:0.09),collect(0.1:0.01:1.0))
des_points = IncompressibleAdjoint.initialize_control_points(spoints;initfun=NACA00)

modelname =create_msh(des_points; AoA=2.5, mesh_ref=4.0)
model = GmshDiscreteModel(modelname)
CLvec= Float64[]
Jtotvec= Vector[]


function compute_lift(uh,ph,nΓ,dΓ,params)
    _, CL = compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)
    return 0.5*(CL-0.75)^2, CL
end


# for iter = 0:1:20
iter = 0.0
writevtk(model, "model$iter")

(uh,duhdt, UH, DUHDT), (ph, PH)=solve_inc_primal_u(model, params; filename="res-unsteady")




CDavg,CLavg = average_CD_CL(model, params, (uh,UH),(ph,PH))
params[:d_boundary] = VectorValue(-1.0,0.0)


uhadj,phadj = solve_inc_adj_s(model, uh, ph,params; filename="res-adj-steady")

Uhadj,Phadj = solve_inc_adj_u(model, (uh,UH),(ph,PH), (uhadj,phadj),params; filename="res-adj-unsteady")

UH

IDXAVG = collect([200:500]...)
function computeJ3(model::DiscreteModel)
    Ω = Triangulation(model)
    dΩ = Measure(Ω, 4)
    unst = Float64[]
    stead_qdm = Float64[]
    stead_cont = Float64[]
for (Δu,p,u, uadj,padj) in zip(DUHDT[IDXAVG],PH[IDXAVG],UH[IDXAVG], reverse(Uhadj)[IDXAVG],reverse(Phadj)[IDXAVG])
    copyto!(duhdt.free_values,Δu)
    copyto!(uh.free_values,u)
    copyto!(ph.free_values,p)
    copyto!(uhadj.free_values,uadj)
    copyto!(phadj.free_values,padj)

    push!(unst,copy(sum(∫(uhadj⋅duhdt)dΩ)))
    push!(stead_qdm,copy(sum(∫(uhadj⋅((∇(uh))'⋅uh+∇(ph)))dΩ)))
    push!(stead_cont,copy(sum(∫(phadj⋅(∇⋅(uh)))dΩ)))

end

tcontrib = unst+stead_qdm+stead_cont

return Statistics.mean(tcontrib)
end

using Plots


J3ref = computeJ3(model)



modelgrid0 = copy(model.grid.node_coordinates)
Nc = IncompressibleAdjoint.get_designparameters_number(des_points)
tag = IncompressibleAdjoint.get_designparameters_tags(des_points)

Γ = BoundaryTriangulation(model; tags=params[:tagname])
dΓ = Measure(Γ, 4)
nΓ = -1 .* get_normal_vector(Γ)
Vector

Γ
ShiftVec = Vector[]

for (x,y) in zip(des_points.xcoordinate,des_points.ycoordinate)
    try vector_value = nΓ(Point(x,y))
    push!(ShiftVec,   [vector_value[1],vector_value[2]])
    catch
        push!(ShiftVec,   ShiftVec[end])

    end
    println(x,y)
end

ShiftVec

vv = ones(Nc) #(tag .== "top") .- (tag .== "bottom")

J1=zeros(Nc)
J2=zeros(Nc)

J3=zeros(Nc)

δ=0.001
shift = vv.*δ

for (i,ss) in enumerate(shift)
    des_tmp = IncompressibleAdjoint.perturb_DesignParameters(des_points, i, ss)
    modelname = create_msh(des_tmp; iter = i+100, mesh_ref=4,AoA=4.0)
    model_tmp = GmshDiscreteModel(modelname)
    model.grid.node_coordinates .= model_tmp.grid.node_coordinates 
    J3tmp = computeJ3(model)
    J3[i] = copy(J3tmp)
    model.grid.node_coordinates .= modelgrid0
end


average_field!(uh,UH[IDXAVG])
average_field!(ph,PH[IDXAVG])

average_field!(uhadj,reverse(Uhadj)[IDXAVG])
average_field!(phadj,reverse(Phadj)[IDXAVG])

fitnessval, S = IncompressibleAdjoint.obj_fun(model,params,uh, ph,compute_lift)

J1ref,J2ref= IncompressibleAdjoint.compute_sensitivity(model,params, uh,ph,uhadj, phadj; objective_function=compute_drag)
J2ref
J3ref

push!(CLvec,S)

average_field!(duhdt,DUHDT[IDXAVG])

writevtk(params[:Ω], "Avg", cellfields=["uh"=>uh,"ph"=>ph,"uhadj"=>uhadj,"phadj"=>phadj,"duhdt"=>duhdt])


for (i,ss) in enumerate(shift)
        des_tmp = IncompressibleAdjoint.perturb_DesignParameters(des_points, i, ss)
        modelname = create_msh(des_tmp; iter = i+100, mesh_ref=4,AoA=2.5)
        model_tmp = GmshDiscreteModel(modelname)
        model.grid.node_coordinates .= model_tmp.grid.node_coordinates 
       J1tmp,J2tmp= IncompressibleAdjoint.compute_sensitivity(model,params, uh,ph,uhadj,phadj; objective_function=compute_drag)
       writevtk(model, "model-$(iter+100+i)")

       J1[i] = J1tmp
       J2[i] = J2tmp

       model.grid.node_coordinates .= modelgrid0
end


J1s = (J1 .- J1ref)./ (δ) #Geometric gradient
J2s =  (J2 .- J2ref)./(δ)

J3s =  (J3 .- J3ref)./(δ)

Jtot = J1s +J2s #Total Gradient
plot(Jtot)

Jtotu =  J1s +J3s

plot(Jtot)

plot(vcat(reverse(Jtotu[7:end]),Jtotu[1:6]))
plot!(xticks = (1:12, string.(1:12)))


α=1.5/(iter+1)
#Step, -1 because opposite direction of the gradient
ΔD = IncompressibleAdjoint.ShiftUpdate(-vv.*α.*Jtot, tag)
iter = iter+1
contr_new = IncompressibleAdjoint.perturb_DesignParameters(des_points, ΔD)
modelname = create_msh(contr_new; iter = iter+1, mesh_ref=4, AoA=4.0)
model = GmshDiscreteModel(modelname)
des_points = contr_new

    # end



CLvec

using Plots
plot(CLvec)



# using IncompressibleAdjoint

function compute_drag(uh,ph,nΓ,dΓ,params)
    CD, _ = compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)
    return CD, CD
end

function compute_efficiency(uh,ph,nΓ,dΓ,params)
    CD, CL = compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)
    return Drag - 0.1*Lift, CL/CD
end


using IncompressibleAdjoint

CDv,CLv =istantaneus_CD_CL(model, params, (uh0,UH),(ph0,PH))
