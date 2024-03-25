using IncompressibleAdjoint
using Gridap,GridapGmsh
using Parameters
using Statistics
using Revise
using Plots
using NLopt

using JLD2,Statistics

"""
NACA0012 - Reynolds 1000, AoA 2.5, CL adjoint optimization, aiming a targhet CL=0.35.
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
Optimization function, it aims for a CL=0.35.
Also the CL is obtained, to check its evolution
"""
function compute_CL(uh,ph,nΓ,dΓ,params)
    _, CL = compute_airfoil_coefficients(uh,ph,nΓ,dΓ,params)
    return 0.5*(CL-0.35)^2, CL
end



βv= [Vector{Float64}(undef, Ndes) for _ = 1:20]

nodes_top = [Vector{Float64}(undef, length(airfoil_nodes_top)) for _ = 1:20] 
nodes_bottom = [Vector{Float64}(undef, length(airfoil_nodes_bottom)) for _ = 1:20] 
cp_top = [Vector{Float64}(undef, length(airfoil_nodes_top)) for _ = 1:20] 
cp_bottom = [Vector{Float64}(undef, length(airfoil_nodes_bottom)) for _ = 1:20] 




merge!(params,Dict(:iter=>0,:i=>0, :fitnessval=>zeros(100), :CL=>zeros(100),:βv=>copy(βv),:fgrad=>copy(βv),
:des_points=>des_points,:Ndes=>Ndes,:shift=>shift,:obj_fun=>compute_CL,
:nodes_top=>nodes_top,:nodes_bottom=>nodes_bottom,:cp_top=>cp_top,:cp_bottom=>cp_bottom))


### Solve the primal unsteady flow

(uh,duhdt, UH, DUHDT), (ph, PH) = solve_inc_primal_u(model, params; filename="res-unsteady")


CL_vec = zeros(length(PH))

for (i,(uhval,phval)) in enumerate(zip(UH,PH))
copyto!(uh.free_values, uhval)
copyto!(ph.free_values, phval)
    fitnessval, CL  = IncompressibleAdjoint.obj_fun(model, params, uh,ph, compute_CL)
    CL_vec[i] = copy(CL)
end

CL_history = [Vector{Float64}(undef, length(PH)) for _ = 1:20] 
CL_history[1] = copy(CL_vec)



fitnessval, CL  = IncompressibleAdjoint.obj_fun(model, params, uh,ph, obj_fun)




### using uh, ph as initial solution for the new unsteady simulations
merge!(params,Dict(:uh00=>uh,:ph00=>ph, :UH=>UH,:PH=>PH, :CL_history=>CL_history))



"""
    optimization_loop(w::Vector,grad::Vector)

It receives:
-w: CST weights 
-grad: gradient dI/dw_i

It is using the solution at the previous iteration on the new geometry in order to reduce the initial transient
"""
function optimization_loop(w::Vector,grad::Vector)
 @unpack des_points,Ndes,shift,AoA, obj_fun, iter, 
            mesh_ref,uh00,ph00,UH,PH=params
    des_points = update_CST_weights(w,des_points)


    modelname =create_msh(des_points; AoA=AoA, mesh_ref=mesh_ref)
    model = GmshDiscreteModel(modelname)
    writevtk(model, "model_$iter")
    uh = uh00
    ph=ph00
  
    (uh,duhdt, UH, DUHDT), (ph, PH) = solve_inc_primal_u(model, params; filename="inc-results-$iter", uh00=uh00,ph00=ph00)
    
    
    
    #### Compute Unsteady CL
    CL_vec = zeros(length(PH))

    for (i,(uhval,phval)) in enumerate(zip(UH,PH))
    copyto!(uh.free_values, uhval)
    copyto!(ph.free_values, phval)
        fitnessval, CL  = IncompressibleAdjoint.obj_fun(model, params, uh,ph, compute_CL)
        CL_vec[i] = copy(CL)
    end

    params[:CL_history][iter+1] = copy(CL_vec)

    
    CD,CL = average_CD_CL(model,params,(uh,UH),(ph,PH))




    #### Adjoint Boundary Conditions
    params[:d_boundary] = VectorValue(0.0, (0.35-CL))
    

    ### run a steady adjoint on the average field
    average_field!(uh,UH[end-50:end])
    average_field!(ph,PH[end-50:end])

    #### average Cp
    airfoil_nodes_top, airfoil_nodes_bottom,_,_ = get_airfoil_characteristics(model, params; tag="airfoil")
    cp_top,cp_bottom = get_aerodynamic_features(params,model, uh,ph)
    params[:nodes_top][iter+1] = copy(getindex.(airfoil_nodes_top,1))
    params[:nodes_bottom][iter+1] =  copy(getindex.(airfoil_nodes_bottom,1))
    params[:cp_top][iter+1] =  copy(cp_top)
    params[:cp_bottom][iter+1] =  copy(cp_bottom)




    uhadj,phadj = solve_inc_adj_s(model, uh, ph,params; filename="res-adj-steady")

    J1ref,(J2_1ref,J2_2ref,J2_3ref,J2_4ref,J2_5ref)= IncompressibleAdjoint.compute_sensitivity(model,params, uh,ph,uhadj,phadj; objective_function=obj_fun)

    fitnessval, CL  = IncompressibleAdjoint.obj_fun(model, params, uh,ph, obj_fun)


    params[:iter]= iter +1 
    params[:βv][iter+1] = copy(w)

    params[:CL][iter+1] = copy(CL)
    params[:fitnessval][iter+1]= copy(fitnessval)
    J1=zeros(Ndes)
    J2_1=zeros(Ndes)
    J2_2=zeros(Ndes)
    J2_3=zeros(Ndes)
    J2_4=zeros(Ndes)
    J2_5=zeros(Ndes)

    copyto!(params[:uh00].free_values, uh.free_values)
    copyto!(params[:ph00].free_values, ph.free_values)


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
    

    grad .= J2s.*des_bool

    params[:fgrad][iter+1] = copy(grad)
    println("CL = $CL")
    jldsave("results/NACA0012_1000_2p5/parameters_2.jld2"; params)
    return J1ref

end



## Set the optimization algorithm
opt = Opt(:LD_LBFGS, Ndes)
## Set the bounds
opt.lower_bounds = [-0.1 .* ones(Int(Ndes/2));-1 .* ones(Int(Ndes/2))]
opt.upper_bounds = [ones(Int(Ndes/2)); 0.1 .* ones(Int(Ndes/2))]

## Set the tolerance
opt.xtol_rel = 2.5e-2

opt.min_objective = optimization_loop

#Initialization
w_init,_ = get_CST_values(des_points) 

###### Resolution
(minf,minx,ret) = optimize(opt,w_init)
#######

######################################################################
# Post processing
######################################################################
@unpack iter = params
iter
plot(params[:fitnessval][1:iter] .*10)
plot(params[:CL][1:iter])



plot(collect(1:1:iter) .-1,params[:CL][1:iter] ,linewidth=2.0,linecolor=:teal,label=false)
plot!(collect(1:1:iter) .-1,params[:CL][1:iter],seriestype=:scatter,markercolor=:teal,label=false)
plot!(collect(1:1:iter) .-1, 0.35.*ones(iter), linestyle=:dash,linecolor=:black, linewidth=1.2,label="target")
plot!(xlabel="iterations",ylabel="CL")

savefig("CL_0p35_conv.pdf")



#Plot series of airfoil


params[:βv][1]
@unpack iter=params

for i=1:1:iter
p_airfoil = update_CST_weights(params[:βv][i],des_points)
plot(p_airfoil.apoints.xu,p_airfoil.apoints.yu,linewidth=2,linecolor=:black,label=false)
plot!(p_airfoil.apoints.xl,p_airfoil.apoints.yl,linewidth=2,linecolor=:black,label=false)
plot!(aspect_ratio = :equal)
savefig("airfoil$i.pdf")
end




for i=1:1:iter

    plot(collect(1:1:i) .-1,params[:CL][1:i] ,linewidth=2.0,linecolor=:teal,label=false)
    plot!(collect(1:1:i) .-1,params[:CL][1:i],seriestype=:scatter,markercolor=:teal,label=false)
    plot!(collect(1:1:iter) .-1, 0.35.*ones(iter), linestyle=:dash,linecolor=:black, linewidth=1.2,label="target")
    plot!(xlabel="iterations",ylabel="CL",xlims=[0.0,8.1],ylims=[0.10,0.40])
    
    savefig("results/CL_$i.png")
end




#############################################
#Unsteady CL
##############################################
using Plots, LaTeXStrings, KissSmoothing
using JLD2
params = load("results/NACA0012_1000_2p5/parameters_2.jld2")["params"]
@unpack iter,dt,tf = params
ttime = collect(0.0:dt:tf)

plot(ttime,params[:CL_history][1], linewidth=1.5, linestyle=:dash, linecolor=:deepskyblue3 ,label="original")
plot!(ttime,params[:CL_history][iter],linewidth=1.5, linestyle=:solid, linecolor=:red ,label="final")
plot!(ylabel=L"C_L",xlabel="t [s]",legend=:topleft)
savefig("CL_035_unsteady.pdf")



plot(getindex.(params[:nodes_top][1],1),params[:cp_top][1], linewidth=1.5, linestyle=:dash, linecolor=:deepskyblue3 ,label="original")
plot!(getindex.(params[:nodes_bottom][1],1),params[:cp_bottom][1], linewidth=1.5, linestyle=:dash, linecolor=:deepskyblue3 ,label=false)

plot!(getindex.(params[:nodes_top][iter],1),denoise(params[:cp_top][iter];factor=0.5), linewidth=1.5, linestyle=:solid, linecolor=:red ,label="final")
plot!(getindex.(params[:nodes_bottom][iter],1),denoise(params[:cp_bottom][iter];factor=0.5), linewidth=1.5, linestyle=:solid, linecolor=:red ,label=false)

plot!(ylabel=L"\overline{C_P}",xlabel="x/c",legend=:bottomright)
yflip!()
savefig("Cp_035_initial-final.pdf")
