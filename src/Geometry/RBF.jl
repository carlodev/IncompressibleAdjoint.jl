using Plots
using IncompressibleAdjoint
using Gridap, GridapGmsh

function fRBF(x::Real,r::Real, fname::Function)
    x_norm = x/r
    if x_norm>1
        return 0
    else
        return fname(x/r)
    end
end

function fRBF(x::AbstractVector,r::Real, fname::Function)
    map(v-> fRBF(v,r, fname), x)
end


function fRBF(x::Real, fname::Function)
    return fname(x)
    
end

function fRBF(x::AbstractVector, fname::Function)
    map(v-> fRBF(v,fname), x)
end




function RBF_CP0(x)
    return (1-x)^2
end

function RBF_CP2(x)
    return (1-x)^4 * (4*x+1)
end

function RBF_CP4(x)
    return (1-x)^6 * (35/3 * x^2 + 6 * x +1)
end

function RBF_CP6(x)
    return (1-x)^8 * (32 * x^3 + 25*x^2 + 8 * x +1)
end

function RBF_CTPS0(x)
    return (1-x)^5
end

function RBF_CTPS1(x)
    return 1 + 80/3*x^2 - 40*x^3 + 15* x^4 - 8/3 * x^5 + 20 * x ^2 * log(x)
end

function RBF_IQBM(x)
    return 1 /(1+x^2)
end

function RBF_GAUSS(x)
    return  exp(-x^2)
end


#### Interpolation

"""
    SPF_Morph()

xb:: boundary/control points coordinates
xc:: nodes coordinates of the model
disp:: prescribed displacement of each boundary point
"""
function SPF_Morph(xb::Vector{VectorValue}, xc::Vector{VectorValue}, disp::Vector{VectorValue})


end


params=Dict(
    :chord => 1.0,
    :u_in=>1.0,
    :D=>2,
    :dt => 0.25,
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
airfoil_cst = AirfoilCST(CST_NACA0012(;N=15,t= 0.0), 0.0)
xx = collect(0:0.0001:1)
des_points = AirfoilCSTDesign(airfoil_cst,xx)

modelname =create_msh(des_points; AoA=4.0, mesh_ref=4.0)
model = GmshDiscreteModel(modelname)

#Store the node_coordinates of the original model
modelgrid0 = copy(model.grid.node_coordinates)



control_points_position = [collect(0.001:0.002:0.01);collect(0.02:0.01:0.1);collect(0.1:0.1:0.9) ]
bnodes_top,bnodes_bottom,bnormals_top,bnormals_bottom = get_control_boundary(control_points_position,model,params)





#### RBF dimension for dimension


using LinearAlgebra, IterativeSolvers, Preconditioners,  IncompleteLU

Nc = length(modelgrid0)
Nb = length([bnodes_top;bnodes_bottom;binlet;boutlet;blimits])
ND = length(bnodes_top[1])
Mbb = ones(Nb,Nb)
# r = sqrt(2) * δ * ND

binlet=IncompressibleAdjoint.get_boundary_nodes(model, "inlet")
boutlet=IncompressibleAdjoint.get_boundary_nodes(model, "outlet")
blimits=IncompressibleAdjoint.get_boundary_nodes(model, "limits")


bnodes = [bnodes_top;bnodes_bottom;binlet;boutlet;blimits]
bnormals = [bnormals_top;bnormals_bottom]
DistM = zeros(Nc,Nb,ND)


s = typeof(modelgrid0)(undef, length(modelgrid0))
for (j,p) in enumerate(modelgrid0)
    for dd in 1:1:ND
        dist = norm.(p[dd] .- getindex.(bnodes,dd))
        DistM[j,:,dd]= fRBF(dist,RBF_IQBM)
    end
end

Shift = typeof(modelgrid0)(undef, length(modelgrid0))

δ = 0.01
support_radius = 0.25
ss = zeros(Nc,ND)

for direction = 1:1:ND


for (i,bn) in enumerate(bnodes)
    dist = norm.(bn[direction] .- getindex.(bnodes,[direction]))
    
    Mbb[i,:]= fRBF(dist,RBF_IQBM)
end


Pb = ones(Nb, ND)
for (i,bn) in enumerate(bnodes)
    Pb[i,:]= [1, bn[direction]]
end

D = hcat(vcat(Mbb,Pb'),vcat(Pb,zeros(ND,ND)))

db=zeros(Nb)
vvd= (direction == 2) ? ones(length(bnormals_top)).* δ : zeros(length(bnormals_top))
db[1:length(bnormals_top)] = vvd
d0 = [db;zeros(ND)]

println(cond(Mbb),cond(D))
#αβ = D\d0
# αβ = gmres(D, d0)
αβ = minres(D, d0)


α = αβ[1:Nb]
β=αβ[end-ND+1:end]


for (j,p) in enumerate(getindex.(modelgrid0,direction))
    ss[j,direction]=sum(DistM[j,:,direction].*α) + β[1]+ β[2] *p
end
end

for j= 1:1:Nc
    Shift[j] = VectorValue(ss[j,:]...)
end


model.grid.node_coordinates .= model.grid.node_coordinates .+ Shift
writevtk(model, "model_rbf")
model.grid.node_coordinates .=modelgrid0


