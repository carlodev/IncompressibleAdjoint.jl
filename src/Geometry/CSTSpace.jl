using Optim
using Optimization, OptimizationBBO
using LinearAlgebra

"""
CST

It stores the weight for the upper and lower side of the airfoil
"""
struct CST
    wu::Vector{Float64}
    wl::Vector{Float64}
end

struct AirfoilCST
    cst::CST
    dz::Float64
    N1::Real
    N2::Real

end

struct AirfoilPoints
xu::Vector{Float64}
xl::Vector{Float64}
yu::Vector{Float64}
yl::Vector{Float64}
end


function AirfoilCST(cstfun::Function)
    cst = cstfun
    return AirfoilCST(cst,0.0)
end

function AirfoilCST(cst::CST, dz::Float64)
    AirfoilCST(cst,dz,0.5,1.0)
end


function get_cst_weight(c::CST)
    c.wu,c.wl
end


function get_cst(a::AirfoilCST)
    return a.cst
end

function get_cst_weight(a::AirfoilCST)
    c = get_cst(a)
    get_cst_weight(c)
end

function AirfoilPoints(v::Vector)
return AirfoilPoints(v,v,v,v)
end

function AirfoilPoints(xu::Vector,xl::Vector)
    return AirfoilPoints(xu,xl,xu,xl)
end


function Base.ones(v::Vector)
    N = length(v)
    return ones(N)
end

function Base.zeros(v::Vector)
    N = length(v)
    return zeros(N)
end

function Base.similar(ap::AirfoilPoints)
nu=length(ap.xu)
nl=length(ap.xl)
AirfoilPoints(ap.xu,ap.xl,  Vector(undef, nu), Vector(undef, nl))
end

"""
    ClassFunction(x::Vector{Float64},N1::Real,N2::Real)

Compute the class function
    ``C = \\phi^N1 \\cdot(1-\\phi)^N2``
"""
function ClassFunction(x::Vector{Float64},N1::Real,N2::Real)
    @assert N1>0
    @assert N2>0 

    C = zeros(length(x))

    for (i,xi) in enumerate(x)
        C[i] = xi^N1*((1-xi)^N2)
    end

    return C    
end


"""
    ShapeFunction(w::Vector,x::Array{Float64})

Compute the shape function

"""
function ShapeFunction(w::Vector,x::Array{Float64})
 
    # Shape function; using Bernstein Polynomials
    n = length(w)-1 # Order of Bernstein polynomials
    
    K = zeros(n+1)
    
    for i = 1:n+1
         K[i] = factorial(n)/(factorial(i-1)*(factorial((n)-(i-1))))
    end
    
    S = zeros(length(x))
    
    for (i,xi) in enumerate(x)
        for j = 1:n+1
            S[i] = S[i] + w[j]*K[j]*xi^(j-1)*((1-xi)^(n-(j-1)))
        end
    end

    return S

end
    
function compute_airfoil_y(w::Vector,x::Array{Float64},N1::Real,N2::Real,dz::Real)

        #Compute Class Function
        C = ClassFunction(x,N1,N2)

        #Compute Shape Function
        S = ShapeFunction(w,x)
        #  Calculate y output
        y = zeros(length(x))
        for (i,xi) in enumerate(x)
           y[i] = C[i]*S[i] + xi*dz;
        end
                
        return y
end


function airfoil_from_CST(a::AirfoilCST,ap::AirfoilPoints)
wu,wl=get_cst_weight(a)
xu = ap.xu
xl = ap.xl
dz = a.dz
N1 = a.N1
N2 = a.N2
yu = compute_airfoil_y(wu,xu,N1,N2,dz)  # Call ClassShape function to determine upper surface y-coordinates
yl = compute_airfoil_y(wl,xl,N1,N2,-dz) # Call ClassShape function to determine lower surface y-coordinates
anew = AirfoilPoints(xu,xl,yu,yl)
return anew
end


function airfoil_from_CST(a::AirfoilCST,v::Vector)
ap = AirfoilPoints(v)
airfoil_from_CST(a,ap)
end


"""
    compute_error(y0,y)
Compute the error of approximation. y0 are the original points, y are the new ones.
"""
function compute_error(y0,y)
    n = length(y0)
    return sqrt(sum((y0-y).^2)/n)
end

"""
    split_w_wlwu(w::Vector{Float64},split_idx::Tuple{Int64})

Creates 2 vector, one for wl and one for wu
"""
function split_w_wlwu(w::Vector{Float64},split_idx::Tuple{Int64,Int64})
    wu = w[1:split_idx[1]]
    wl = w[split_idx[1]+1:split_idx[1]+split_idx[2]]
    return wu,wl
end

"""
    compute_cst_error(w,p)
Compute the error of approximation
"""
function compute_cst_error(w,params)
    split_idx,ap,dz,y0 = params
    wu,wl = split_w_wlwu(w,split_idx)
    println(wu)

    cst_tmp = CST(wu,wl)
    
    acst_tmp = AirfoilCST(cst_tmp,dz)
    ap_tmp = similar(ap)

    ap_tmp = airfoil_from_CST(acst_tmp,ap_tmp)
    ytmp = [ap_tmp.yu...,ap_tmp.yl...]
    err = compute_error(y0,ytmp)

    return err
end


function cst_from_points(acst::AirfoilCST,ap::AirfoilPoints; maxiters = 50.0, maxtime=50.0)
    yu = ap.yu
    yl = ap.yl
    dz = acst.dz
    wu,wl=get_cst_weight(acst)


    y0 = copy([yu;yl])
    split_idx = length(wu),length(wl)
    
    w0 = copy([wu...,wl...])
    ubound = [ones(wu)...,zeros(wl)...]
    lbound = [zeros(wu)...,-1 .*ones(wl)...]

    error_function = OptimizationFunction(compute_cst_error)
        
    params = (split_idx,ap,dz,y0)

    prob = Optimization.OptimizationProblem(error_function, w0, params,  ub = ubound, lb = lbound)
    
    sol = Optimization.solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters = maxiters,  maxtime = maxtime)
    sol = collect(sol)
    
    wl, wu = split_w_wlwu(sol, split_idx)
   
    return CST(wl,wu) 
end


