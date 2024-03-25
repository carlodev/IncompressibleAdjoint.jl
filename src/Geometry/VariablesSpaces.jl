
abstract type DesignParameters end

"""
It provides the arifoil points to Gmsh to create the mesh. In tag are stored the "top" or "bottom" string for each point to identify if it is on the top or bottom surface
"""
struct SplinePoints
    xcoordinate::Vector{Float64}
    ycoordinate::Vector{Float64}
    tag::Vector{String}
end

"""
It stores the coordinate of the spline control points.
"""
struct ControlPoints<:DesignParameters
    xcoordinate::Vector{Float64}
    ycoordinate::Vector{Float64}
    tag::Vector{String}
end

Base.copy(a::ControlPoints) = ControlPoints(a.xcoordinate,a.xcoordinate,a.tag)

function get_designparameters_number(cp::ControlPoints)
    return length(cp.tag)
end

function get_designparameters_tags(cp::ControlPoints)
    return cp.tag
end

"""
It stores the update for each degree of freedom of the the airfoil shape
"""
struct ShiftUpdate
    δ::Vector{Float64}
    tag::Vector{String}
end



function perturb_DesignParameters(a::ControlPoints, i::Int64, δ::Float64)
yc = copy(a.ycoordinate)
yc[i] = yc[i] + δ

return ControlPoints(a.xcoordinate,yc,a.tag)
end


function perturb_DesignParameters(a::ControlPoints, i::Int64, δ::Vector{Float64})
    xc = copy(a.xcoordinate)
    yc = copy(a.ycoordinate)
    xc[i] = xc[i] + δ[1]
    yc[i] = yc[i] + δ[2]
    
    return ControlPoints(xc,yc,a.tag)
end
    

function perturb_DesignParameters(a::ControlPoints, u::ShiftUpdate; tol=0.0025)
    @assert a.tag == u.tag
    δt = get_tag_shift(u,"top")
    δb = get_tag_shift(u,"bottom")

    xt,yt=get_tag_coordinates(a,"top")
    xb,yb=get_tag_coordinates(a,"bottom")
    yt0 = copy(yt)
    yb0 = copy(yb)



    yt_tmp = yt0 + δt
    yb_tmp = yb0 + δb

   
    d = yt_tmp .- yb_tmp
    d0 = yt0 .- yb0
    for i in findall( <(tol), d)
      k = (tol-d0[i])/( δt[i]- δb[i])
      println("d = $(d[i]), dt =$(δt[i]) db =$(δt[i])")
      println("k=$k")
      yt_tmp[i] = yt0[i] + k*δt[i]
      yb_tmp[i]= yb0[i] + k*δb[i]
    end
    yc_new = [yt_tmp...,yb_tmp...]

    return  ControlPoints(a.xcoordinate,yc_new,u.tag)
end


function initialize_control_points(xcontrol::Vector{Float64}; chord = 1.0, initfun=circle)

    ycontrol_top = initfun(xcontrol) #circle.(xcontrol,1,chord)
    ycontrol_bottom = -1 .* initfun.(xcontrol) #circle.(xcontrol,-1,chord)

    toptags = repeat(["top"],length(ycontrol_top))
    bottomtags = repeat(["bottom"],length(ycontrol_bottom))
    tags = vcat(toptags,bottomtags)
    xcoordinate = vcat(xcontrol,xcontrol)
    ycoordinate = vcat(ycontrol_top,ycontrol_bottom)
    return ControlPoints(xcoordinate,ycoordinate,tags)
end

function get_tag_coordinates(cp::ControlPoints,tagname::String)
    idx = findall(a->a==tagname,cp.tag)
    return cp.xcoordinate[idx],cp.ycoordinate[idx]
end

function get_tag_coordinates(cp::SplinePoints,tagname::String)
    idx = findall(a->a==tagname,cp.tag)
    return cp.xcoordinate[idx],cp.ycoordinate[idx]
end

function get_tag_shift(u::ShiftUpdate,tagname::String)
    idx = findall(a->a==tagname,u.tag)
    return u.δ[idx]
end





function add_control_point_(a::ControlPoints,xnew::Float64,tag::String;chord=1.0)
    xc,yc = get_tag_coordinates(a,tag)
    xcc = copy(xc)
    ycc = copy(yc)
    itp = Interpolations.linear_interpolation([0.0,xc...,chord],[0.0,yc...,0.0]) #BSplineKit.interpolate(xc, yc, BSplineOrder(1))
    ynew = copy(itp(xnew))
    
    pos = findfirst(a-> a> 0, xc .- xnew)
    if isnothing(pos)
        pos = length(xcc)+1
    end
    insert!(xcc,pos,xnew)
    insert!(ycc,pos,ynew)
    tag = repeat([tag],length(xcc))
    return xcc,ycc,tag
end


function add_control_point(a::ControlPoints,xnew::Float64)
    xcc_t,ycc_t,tag_t= add_control_point_(a,xnew,"top")
    xcc_b,ycc_b,tag_b= add_control_point_(a,xnew,"bottom")
    anew= ControlPoints([xcc_t...,xcc_b...],[ycc_t...,ycc_b...],[tag_t...,tag_b...])
    return anew
end

####################################################################################################
#CST interface
####################################################################################################

struct CSTWeights
    w::Vector{Float64}
    tag::Vector{String}
end

struct AirfoilCSTDesign<:DesignParameters
    acst::AirfoilCST
    apoints::AirfoilPoints
    acstw::CSTWeights
end

function AirfoilCSTDesign(acst::AirfoilCST,v::Vector)
    apoints = airfoil_from_CST(acst,v)
    acstw = CSTWeights(acst.cst)
    AirfoilCSTDesign(acst,apoints,acstw)
end

function get_designparameters_number(ad::AirfoilCSTDesign)
    return length(ad.acstw.tag)
end

function get_designparameters_tags(ad::AirfoilCSTDesign)
    return ad.acstw.tag
end

function get_CST_values(ad::AirfoilCSTDesign)
    weights = ad.acstw.w
    tags = ad.acstw.tag
    return weights,tags
end

function update_CST_weights(w::Vector, ad::AirfoilCSTDesign)
    @assert length(ad.acstw.w)==length(w)
    _,tags = get_CST_values(ad)
    newCSTW = CSTWeights(w,tags)
    new_ad = update_AirfoilCSTDesign(newCSTW,ad)
    return new_ad
end

function CSTWeights(c::CST)
    wu = c.wu
    wl = c.wl
    nu = length(wu)
    nl=length(wl)
    toptags = repeat(["top"],nu)
    bottomtags = repeat(["bottom"],nl)
    w = [wu;wl]
    tag=[toptags;bottomtags]
    return CSTWeights(w,tag)
end


function CST(c::CSTWeights)
    w = c.w
    tag=c.tag
    topidx = findall(tag .== "top" )
    bottomidx = findall(tag .== "bottom" )
    wu = w[topidx]
    wl = w[bottomidx]
    return CST(wu,wl)
end

function perturb_DesignParameters(ad::AirfoilCSTDesign, i::Int64, δ::Float64)
    n_acstw = perturb_DesignParameters(ad.acstw,i,δ)
    return update_AirfoilCSTDesign(n_acstw,ad)
end

function perturb_DesignParameters(ad::AirfoilCSTDesign, u::ShiftUpdate)
    n_acstw = perturb_DesignParameters(ad.acstw,u)
    return  update_AirfoilCSTDesign(n_acstw,ad)
end

function perturb_DesignParameters(c::CSTWeights, i::Int64, δ::Float64)
    w = copy(c.w)
    w[i] = w[i] + δ
    
    return CSTWeights(w,c.tag)
end
    

function perturb_DesignParameters(c::CSTWeights, u::ShiftUpdate)
    @assert c.tag == u.tag
    w = copy(c.w)
    new_w = w .+ u.δ

    return  CSTWeights(new_w,c.tag)
end


function update_AirfoilCSTDesign(newCSTW::CSTWeights, ad::AirfoilCSTDesign)
acst = ad.acst
apoints = ad.apoints
n_cst = CST(newCSTW)
n_airfoilCST = AirfoilCST(n_cst,acst.dz,acst.N1,acst.N2)
n_apoints = airfoil_from_CST(n_airfoilCST,apoints)
return AirfoilCSTDesign(n_airfoilCST,n_apoints,newCSTW)
end

################################################################
#Get SPline Points
################################################################
function SplinePoints(c::ControlPoints)
    SplinePoints(c.xcoordinate,c.ycoordinate,c.tag)
end

function SplinePoints(ap::AirfoilPoints)
    xu = ap.xu
    xl = ap.xl
    yu = ap.yu
    yl = ap.yl
    nu = length(xu)
    nl = length(xl)
    toptags = repeat(["top"],nu)
    bottomtags = repeat(["bottom"],nl)
    xcoordinate = [xu;xl]
    ycoordinate = [yu;yl]
    tag=[toptags;bottomtags]

    return SplinePoints(xcoordinate,ycoordinate,tag)
end

function SplinePoints(a::AirfoilCSTDesign)
    new_apoints = airfoil_from_CST(a.acst,a.apoints)
    SplinePoints(new_apoints)
end