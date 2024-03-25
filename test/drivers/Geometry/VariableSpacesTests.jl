module VariableSpacesTests
using IncompressibleAdjoint
using IncompressibleAdjoint.Geometry
using Test

xx = collect(0.0:0.001:1.0)

contr_point = initialize_control_points(xx;initfun=NACA00)
@test typeof(contr_point) == ControlPoints

NACA0012 = CST_NACA0012(;N=20, t=0.0)

airfoil_cst = AirfoilCST(NACA0012, 0.0)
des_points = AirfoilCSTDesign(airfoil_cst,xx)

@test typeof(des_points) == AirfoilCSTDesign
@test typeof(des_points) <: DesignParameters

@test typeof(SplinePoints(contr_point)) == SplinePoints
@test typeof(SplinePoints(des_points))==SplinePoints

### CST Interfaces
new_w = ones(des_points.acstw.w)

get_CST_values(des_points)
update_des_points = update_CST_weights(new_w,des_points)


@test update_des_points.acstw.w == new_w

w,tag=get_CST_values(update_des_points)

@test typeof(w) <: Vector
@test typeof(tag) <: Vector{String}

end