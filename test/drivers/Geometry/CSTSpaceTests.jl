module CSTSpace
using Test
using IncompressibleAdjoint
using IncompressibleAdjoint.Geometry


xx = collect(0.0:0.001:1.0)


### Test Geometry Creation
CIRCLE = circle(xx)
@test typeof(CIRCLE) == Vector{Float64}

NACA0015 = NACA00(xx;t=0.15)
@test typeof(NACA0015) == Vector{Float64}

NACA0012 = CST_NACA0012(;N=20, t=0.0)
@test typeof(NACA0012) == CST


airfoil_cst = AirfoilCST(NACA0012, 0.0)
@test typeof(airfoil_cst) == AirfoilCST


des_points = AirfoilCSTDesign(airfoil_cst,xx)
@test typeof(des_points) == AirfoilCSTDesign

ap=AirfoilPoints(xx)
ap_new = IncompressibleAdjoint.Geometry.airfoil_from_CST(airfoil_cst,ap)

@test typeof(ap_new) == AirfoilPoints


ww = 0.1 .*ones(6)
dz = 0.0

ccst = CST(ww,-ww)
acst = AirfoilCST(ccst,dz)

apNACA0015 = AirfoilPoints(xx,xx,NACA0015, -NACA0015)

naca0015cst = IncompressibleAdjoint.Geometry.cst_from_points(acst,ap; maxiters = 20.0, maxtime=20.0)

@test typeof(naca0015cst) == CST


end