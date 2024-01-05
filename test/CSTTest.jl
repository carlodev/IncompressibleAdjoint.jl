using IncompressibleAdjoint
xx = collect(LinRange(0.0,1.0,1000))
yy = NACA00.(xx)

ww = 0.1 .*ones(3)
dz = 0.0

ccst = CST(ww,-ww)
acst = AirfoilCST(ccst,dz)

ap = AirfoilPoints(xx,xx,yy, -yy)

naca0012cst = IncompressibleAdjoint.cst_from_points(acst,ap; maxiters = 50.0, maxtime=50.0)