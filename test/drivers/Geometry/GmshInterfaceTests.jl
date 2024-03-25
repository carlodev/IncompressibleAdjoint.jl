module GmshInterfaceTests
using IncompressibleAdjoint
using IncompressibleAdjoint.Geometry
using Test
using Gmsh


airfoil_cst = AirfoilCST(CST_NACA0012(;N=15,t= 0.0), 0.0)
xx = collect(0:0.0001:1)
des_points = AirfoilCSTDesign(airfoil_cst,xx)

mesh_folder = joinpath("test", "TestMesh")

modelname =create_msh(des_points; AoA=5.0, mesh_ref=4.0, folder=mesh_folder)


@test typeof(modelname) == String

rm(mesh_folder; force=true, recursive=true)


end