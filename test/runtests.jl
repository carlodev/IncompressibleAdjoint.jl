using IncompressibleAdjoint
using Test

@testset "GeometryTests" begin include(joinpath("drivers","Geometry", "GeometryTests.jl")) end

@testset "Cases.jl" begin 
    include(joinpath("drivers", "NACA0012_CST_Opt.jl")) 
end

