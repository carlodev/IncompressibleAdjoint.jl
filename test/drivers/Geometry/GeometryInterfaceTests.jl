module GeometryInterfaceTests
using Test
using IncompressibleAdjoint

@testset "GmshInterfaceTests" begin include("GmshInterfaceTests.jl") end
@testset "VariableSpacesTests" begin include("VariableSpacesTests.jl") end
@testset "CSTSpaceTests" begin include("CSTSpaceTests.jl") end


end