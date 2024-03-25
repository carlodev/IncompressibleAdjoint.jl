module InitializeParametersTests
using Test
using IncompressibleAdjoint

params = Dict(:a=>rand(5), :b=>"efg")
IncompressibleAdjoint.verifykey(params,:c)

@test  params[:c]==false

updatekey(params, :c, 70)

@test  params[:c]==70

end