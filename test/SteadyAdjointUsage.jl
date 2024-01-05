using IncompressibleAdjoint
using Gridap,GridapGmsh
using Parameters

Reynolds = 10_000

params=Dict(
    :chord => 1.0,
    :D=>2,
    :Re => Reynolds,
    :dt => 0.05,
    :ν =>1 / Reynolds,
    :D =>2,
    :order => 1,
    :method=>:SUPG,
    :tagname=>"airfoil",
    :t0=>0.0,
    :tf=>20.0,
    :t_endramp=>0.0,
    :θ=>1.0,
)
