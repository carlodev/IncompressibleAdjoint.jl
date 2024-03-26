module IncompressibleSolvers

using IncompressibleAdjoint
using IncompressibleAdjoint.Equations
using Gridap
using GridapDistributed
using Parameters

using Gridap.Algebra
using Gridap.FESpaces

export solve_inc_primal_u
export solve_inc_primal_s
include("SolvePrimal.jl")

export solve_inc_adj_u
export solve_inc_adj_s
include("SolveAdjoint.jl")

export solve_inc_direct_differentiation_u
export solve_inc_direct_differentiation_s
include("SolveDirectDifferentiation.jl")

end