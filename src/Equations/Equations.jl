module Equations
using Gridap
using GridapDistributed
using Parameters

import IncompressibleAdjoint: updatekey

export eq_primal_steady
export eq_primal_unsteady
export eq_adjoint_steady
export eq_adjoint_unsteady

export eq_direct_differentiation_steady
export eq_direct_differentiation_unsteady

include("EquationsOperations.jl")
include("StabParams.jl")
include("PrimalEquations.jl")
include("AdjointEquations.jl")
include("DirectDifferentiation.jl")

end