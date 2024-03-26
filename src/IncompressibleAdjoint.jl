module IncompressibleAdjoint

using Gridap
using GridapGmsh
using Statistics
using LinearAlgebra
using Gmsh
using SparseArrays
using Parameters
using FillArrays

using Gridap.Algebra
import Gridap.Algebra:AffineOperator

using GridapPETSc
### Create Figure Pkgs
# using Plots, XLSX, DataFrames, Images, ImageView
using OffsetArrays 

#CST support
using Optim
using Optimization, OptimizationBBO




export updatekey
include("InitializeParams.jl")

export ThetaMethodBackw
include("ThetaMethodBackward.jl")


include(joinpath("Equations","Equations.jl"))
include(joinpath("Geometry","Geometry.jl"))


export create_ũ_vector
export update_ũ
export update_ũ_vector!
include("LinearUtilities.jl")


include(joinpath("IncompressibleSolvers","IncompressibleSolvers.jl"))

export iterate_optimization
export compute_airfoil_forces
export compute_airfoil_coefficients
export rotation
export compute_drag
export compute_lift
include("AdjointIterators.jl")

export average_field!
export istantaneus_CD_CL
export average_CD_CL
include("TimeAverage.jl")

include(joinpath("Utils","Utils.jl"))
end
