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


export create_ũ_vector
export update_ũ
export update_ũ_vector!
include("LinearUtilities.jl")

export CST
export AirfoilCST
export AirfoilPoints
include("CSTSpace.jl")

export AirfoilCSTDesign
export DesignParameters
export SplinePoints
export ControlPoints
include("VariablesSpaces.jl")

export create_msh
include("GmshInterface.jl")

export circle
export NACA00
export CST_NACA0012
include("GeometryShapes.jl")

export solve_inc_primal_u
export solve_inc_primal_s
export solve_inc_adj_u
export solve_inc_adj_s
include("Solve_primal_and_adjoint_mix.jl")


export get_control_boundary
export morph_kernel
export get_radius_shift
export compute_point_dist
export get_shift_vec
include("Morph.jl")

export solve_inc_direct_differentiation_s
include("Solve_direct_differentiation.jl")


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

end
