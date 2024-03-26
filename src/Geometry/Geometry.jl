module Geometry

using Gridap, GridapGmsh
using Gmsh
using LinearAlgebra
using Parameters

export CST
export AirfoilCST
export AirfoilPoints
include("CSTSpace.jl")

export AirfoilCSTDesign
export DesignParameters
export SplinePoints
export ControlPoints
export get_CST_values
export update_CST_weights
export initialize_control_points
export get_designparameters_number
export get_designparameters_tags
export perturb_DesignParameters
include("VariablesSpaces.jl")

export create_msh
include("GmshInterface.jl")

export circle
export NACA00
export CST_NACA0012
include("GeometryShapes.jl")

export get_control_boundary
export morph_kernel
export get_radius_shift
export compute_point_dist
export get_shift_vec
include("Morph.jl")

end