using Plots
using IncompressibleAdjoint
using IncompressibleAdjoint.Utils
using IncompressibleAdjoint.Geometry
using Gridap, GridapGmsh
using IterativeSolvers

modelname = "MeshFiles/Rectangle.msh"
model = GmshDiscreteModel(modelname)
writevtk(model, "Model")


#Store the node_coordinates of the original model
modelgrid0 = copy(model.grid.node_coordinates)

displacement =  [VectorValue(0.75,0.25) for _ in 1:56]
model.grid.node_coordinates .= MorphRBF(model, "rectangle", ["limits"], displacement)
writevtk(model, "model_rbf")
model.grid.node_coordinates .= modelgrid0
