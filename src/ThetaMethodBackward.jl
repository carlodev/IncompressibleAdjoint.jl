using Gridap.FESpaces: assemble_matrix_add!
using Gridap.ODEs
using Gridap.ODEs.TransientFETools
using Gridap.ODEs.ODETools
using Gridap.Algebra
using Gridap.FESpaces
using Gridap.Helpers

using Gridap.ODEs.TransientFETools: ODEOpFromFEOp
using Gridap.ODEs.TransientFETools: TransientFEOperatorFromWeakForm
# using Gridap.ODEs.TransientFETools: TransientFEOperatorFromWeakForm
using Gridap.ODEs.TransientFETools: jacobians!
using Gridap.ODEs.TransientFETools: fill_jacobians
using Gridap.ODEs.TransientFETools: _matdata_jacobian
using Gridap.ODEs.TransientFETools: _vcat_matdata
using Gridap.ODEs.TransientFETools: fill_initial_jacobians
using Gridap.FESpaces: assemble_matrix_add!
using Gridap.FESpaces: collect_cell_matrix
using Gridap.FESpaces: EvaluationFunction
using Gridap.FESpaces: get_test
using Gridap.FESpaces: get_trial
using Gridap.FESpaces: collect_cell_vector
using Gridap.FESpaces: assemble_vector!
using LinearAlgebra: fillstored!



struct ThetaMethodBackw <: ODESolver
    nls::NonlinearSolver
    dt::Float64
    θ::Real
end

"""
    Gridap.ODEs.ODETools.solve_step!(uf::AbstractVector,
    solver::ThetaMethodBackw,
    op::AffineODEOperator,
    u0::AbstractVector,
    t0::Real,
    cache)

Support the new type `ThetaMethodBackw` which allows integration back in time, useful for unsteady adjoint problems
"""
function Gridap.ODEs.ODETools.solve_step!(uf::AbstractVector,
    solver::ThetaMethodBackw,
    op::AffineODEOperator,
    u0::AbstractVector,
    t0::Real,
    cache) # -> (uF,tF)

dt = solver.dt
α = 1-solver.θ
dtα = dt*α

tα = t0+dtα


if cache === nothing
ode_cache = allocate_cache(op)
vθ = similar(u0)
vθ .= 0.0
l_cache = nothing
A, b = Gridap.ODEs.ODETools._allocate_matrix_and_vector(op,t0,u0,ode_cache)
else
ode_cache, vθ, A, b, l_cache = cache
end

ode_cache = update_cache!(ode_cache,op,tα)

_matrix_and_vector_back!(A,b,op,tα,dtα,u0,ode_cache,vθ)

## Useful to check if simulation is diverging
# println(Statistics.maximum(abs.(b)))
# println(norm(diag(A)))

afop = Gridap.Algebra.AffineOperator(A,b)


newmatrix = true
l_cache = Gridap.Algebra.solve!(uf,solver.nls,afop,l_cache,newmatrix)

uf = u0 - uf/α


cache = (ode_cache, vθ, A, b, l_cache)

tf = t0+dt
return (uf,tf,cache)

end

function _matrix_and_vector_back!(A,b,odeop,tθ,dtθ,u0,ode_cache,vθ)
    _matrix_back!(A,odeop,tθ,dtθ,u0,ode_cache,vθ)
    Gridap.ODEs.ODETools._vector!(b,odeop,tθ,dtθ,u0,ode_cache,vθ)
end

function _matrix_back!(A,odeop,tθ,dtθ,u0,ode_cache,vθ)
    z = zero(eltype(A))
    fillstored!(A,z)
    jacobians_back!(A,odeop,tθ,(vθ,vθ),(-1.0,1/dtθ),ode_cache)
end

function jacobians_back!(
    J::AbstractMatrix,
    op::ODEOpFromFEOp,
    t::Real,
    xhF::Tuple{Vararg{AbstractVector}},
    γ::Tuple{Vararg{Real}},
    ode_cache)
    Xh, = ode_cache
    dxh = ()
    for i in 2:Gridap.ODEs.TransientFETools.get_order(op)+1
      dxh = (dxh...,EvaluationFunction(Xh[i],xhF[i]))
    end
    xh=TransientCellField(EvaluationFunction(Xh[1],xhF[1]),dxh)
    jacobians_back!(J,op.feop,t,xh,γ,ode_cache)
 end

  function jacobians_back!(
    A::AbstractMatrix,
    op::TransientFEOperatorFromWeakForm,
    t::Real,
    xh::TransientCellField,
    γ::Tuple{Vararg{Real}},
    cache)
    _matdata_jacobians = fill_jacobians_back(op,t,xh,γ)
    matdata = Gridap.ODEs.TransientFETools._vcat_matdata(_matdata_jacobians)
    assemble_matrix_add!(A,op.assem_t, matdata)
    A
  end
  
  
  function fill_jacobians_back(
    op::TransientFEOperatorFromWeakForm,
    t::Real,
    xh::T,
    γ::Tuple{Vararg{Real}}) where T
    _matdata = ()
    for i in 1:Gridap.ODEs.TransientFETools.get_order(op)+1
      # if (γ[i] > 0.0)
        _matdata = (_matdata...,Gridap.ODEs.TransientFETools._matdata_jacobian(op,t,xh,i,γ[i]))
      # end
    end
    return _matdata
  end