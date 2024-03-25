using NLopt, Plots

###Store the values at each optimization cicle
fvals_vec = Float64[]

#Store the value of each variable at each cicle
fvar_vec = typeof(rand(3))[]


#Function (also black box). It receives the unknown x
# It computes the value (f(x)) and gradient (vector df/dx_i). Here it is computed analytically

function myfunc(x::Vector, grad::Vector)
    if length(grad) > 0
        grad[1] = 4*(x[1]-3)^3
        grad[2] = 4*(x[2]-5)^3
        grad[3] = 2*(x[3]+4)
    end
    val = (x[1]-3)^4+(x[2]-5)^4  +(x[3]+4)^2+ 2
    println("f(x) = $val, x =$x")
    push!(fvals_vec,val)
    push!(fvar_vec,copy(x))

    return val
end

## Set the optimization algorithm
opt = Opt(:LD_LBFGS, 3)
## Set the bounds
opt.lower_bounds = [-10, -10,-10 ]
opt.upper_bounds = [10, 10,10 ]

## Set the tolerance
opt.xtol_rel = 1e-3

opt.min_objective = myfunc


fvals_vec = Float64[]
### return the min value, the solution vector, and stop criteria reached
(minf,minx,ret) = optimize(opt, [0.0,0.0,0.0])
numevals = opt.numevals # the number of function evaluations
println("got $minf at $minx after $numevals iterations (returned $ret)")

plot(fvals_vec[3:end])
plot(getindex.(fvar_vec[3:end],1))
plot!(getindex.(fvar_vec[3:end],2))
plot!(getindex.(fvar_vec[3:end],3))

# methodswith(typeof(opt))
# propertynames(opt)
