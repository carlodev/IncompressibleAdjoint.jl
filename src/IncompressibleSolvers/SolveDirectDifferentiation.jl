function create_direct_differentiation_spaces(model, params::Dict{Symbol,Any})
    @unpack tagname, order, D = params
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{D,Float64}, order )
    V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=["inlet", tagname])
    reffeₚ = ReferenceFE(lagrangian, Float64, order)
    Q = TestFESpace(model, reffeₚ, conformity=:H1, dirichlet_tags=["outlet"])
    return V,Q
end



function inc_direct_differentiation(model,uh::SingleFieldFEFunction, params::Dict{Symbol,Any}, point::VectorValue, nb::VectorValue)
    @unpack D,order,tagname,δ = params

    V,Q = create_direct_differentiation_spaces(model,params)

    function node_filter(x)
        R = 0.001 #get_radius_shift(nb.*δ)
        dist = compute_point_dist(x, point)
        morph_kernel(dist,R)
    end

    u_walls = VectorValue(zeros(D)...)

    P = TrialFESpace(Q, 0.0)
    nf = interpolate_everywhere(node_filter,P)


    utagname = -transpose(∇(uh))⋅(nb*nf)

    U = TrialFESpace(V, [u_walls,utagname])

    Y = MultiFieldFESpace([V, Q])
    X = MultiFieldFESpace([U, P])

    degree = 8
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)
    updatekey(params,:Ω,Ω)
    updatekey(params,:dΩ,dΩ)
    

    writevtk(Ω, "DirectDifferentiation/BC$(point[1])", cellfields=["BCU"=>utagname, "node_filter"=>nb*nf])

    uhb0 = interpolate(u_walls, U)
    phb0 = interpolate(0.0, P)
    xhb0 = interpolate([uhb0, phb0], X)

    return params, X,Y,xhb0
end

function solve_inc_direct_differentiation_s(model,uh::SingleFieldFEFunction, params::Dict{Symbol,Any}, point::VectorValue, nb::VectorValue; filename="inc-direct-diff")
    
    params, X,Y,xhb0 = inc_direct_differentiation(model,uh, params, point, nb)
    
    res, rhs = eq_direct_differentiation_steady(params)

    op = AffineFEOperator(res, rhs, X, Y)

    ls = LUSolver()

    solver = LinearFESolver(ls)

    Gridap.solve!(xhb0, solver, op)
    uhb, phb = xhb0


    if !isnothing(filename)
        writevtk(Ω, filename, cellfields=["uhb" => uhb, "phb" => phb])
    end
    

    return uhb,phb
end

function solve_inc_direct_differentiation_u(model,primal_sol_uh::Tuple, primal_sol_ph::Tuple, params::Dict{Symbol,Any}, point::VectorValue, nb::VectorValue; filename="inc-direct-diff")
    
    @unpack t0,tf,dt, θ=params
     
    uh0, UH = primal_sol_uh
    ph0, PH = primal_sol_ph
    
    params, X,Y,xhb0 = inc_direct_differentiation(model,uh, params, point, nb)

    ode_solver = ThetaMethod(ls, dt, θ)
    
    m, res, rhs = eq_direct_differentiation_unsteady(params)

    op = TransientAffineFEOperator(m, res, rhs, X, Y)

    ls = LUSolver()

    ode_solver = ThetaMethod(ls,dt,θ)

    sol = Gridap.solve(ode_solver, op, xhb0, t0, tf)
   
    uhb0,phb0=xhb0

    UbH = [copy(uhb0.free_values)]
    PbH = [copy(phb0.free_values)]

    updatekey(params, :uh,uh0)
    updatekey(params, :ph,ph0)

    res_path = "DirectDifferentiation"
    mkpath(res_path)

    createpvd(filename) do pvd
        for (idx,(xhtn, t)) in enumerate(sol)
            ubh = xhtn[1]
            pbh = xhtn[2]
            pvd[t] = createvtk(Ω, joinpath(res_path, "$(filename)_$t" * ".vtu"), cellfields=["ubh" => ubh, "pbh" => pbh])
            copyto!(params[:uh].free_values, UH[idx])
            copyto!(params[:ph].free_values, PH[idx])

            println("DirectDifferentiation solved at time step $t")

            push!(UbH, copy(ubh.free_values))
            push!(PbH, copy(pbh.free_values))


        end
    end

    return (uh,UH), (ph, PH)
end