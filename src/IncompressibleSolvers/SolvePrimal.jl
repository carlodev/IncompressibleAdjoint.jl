function create_primal_spaces(model, params::Dict{Symbol,Any})
    @unpack tagname, order, D = params
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{D,Float64}, order )
    V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=["inlet", "limits", tagname],  dirichlet_masks=[(true,true), (false,true),(true,true) ])
    reffeₚ = ReferenceFE(lagrangian, Float64, order)
    Q = TestFESpace(model, reffeₚ, conformity=:H1, dirichlet_tags=["outlet"])

    return V,Q
end


function solve_inc_primal_u(model, params::Dict{Symbol,Any}; filename="inc-results", uh00=nothing,ph00=nothing)
    
    @unpack D,order,t_endramp,t0,tf,θ,dt = params


    V,Q = create_primal_spaces(model,params)

    uin(t) = (t < t_endramp) ? (1.0 - 1 .*(t_endramp-t)/t_endramp) : 1.0

    u0(x,t) = D == 2 ? VectorValue(uin(t), 0.0) :  VectorValue(uin(t), 0.0, 0.0)
    u0(t::Real) = x -> u0(x,t)


    u_walls(x,t) = D == 2 ? VectorValue(0.0, 0.0) :  VectorValue(0.0, 0.0, 0.0)
    u_walls(t::Real) = x -> u_walls(x,t)
    p0(x,t) = 0.0
    p0(t::Real) = x -> p0(x,t)

    U = TransientTrialFESpace(V, [u0, u0, u_walls])
    P = TransientTrialFESpace(Q, p0)

    Y = TransientMultiFieldFESpace([V, Q])
    X = TransientMultiFieldFESpace([U, P])

    degree = order * 2 
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)
    updatekey(params,:Ω,Ω)
    updatekey(params,:dΩ,dΩ)
    
    uh0 = interpolate(u0(0.0), U(0.0))
    
    ph0 = interpolate(p0(0.0), P(0.0))

    if !isnothing(uh00)
       uh0.free_values .=  uh00.free_values
    end
    
    if !isnothing(ph00)
        ph0.free_values .=  ph00.free_values
    end

    dudt = interpolate(u_walls(0.0), U(0.0))

    xh0 = interpolate([uh0, ph0], X(0.0))

    updatekey(params, :uh,uh0)
    m, res, rhs = eq_primal_unsteady(params)

    op = TransientAffineFEOperator(m, res, rhs, X, Y)

    ls = LUSolver()

    ode_solver = ThetaMethod(ls,dt,θ)

    sol = Gridap.solve(ode_solver, op, xh0, t0, tf)

    UH = [copy(uh0.free_values)]
    PH = [copy(ph0.free_values)]
    
    uh = uh0
    ph = ph0

    uvector = create_ũ_vector(uh0.free_values)

    res_path = "Results_primal"
    mkpath(res_path)

    createpvd(filename) do pvd
        for (xhtn, t) in sol
            uh = xhtn[1]
            ph = xhtn[2]
            pvd[t] = createvtk(Ω, joinpath(res_path, "$(filename)_$t" * ".vtu"), cellfields=["uh" => uh, "ph" => ph])
            push!(UH, copy(uh.free_values))
            push!(PH, copy(ph.free_values))
            println("Primal solved at time step $t")

            update_ũ_vector!(uvector,uh.free_values)
            u_new=update_ũ(uvector)
            copyto!(params[:uh].free_values,u_new)


        end
    end

    DUHDT = vcat(0.0,UH[2:end] -UH[1:end-1])./dt
    return (uh,dudt, UH, DUHDT), (ph, PH)

end



function solve_inc_primal_s(model, params::Dict{Symbol,Any}; filename="inc-steady")
    @unpack D,order,t_endramp,t0,tf,θ,dt,u_in = params

    V,Q = create_primal_spaces(model,params)

    u0 = D == 2 ? VectorValue(u_in, 0.0) :  VectorValue(u_in, 0.0, 0.0)
    u_walls = D == 2 ? VectorValue(0.0, 0.0) :  VectorValue(0.0, 0.0, 0.0)
    p0 = 0.0

    U = TrialFESpace(V, [u0, u0, u_walls])
    P = TrialFESpace(Q, [p0])

    Y = MultiFieldFESpace([V, Q])
    X = MultiFieldFESpace([U, P])

    degree = 8
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)
    updatekey(params,:Ω,Ω)
    updatekey(params,:dΩ,dΩ)


    uh = interpolate(u0, U)
    ph = interpolate(-5.0, P)
    xh = interpolate([uh, ph], X)

    updatekey(params, :uh,uh)
    res, rhs = eq_primal_steady(params)

    op = AffineFEOperator(res, rhs, X, Y)

  
    ls = LUSolver()
    solver = LinearFESolver(ls)

    Gridap.solve!(xh, solver, op)
    uh, ph = xh 


    if !isnothing(filename)
        writevtk(Ω, filename, cellfields=["uh" => uh, "ph" => ph])
    end
    

    return uh,ph
end
