using IncompressibleAdjoint.Equations


function create_primal_spaces(model, params::Dict{Symbol,Any})
    @unpack tagname, order, D = params
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{D,Float64}, order )
    V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=["inlet", "limits", tagname])
    reffeₚ = ReferenceFE(lagrangian, Float64, order)
    Q = TestFESpace(model, reffeₚ, conformity=:H1, dirichlet_tags=["outlet"])

    return V,Q
end

function create_adjoint_spaces(model, params::Dict{Symbol,Any})
    @unpack tagname, order, D = params
    reffe_u_adj = ReferenceFE(lagrangian, VectorValue{D,Float64}, order)
    V_adj = TestFESpace(model, reffe_u_adj, conformity=:H1, dirichlet_tags=[tagname, "outlet"])
    reffe_p_adj = ReferenceFE(lagrangian, Float64, order)
    Q_adj = TestFESpace(model, reffe_p_adj, conformity=:H1, dirichlet_tags="inlet")
    
    return V_adj,Q_adj
end



function solve_inc_primal_u(model, params::Dict{Symbol,Any}; filename="inc-results")
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


    uh0 = interpolate_everywhere(u0(0.0), U(0.0))
    ph0 = interpolate_everywhere(p0(0.0), P(0.0))
    xh0 = interpolate_everywhere([uh0, ph0], X(0.0))

    updatekey(params, :uh,uh0)
    m, res, rhs = eq_primal_unsteady(params)

    op = TransientAffineFEOperator(m, res, rhs, X, Y)

    ls = LUSolver()

    ode_solver = ThetaMethod(ls,dt,θ)


    # θ_params = create_theta_params(Δt, θ, 2, U(0.0), P(0.0), X(0.0))

    # ode_solver = ThetaMethodMix(ls, Δt, θ_params)

    ode_solver = ThetaMethod(ls, dt, θ)

    sol = Gridap.solve(ode_solver, op, xh0, t0, tf)

    UH = [copy(uh0.free_values)]
    PH = [copy(ph0.free_values)]
    uh = uh0
    ph = ph0

    uvector = create_ũ_vector(uh.free_values)
    println("uvec created")

    res_path = "Results_primal"
    mkpath(res_path)

    createpvd(filename) do pvd
        for (xhtn, t) in sol
            uh = xhtn[1]
            ph = xhtn[2]
            pvd[t] = createvtk(Ω, joinpath(res_path, "$(filename)_$t" * ".vtu"), cellfields=["uh" => uh, "ph" => ph])
            push!(UH, copy(uh.free_values))
            push!(PH, copy(ph.free_values))

            update_ũ_vector!(uvector,uh.free_values)
            println("udate uvec")

            copyto!(params[:uh].free_values,update_ũ(uvector))
        end
    end

    return (uh, UH), (ph, PH)
end



function solve_inc_primal_s(model, params::Dict{Symbol,Any}; filename="inc-steady")
    @unpack D,order,t_endramp,t0,tf,θ,dt,u_in = params

    V,Q = create_primal_spaces(model,params)

    u0 = D == 2 ? VectorValue(u_in, 0.0) :  VectorValue(u_in, 0.0, 0.0)
    u_walls = D == 2 ? VectorValue(0.0, 0.0) :  VectorValue(0.0, 0.0, 0.0)
    p0 = 0.0

    U = TrialFESpace(V, [u0, u0, u_walls])
    P = TrialFESpace(Q, p0)

    Y = MultiFieldFESpace([V, Q])
    X = MultiFieldFESpace([U, P])

    degree = order * 2 
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)
    updatekey(params,:Ω,Ω)
    updatekey(params,:dΩ,dΩ)


    uh = interpolate_everywhere(u0, U)
    ph = interpolate_everywhere(p0, P)
    xh = interpolate_everywhere([uh, ph], X)

    updatekey(params, :uh,uh)
    res, rhs = eq_primal_steady(params)

    op = AffineFEOperator(res, rhs, X, Y)

  
    ls = LUSolver()
    solver = LinearFESolver(ls)

    Gridap.solve!(xh, solver, op)
    uh, ph = xh 

    #Correct bug Gridap on airfoil boundary
    for (i,val) in enumerate(Gridap.FESpaces.get_dirichlet_dof_tag(uh.fe_space) .== 3)
        if val
            uh.dirichlet_values[i] = 0.0
        end
    end


    if !isnothing(filename)
        writevtk(Ω, filename, cellfields=["uh" => uh, "ph" => ph])
    end
    

    return uh,ph
end


function solve_inc_adj_u(model, primal_sol_uh::Tuple, primal_sol_ph::Tuple, adjstart::Tuple, params::Dict{Symbol,Any}; filename="res-adj-unsteady")
    @unpack D,order,t_endramp,t0,tf,θ,dt,Ω,d_boundary = params

    uh0, UH = primal_sol_uh
    ph0, PH = primal_sol_ph

    ϕuh0, ϕph0 = adjstart
    
    V_adj,Q_adj = create_adjoint_spaces(model, params)

    #To be shifted outside
    # d_lift = VectorValue(0, 0.0)
    # d_drag = VectorValue(1.0, 0)

    # d_boundary = d_drag + d_lift


    uin(t) = (t < t_endramp) ? (1.0 - 1 .*(t_endramp-t)/t_endramp) : 1.0

    d0(x,t) = d_boundary #D == 2 ? VectorValue(-1 .*uin(t), 0.0) :  VectorValue(-1 .*uin(t), 0.0, 0.0)
    d0(t::Real) = x -> d0(x,t)
    u_walls(x,t) = VectorValue(zeros(D)...)
    u_walls(t::Real) = x -> u_walls(x,t)

    U_adj = TransientTrialFESpace(V_adj, [d0, u_walls])
    P_adj = TrialFESpace(Q_adj, 0.0)

    Y_adj = MultiFieldFESpace([V_adj, Q_adj])
    X_adj = TransientMultiFieldFESpace([U_adj, P_adj])

    ϕuh0 = interpolate_everywhere(VectorValue(0, 0), U_adj(0.0))
    ϕph0 = interpolate_everywhere(0.0, P_adj)

    ϕxh0 = interpolate_everywhere([ϕuh0, ϕph0], X_adj(0.0))
    
    Nfields = length(UH)

    ΦUH = [copy(UH[Nfields])]
    ΦPH = [copy(PH[Nfields])]

    copyto!(params[:uh].free_values, UH[Nfields])

    uh = uh0
    ph = ph0

    updatekey(params, :uh,uh)
    updatekey(params, :ph,ph)


    m_adj, res_adj, rhs = eq_adjoint_unsteady(params)
    op = TransientAffineFEOperator(m_adj, res_adj, rhs, X_adj, Y_adj)



    ls = LUSolver()

    θ_adj = 0.0
    ode_solver = ThetaMethodBackw(ls, dt, θ_adj)
  

    sol_adj = Gridap.solve(ode_solver, op, ϕxh0, t0, tf)

    res_path = "Results_adj"
    mkpath(res_path)

    #Adjoint going backwards
    createpvd(filename) do pvd
        for (idx, (ϕxh, t)) in enumerate(sol_adj)
            ϕuh = ϕxh[1]
            ϕph = ϕxh[2]


            pvd[t] = createvtk(Ω, joinpath(res_path, "$(filename)_$t" * ".vtu"), cellfields=["phi-uh" => ϕuh, "phi-ph" => ϕph,
                "uh" => uh, "ph" => ph])
            push!(ΦUH, copy(ϕuh.free_values))
            push!(ΦPH, copy(ϕph.free_values))

            IDX = Nfields - idx + 1

            println(IDX)
            copyto!(params[:uh].free_values, UH[IDX])
            copyto!(params[:ph].free_values, PH[IDX])

        end
    end

    return ΦUH, ΦPH
end



function solve_inc_adj_s(model, uh, ph, params::Dict{Symbol,Any}; filename="res-adj-steady")
    
    @unpack Ω, d_boundary = params
    copyto!(params[:uh].free_values, uh.free_values)

  
    # #To be shifted outside
    # d_lift = VectorValue(0, 0.0)
    # d_drag = VectorValue(1.0, 0)

    # d_boundary = d_drag + d_lift


    V_adj,Q_adj = create_adjoint_spaces(model, params)

    U_adj = TrialFESpace(V_adj, [d_boundary, VectorValue(0, 0)])
    P_adj = TrialFESpace(Q_adj, 0.0)

    Y_adj = MultiFieldFESpace([V_adj, Q_adj])
    X_adj = MultiFieldFESpace([U_adj, P_adj])

    res_adj, rhs = eq_adjoint_steady(params)


    op_adj = AffineFEOperator(res_adj, rhs, X_adj, Y_adj)

    ls = LUSolver()
    solver = LinearFESolver(ls)

    ϕu, ϕp = Gridap.solve(solver, op_adj)
    
    res_path = "Results_adj"
    mkpath(res_path)

    if !isnothing(filename)
        writevtk(Ω, joinpath(res_path, "$(filename)" * ".vtu"), cellfields=["phi-u" => ϕu, "phi-p" => ϕp, "uh0" => params[:uh]])
    end

    return ϕu, ϕp
end
