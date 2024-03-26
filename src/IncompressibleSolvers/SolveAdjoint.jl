function create_adjoint_spaces(model, params::Dict{Symbol,Any})
    @unpack tagname, order, D = params
    reffe_u_adj = ReferenceFE(lagrangian, VectorValue{D,Float64}, order)
    V_adj = TestFESpace(model, reffe_u_adj, conformity=:H1, dirichlet_tags=[tagname, "outlet","limits"], dirichlet_masks=[(true,true), (true,true),(false,true) ])
    reffe_p_adj = ReferenceFE(lagrangian, Float64, order)
    Q_adj = TestFESpace(model, reffe_p_adj, conformity=:H1, dirichlet_tags="inlet")
    
    return V_adj,Q_adj
end


function solve_inc_adj_u(model, primal_sol_uh::Tuple, primal_sol_ph::Tuple, adjstart::Tuple, params::Dict{Symbol,Any}; filename="res-adj-unsteady")
    @unpack D,order,t_endramp,t0,tf,θ,dt,Ω,d_boundary = params

    uh0, UH = primal_sol_uh
    ph0, PH = primal_sol_ph

    ϕuh0, ϕph0 = adjstart
    
    V_adj,Q_adj = create_adjoint_spaces(model, params)

  

    uin(t) = (t < t_endramp) ? (1.0 - 1 .*(t_endramp-t)/t_endramp) : 1.0

    d0(x,t) = d_boundary #D == 2 ? VectorValue(-1 .*uin(t), 0.0) :  VectorValue(-1 .*uin(t), 0.0, 0.0)
    d0(t::Real) = x -> d0(x,t)
    u_walls(x,t) = VectorValue(zeros(D)...)
    u_walls(t::Real) = x -> u_walls(x,t)

    U_adj = TransientTrialFESpace(V_adj, [d0, u_walls,u_walls])
    P_adj = TrialFESpace(Q_adj, 0.0)

    Y_adj = MultiFieldFESpace([V_adj, Q_adj])
    X_adj = TransientMultiFieldFESpace([U_adj, P_adj])

    ϕuh0 = interpolate(VectorValue(0, 0), U_adj(0.0),)
    ϕph0 = interpolate(0.0, P_adj)

    ϕxh0 = interpolate([ϕuh0, ϕph0], X_adj(0.0))
    
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
            println("Adjoint solved at time step $t")
            copyto!(params[:uh].free_values, UH[IDX])
            copyto!(params[:ph].free_values, PH[IDX])

        end
    end

    return ΦUH, ΦPH
end



function solve_inc_adj_s(model, uh, ph, params::Dict{Symbol,Any}; filename="res-adj-steady")
    
    @unpack Ω, d_boundary = params
    copyto!(params[:uh].free_values, uh.free_values)


    Γout = BoundaryTriangulation(model; tags="outlet")
    dΓout = Measure(Γout, 4)
    nΓout = -1 .* get_normal_vector(Γout)

    Γlim = BoundaryTriangulation(model; tags="limits")
    dΓlim = Measure(Γlim, 4)
    nΓlim = -1 .* get_normal_vector(Γlim)

    Γairfoil = BoundaryTriangulation(model; tags="airfoil")
    dΓairfoil = Measure(Γairfoil, 4)
    nΓairfoil = -1 .* get_normal_vector(Γairfoil)


    updatekey(params, :dΓout,dΓout)
    updatekey(params, :nΓout,nΓout)

    updatekey(params, :dΓlim,dΓlim)
    updatekey(params, :nΓlim,nΓlim)
    
    updatekey(params, :dΓairfoil,dΓairfoil)
    updatekey(params, :nΓairfoil,nΓairfoil)

    V_adj,Q_adj = create_adjoint_spaces(model, params)
    println(d_boundary)
    U_adj = TrialFESpace(V_adj, [d_boundary, VectorValue(0, 0),VectorValue(0, 0)])
    P_adj = TrialFESpace(Q_adj,0.0)

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
