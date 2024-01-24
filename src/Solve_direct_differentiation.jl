using IncompressibleAdjoint.Equations
include("Morph.jl")

function create_direct_differentiation_spaces(model, params::Dict{Symbol,Any})
    @unpack tagname, order, D = params
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{D,Float64}, order )
    V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=["inlet", tagname])
    reffeₚ = ReferenceFE(lagrangian, Float64, order)
    Q = TestFESpace(model, reffeₚ, conformity=:H1, dirichlet_tags=["outlet"])

    return V,Q
end




function solve_inc_direct_differentiation_s(model,uh::SingleFieldFEFunction, params::Dict{Symbol,Any}, point::VectorValue, nb::VectorValue; filename="inc-direct-diff")
    @unpack D,order,tagname = params

    V,Q = create_direct_differentiation_spaces(model,params)

    function morph_kernel(x)
        R = get_radius_shift(δ)
        println(point)
    
        dist = compute_point_dist(x, point)
        δ.*morph_kernel(dist,R)
    end

    u_walls = VectorValue(zeros(D)...)

    utagname = -∇(uh)⋅nb*morph_kernel(x)

    U = TrialFESpace(V, [u_walls,utagname])
    P = TrialFESpace(Q, 0.0)

    Y = MultiFieldFESpace([V, Q])
    X = MultiFieldFESpace([U, P])

    degree = order * 2 
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)
    updatekey(params,:Ω,Ω)
    updatekey(params,:dΩ,dΩ)


    uhb0 = interpolate(u_walls, U)
    phb0 = interpolate(0.0, P)
    xhb0 = interpolate([uhb0, phb0], X)

    res, rhs = eq_direct_differentiation_steady(params)

    op = AffineFEOperator(res, rhs, X, Y)

    ls = LUSolver()

    solver = LinearFESolver(ls)

    # for i = 1:1:5
        Gridap.solve!(xhb0, solver, op)
        uhb, phb = xhb0
    # end


    if !isnothing(filename)
        writevtk(Ω, filename, cellfields=["uhb" => uhb, "phb" => phb])
    end
    

    return uhb,phb
end
