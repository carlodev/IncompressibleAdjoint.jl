using IncompressibleAdjoint.Equations
using IncompressibleAdjoint
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
