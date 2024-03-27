abstract type RBFFunction end

struct RBFFunctionGlobalSupport <: RBFFunction
    fun::Function
end

struct RBFFunctionLocalSupport <: RBFFunction
    fun::Function
    support_radius::Float64
end

function get_rbf_function(rbffun::RBFFunction)
    return rbffun.fun
end





function fRBF(x::Real, r::Real, fname::Function)
    x_norm = x / r
    if x_norm > 1
        return 0
    else
        return fname(x / r)
    end
end

function fRBF(x::AbstractVector, r::Real, fname::Function)
    map(v -> fRBF(v, r, fname), x)
end

function fRBF(x::AbstractVector, rbffun::RBFFunctionLocalSupport)
    fname = get_rbf_function(rbffun)
    r = rbffun.support_radius
    return fRBF(x, r, fname)
end

function fRBF(x::AbstractVector, rbffun::RBFFunctionGlobalSupport)
    fname = get_rbf_function(rbffun)
    return fRBF(x, fname)
end

function fRBF(x::Real, fname::Function)
    return fname(x)

end

function fRBF(x::AbstractVector, fname::Function)
    map(v -> fRBF(v, fname), x)
end




#Local Support
function RBF_CP0(x)
    return (1 - x)^2
end

function RBF_CP2(x)
    return (1 - x)^4 * (4 * x + 1)
end

function RBF_CP4(x)
    return (1 - x)^6 * (35 / 3 * x^2 + 6 * x + 1)
end

function RBF_CP6(x)
    return (1 - x)^8 * (32 * x^3 + 25 * x^2 + 8 * x + 1)
end

function RBF_CTPS0(x)
    return (1 - x)^5
end

function RBF_CTPS1(x)
    return 1 + 80 / 3 * x^2 - 40 * x^3 + 15 * x^4 - 8 / 3 * x^5 + 20 * x^2 * log(x)
end


#### Global Support
function RBF_IQB(x)
    return 1 / (1 + x^2)
end

function RBF_GAUSS(x)
    return exp(-x^2)
end


#### RBF Algorithm
function MorphRBF(model::Gridap.Geometry.UnstructuredDiscreteModel, movetag::String, fixtag::Vector{String}, displacement::Vector)
    RBFfun = RBFFunctionLocalSupport(RBF_CP4,4.0)
    MorphRBF(model, movetag, fixtag, displacement,RBFfun)
end


function MorphRBF(model::Gridap.Geometry.UnstructuredDiscreteModel, movetag::String, fixtag::Vector{String}, displacement::Vector,RBFfun::RBFFunction)
    modelgrid0 = copy(model.grid.node_coordinates)
    movenodes = get_boundary_nodes(model, movetag)
    fixnodes = get_boundary_nodes(model, fixtag)
    MorphRBF(modelgrid0, movenodes, fixnodes, displacement, RBFfun)
end

function MorphRBF(modelgrid0::Vector, movenodes::Vector, fixnodes::Vector, displacement::Vector, RBFfun::RBFFunction)
    bnodes = [movenodes; fixnodes]

    D,ΦM = allocate_RBF_matrix(modelgrid0, bnodes, RBFfun)
    rhs = allocate_RBF_rhs(movenodes, bnodes, displacement)
    Shift = solveRBF(modelgrid0, D, rhs,ΦM)

    modelgrid0 .= modelgrid0 .+ Shift


end



function allocate_RBF_matrix(modelgrid0::Vector, bnodes::Vector, RBFfun::RBFFunction)
    #Number of nodes in the domain
    Nc = length(modelgrid0)
    #Number of boundary nodes
    Nb = length(bnodes)
    #Problem dimension 2 or 3
    ND = length(bnodes[1])
    Mbb = ones(Nb, Nb)

    ΦM = zeros(Nc, Nb)

    for (j, p) in enumerate(modelgrid0)
        dist = norm.(p .- bnodes)
        ΦM[j, :] = fRBF(dist, RBFfun)
    end

    for (i, bn) in enumerate(bnodes)
        dist = norm.(bn .- bnodes)

        Mbb[i, :] = fRBF(dist, RBFfun)
    end

    Pb = ones(Nb, ND + 1)

    for (i, bn) in enumerate(bnodes)
        Pb[i, :] = [1, bn...]
    end

    D = hcat(vcat(Mbb, Pb'), vcat(Pb, zeros(ND + 1, ND + 1)))
    return D,ΦM
end

function allocate_RBF_rhs(movenodes::Vector, bnodes::Vector, displacement::Vector)
    ND = length(movenodes[1])
    Nb = length(bnodes)
    Nmove = length(movenodes)
    rhs = [zeros(Nb + ND + 1) for _ in 1:ND]
    for direction in 1:ND
        rhs[direction][1:Nmove] .= getindex.(displacement, direction)
    end

    return rhs
end


function solveRBF(modelgrid0::Vector, D::Matrix{Float64}, rhs::Vector,ΦM::Matrix{Float64})
    ND = length(modelgrid0[1])
    Nc = length(modelgrid0)
    ss = zeros(Nc, ND)
    Nb = size(D)[1]-ND-1

    for (direction, d0) in enumerate(rhs)
        αβ = gmres(D, d0)
        # αβ = minres(D, d0)


        α = αβ[1:Nb]
        β = αβ[end-ND:end]


        #β = β[1] + β[2] x + β[3] y + β[4] z

        βv = VectorValue(β[2:end]...)

        for (j, p) in enumerate(modelgrid0)
            ss[j, direction] = sum(ΦM[j, :] .* α) + β[1] + βv ⋅ p
        end
    end
    Shift = typeof(modelgrid0)(undef, Nc)

    for j = 1:1:Nc
        Shift[j] = VectorValue(ss[j, :]...)
    end

    return Shift
end
