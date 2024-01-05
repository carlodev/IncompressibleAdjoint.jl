"""
    verifykey(params::Dict{Symbol,Any},keyname; val = false)

It check if the dictionary params has the entry `keyname`. If not it adds the new entry with the value val. It is used to add default values
"""
function verifykey(params::Dict{Symbol,Any},keyname; val = false)
    if !haskey(params, keyname)
        merge!(params,Dict(keyname=>val))
    end
end

function updatekey(params::Dict{Symbol,Any},keyname::Symbol, val)
    if !haskey(params, keyname)
        merge!(params,Dict(keyname=>val))
    else
        params[keyname] = val
    end
end


function create_theta_params(dt::Real, θ::Vector, D::Int64, U0, P0, X0)
    dtθ = dt .* θ

    if D == 2
        uhθt = interpolate_everywhere(VectorValue(θ[1], θ[1]), U0)
        uhθ = interpolate_everywhere(VectorValue(dtθ[1], dtθ[1]), U0)

    elseif D == 3
        uhθt = interpolate_everywhere(VectorValue(θ[1], θ[1], θ[1]), U0)
        uhθ = interpolate_everywhere(VectorValue(dtθ[1], dtθ[1], dtθ[1]), U0)
    end
    phθt = interpolate_everywhere(θ[2], P0)
    xhθt = interpolate_everywhere([uhθt, phθt], X0)

    θ_vec = get_free_dof_values(xhθt)

    phθ = interpolate_everywhere(dtθ[2], P0)
    xhθ = interpolate_everywhere([uhθ, phθ], X0)

    dtθ_vec = get_free_dof_values(xhθ)

    θ_params = [θ, θ_vec, dtθ, dtθ_vec]
    return θ_params
end