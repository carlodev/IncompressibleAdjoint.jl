function Base.abs(a::VectorValue)
    VectorValue(map(b-> abs(b), [a...]))
end

function LinearAlgebra.norm(v::VectorValue)
    return sqrt(sum(map(a-> v[a].^2, eachindex(v))))
 
 end

"""
    get_radius_shift(shift)

Given the shift, it computes the radius-influence of it.
"""
function get_radius_shift(shift::VectorValue)
    return sqrt(2)*length(shift)^2 * norm(shift) #sum(abs(shift))
end


function compute_point_dist(coordinate, point::VectorValue)
    return norm.(coordinate .- point)
end


#### This is the morphing function
function morph_kernel(x::Real,r::Real)
    @assert x >= 0 
    xr = x/r
    if xr<1
        return (1-xr)^6 * (35/3 * xr^2 + 6 * xr +1)
        #return (cos(pi/2 * (x / r)^2 ))^2
    else 
        return 0.0
    end
end

function scale_shift(shift::VectorValue, x::Vector,r::Real)
    scale = map(a-> morph_kernel(a,r),x)
    return map(s->s*shift ,scale)
end

function get_shift_vec(model::DiscreteModel, point::VectorValue, r::Real, shift::VectorValue)
    dist = compute_point_dist(model.grid.node_coordinates, point)
    shift_vec = scale_shift(shift, dist, r)
    return shift_vec
end

function get_shift_vec(model::DiscreteModel, point::VectorValue, shift::VectorValue)
    r = get_radius_shift(shift)
    @assert norm(shift) < r
    shift_vec = get_shift_vec(model, point, r, shift)
    return shift_vec
end



function compute_radius(scale::Vector,radius::Float64)
    return radius .*(scale.>0)
end

### Geometry Morph ###
function conv_VectorValue(v::VectorValue)
    [v...]
end

function get_visgrid(model::DiscreteModel, tagname::String)
    Γ = BoundaryTriangulation(model; tags=tagname)

    f = (reffe) -> Gridap.Geometry.UnstructuredGrid(reffe)
    ref_grids = map(f, Gridap.Geometry.get_reffes(Γ))
    visgrid = Gridap.Visualization.VisualizationGrid(Γ, ref_grids)
    return visgrid
end

function get_boundary_nodes(model::DiscreteModel, tagname::String)
    visgrid = get_visgrid(model, tagname)
    idx = get_nodes_idx(model, tagname)

    boundary_nodes=visgrid.sub_grid.node_coordinates[idx]
    return boundary_nodes
end

function get_nodes_idx(model::DiscreteModel, tagname::String)
    visgrid = get_visgrid(model, tagname)
    visgrid_ = conv_VectorValue.(visgrid.sub_grid.node_coordinates)

    unique_idx = unique(i -> round.(visgrid_[i], digits=5), eachindex(visgrid_))
    #Simple constrain added
    idx = unique_idx[findall(b -> b == true, map(a -> getindex(a, 1) > -10.04, visgrid_[unique_idx]))]

    return idx
end

function get_normals(model::DiscreteModel, tagname::String)
    Γ = BoundaryTriangulation(model; tags=tagname)
    nΓ = -1 .* get_normal_vector(Γ)

    visgrid = get_visgrid(model, tagname)

    idx = get_nodes_idx(model, tagname)
    pdat = Gridap.Visualization._prepare_pdata(Γ, Dict("nΓ" => nΓ), visgrid.cell_to_refpoints)
    nnΓ = pdat["nΓ"]
    return nnΓ[idx]
end

# function boundary_nodes_index(model::DiscreteModel, tag::String)
#     labels = get_face_labeling(model)
#     idx_tag_to_name = findfirst(a -> a == tag, labels.tag_to_name)

#     idx_tmp = Int64[]
#     for (i, val) in enumerate(labels.d_to_dface_to_entity[1])
#         if sum(val .== labels.tag_to_entities[idx_tag_to_name]) > 0
#             push!(idx_tmp, i)
#         end
#     end
#     return idx_tmp
# end

# function map_model_to_boundary(model::DiscreteModel, tagname::String)
#     bni = boundary_nodes_index(model, tagname)
#     aa = model.grid.node_coordinates[bni]
#     visgrid = get_visgrid(model, tagname)
#     visgrid_ = conv_VectorValue.(visgrid.sub_grid.node_coordinates)
#     idx = get_nodes_idx(model, tagname)

#     mv = Int64[]
#     for vg in visgrid_[idx]
#         for (i, aai) in enumerate(aa)
#             if norm(conv_VectorValue(aai) - vg) < 1e-5
#                 push!(mv, i)
#             end
#         end

#     end

#     @assert length(mv) == length(idx)
#     return bni[mv]
# end


# function move_node!(model::DiscreteModel,tagname::String, I::Int64, δ::VectorValue)
#     nng = get_normals(model, tagname)[I]
#     mvv = map_model_to_boundary(model,tagname)[I]
#     sv = get_shift_vec(model, model.grid.node_coordinates[mvv], - nng .*δ)
#     model.grid.node_coordinates .= model.grid.node_coordinates .+ sv
# end

# function move_node!(model::DiscreteModel,nng::VectorValue,mvv::Int64, δ::VectorValue)
#     sv = get_shift_vec(model, model.grid.node_coordinates[mvv], - nng .*δ)
#     model.grid.node_coordinates .= model.grid.node_coordinates .+ sv
# end

# function move_node!(model::DiscreteModel,tagname::String, I::Vector{Int64}, δ::VectorValue)
#     nng = get_normals(model, tagname)[I]
#     mvv = map_model_to_boundary(model,tagname)[I]
#     map((a,b)->move_node!(model,a,b, δ),nng,mvv)
# end


# function update_model!(model::DiscreteModel, sv::Vector, δ::VectorValue)
#     shifts_u = limiter(δ,sv)
#     model.grid.node_coordinates .= model.grid.node_coordinates .+ shifts_u
# end

# function limiter(δ::VectorValue,sv::Vector)
#     sv1 = getindex.(sv,1)
#     sv2 = getindex.(sv,2)
    
#     i1 = findall(a-> abs(a).>δ[1], sv1)
#     i2 = findall(a->abs(a).>δ[2], sv2)
    
    
#     for i in i1
#     sv1[i] = δ[1] *sign(sv1[i])
#     end
    
#     for i in i2
#         sv2[i] =δ[2]*sign(sv2[i])
#     end
        
    
#     shifts_u = map((a,b)->VectorValue(a,b), sv1,sv2)
#     return shifts_u
# end



###
#MorphUtilities
###

"""
Given a vector of the position of control points in `control_points_position` it provides the indexes of the points on the boundary
"""
function get_idx_control_points(control_points_position::Vector{Float64},model,params)
    @unpack AoA,tagname = params
    control_points_position = control_points_position.*cosd(AoA)
    bnodes = get_boundary_nodes(model, tagname)
    IDX_CONTR = Int64[]
    for cpp in control_points_position
        idx_control=sortperm(abs.(getindex.(bnodes,1) .- cpp))
        push!(IDX_CONTR, idx_control[1])
        cond = true
        i = 2

        while cond
        if abs.(bnodes[idx_control[1]][2] - bnodes[idx_control[i]][2]) >1e-3
            push!(IDX_CONTR, copy(idx_control[i]))
            cond = false
        end
        i = i+1
    end
    end
    return IDX_CONTR
end

function get_control_boundary(control_points_position, model, params)
    @unpack tagname = params
    idxs = get_idx_control_points(control_points_position,model,params)
    bnodes = get_boundary_nodes(model, tagname)[idxs]
    bnormals = get_normals(model, tagname)[idxs]
    bnodes_top,bnodes_bottom,bnormals_top,bnormals_bottom= split_top_bottom_boundaries(bnodes,bnormals,params)
    return bnodes_top,bnodes_bottom,bnormals_top,bnormals_bottom
end
 

function is_top(bnodes,params)
    @unpack AoA = params
    topnodes = typeof(bnodes[1])[]
    bottomnodes = typeof(bnodes[1])[]

    m = tand(AoA)
    top=Bool.(zeros(length(bnodes)))
    for (i,n) in enumerate(bnodes)
        if -n[1]*m < n[2]
            top[i]=true
            push!(topnodes,n)
        else
            push!(bottomnodes,n)
        end
    end
    return top,topnodes,bottomnodes
end

function normal_correction(bnormals,top)
    topnormals = typeof(bnormals[1])[]
    bottomnormals = typeof(bnormals[1])[]
    for (i,t) in enumerate(top)
        if t && bnormals[i][2]<0
            bnormals[i] = -1 .*bnormals[i]
        elseif t==false && bnormals[2]>0
            bnormals[i] = -1 .*bnormals[i]
        end
    end

    for (t,bn) in zip(top,bnormals)
        if t 
            push!(topnormals,bn)
        else
            push!(bottomnormals,bn)     
        end
    end

    return topnormals,bottomnormals
end

function split_top_bottom_boundaries(bnodes,bnormals,params)
 
    istop,topnodes,bottomnodes=is_top(bnodes,params)
    topnormals,bottomnormals =  normal_correction(bnormals,istop)
    
    idxtop=sortperm(map(tn->tn[1],topnodes))
    idxbottom=sortperm(map(bn->bn[1],bottomnodes))
  return topnodes[idxtop],bottomnodes[idxbottom],topnormals[idxtop],bottomnormals[idxbottom]
    

end