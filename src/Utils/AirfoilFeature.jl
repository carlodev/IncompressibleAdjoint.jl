uniqueidx(v) = unique(i -> v[i], eachindex(v))


function get_nodes_idx(model, params)
    @unpack AoA = params


    f = (reffe) -> Gridap.Geometry.UnstructuredGrid(reffe)
    Γ = BoundaryTriangulation(model; tags="airfoil")
    trian_tag = Γ.trian
    ref_grids0 = map(f, Gridap.Geometry.get_reffes(trian_tag))
    visgrid0 = Gridap.Visualization.VisualizationGrid(trian_tag, ref_grids0)
    airfoil_points0 = visgrid0.sub_grid.node_coordinates
    idx_uniques = uniqueidx(airfoil_points0)
    airfoil_points = airfoil_points0[idx_uniques]
    idx_above = is_above.(airfoil_points;AoA)
    idx_top = findall(idx_above.> 0)
    idx_bottom = findall(idx_above .< 0)
    perm_top = sortperm(getindex.(airfoil_points[idx_top],1))
    perm_bottom = sortperm(getindex.(airfoil_points[idx_bottom],1))
    IDX_TOP = idx_top[perm_top]
    IDX_BOTTOM = idx_bottom[perm_bottom]
    
    merge!(params, Dict(:IDX_TOP=>idx_uniques[IDX_TOP],:IDX_BOTTOM=>idx_uniques[IDX_BOTTOM]))

    return airfoil_points[IDX_TOP],airfoil_points[IDX_BOTTOM], idx_uniques[IDX_TOP], idx_uniques[IDX_BOTTOM]
end

function get_normals(model, params)

    Γ = BoundaryTriangulation(model; tags="airfoil")
    n_Γ = -1 .*get_normal_vector(Γ)
    n_Γ_dict = Dict("n_Γ" => n_Γ)
    f = (reffe) -> Gridap.Geometry.UnstructuredGrid(reffe)

    n_trian_tag = n_Γ.trian
    n_ref_grids0 = map(f, Gridap.Geometry.get_reffes(n_trian_tag))
    n_visgrid0 = Gridap.Visualization.VisualizationGrid(n_trian_tag, n_ref_grids0)
    n_Γ_data = Gridap.Visualization._prepare_pdata(n_trian_tag, n_Γ_dict, n_visgrid0.cell_to_refpoints)
    airfoil_normals = n_Γ_data["n_Γ"]
    

    _,_,IDX_TOP, IDX_BOTTOM =  get_nodes_idx(model, params)

    airfoil_normals_top = airfoil_normals[IDX_TOP]
    airfoil_normals_bottom = airfoil_normals[IDX_BOTTOM]
    return airfoil_normals_top,airfoil_normals_bottom
end



function get_airfoil_characteristics(model, params; tag="airfoil")
    airfoil_nodes_top, airfoil_nodes_bottom = get_nodes_idx(model, params)
    airfoil_normals_top,airfoil_normals_bottom = get_normals(model, params)

    return airfoil_nodes_top, airfoil_nodes_bottom,airfoil_normals_top,airfoil_normals_bottom
end


function get_aerodynamic_features(params::Dict{Symbol,Any},model, uh,ph)
    @unpack IDX_TOP,IDX_BOTTOM,u_in=params

    q = 0.5 .* u_in^2

    f = (reffe) -> Gridap.Geometry.UnstructuredGrid(reffe)

    Γ = BoundaryTriangulation(model; tags="airfoil")
    cf = Dict("uh"=>uh,"ph"=>ph)
    trian_tag = Γ.trian
    
    ref_grids0 = map(f, Gridap.Geometry.get_reffes(trian_tag))
    visgrid0 = Gridap.Visualization.VisualizationGrid(trian_tag, ref_grids0)
    

    
    pdata = Gridap.Visualization._prepare_pdata(Γ,cf,visgrid0.cell_to_refpoints)

    
    pressure_top = pdata["ph"][IDX_TOP]
    pressure_bottom  = pdata["ph"][IDX_BOTTOM]
    cp_top = pressure_top ./ q
    cp_bottom = pressure_bottom./ q

     return cp_top,cp_bottom
end




function is_above(p; AoA)

    p2 = Point(cosd(AoA),-sind(AoA))
    return is_above(p, p2)
end

function is_above(p, p2; der_slope=-1)
    der = p2[2]/p2[1] .* -der_slope # -0.20 #derivative in trailing edge Manage the (-0.20) factor
    c = -der/8 + der #derivative in leading edge
    
    a = (2*p2[2]-p2[1]*der - c*p2[1])/(-p2[1]^3)
    
    b =(der-c-3*a*p2[1]^2)/(2 * p2[1])
    
    treshold_fun(x) = a*x^3 +b*x^2 +c*x
    
        res = p[2] - treshold_fun(p[1])
        if res > 0
            return 1
        else
            return -1
        end
end

