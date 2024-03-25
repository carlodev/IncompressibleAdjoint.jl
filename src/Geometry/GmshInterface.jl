using Gmsh
import Gmsh: gmsh

function create_msh(des_tmp::DesignParameters;AoA=0.0, iter =0,chord = 1.0, mesh_ref=1.0, folder="MeshFiles")
    spline_points = SplinePoints(des_tmp)
    create_msh(spline_points;AoA=AoA, iter = iter, chord = chord, mesh_ref=mesh_ref, folder=folder)
end


function create_msh(xcontrol::Vector{Float64}; AoA=0.0, iter = 0, chord = 1.0, initfun=circle, folder="MeshFiles")
control_points = initialize_control_points(xcontrol; chord = chord, initfun=initfun)
create_msh(control_points; AoA=AoA, iter = iter, chord = chord, folder=folder)
end


function split_splines_points(spline_points::SplinePoints, AoA::Float64; pos=0.065, chord = 1.0)
         

    tgx_coord_top,tgy_coord_top = get_tag_coordinates(spline_points,"top")
    tgx_coord_bottom,tgy_coord_bottom = get_tag_coordinates(spline_points,"bottom")

    threshold = pos * chord
    top_LE_point = findmin(abs.(tgx_coord_top .- threshold))[2]
    bottom_LE_point = findmin(abs.(tgx_coord_bottom .- threshold))[2]
   

    Mtop = rotate_points([tgx_coord_top[top_LE_point+1:end],tgy_coord_top[top_LE_point+1:end]], AoA)
    Mbottom = rotate_points([tgx_coord_bottom[bottom_LE_point+1:end],tgy_coord_bottom[bottom_LE_point+1:end]], AoA)
    Mle = rotate_points([ [reverse(tgx_coord_bottom[1:bottom_LE_point]);tgx_coord_top[1:top_LE_point]],
    [reverse(tgy_coord_bottom[1:bottom_LE_point]);tgy_coord_top[1:top_LE_point]]  ], AoA)


     return Mtop, Mbottom, Mle
end


function rotate_points(Mpoints::Vector{Vector{Float64}}, AoA::Float64)
    xr = Float64[]
    yr = Float64[]
    for (x,y) in zip(Mpoints...)
        xrt= x * cosd(AoA) + y*sind(AoA)
        yrt = -1*x * sind(AoA) + y *cosd(AoA)
        push!(xr,xrt)
        push!(yr,yrt)
    end
    return [xr,yr]
end

function rotate_points(v::Vector{Float64}, AoA::Float64)
    x = v[1] 
    y = v[2]
    
    xrt= x * cosd(AoA) + y*sind(AoA)
    yrt = -1*x * sind(AoA) + y *cosd(AoA)

    return xrt,yrt
end

function find_origin_idx(leading_edge_points::Vector)
    _,idx = findmin(norm.(leading_edge_points))
    return idx
end

"""
    create_msh(spline_points::SplinePoints; AoA=0.0, iter = 0, chord= 1.0, mesh_ref=1.0)

From a set of `spline_points` it creates the .msh file. Incresing `mesh_ref` is increasing the mesh density.
"""
function create_msh(spline_points::SplinePoints; AoA=0.0, iter = 0, chord= 1.0, mesh_ref=1.0, folder="MeshFiles")


    gmsh.initialize()
    
    gmsh.model.add("Model1")
    Lback = 16*chord
    H=6*chord
    offset =2.35
    
    
    gmsh.model.geo.addPoint(Lback, -H, 0)
    gmsh.model.geo.addPoint(Lback, H, 0)
    
    gmsh.model.geo.addPoint(0.0, -H, 0)
    gmsh.model.geo.addPoint(chord, -H, 0)
    
    pback = rotate_points([Lback,0.0], AoA)
    gmsh.model.geo.addPoint(Lback, pback[2], 0)

    
    gmsh.model.geo.addPoint(chord, H, 0)
    gmsh.model.geo.addPoint(0.0, H, 0)
    
    # gmsh.model.geo.addPoint(Lfront, 0.0, 0)
    
    
    top_points = Int32[]
    bottom_points = Int32[]
    leading_edge_points = Int32[]
    leading_edge_points_coordinates = Vector[]


    Mtop, Mbottom, Mle = split_splines_points(spline_points, AoA)
    
    

    for (xp,yp) in zip(Mtop...)
            idx = gmsh.model.geo.addPoint(xp, yp, 0)
            push!(top_points,idx)
    end
    
    for  (xp,yp) in zip(Mbottom...)
            idx = gmsh.model.geo.addPoint(xp, yp, 0)
            push!(bottom_points,idx)
    end
    
    for  (xp,yp) in zip(Mle...)
        idx = gmsh.model.geo.addPoint(xp, yp, 0)
        push!(leading_edge_points,idx)
        push!(leading_edge_points_coordinates, [xp,yp])
    end
    
    trailing_coordinate =rotate_points([1.0,0.0],AoA)

    trailing = gmsh.model.geo.addPoint(trailing_coordinate[1],trailing_coordinate[2], 0)

    trailing_e_top = gmsh.model.geo.addPoint(trailing_coordinate[1],trailing_coordinate[2]+offset, 0)
    trailing_e_bottom = gmsh.model.geo.addPoint(trailing_coordinate[1],trailing_coordinate[2]-offset, 0)


    top_le_point = leading_edge_points[end]
    bottom_le_point = leading_edge_points[1]

    top_e_le_point = gmsh.model.geo.addPoint(0, offset, 0)
    bottom_e_le_point = gmsh.model.geo.addPoint(0,-offset, 0)

    top_back_point = gmsh.model.geo.addPoint(Lback, pback[2]+offset, 0)
    bottom_back_point = gmsh.model.geo.addPoint(Lback, pback[2]-offset, 0)

    origin_point0 = find_origin_idx(leading_edge_points_coordinates)
    origin_point = gmsh.model.geo.addPoint(0.0, 0.0, 0)

    limits_lines = zeros(Int32,4)
    inlet_lines = zeros(Int32,1)
    outlet_lines = zeros(Int32,4)
    


    #External Boundary Lines
    limits_lines[1] =  gmsh.model.geo.addLine(3,4)
    limits_lines[2] =  gmsh.model.geo.addLine(4,1)
    limits_lines[3] =  gmsh.model.geo.addLine(7,6)
    limits_lines[4] =  gmsh.model.geo.addLine(6,2)


    outlet_lines[1] = gmsh.model.geo.addLine(bottom_back_point,1)
    outlet_lines[2] = gmsh.model.geo.addLine(5,bottom_back_point)
    outlet_lines[3] = gmsh.model.geo.addLine(5,top_back_point)
    outlet_lines[4] = gmsh.model.geo.addLine(top_back_point,2)
    inlet_lines[1] = gmsh.model.geo.addCircleArc(7,origin_point,3)

    
    #Internal Lines
    gmsh.model.geo.addLine(trailing,5)
    gmsh.model.geo.addLine(trailing,trailing_e_top)
    gmsh.model.geo.addLine(trailing,trailing_e_bottom)

    gmsh.model.geo.addLine(trailing_e_top,6)
    gmsh.model.geo.addLine(trailing_e_top,top_back_point)

    gmsh.model.geo.addLine(trailing_e_bottom,4)
    gmsh.model.geo.addLine(trailing_e_bottom,bottom_back_point)

    gmsh.model.geo.addLine(top_e_le_point,trailing_e_top)
    gmsh.model.geo.addLine(bottom_e_le_point,trailing_e_bottom)

    gmsh.model.geo.addCircleArc(top_e_le_point,origin_point,bottom_e_le_point)

    gmsh.model.geo.addLine(top_e_le_point,7)
    gmsh.model.geo.addLine(bottom_e_le_point,3)

    gmsh.model.geo.addLine(top_le_point,top_e_le_point)
    gmsh.model.geo.addLine(bottom_le_point,bottom_e_le_point)

    #Airfoil Splines
   
    top_spline = gmsh.model.geo.addSpline([top_le_point, top_points...,trailing])
    bottom_spline = gmsh.model.geo.addSpline([bottom_le_point,bottom_points...,trailing])
    gmsh.model.geo.addSpline(leading_edge_points)
    
    

    #Curve Loops
    gmsh.model.geo.addCurveLoop([9,-21,-19,20])
    gmsh.model.geo.addCurveLoop([19,-23,26,22])
    gmsh.model.geo.addCurveLoop([21,1,-15,-18])
    gmsh.model.geo.addCurveLoop([-17,20,3,-13])
    gmsh.model.geo.addCurveLoop([-18,-23,25,12])
    gmsh.model.geo.addCurveLoop([-22,24,11,-17])
    gmsh.model.geo.addCurveLoop([15,2,-5,-16])
    gmsh.model.geo.addCurveLoop([12,16,-6,-10])
    gmsh.model.geo.addCurveLoop([10,7,-14,-11])
    gmsh.model.geo.addCurveLoop([13,4,-8,-14])
    
    
    for i = 1:1:10
        gmsh.model.geo.addPlaneSurface([i])
    end
    
    
    #vertical outer lines
    for i in [20,13,8,21,15,5]
        gmsh.model.geo.mesh.setTransfiniteCurve(i, 40, "Progression", 1.0)
    end

        
    #vertical inner lines
    for i in [7,11,22,23,12,6]
        gmsh.model.geo.mesh.setTransfiniteCurve(i, maximum([Int32(40), Int32(round(80*mesh_ref))]), "Progression", 1.12) # 1.02
    end
    
    
    #inlet and leading edge
    for i in [9,19,26]
        gmsh.model.geo.mesh.setTransfiniteCurve(i, maximum([Int32(101), Int32(round(20*mesh_ref))]), "Progression", 1.0)
    end

    #top airfoil
    for i in [24,17,3]
        gmsh.model.geo.mesh.setTransfiniteCurve(i, maximum([Int32(100), Int32(round(150*mesh_ref))]), "Progression", 1.0)
    end
    
    #bottom airfoil
    for i in [25,18,1]
        gmsh.model.geo.mesh.setTransfiniteCurve(i, maximum([Int32(50), Int32(round(50*mesh_ref))]), "Progression", 1.0)
    end
    
 
    
    #Shear Curves
    for i in [4,14,10,16,2]
        gmsh.model.geo.mesh.setTransfiniteCurve(i, maximum([Int32(40), Int32(round(20*mesh_ref))]), "Progression", 1.15)
    end
    
    
    for i = 1:1:10
        gmsh.model.geo.mesh.setTransfiniteSurface(i)
    end
    
    
    
    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    
    #Points
    gmsh.model.addPhysicalGroup(0, [trailing,top_le_point,bottom_le_point], -1, "airfoil")
    gmsh.model.addPhysicalGroup(0, [trailing], -1, "trailing")


    gmsh.model.addPhysicalGroup(0, [top_e_le_point], -1, "top_e_le_point")

    gmsh.model.addPhysicalGroup(0, [5,top_back_point,bottom_back_point],-1, "outlet")
    gmsh.model.addPhysicalGroup(0, [3,4,1,7,6,2],-1,"limits")
    
    #Lines
    gmsh.model.addPhysicalGroup(1, [24,25,26],-1, "airfoil")
    gmsh.model.addPhysicalGroup(1, [1,2,3,4],-1, "limits")
    gmsh.model.addPhysicalGroup(1, [5,6,7,8],-1, "outlet")
    gmsh.model.addPhysicalGroup(1, [9],-1, "inlet")
   

    #Surfaces
    gmsh.model.addPhysicalGroup(2, collect(1:10),-1, "fluid")
    
    mkpath(folder)

    mesh_filename = joinpath(folder,"Mesh$iter.msh")

    gmsh.model.mesh.generate(2)
    gmsh.write(mesh_filename)
    gmsh.finalize()
    return mesh_filename
end
