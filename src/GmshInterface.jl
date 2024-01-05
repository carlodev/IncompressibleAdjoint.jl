using Gmsh
import Gmsh: gmsh

function create_msh(des_tmp::DesignParameters;AoA=0.0, iter =0,chord = 1.0, mesh_ref=1.0)
    spline_points = SplinePoints(des_tmp)
    create_msh(spline_points;AoA=AoA, iter = iter, chord = chord, mesh_ref=mesh_ref)
end


function create_msh(xcontrol::Vector{Float64}; AoA=0.0, iter = 0, chord = 1.0, initfun=circle)
control_points = initialize_control_points(xcontrol; chord = chord, initfun=initfun)
create_msh(control_points; AoA=AoA, iter = iter, chord = chord)
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

function create_msh(spline_points::SplinePoints; AoA=0.0, iter = 0, chord = 1.0, mesh_ref=1.0)


    gmsh.initialize()
    
    gmsh.model.add("Model1")
    Lfront = -2*chord
    Lback = 6*chord
    H=2*chord
    
    
    
    gmsh.model.geo.addPoint(Lback, -H, 0)
    gmsh.model.geo.addPoint(Lback, H, 0)
    
    gmsh.model.geo.addPoint(0.0, -H, 0)
    gmsh.model.geo.addPoint(chord, -H, 0)
    
    bpack = rotate_points([Lback,0.0], AoA)
    gmsh.model.geo.addPoint(bpack[1], bpack[2], 0)

    
    gmsh.model.geo.addPoint(chord, H, 0)
    gmsh.model.geo.addPoint(0.0, H, 0)
    
    gmsh.model.geo.addPoint(Lfront, 0.0, 0)
    
    
    top_points = Int32[]
    bottom_points = Int32[]
    leading_edge_points = Int32[]
   

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
    end
    
    trailing_coordinate =rotate_points([1.0,0.0],AoA)

    trailing = gmsh.model.geo.addPoint(trailing_coordinate[1],trailing_coordinate[2], 0)

    top_le_point = leading_edge_points[end]
    bottom_le_point = leading_edge_points[1]

    origin_point = gmsh.model.geo.addPoint(0.0, 0.0, 0)

    limits_lines = zeros(Int32,4)
    inlet_lines = zeros(Int32,2)
    outlet_lines = zeros(Int32,2)
    


    #External Boundary Lines
    limits_lines[1] =  gmsh.model.geo.addLine(3,4)
    limits_lines[2] =  gmsh.model.geo.addLine(4,1)
    limits_lines[3] =  gmsh.model.geo.addLine(7,6)
    limits_lines[4] =  gmsh.model.geo.addLine(6,2)

    outlet_lines[1] = gmsh.model.geo.addLine(5,2)
    outlet_lines[2] = gmsh.model.geo.addLine(5,1)
 
    inlet_lines[1] = gmsh.model.geo.addCircleArc(7,origin_point,3)
    
    #Internal Lines
    gmsh.model.geo.addLine(trailing,4)
    gmsh.model.geo.addLine(trailing,5)
    gmsh.model.geo.addLine(trailing,6)
    gmsh.model.geo.addLine(top_le_point,7)
    gmsh.model.geo.addLine(bottom_le_point,3)

    #Airfoil Splines
   
    top_spline = gmsh.model.geo.addSpline([top_le_point, top_points...,trailing])
    bottom_spline = gmsh.model.geo.addSpline([bottom_le_point,bottom_points...,trailing])
    inlet_lines[2] = gmsh.model.geo.addSpline(leading_edge_points)
    
    
    gmsh.model.geo.addCurveLoop([7,-12,inlet_lines[2], 11])
    gmsh.model.geo.addCurveLoop([-12, bottom_spline, 8,-1])
    gmsh.model.geo.addCurveLoop([2,-6,-9,8])
    gmsh.model.geo.addCurveLoop([9,5,-4,-10])
    gmsh.model.geo.addCurveLoop([top_spline, 10,-3,-11])

    
    for i = 1:1:5
        gmsh.model.geo.addPlaneSurface([i])
    end
    
    
    #vertical lines
    for i in [outlet_lines...,11,10,12,8]
        gmsh.model.geo.mesh.setTransfiniteCurve(i, minimum([Int32(100), Int32(round(25*mesh_ref))]), "Progression", 1.1)
    end
    
    #inlet and leading edge
    for i in [inlet_lines...,]
        gmsh.model.geo.mesh.setTransfiniteCurve(i, minimum([Int32(50), Int32(round(20*mesh_ref))]), "Progression", 1.0)
    end

    for i in [top_spline,3]
        gmsh.model.geo.mesh.setTransfiniteCurve(i, minimum([Int32(300), Int32(round(25*mesh_ref))]), "Progression", 1.0)
    end
    
    
    for i in [bottom_spline,1]
        gmsh.model.geo.mesh.setTransfiniteCurve(i, minimum([Int32(60), Int32(round(25*mesh_ref))]), "Progression", 1.0)
    end
    
 
    
    #Shear Curves
    for i in [2,9,4]
        gmsh.model.geo.mesh.setTransfiniteCurve(i, minimum([Int32(40), Int32(round(10*mesh_ref))]), "Progression", 1.1)
    end
    
    
    for i = 1:1:5
        gmsh.model.geo.mesh.setTransfiniteSurface(i)
    end
    
    
    
    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    
    #Points
    gmsh.model.addPhysicalGroup(0, [trailing,top_le_point,bottom_le_point], -1, "airfoil")
    gmsh.model.addPhysicalGroup(0, [trailing], -1, "trailing")

    # gmsh.model.addPhysicalGroup(0, [],-1, "inlet")
    gmsh.model.addPhysicalGroup(0, [5],-1, "outlet")
    gmsh.model.addPhysicalGroup(0, [3,4,1,7,6,2],-1,"limits")
    
    #Lines
    gmsh.model.addPhysicalGroup(1, [top_spline,bottom_spline,  inlet_lines[2]],-1, "airfoil")
    gmsh.model.addPhysicalGroup(1, limits_lines,-1, "limits")
    gmsh.model.addPhysicalGroup(1, outlet_lines,-1, "outlet")
    gmsh.model.addPhysicalGroup(1, [inlet_lines[1]], -1,"inlet")
    
    #Surfaces
    gmsh.model.addPhysicalGroup(2, collect(1:5),-1, "fluid")
    
    mkpath("MeshFiles")

    mesh_filename = joinpath("MeshFiles","Mesh$iter.msh")
    gmsh.model.mesh.generate(2)
    gmsh.write(mesh_filename)
    gmsh.finalize()
    return mesh_filename
end
