
function compute_fd_gradient(model::DiscreteModel, contr_points::ControlPoints,tagname::String; δ=0.02, objective_function=compute_drag)
    
    uh0,ph0 = solve_inc_primal(model, tagname; filename=joinpath("FD","FD-0"))
    Nc = length(contr_points.xcoordinate)
    fvref,Sref = obj_fun(model,tagname, ph0,objective_function)

    fv = zeros(Nc)
    S = zeros(Nc)

    vv = (contr_points.tag .== "top") .- (contr_points.tag .== "bottom")
    shift = vv.*δ
    for i = 1:1:Nc
        contr_tmp = perturb_ControlPoint(contr_points, i, shift[i])
        modelname = create_msh(contr_tmp; iter = i, mesh_ref=4)
        model_tmp = GmshDiscreteModel(modelname)
        uh,ph = solve_inc_primal(model_tmp, tagname; filename=joinpath("FD","FD-$i"))
        fv_tmp,Stmp = obj_fun(model_tmp,tagname, ph,objective_function)
        fv[i] = fv_tmp
        S[i] = Stmp

    end
    fv_grad = (fv .- fvref) ./δ
    Sgrad = (S .- Sref) ./δ

    return fv_grad,Sgrad
end