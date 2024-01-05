function istantaneus_CL_CD(model::DiscreteModel, params::Dict{Symbol,Any}, (uh0,UH),(ph0,PH))
    @unpack tagname= params
    Γ = BoundaryTriangulation(model; tags=tagname)
    dΓ = Measure(Γ, 4)
    nΓ = -1 .* get_normal_vector(Γ)
    CD_vec = Float64[]
    CL_vec=Float64[]
    for (uhval,phval) in zip(UH,PH)
        copyto!(uh0.free_values, uhval)
        copyto!(ph0.free_values, phval)
        CDtmp,CLtmp = compute_airfoil_coefficients(uh0,ph0,nΓ,dΓ,params)
        push!(CD_vec,CDtmp)
        push!(CL_vec,CLtmp)
    end
return CD_vec,CL_vec
end

function average_CL_CD(model::DiscreteModel, params::Dict{Symbol,Any}, (uh0,UH),(ph0,PH))
    @unpack t0,tf,dtj,t_primal_start_avg = params
    CD_vec,CL_vec = istantaneus_CL_CD(model, params, (uh0,UH),(ph0,PH))
    idx_start = findfirst(x->x==t_primal_start_avg,collect(t0:dt:tf))
    CD_avg = Statistics.mean(CD_vec[idx_start:end])
    CL_avg = Statistics.mean(CL_vec[idx_start:end])
return CD_avg,CL_avg
end



function average_CL_CD(CD_vec::AbstractVector,CL_vec::AbstractVector, params::Dict{Symbol,Any})
    @unpack t0,tf,dt = params
    idx_start = findfirst(x->x==5.0,collect(t0:dt:tf))
    CD_avg = Statistics.mean(CD_vec[idx_start:end])
    CL_avg = Statistics.mean(CL_vec[idx_start:end])
return CD_avg,CL_avg
end