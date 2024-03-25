function istantaneus_CD_CL(model::DiscreteModel, params::Dict{Symbol,Any}, (uh0,UH),(ph0,PH))
    @unpack tagname= params
    Γ = BoundaryTriangulation(model; tags=tagname)
    dΓ = Measure(Γ, 4)
    nΓ = -1 .* get_normal_vector(Γ)
    CD_vec = Float64[]
    CL_vec = Float64[]
    for (uhval,phval) in zip(UH,PH)
        copyto!(uh0.free_values, uhval)
        copyto!(ph0.free_values, phval)
        CDtmp,CLtmp = compute_airfoil_coefficients(uh0,ph0,nΓ,dΓ,params)
        push!(CD_vec,CDtmp)
        push!(CL_vec,CLtmp)
    end
return CD_vec,CL_vec
end

function average_CD_CL(model::DiscreteModel, params::Dict{Symbol,Any}, (uh0,UH),(ph0,PH))
    CD_vec,CL_vec = istantaneus_CD_CL(model, params, (uh0,UH),(ph0,PH))
    CD_avg,CL_avg = average_CD_CL(CD_vec,CL_vec, params)
return CD_avg,CL_avg
end



function average_CD_CL(CD_vec::AbstractVector,CL_vec::AbstractVector, params::Dict{Symbol,Any})
    @unpack t0,tf,dt = params
    idx_start = findfirst(x->x==tf*0.75,collect(t0:dt:tf))
    CD_avg = Statistics.mean(CD_vec[idx_start:end])
    CL_avg = Statistics.mean(CL_vec[idx_start:end])
return CD_avg,CL_avg
end


function average_field!(fh::SingleFieldFEFunction,fv::AbstractArray)
    Np = length(fv[1])
    fvavg = zeros(Np)
    for i = 1:1:Np
        vals = getindex.(fv,i)
        fvavg[i] = Statistics.mean(vals)
    end
    copyto!(fh.free_values,fvavg)
end