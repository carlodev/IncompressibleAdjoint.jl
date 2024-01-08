function create_ũ_vector(zfv1::AbstractVector)
    return [deepcopy(zfv1), deepcopy(zfv1), deepcopy(zfv1), deepcopy(zfv1)]
end


function update_ũ(ũ_vec::Vector)
    coeff = [2.1875, -2.1875, 1.3125, -0.3125]
    updt_ũ = ũ_vec[1]*coeff[1] + ũ_vec[2] *coeff[2] + ũ_vec[3] *coeff[3] + ũ_vec[4]*coeff[4]
      return updt_ũ
end

function  update_ũ_vector!(ũ_vec::Vector, uh_new::AbstractVector)
    circshift!(ũ_vec,-1)
    ũ_vec[1] = deepcopy(uh_new)
  end