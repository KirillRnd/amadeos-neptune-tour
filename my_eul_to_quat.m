function q=my_eul_to_quat(alpha_rad, beta_rad, gamma_rad, order)

qtmp=eul_to_quat(alpha_rad, beta_rad, gamma_rad, order);
q = [qtmp.w qtmp.x qtmp.y qtmp.z];
end
