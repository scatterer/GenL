function rho_e = generate_density(z_ss,rho_e,sorted_idx,o,N,lat_par,area,Q0,ff,z,strain)
  ll = 1;
  for l = 1:N
    for s = 1:length(z_ss{o})
      pos_z                   = strain(ll);
      idx = z >= pos_z - 1*lat_par(o) & z <=  pos_z + 1*lat_par(o);

      rho_e(idx) = rho_e(idx) + 1/area(o)*trapz(Q0,ff{o}(:,sorted_idx(s)).*exp(1i*Q0*(z(idx) - pos_z)));
      ll = ll + 1;
    end
  end


end

