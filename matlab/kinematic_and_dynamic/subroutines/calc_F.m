function [F_tot,strain] = calc_F(Q,z_s,lat_par,N,i,startz,stack,f,v_uc)
  r0 = 2.814042735053330e-05;
  [pos_vector,sorted_idx,f_index] = generate_strain(z_s,lat_par,N,i,startz,stack);

  strain{i} = pos_vector;

  F_tot  = zeros(size(Q));
  F_tot_test  = zeros(size(Q));

  m = length(Q);

  slices = min(5000,m);

  iter = floor(m/slices);

  pre_factor = -1i*4*pi*r0*lat_par(i)/v_uc;
 % f0 = f;
  f = f*pre_factor;
  iQ = 1i*Q;
  for l = 1:iter
    idx        = slices*(l-1)+1:slices*l;
    % A = f(idx,f_index);
    % C = iQ(idx) * pos_vector';
    % B = exp(C);
    % F_tot_test(idx) = sum(A.*B,2);
    F_tot(idx) = sum(f(idx,f_index).*exp(iQ(idx) * pos_vector'),2);
  end
  % 
  % F_tota = zeros(size(Q));
  % for l = 1:length(pos_vector)
  %   %g      = (4*pi*r0*lat_par(i)/v_uc)*f(:,f_index(l))./Q;
  %   %g0     = (4*pi*r0*lat_par(i)/v_uc)*f0(:,f_index(l))./Q;
  %   %F_tot  = F_tot - 1i*g.*exp(1i.*(Q.*pos_vector(l)));
  % 
  %  F_tota  = F_tota - f0(:,f_index(l)).*exp(1i.*(Q.*pos_vector(l)));
  % end
  % F_tota = F_tota*(1i*4*pi*r0*lat_par(i)/v_uc);

end