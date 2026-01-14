function [strained_position,sorted_idx,f_index] = generate_strain(z_s,lat_par,N,o,startz,stack)


  pos_vector   = zeros(N*length(z_s),1);
  f_index      = zeros(N*length(z_s),1);
  [z_s,sorted_idx] = sort(z_s);

  ll = 1;
  for l = 1:N
    for s = 1:length(z_s)
      pos_vector(ll) = z_s(s) + lat_par(o)*(l-1);
      f_index(ll)    = sorted_idx(s); % pointer to the right form factor at this position

      ll = ll + 1;
    end
  end


  pos_vector = pos_vector - pos_vector(1) + startz;

  % The new strained position of atom s in unit cell N in stack o
  % The atom closest to the interface is unaffected by strain and it is
  % modified by dinterface

  strained_position = zeros(size(pos_vector));
  displacement          = zeros(size(strained_position));

  % pos_vector is sorted but may have more than one atom at the lowest z
  strained_position(1) = pos_vector(1);

  bottom_strain_boundary = round(stack{o}.bottom_strain_end);
  top_strain_boundary    = length(pos_vector)-round(stack{o}.top_strain_end);

  for ll = 2:length(pos_vector)
    if ll > stack{o}.bottom_strain_end % add no extra strain 
      displacement(ll) = pos_vector(ll) - pos_vector(ll-1);
    else
      if (pos_vector(ll) - pos_vector(ll-1)) == 0 % if two atoms are at the same coordinate in z
        displacement(ll) = 0;
      else
        displacement(ll) = (pos_vector(ll) - pos_vector(ll-1)) + stack{o}.bottom_strain_amplitude*(pos_vector(bottom_strain_boundary) - pos_vector(ll));
      end
    end
    if ll > length(pos_vector) - stack{o}.top_strain_end
      if (pos_vector(ll) - pos_vector(ll-1)) == 0
        displacement(ll) = 0;
      else
        displacement(ll) = displacement(ll) + stack{o}.top_strain_amplitude*( pos_vector(ll) - pos_vector(top_strain_boundary));
      end
    end
    strained_position(ll) = strained_position(ll-1) + displacement(ll);
  end
end