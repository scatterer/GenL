function output = parratt_matrix_repeated_rhobuiltin(Q,lambda,stack,control,instrument)

output.trans      = [];
output.absorption = [];
output.I          = [];
output.refl       = [];


% Assuming the z axis of the unit cell is aligned with Q
Z = (1:118);
Zt = {'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr',...
  'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',...
  'In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',...
  'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra',...
  'Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'};

% Read poscar and generate lattice vectors, area and volume.
% Generate strain profiles

Q0 = (-control.maxQ0:control.stepQ0:control.maxQ0)';
strain = cell(length(stack),1);
for o = 1:length(stack)

    layer   = stack{o};
    N       = round(layer.N);
    if isempty(stack{o}.pre_calc_f)
      calculate_from_scratch = true;
      pre_calc_f = [];
      [type,type_nr,r,a1,a2,a3] = read_poscar(stack{o}.filename);
      pre_calc_f(1).a1 = a1;
      pre_calc_f(1).a2 = a2;
      pre_calc_f(1).a3 = a3;
      pre_calc_f(1).type = type;
      pre_calc_f(1).type_nr = type_nr;
      pre_calc_f(1).r = r;
    else
      calculate_from_scratch = false;
      pre_calc_f = stack{o}.pre_calc_f;
    end

    as = sqrt(stack{o}.area_scale);

   switch stack{o}.direction
     case 1
       scaling = pre_calc_f(1).a1*stack{o}.scale;
       area(o)    = abs(norm(cross(pre_calc_f(1).a2*as,pre_calc_f(1).a3*as)));
       v_uc    = abs(dot(scaling,cross(pre_calc_f(1).a2*as,pre_calc_f(1).a3*as)));
     case 2
       scaling = pre_calc_f(1).a2*stack{o}.scale;
       area(o)    = abs(norm(cross(pre_calc_f(1).a1*as,pre_calc_f(1).a3*as)));
       v_uc    = abs(dot(scaling,cross(pre_calc_f(1).a1*as,pre_calc_f(1).a3*as)));
     case 3
       scaling = pre_calc_f(1).a3*stack{o}.scale;
       area(o)    = abs(norm(cross(pre_calc_f(1).a1*as,pre_calc_f(1).a2*as)));
       v_uc    = abs(dot(scaling,cross(pre_calc_f(1).a1*as,pre_calc_f(1).a2*as)));
   end

    lat_par(o) = norm(scaling);
    z_s = zeros(size(pre_calc_f(1).r,1),1);
    for l = 1:size(pre_calc_f(1).r,1)
      z_s(l) = pre_calc_f(1).r(l,:)*scaling';
    end
    z_ss{o} = z_s;
    r_ss{o} = pre_calc_f(1).r;
    


    k = 0;
    l = 0;
    f = zeros(length(Q0),length(z_ss{o}));

    for i = 1:size(pre_calc_f(1).type,1) %
      if ~isempty(pre_calc_f(1).type{i})
        k = k + 1;
        cur_Z = Z(strcmp(Zt,pre_calc_f(1).type{i}));
        if calculate_from_scratch
          [g , ~, ] = Form_factors(Q0(:),  Read_form_factor_coefficients(cur_Z,lambda), 100/100);
          B = DB_prefactor(cur_Z);
          pre_calc_f(k).f = g.*exp(-B*(Q0/4/pi).^2);
        end
        for j = 1:pre_calc_f(1).type_nr(k)
          l = l + 1;
                
          f(:,l)  = pre_calc_f(k).f;
          coeffs(l) =  Read_form_factor_coefficients(cur_Z,lambda);

        end
      end
    end
    ff{o} = f;

    full_stack = false;
    % Totthick when doing the repeated uc trick
    if o == 1
      if full_stack
        totthick = lat_par(o)*N;
        last_atom_z = lat_par(1)*(N-1)+max(z_ss{1});
        strain{o} = last_atom_z;
      else
        totthick = lat_par(o)*3;
        last_atom_z = lat_par(1)*(3-1)+max(z_ss{1});
        strain{o} = last_atom_z;
      end
    end
    if o > 1
      if o > 2
        last_atom_z = strain{o-1}(end);
      end

      startz = last_atom_z + stack{o}.dinterface;
      if stack{o}.roughness
        sorted_idx = {};
        strained_position = {};
        sigma  = stack{o}.sigma;
        nvector = max(1,round(N-3*sigma)):1:max(1,round(N+3*sigma));

        B   = zeros(size(nvector));

        if sigma == 0
          A = 1;
        else
          A = 1/(sqrt(2*pi)*sigma);
        end

        for NN = nvector
          if sigma ==0
            B(NN) = 1;
          else
            B(NN)         = A*exp(-(NN-N)^2/(2*sigma^2));
          end
        end
        for m = nvector
          [strained_position{m},sorted_idx{m}] = generate_strain(z_ss{o},lat_par,m,o,startz,stack);
        end
      else
        [strained_position,sorted_idx] = generate_strain(z_ss{o},lat_par,N,o,startz,stack);
      end

      strain{o} = strained_position;
      sorted{o} = sorted_idx;
    end
end

 if stack{o}.roughness
   s = strain{end};
   for lll = 1:length(s)
 
    tt(lll) = max([max(s{lll}),0]);
   end
   totthick = max(tt);

else
  totthick = max(strain{end});
 end
 tic
disp('Generating density...')
for o = 1:length(stack)
    N = round(stack{o}.N);
    if o == 1 % substrate
      slices = control.slices;
      dz     = lat_par(o)/slices; % A
      vacuum_slices = round(control.vacuum_thick/dz);
      vacuum_thick  = dz*vacuum_slices;

      z = -vacuum_thick:dz:totthick+vacuum_thick;
      rho_e = zeros(size(z));

      if ~full_stack
        for l = 1:3 
            for s = 1:length(z_ss{o})
                idx = z >= z_ss{o}(s) + lat_par(o)*(l-1) - 1*lat_par(o) & z <= z_ss{o}(s) + lat_par(o)*(l-1) + 1*lat_par(o);

                rho_e(idx) = rho_e(idx) + 1/area(o)*trapz(Q0,ff{o}(:,s).*exp(1i*Q0*(z(idx) - z_ss{o}(s) - lat_par(o)*(l-1))));
           end
        end
        last_atom_z = lat_par(o)*(3-1)+max(z_ss{o});
        substrate_end_z = last_atom_z;
      else
        %for calculating the full stack
        for l = 1:N  
          for s = 1:length(z_ss{o})
            idx = z >= z_ss{o}(s) + lat_par(o)*(l-1) - 1*lat_par(o) & z <= z_ss{o}(s) + lat_par(o)*(l-1) + 1*lat_par(o);
            rho_e(idx) = rho_e(idx) + 1/area(o)*trapz(Q0,ff{o}(:,s).*exp(1i*Q0*(z(idx) - z_ss{o}(s) - lat_par(o)*(l-1))));
          end
        end
        last_atom_z = lat_par(o)*(N-1)+max(z_ss{o});
        substrate_end_z = lat_atom_z;
      end

        start_idx = slices;
        z_0   = z(1:vacuum_slices+start_idx);
        z_1   = z(vacuum_slices+start_idx:vacuum_slices+start_idx*2);
        z_2   = z(vacuum_slices+start_idx*2:vacuum_slices+start_idx*3); % was 5
        rho_0 = rho_e(1:vacuum_slices + start_idx);
        rho_1 = rho_e(vacuum_slices + start_idx:vacuum_slices+start_idx*2);
        rho_2 = rho_e(vacuum_slices + start_idx*2:vacuum_slices + start_idx*3);
 
    else % film

      if stack{o}.roughness
        startz = last_atom_z + stack{o}.dinterface;

        sigma  = stack{o}.sigma;
        nvector = max(1,round(N-3*sigma)):1:max(1,round(N+3*sigma));

        B   = zeros(size(nvector));
        
        if sigma == 0
          A = 1;
        else
          A = 1/(sqrt(2*pi)*sigma);
        end
        
        

        for NN = nvector

          % including roughness
          if sigma ==0
            B(NN) = 1;
          else
            B(NN)         = A*exp(-(NN-N)^2/(2*sigma^2));
          end
        end

        % Assuming rho_e is large enough to cover all the unit cells
        rho_e_rough = cell(length(nvector),1);
        for m = nvector
         rho_e_rough{m} = generate_density(z_ss,rho_e,sorted_idx{m},o,m,lat_par,area,Q0,ff,z,strained_position{m});
        end

      else % no roughness
        % We put the next unit cell ontop of the last atom in the previous unit cell plus some interface distance.
        % Beware if none of the coordinates in the new unit cell are zero, you get an additional distance.
        % rho_e = generate_density(z_ss,rho_e,sorted{o},o,N,lat_par,area,Q0,ff,z,strain{o});
         ll = 1;
         sorted_idx = sorted{o};
         for l = 1:N
           for s = 1:length(z_ss{o})
             pos_z                   = strain{o}(ll);
             idx = z >= pos_z - 1*lat_par(o) & z <=  pos_z + 1*lat_par(o);

             rho_e(idx) = rho_e(idx) + 1/area(o)*trapz(Q0,ff{o}(:,sorted_idx(s)).*exp(1i*Q0*(z(idx) - pos_z)));
             ll = ll + 1;
           end
         end
        last_atom_z = startz + lat_par(o)*(N-1) + max(z_ss{o});

      end
    end

end
disp('done')
toc

% We have the density vs z. Now lets propagate

if ~full_stack && ~control.dotransmission
  rho_0r = rho_0(end:-1:1);
  rho_1r = rho_1(end:-1:1);
  rho_rest = rho_e(vacuum_slices+start_idx*2:end);
  
  if control.plot_density
    figure(2)
    clf
    plot(z_0,real(rho_0),'-r')
    hold on;
    plot(z_1,real(rho_1),'-k')
    plot(z_2,real(rho_2),'-b')
    z_rest   = z(vacuum_slices+start_idx*2:end);

    plot(z_rest,real(rho_rest),'-m');
  end
  disp('Propagating...')
  tic
  if stack{o}.roughness
    substrate_end = vacuum_slices+start_idx*2;
    if control.pol == 0 || control.pol == 1
      output.refl   = propagate_vectorized_chunks_roughness(Q,lambda,rho_0r,rho_1r,round(stack{1}.N),rho_e_rough,dz,substrate_end,nvector,B,control.pol);
    elseif control.pol == 2
      refl_p   = propagate_vectorized_chunks_roughness(Q,lambda,rho_0r,rho_1r,round(stack{1}.N),rho_e_rough,dz,substrate_end,nvector,B,0);
      refl_s   = propagate_vectorized_chunks_roughness(Q,lambda,rho_0r,rho_1r,round(stack{1}.N),rho_e_rough,dz,substrate_end,nvector,B,1);

      output.refl = (refl_s + cosd(instrument.theta_m*2)^2.*refl_p)./(1 + cosd(instrument.theta_m*2).^2);
    end
  else % no roughness

    rho_restr = rho_rest(end:-1:1);

    if control.pol == 0 || control.pol == 1
      output.refl = propagate_vectorized_chunks(Q,lambda,rho_0r,rho_1r,round(stack{1}.N),rho_restr,dz,stack{o}.sigma,control.pol);

    elseif control.pol == 2
      refl_p = propagate_vectorized_chunks(Q,lambda,rho_0r,rho_1r,round(stack{1}.N),rho_restr,dz,stack{o}.sigma,0);
      refl_s = propagate_vectorized_chunks(Q,lambda,rho_0r,rho_1r,round(stack{1}.N),rho_restr,dz,stack{o}.sigma,1);

      % Parratt's formalism includes the cosd(theta*2)^2 from the sample
      % Here we correct the pi polarization for the polarization picked up by the monochromator
      % and then we divide by the intensity just before the sample, which has polarization 1+cosd(theta_m*2)^2
      output.refl = (refl_s + cosd(instrument.theta_m*2)^2.*refl_p)./(1 + cosd(instrument.theta_m*2).^2);
    end
  end
toc
disp('done')
else % for debugging
end
end