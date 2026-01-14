function [output,stack] = calc_kinematic_full(Q,lambda,stack,control,instrument)
  % Assuming the z axis of the unit cell is aligned with Q
  Z = (1:118);
  Zt = {'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr',...
    'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',...
    'In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',...
    'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra',...
    'Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'};

  r0 = 2.814042735053330e-05;
  strain = cell(length(stack),1);
  R = 0;

  for i = 1:length(stack)

    if isempty(stack{i}.pre_calc_f)
      calculate_from_scratch = true;
      pre_calc_f = [];
      [type,type_nr,r,a1,a2,a3] = read_poscar(stack{i}.filename);
      pre_calc_f(1).a1 = a1;
      pre_calc_f(1).a2 = a2;
      pre_calc_f(1).a3 = a3;
      pre_calc_f(1).type = type;
      pre_calc_f(1).type_nr = type_nr;
      pre_calc_f(1).r = r;
    else
      calculate_from_scratch = false;
      pre_calc_f = stack{i}.pre_calc_f;
    end

    as = sqrt(stack{i}.area_scale);

    switch stack{i}.direction
      case 1
        scaling = pre_calc_f(1).a1*stack{i}.scale;
        area    = abs(norm(cross(pre_calc_f(1).a2*as,pre_calc_f(1).a3*as)));
        v_uc    = abs(dot(scaling,cross(pre_calc_f(1).a2*as,pre_calc_f(1).a3*as)));
      case 2
        scaling = pre_calc_f.a2(1)*stack{i}.scale;
        area    = abs(norm(cross(pre_calc_f(1).a1*as,pre_calc_f(1).a3*as)));
        v_uc    = abs(dot(scaling,cross(pre_calc_f(1).a1*as,pre_calc_f(1).a3*as)));
      case 3
        scaling = pre_calc_f(1).a3*stack{i}.scale;
        area    = abs(norm(cross(pre_calc_f(1).a1*as,pre_calc_f(1).a2*as)));
        v_uc    = abs(dot(scaling,cross(pre_calc_f(1).a1*as,pre_calc_f(1).a2*as)));
    end

    lat_par(i) = norm(scaling);

    z_s = zeros(size(pre_calc_f(1).r,1),1);
    for l = 1:size(pre_calc_f(1).r,1)
      z_s(l) = pre_calc_f(1).r(l,:)*scaling';
    end

    k = 0;
    l = 0;

    f  = zeros(length(Q),length(z_s));
    %f0 = zeros(1        ,length(z_s));
    sub_F = zeros(size(Q));
    sub_F0 = zeros(size(Q));
    for j = 1:size(pre_calc_f(1).type,1) %
      if ~isempty(pre_calc_f(1).type{j})
        k = k + 1;
        cur_Z = Z(strcmp(Zt,pre_calc_f(1).type{j}));
        if calculate_from_scratch
          [g , ~, ]  = Form_factors(Q(:),  Read_form_factor_coefficients(cur_Z,lambda), 100/100);

          B = DB_prefactor(cur_Z);
          pre_calc_f(k).f  = g.*exp(-B*(Q/4/pi).^2)./Q;
          if i == 1
            [g0 , ~, ] = Form_factors(0   ,  Read_form_factor_coefficients(cur_Z,lambda), 100/100);
            pre_calc_f(k).f0 = g0;
          end
        end
        for m = 1:pre_calc_f(1).type_nr(k)
          l = l + 1;

          f(:,l)  = pre_calc_f(k).f;
          if i == 1
            sub_F  = sub_F  + pre_calc_f(k).f.*exp(1i*Q.*pre_calc_f(1).r(l,:)*scaling'); %r(l,direction)*scaling %f.*exp(1i*Q*r(l,direction)*scaling);
            sub_F0 = sub_F0 + pre_calc_f(k).f0;
          end
        end
      end
    end

    N  = round(stack{i}.N);


    if i == 1
      last_atom_z = 0;

      g  = 4*pi.*sub_F./v_uc.*lat_par(i).*r0;
      g0 = 4*pi.*sub_F0./v_uc.*lat_par(i).*r0./Q;
      N  = round(stack{i}.N);

      if stack{i}.roughness
        sigma   = stack{i}.sigma;
        nvector = max(1,round(N-3*sigma)):1:max(1,round(N+3*sigma));
        B       = zeros(size(nvector));

        if sigma == 0
          A = 1;
        else
          A   = 1/(sqrt(2*pi)*sigma);
        end

        for NN = nvector
          if sigma == 0
            B(NN) = 1;
          else
            B(NN) = A*exp(-(NN-N)^2/(2*sigma^2));
          end
        end


        F = zeros(length(Q),length(nvector));
        r = 1;
        tally = 0;
        for m = nvector
          F(:,r) = R; % R is just a complex number here at a given angle?
          F(:,r) = F(:,r) + exp(1i*Q*(stack{i}.dinterface)).*-1i.*g.*(1 - exp(1i*(Q*lat_par(i)-2*g0)*m) )./(1 - exp(1i*(Q*lat_par(i)-2*g0)));
          F(:,r) = F(:,r)*sqrt(B(m));
          tally = tally + sqrt(B(m));
          r = r + 1;
        end
        F_tot = sum(F/tally,2); % Obs! Intensity, not amplitude squared!


      else
        %F_tot = -1i*g.*(1 - exp(1i*(Q*lat_par-2*g0)*N) )./(1 - exp(1i*(Q*lat_par-2*g0)));
        F_tot = -1i*g.*(1 - exp(1i*(Q*lat_par(i)-0*2*g0)*N) )./(1 - exp(1i*(Q*lat_par(i)-0*2*g0)));

      end
      R = R + F_tot;
      last_atom_z = N*lat_par(i) + max(z_s);
    else

      if i > 2
        last_atom_z = strain{i-1}(end);
      end

      startz = last_atom_z + stack{i}.dinterface;

      if stack{i}.roughness
        sigma   = stack{i}.sigma;
        nvector = max(1,round(N-3*sigma)):1:max(1,round(N+3*sigma));
        B       = zeros(size(nvector));

        if sigma == 0
          A = 1;
        else
          A   = 1/(sqrt(2*pi)*sigma);
        end

        for NN = nvector
          if sigma == 0
            B(NN) = 1;
          else
            B(NN) = A*exp(-(NN-N)^2/(2*sigma^2));
          end
        end


        F = zeros(length(Q),length(nvector));
        r = 1;
        tally = 0;
        for m = nvector
          F(:,r) = R; % R is just a complex number here at a given angle?
          [F_tot,strain] = calc_F(Q,z_s,lat_par,m,i,startz,stack,f,v_uc);
          F(:,r) = F(:,r) + exp(1i*Q*(stack{i}.dinterface)).*F_tot;
          F(:,r) = F(:,r)*sqrt(B(m));
          tally = tally + sqrt(B(m));
          r = r + 1;
        end
        R = sum(F/tally,2); % Obs! Intensity, not amplitude squared!


      else

        [F_tot,strain] = calc_F(Q,z_s,lat_par,N,i,startz,stack,f,v_uc);
        R = R + exp(1i*Q*(stack{i}.dinterface)).*F_tot;
      end
      
    end
    stack{i}.pre_calc_f = pre_calc_f;
  end



  I = abs(R).^2;
  % optics on 1: incidence side, -1: detector side
  if instrument.theta_mPath == 1
    if control.pol == 0     % sigma polarization (in the scattering plane)
      output.refl = I;
    elseif control.pol == 1 % pi polarization    (out of scattering plane)
      output.refl = I.*cosd(instrument.theta*2).^2;
    elseif control.pol == 2
      P = (1 + cosd(instrument.theta_m*2)^2.*cosd(instrument.theta*2).^2)./(1 + cosd(instrument.theta_m*2).^2);
      output.refl = I.*P;
    end
  elseif instrument.theta_mPath == 2
    output.refl = I;
    % not implemented, detector side
  elseif instrument.theta_mPath == 0
    output.refl = I;
  end
end






% Code for combining atoms at the same depth
%  [z_s,sorted_idx] = sort(z_s);
% f = f(:,sorted_idx);
% f0 = f(:,sorted_idx);
% k = 1;
% idx = [];
% tally = ones(size(z_s));
% for l = 1:length(z_s)
%   curz = z_s(l);
%
%   for ll = l+1:length(z_s)
%     if z_s(ll) == curz
%       f(:,l)  = f(:,l) + f(:,ll);
%       f0(:,l) = f0(:,l) + f0(:,ll);
%       idx(end+1) = l;
%       tally(l) = tally(l) + 1;
%     end
%   end
% end
% z_s = unique(z_s);
% if isempty(idx)
%   idx = 1:length(z_s);
%   tally = 1;
% end
% f = f(:,idx);
% f0 = f0(:,idx);
%
% pos_vector = [];
% ll = 1;
% for l = 1:N
%   for s = 1:length(z_s)
%     pos_vector(ll) = z_s(s) + lat_par*(l-1);
%     ff(:,ll)       = f(:,s);
%     ff0(:,ll)      = f0(:,s);
%     ll             = ll + 1;
%   end
% end
%