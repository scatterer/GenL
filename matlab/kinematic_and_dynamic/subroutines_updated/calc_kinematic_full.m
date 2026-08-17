%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:           Gunnar Palsson / roughness corrections by OpenAI
% Timestamp:        2026-08-14 15:22 +02:00
% Description:      calculate kinematic approximation
%
% Crystal-thickness roughness:
%   control.roughness_averaging = 'coherent' or 'incoherent'
%
% Kinematic refraction mode:
%   control.kinematic_refraction = 'vanilla'  (default)
%   control.kinematic_refraction = 'full'
%
% Kinematic absorption:
%   control.kinematic_absorption = 'full'      (default)
%   control.kinematic_absorption = 'none'
%
% In 'full' absorption mode, each material gets its own linear intensity
% attenuation coefficient from the imaginary Q=0 unit-cell form factor:
%
%   mu = 2*r_e*lambda*abs(Im(F0))/V_uc        [1/Angstrom]
%
% Absorption is applied at the COMPLEX AMPLITUDE level.  For symmetric
% reflection, a contribution at depth d below the current surface is
% multiplied by exp(-mu*d/sin(theta)).  Lower layers are attenuated by the
% full thickness of every layer above them before amplitudes are combined.
%
% 'full' replaces the old high-Q g0 ~ 1/Q REFRACTION phase correction
% by its exact square-root counterpart for the substrate cell-to-cell phase.
% Only the real forward-scattering density is used in this refraction shift;
% absorption is handled separately at amplitude level as described above. The
% intra-unit-cell structure factor remains evaluated at the external Q, just
% as in the original g0 implementation.  No Fresnel multiple reflections are
% introduced.
%
% calc_F is included as a local function at the end of this file.
% Its old >5000-point remainder bug is fixed here.
%
% Current limitation (intentional): crystal-thickness roughness is
% supported only for the TOPMOST layer.  This matches the current
% dynamical implementation and avoids silently averaging an intermediate
% layer before the interference with layers above it is known.
%
% To do as user:    Nothing.
%
% Note:
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output,stack] = calc_kinematic_full(Q,lambda,stack,control,instrument)

  % The algebra in this routine assumes Q is a column vector.  Force that
  % convention explicitly; the old code could create an M-by-M array by
  % implicit expansion when a row-vector Q was supplied.
  Q = Q(:);

  if isempty(Q)
    output.refl = zeros(size(Q));
    return;
  end

  if any(Q == 0)
    error(['calc_kinematic_full is undefined at Q = 0 because the ' ...
           'kinematic amplitude contains a 1/Q factor.']);
  end

  % Parse control modes once, not inside any hot loop.
  [coherent_roughness,~] = get_roughness_mode(control);
  [full_refraction,~] = get_kinematic_refraction_mode(control);
  [full_absorption,~] = get_kinematic_absorption_mode(control);

  % Symmetric geometry: Q = 4*pi*sin(theta)/lambda.  Use Q itself rather
  % than a separately stored theta vector so attenuation is exactly aligned
  % with every calculated Q point.
  sin_theta = abs(Q) * lambda / (4*pi);
  sin_theta = min(sin_theta,1);
  sin_theta(sin_theta < realmin('double')) = realmin('double');

  % The current thickness-distribution model is only correct for a rough
  % top layer.  Averaging an intermediate layer and then adding a layer on
  % top would destroy the realization-by-realization interference.
  nstack = numel(stack);
  for ii = 1:max(0,nstack-1)
    if isfield(stack{ii},'roughness') && ~isempty(stack{ii}.roughness) && stack{ii}.roughness
      error(['Crystal-thickness roughness is currently supported only ' ...
             'for the topmost layer. Layer %d of %d is marked rough.'], ...
             ii,nstack);
    end
  end

  % Assuming the z axis of the unit cell is aligned with Q
  Z = (1:118);
  Zt = {'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr',...
    'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',...
    'In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',...
    'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra',...
    'Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'};

  r0 = 2.814042735053330e-05;
  strain = cell(length(stack),1);
  R = complex(zeros(size(Q)));

  % Only used when the final/top layer is averaged incoherently.  In that
  % case there is no single complex amplitude representing the ensemble.
  have_incoherent_intensity = false;
  I_incoherent = [];

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
        % FIX: the old code used pre_calc_f.a2(1), which is not the a2
        % lattice vector and can also misbehave for a struct array.
        scaling = pre_calc_f(1).a2*stack{i}.scale;
        area    = abs(norm(cross(pre_calc_f(1).a1*as,pre_calc_f(1).a3*as)));
        v_uc    = abs(dot(scaling,cross(pre_calc_f(1).a1*as,pre_calc_f(1).a3*as)));
      case 3
        scaling = pre_calc_f(1).a3*stack{i}.scale;
        area    = abs(norm(cross(pre_calc_f(1).a1*as,pre_calc_f(1).a2*as)));
        v_uc    = abs(dot(scaling,cross(pre_calc_f(1).a1*as,pre_calc_f(1).a2*as)));
      otherwise
        error('stack{%d}.direction must be 1, 2, or 3.',i);
    end

    lat_par(i) = norm(scaling); %#ok<AGROW>

    z_s = zeros(size(pre_calc_f(1).r,1),1);
    for l = 1:size(pre_calc_f(1).r,1)
      z_s(l) = pre_calc_f(1).r(l,:)*scaling';
    end

    % ------------------------------------------------------------------
    % Safe kinematic form-factor cache.
    %
    % Q is fixed during one fit/simulation, so avoid storing/comparing the
    % complete Q vector on every model evaluation.  A short signature still
    % invalidates the cache if the same stack object is reused for a different
    % dataset.  Use a separate f_kin field so this routine never overwrites
    % the Q0-dependent .f cache used by the dynamical density calculation.
    % ------------------------------------------------------------------
    q_signature = [numel(Q), Q(1), Q(end), sum(Q), sum(Q.^2)];
    cache_matches_Q = false;
    if ~calculate_from_scratch && ...
       isfield(pre_calc_f,'kinematic_Q_signature') && ...
       isfield(pre_calc_f,'kinematic_lambda') && ...
       ~isempty(pre_calc_f(1).kinematic_Q_signature) && ...
       ~isempty(pre_calc_f(1).kinematic_lambda)
      cache_matches_Q = isequal(pre_calc_f(1).kinematic_Q_signature,q_signature) && ...
                        isequal(pre_calc_f(1).kinematic_lambda,lambda);
    end

    k = 0;
    l = 0;

    f = complex(zeros(length(Q),length(z_s)));

    % Q=0 unit-cell forward-scattering sum.  It is needed for
    % layer-specific absorption and, for the substrate, for refraction.
    % Keep it scalar; the old code expanded it to the size of Q.
    layer_F0 = 0;
    need_forward_scattering = full_absorption || (i == 1 && full_refraction);

    for j = 1:size(pre_calc_f(1).type,1)
      if ~isempty(pre_calc_f(1).type{j})
        k = k + 1;
        cur_Z = Z(strcmp(Zt,pre_calc_f(1).type{j}));
        if isempty(cur_Z)
          error('Unknown element symbol ''%s''.',pre_calc_f(1).type{j});
        end

        % Cache the wavelength-dependent form-factor coefficients and the
        % Q-independent Debye-Waller prefactor separately.
        need_coeff = calculate_from_scratch || ...
                     numel(pre_calc_f) < k || ...
                     ~isfield(pre_calc_f,'kin_coeff') || ...
                     isempty(pre_calc_f(k).kin_coeff) || ...
                     ~isfield(pre_calc_f,'kin_coeff_lambda') || ...
                     isempty(pre_calc_f(k).kin_coeff_lambda) || ...
                     ~isequal(pre_calc_f(k).kin_coeff_lambda,lambda);

        if need_coeff
          pre_calc_f(k).kin_coeff = Read_form_factor_coefficients(cur_Z,lambda);
          pre_calc_f(k).kin_coeff_lambda = lambda;
        end

        if calculate_from_scratch || ...
           numel(pre_calc_f) < k || ...
           ~isfield(pre_calc_f,'kin_DW') || ...
           isempty(pre_calc_f(k).kin_DW)
          pre_calc_f(k).kin_DW = DB_prefactor(cur_Z);
        end

        need_f = need_coeff || ~cache_matches_Q || ...
                 numel(pre_calc_f) < k || ...
                 ~isfield(pre_calc_f,'f_kin') || ...
                 isempty(pre_calc_f(k).f_kin) || ...
                 numel(pre_calc_f(k).f_kin) ~= numel(Q);

        if need_f
          [g,~,~] = Form_factors(Q,pre_calc_f(k).kin_coeff,100/100);
          pre_calc_f(k).f_kin = g .* ...
              exp(-pre_calc_f(k).kin_DW*(Q/4/pi).^2) ./ Q;
        end

        % Q=0 forward scattering gives both the real refraction density and
        % the imaginary absorption density.  Cache it by wavelength.
        if need_forward_scattering
          need_f0 = need_coeff || ...
                    numel(pre_calc_f) < k || ...
                    ~isfield(pre_calc_f,'f0_kin') || ...
                    isempty(pre_calc_f(k).f0_kin);
          if need_f0
            [g_forward,~,~] = Form_factors(0,pre_calc_f(k).kin_coeff,100/100);
            pre_calc_f(k).f0_kin = g_forward;
          end
        end

        for m = 1:pre_calc_f(1).type_nr(k)
          l = l + 1;
          f(:,l) = pre_calc_f(k).f_kin;

          if need_forward_scattering
            layer_F0 = layer_F0 + pre_calc_f(k).f0_kin;
          end
        end
      end
    end

    pre_calc_f(1).kinematic_Q_signature = q_signature;
    pre_calc_f(1).kinematic_lambda = lambda;

    % Layer-specific linear INTENSITY attenuation coefficient.  With
    % n = 1-delta+i*beta and beta = r_e*lambda^2*f2/(2*pi*V),
    % mu = 4*pi*beta/lambda = 2*r_e*lambda*f2/V.  Form-factor libraries
    % differ in the sign convention for Im(f), so use its magnitude to
    % enforce the physical mu >= 0.
    if full_absorption
      mu_layer = 2*r0*lambda*abs(imag(layer_F0))/v_uc;  % 1/Angstrom
    else
      mu_layer = 0;
    end
    alpha_layer = mu_layer ./ sin_theta;  % double-pass AMPLITUDE decay / Angstrom

    % Expose the currently calculated material absorption for diagnostics.
    stack{i}.kinematic_mu = mu_layer;
    if mu_layer > 0
      stack{i}.kinematic_absorption_length = 1/mu_layer;
    else
      stack{i}.kinematic_absorption_length = Inf;
    end

    N = round(stack{i}.N);
    if N < 1
      error('stack{%d}.N must round to at least one unit cell.',i);
    end

    layer_is_rough = isfield(stack{i},'roughness') && ...
                     ~isempty(stack{i}.roughness) && stack{i}.roughness;

    if i == 1
      % ---------------------------------------------------------------
      % SUBSTRATE
      % ---------------------------------------------------------------
      % Re-anchor the coordinate system at the TOP substrate unit cell.
      % This is only a global phase transformation in vanilla mode, so it
      % leaves the calculated intensity unchanged, but it lets stack{1}.N
      % represent a physically thick substrate without creating enormous
      % absolute z coordinates for every film above it.
      last_atom_z = max(z_s);

      % Refraction modifies only the substrate CELL-TO-CELL phase, exactly
      % like the old g0 correction intended.  The unit-cell structure factor
      % itself remains evaluated at the measured external Q.
      if full_refraction
        rho_e0_real = real(layer_F0) / v_uc;
        Qc2 = 16*pi*r0*rho_e0_real;

        % The substrate extends toward decreasing z in the surface-anchored
        % representation.  Choose the branch with Im(Qphase)<=0 below the
        % critical condition so exp(-i*Qphase*d) decays with depth d>0.
        Qphase = conj(sqrt(complex(Q.^2 - Qc2,0)));
      else
        Qphase = Q;
      end

      substrate_top_z = max(z_s);
      sub_F = complex(zeros(size(Q)));
      for atom_idx = 1:numel(z_s)
        atom_term = exp(1i*Q.*z_s(atom_idx));
        if full_absorption
          atom_depth = max(0,substrate_top_z - z_s(atom_idx));
          atom_term = atom_term .* exp(-alpha_layer .* atom_depth);
        end
        sub_F = sub_F + f(:,atom_idx).*atom_term;
      end

      g = 4*pi.*sub_F./v_uc.*lat_par(i).*r0;
      substrate_phi = Qphase*lat_par(i);
      gamma_cell = alpha_layer*lat_par(i);

      if layer_is_rough
        sigma_N = stack{i}.sigma;
        [nvector,weights] = thickness_distribution(N,sigma_N);

        if coherent_roughness
          Ravg = complex(zeros(size(Q)));
          for kk = 1:numel(nvector)
            NN = nvector(kk);
            Ssub = damped_lattice_sum_down(substrate_phi,NN,gamma_cell);
            Fm = -1i*g .* Ssub;
            Ravg = Ravg + weights(kk).*Fm;
          end
          R = Ravg;
        else
          Iavg = zeros(size(Q));
          for kk = 1:numel(nvector)
            NN = nvector(kk);
            Ssub = damped_lattice_sum_down(substrate_phi,NN,gamma_cell);
            Fm = -1i*g .* Ssub;
            Iavg = Iavg + weights(kk).*abs(Fm).^2;
          end
          I_incoherent = Iavg;
          have_incoherent_intensity = true;
        end

      else
        Ssub = damped_lattice_sum_down(substrate_phi,N,gamma_cell);
        R = R - 1i*g .* Ssub;
      end

    else

      if i > 2
        last_atom_z = strain{i-1}(end);
      end

      % generate_strain receives startz and returns absolute atomic z
      % positions.  calc_F therefore already contains the dinterface phase.
      % Vacuum/interface gaps themselves do not absorb.
      startz = last_atom_z + stack{i}.dinterface;

      if layer_is_rough
        sigma_N = stack{i}.sigma;
        [nvector,weights] = thickness_distribution(N,sigma_N);

        % Top-layer thickness ensemble.  Absorption changes BOTH pieces:
        %   1) each realization has a self-absorbed layer amplitude F_m;
        %   2) the complete lower-stack amplitude is attenuated by the full
        %      thickness of that realization before it interferes with F_m.
        if coherent_roughness
          Ravg = complex(zeros(size(Q)));
          for kk = 1:numel(nvector)
            NN = nvector(kk);
            [F_m,~,thickness_m] = calc_F(Q,z_s,lat_par,NN,i,startz, ...
                stack,f,v_uc,alpha_layer,full_absorption);

            if full_absorption
              Tm = exp(-alpha_layer .* thickness_m);
              Rm = Tm .* R + F_m;
            else
              Rm = R + F_m;
            end
            Ravg = Ravg + weights(kk).*Rm;
          end
          R = Ravg;

        else
          Iavg = zeros(size(Q));
          for kk = 1:numel(nvector)
            NN = nvector(kk);
            [F_m,~,thickness_m] = calc_F(Q,z_s,lat_par,NN,i,startz, ...
                stack,f,v_uc,alpha_layer,full_absorption);

            if full_absorption
              Tm = exp(-alpha_layer .* thickness_m);
              Rm = Tm .* R + F_m;
            else
              Rm = R + F_m;
            end
            Iavg = Iavg + weights(kk).*abs(Rm).^2;
          end

          I_incoherent = Iavg;
          have_incoherent_intensity = true;
        end

      else
        [F_tot,strain,layer_thickness] = calc_F(Q,z_s,lat_par,N,i,startz, ...
            stack,f,v_uc,alpha_layer,full_absorption);

        if full_absorption
          Tlayer = exp(-alpha_layer .* layer_thickness);
          R = Tlayer .* R + F_tot;
        else
          R = R + F_tot;
        end
      end
    end

    stack{i}.pre_calc_f = pre_calc_f;
  end

  if have_incoherent_intensity
    % Only possible for the topmost layer (enforced above).
    I = I_incoherent;
  else
    I = abs(R).^2;
  end

  % optics on 1: incidence side, 2: detector side
  if instrument.theta_mPath == 1
    if control.pol == 0
      % sigma polarization: E perpendicular to the scattering plane.
      output.refl = I;
    elseif control.pol == 1
      % pi polarization: E parallel to the scattering plane.
      output.refl = I.*cosd(instrument.theta*2).^2;
    elseif control.pol == 2
      P = (1 + cosd(instrument.theta_m*2)^2.*cosd(instrument.theta*2).^2) ./ ...
          (1 + cosd(instrument.theta_m*2).^2);
      output.refl = I.*P;
    else
      error('control.pol must be 0 (sigma), 1 (pi), or 2 (mixed).');
    end
  elseif instrument.theta_mPath == 2
    output.refl = I;
    % not implemented, detector side
  elseif instrument.theta_mPath == 0
    output.refl = I;
  else
    error('instrument.theta_mPath must be 0, 1, or 2.');
  end

end


function [F_tot,strain,layer_thickness] = calc_F( ...
    Q,z_s,lat_par,N,i,startz,stack,f,v_uc,alpha_layer,full_absorption)
%CALC_F Kinematic amplitude of one finite layer at its absolute z positions.
%
% alpha_layer(Q) = mu_layer/sin(theta), where mu_layer is the ordinary
% linear INTENSITY attenuation coefficient.  A scattering contribution at
% depth d below the layer surface therefore receives the symmetric-geometry
% double-pass AMPLITUDE factor exp(-alpha_layer*d).
%
% generate_strain uses startz when constructing pos_vector.  Consequently
% exp(i*Q*pos_vector) already includes stack{i}.dinterface through startz.
% The caller must not apply a second exp(i*Q*dinterface) phase factor.

  r0 = 2.814042735053330e-05;
  [pos_vector,~,f_index] = generate_strain(z_s,lat_par,N,i,startz,stack);

  strain = cell(length(stack),1);
  strain{i} = pos_vector;

  m = numel(Q);
  F_tot = complex(zeros(size(Q)));
  layer_thickness = 0;
  if m == 0 || isempty(pos_vector)
    return;
  end

  % The parent geometry defines the next interface from the last/topmost
  % atomic position, so use the same surface definition for absorption.
  top_z = max(pos_vector);
  layer_thickness = max(0,top_z - startz);
  atom_depth = max(0,top_z - pos_vector(:));

  pre_factor = -1i*4*pi*r0*lat_par(i)/v_uc;
  f_scaled = f .* pre_factor;
  iQ = 1i*Q;

  chunk_size = min(5000,m);
  for first = 1:chunk_size:m
    last = min(first + chunk_size - 1,m);
    idx = first:last;

    phase = exp(iQ(idx) * pos_vector.');
    if full_absorption
      attenuation = exp(-alpha_layer(idx) * atom_depth.');
      phase = phase .* attenuation;
    end

    F_tot(idx) = sum(f_scaled(idx,f_index).*phase,2);
  end
end


function [full_refraction,mode] = get_kinematic_refraction_mode(control)
%GET_KINEMATIC_REFRACTION_MODE Parse control.kinematic_refraction once.
%
%   'vanilla' : ordinary kinematic phase, Qphase = Q.  This is the default.
%   'full'    : exact square-root substrate refraction phase.  It replaces
%               the old first-order g0 ~ 1/Q correction, using real F(0)
%               only and leaving the intra-unit-cell structure factor at Q.
%
% This option changes the kinematic phase only; it does not turn the model
% into a Parratt/DWBA calculation and does not add Fresnel multiple
% reflections.

  if isfield(control,'kinematic_refraction') && ...
          ~isempty(control.kinematic_refraction)
    mode = control.kinematic_refraction;
  else
    mode = 'vanilla';
  end

  if isstring(mode)
    if ~isscalar(mode)
      error('control.kinematic_refraction must be a scalar string or char vector.');
    end
    mode = char(mode);
  end

  if ~ischar(mode)
    error('control.kinematic_refraction must be ''vanilla'' or ''full''.');
  end

  mode = lower(strtrim(mode));
  if strcmp(mode,'vanilla')
    full_refraction = false;
  elseif strcmp(mode,'full')
    full_refraction = true;
  else
    error('control.kinematic_refraction must be ''vanilla'' or ''full''.');
  end
end


function [full_absorption,mode] = get_kinematic_absorption_mode(control)
%GET_KINEMATIC_ABSORPTION_MODE Parse control.kinematic_absorption once.
%
%   'full' : layer-specific absorption from Im(F0), applied at amplitude
%            level to self-scattering and to all lower layers.  Default.
%   'none' : reproduce the undamped kinematic model.

  if isfield(control,'kinematic_absorption') && ...
          ~isempty(control.kinematic_absorption)
    mode = control.kinematic_absorption;
  else
    mode = 'full';
  end

  if isstring(mode)
    if ~isscalar(mode)
      error('control.kinematic_absorption must be a scalar string or char vector.');
    end
    mode = char(mode);
  end

  if ~ischar(mode)
    error('control.kinematic_absorption must be ''full'' or ''none''.');
  end

  mode = lower(strtrim(mode));
  if strcmp(mode,'full')
    full_absorption = true;
  elseif strcmp(mode,'none')
    full_absorption = false;
  else
    error('control.kinematic_absorption must be ''full'' or ''none''.');
  end
end


function [coherent_mode,mode] = get_roughness_mode(control)
%GET_ROUGHNESS_MODE Read and validate control.roughness_averaging once.

  if isfield(control,'roughness_averaging') && ...
          ~isempty(control.roughness_averaging)
    mode = control.roughness_averaging;
  else
    % Backward-compatible default.  Set the field explicitly when testing
    % coherent versus incoherent thickness averaging.
    mode = 'coherent';
  end

  if isstring(mode)
    if ~isscalar(mode)
      error('control.roughness_averaging must be a scalar string or char vector.');
    end
    mode = char(mode);
  end

  if ~ischar(mode)
    error('control.roughness_averaging must be ''coherent'' or ''incoherent''.');
  end

  mode = lower(strtrim(mode));
  if strcmp(mode,'coherent')
    coherent_mode = true;
  elseif strcmp(mode,'incoherent')
    coherent_mode = false;
  else
    error('control.roughness_averaging must be ''coherent'' or ''incoherent''.');
  end
end


function [nvector,weights] = thickness_distribution(N,sigma_N)
%THICKNESS_DISTRIBUTION Discrete, truncated Gaussian in unit-cell count.
%
% sigma_N is in UNIT CELLS, not Angstrom and not interface microroughness.

  if ~isscalar(sigma_N) || ~isfinite(sigma_N) || sigma_N < 0
    error('Crystal-thickness sigma must be a finite nonnegative scalar.');
  end

  if sigma_N == 0
    nvector = N;
    weights = 1;
    return;
  end

  nmin = max(1,round(N - 3*sigma_N));
  nmax = max(1,round(N + 3*sigma_N));
  nvector = nmin:nmax;

  % The continuous Gaussian prefactor cancels when the discrete truncated
  % distribution is normalized, so do not compute it.
  weights = exp(-0.5*((nvector - N)/sigma_N).^2);
  weights = weights ./ sum(weights);
end


function S = damped_lattice_sum_down(phi,N,gamma)
%DAMPED_LATTICE_SUM_DOWN Surface-anchored substrate lattice sum.
%
% The top substrate unit cell is d=0 and deeper cells are d=1,2,... .
% Their relative amplitude ratio is
%
%   r = exp(-gamma - i*phi),
%
% where gamma = mu*a/sin(theta) is the symmetric double-pass AMPLITUDE
% attenuation per unit cell.  Thus
%
%   S = 1 + r + ... + r^(N-1) = (1-r^N)/(1-r).
%
% Using expm1 keeps the expression stable both at Bragg conditions and for
% very large physical substrate N.  When absorption is appreciable, r^N
% naturally vanishes and the finite substrate becomes effectively
% semi-infinite without any ad-hoc smoothing.

  % Reduce the REAL phase modulo 2*pi before using expm1.  This makes the
  % expression accurate close to any Bragg order, not only near phi=0.
  phi_reduced = mod(real(phi) + pi,2*pi) - pi + 1i*imag(phi);
  log_r = -gamma - 1i*phi_reduced;

  S = complex(zeros(size(phi)));
  exactly_one = (log_r == 0);
  other = ~exactly_one;

  if any(other)
    S(other) = expm1(N .* log_r(other)) ./ expm1(log_r(other));
  end

  if any(exactly_one)
    S(exactly_one) = N;
  end
end
