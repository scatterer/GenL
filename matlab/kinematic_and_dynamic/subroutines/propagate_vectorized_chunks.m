%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:           Gunnar Palsson
% Description:      propagate the matrix.
%
% To do as user:    Nothing.
%
% Note:
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refl = propagate_vectorized_chunks(Q,lambda,rho_0r,rho_1r,N,rho_restr,dz,sigma,pol)

  m = length(Q);

  slices = min(5000,m);

  iter = floor(m/slices);
  refl = zeros(size(Q));

  for i = 1:iter
    idx      = slices*(i-1)+1:slices*i;
    mm       = length(idx);
    currentQ = Q(idx);

    [A1,A2,A3,A4]  =  make_A_matrix_vec(currentQ,lambda,rho_1r,dz,0,pol);

    B1 = do_matrix_propagation_optimzed_vec_winograd(A1,A2,A3,A4);

    B2  = zeros(2,2,length(idx));

    for k = 1:length(idx)
      B2(:,:,k)  = fastMatrixExponentiation_single(B1(:,:,k),N-2);
    end

    [A1,A2,A3,A4]   = make_A_matrix_vec(currentQ,lambda,rho_0r,dz,0,pol);

    B0  = do_matrix_propagation_optimzed_vec_winograd(A1,A2,A3,A4);


    
    for k = 1:mm
      B2(:,:,k) = B2(:,:,k)*B0(:,:,k);
    end

    [A1,A2,A3,A4] = make_A_matrix_vec(currentQ,lambda,rho_restr,dz,sigma,pol);

    B3 = do_matrix_propagation_optimzed_vec_winograd(A1,A2,A3,A4);

    for k = 1:mm
      B2(:,:,k) = B3(:,:,k)*B2(:,:,k);
    end

    refl(idx) = abs(reshape(B2(2,1,:)./B2(1,1,:),[1,mm])).^2;
  end

  idx      = slices*i+1:(slices*i + length(Q)-slices*i);
  if ~isempty(idx)
    mm       = length(idx);
    currentQ = Q(idx);

    [A1,A2,A3,A4]  =  make_A_matrix_vec(currentQ,lambda,rho_1r,dz,0,pol);

    B1 = do_matrix_propagation_optimzed_vec_winograd(A1,A2,A3,A4);

    B2  = zeros(2,2,length(idx));

    for k = 1:length(idx)
      B2(:,:,k)  = fastMatrixExponentiation_single(B1(:,:,k),N-2);
    end

    [A1,A2,A3,A4]   = make_A_matrix_vec(currentQ,lambda,rho_0r,dz,0,pol);

    B0  = do_matrix_propagation_optimzed_vec_winograd(A1,A2,A3,A4);
    for k = 1:mm
      B2(:,:,k) = B2(:,:,k)*B0(:,:,k);
    end

    [A1,A2,A3,A4] = make_A_matrix_vec(currentQ,lambda,rho_restr,dz,sigma,pol);

    B3 = do_matrix_propagation_optimzed_vec_winograd(A1,A2,A3,A4);

    for k = 1:mm
      B2(:,:,k) = B3(:,:,k)*B2(:,:,k);
    end

    refl(idx) = abs(reshape(B2(2,1,:)./B2(1,1,:),[1,mm])).^2;
  end


end


function [A1,A2,A3,A4] = make_A_matrix_vec(Q,lambda,rho_e,dz,sigma,pol)
    % Constants
    re = 2.814042735053330e-05;
    k0 = (2*pi/lambda);

    factor      =      8 * k0^2 * lambda^2 / (2 * pi);

    % Precompute layer properties
    sld   = rho_e' * re / (2 * pi ); 
    delta = factor * real(sld);      % This is delta*8*k0^2
    beta  = 1i*factor * imag(sld); % This is beta*8*k0^2

    nn = length(rho_e);
    m  = length(Q);

    Qsqr = (Q(:)').^2;

    Qsqr_rep   = repmat(Qsqr,[nn,1]);

    kz_all     = (sqrt(Qsqr_rep - repmat(delta,[1,m]) + repmat(beta,[1,m]))').';

    clear Qsqr_rep Qsqr;

    kz_shifted = circshift(kz_all,1,1);

    exp_pos    = exp(kz_all * (1i*dz/2 ));


    % Kao's matrix for roughness
    %ehp = exp( ( kz_all + kz_shifted).^2*(-sigma^2/2));
    %ehm = exp( ( kz_all - kz_shifted).^2*(-sigma^2/2));

    % A1_s = ehm.*invtp;
    % A2_s = r.*ehp.*invtm;
    % A3_s = r.*ehp.*invtp;
    % A4_s = ehm.*invtm;

    if pol == 0 % pol=0 for sigma polarization (perpendicular)
      % sigma polarization
      denom      = 1./(kz_shifted + kz_all);

      r          = (kz_shifted - kz_all).*denom;
      t          = 2*kz_shifted.*denom;

      invtp      = exp_pos./t;
      invtm      = 1./(t.*exp_pos);

      A1 = invtp;
      A2 = invtm .* r;
      A3 = invtp .* r;
      A4 = invtm;

    elseif pol == 1 % pol=1 for pi polarization (parallel to the scattering plane)
      % pi polarization
      %deltan = lambda^2*real(sld)/2/pi;    % This is delta
      %betan  = 1i*lambda/(4*pi)*imag(sld); % This is beta

      deltan = lambda^2*real(sld)/2/pi;    % This is delta
      betan  = 1i*lambda^2/(2*pi)*imag(sld); % This is beta

   % beta_factor = 1i * 8 * k0^2 * lambda^2 / (2 * pi);

      deltan_rep  = repmat(deltan,[1,m]);
      betan_rep   = repmat(betan,[1,m]);

      n          = 1 - deltan_rep - betan_rep;
      n_shifted  = circshift(n,1,1);

      fact1 = n_shifted.^2.*kz_all;
      fact2 = n.^2.*kz_shifted;
      denom      = (fact1 + fact2);
      r_pi       = (fact1 - fact2)./denom;
      t_pi       = 2*fact1./denom;

      invtp_pi   = exp_pos./t_pi;
      invtm_pi   = 1./(t_pi.*exp_pos);

      A1 = invtp_pi;
      A2 = r_pi.*invtm_pi;
      A3 = r_pi.*invtp_pi;
      A4 = invtm_pi;
    end
end


function B = do_matrix_propagation_optimzed_vec_winograd(A1,A2,A3,A4) %,z)
  % [length(rho) x length(Q) ]
  [o,l] = size(A1);
  
  % [ A1 A2   * [ B1 B2
  %   A3 A4 ]     B3 B4 ]
  
   B1 = ones(1,l);
   B2 = zeros(1,l);
   B3 = B2;
   B4 = B1;



   a1 = A3 - A1;
   a2 = A3 + A4;
   a3 = a1 + A4; % A3 - A1 + A4
   a4 = A2 - a3; % A1 + A2 - A3 - A4;

   for j = o:-1:2
     t  = A1(j,:).*B1;
     b1 = B2 - B4;
     u  = a1(j,:).*b1;
     v  = a2(j,:).*(B2 - B1);
     w  = t + a3(j,:).*(B1 - b1);

     w1  = w + u;
     B3o = B3;

     B3 = w1 + A4(j,:).*(B3o - B1 + b1);
     B1 = t + A2(j,:).*B3o;

     B2 = w + v + a4(j,:).*B4;
     B4 = w1 + v;

   end

   B = zeros(2, 2, l, 'like', A1);
   B(1,1,:) = B1;
   B(1,2,:) = B2;
   B(2,1,:) = B3;
   B(2,2,:) = B4;

end


