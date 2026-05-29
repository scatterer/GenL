function refl = propagate_vectorized_chunks_roughness(Q,lambda,rho_0r,rho_1r,N,rho_e_rough,dz,substrate_end,nvector,BB,pol)
    
  F_film = zeros(length(nvector),length(Q));
  
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
    tally = 0;

    r = 1;
    for n = nvector
      rho_rest  = rho_e_rough{n}(substrate_end:end);
      rho_restr = rho_rest(end:-1:1);
      [A1,A2,A3,A4] = make_A_matrix_vec(currentQ,lambda,rho_restr,dz,0,pol);
      B3 = do_matrix_propagation_optimzed_vec_winograd(A1,A2,A3,A4);

      % for k = 1:mm
      %   B2(:,:,k) = B3(:,:,k)*B2(:,:,k);
      % end

      Btemp = B2;
      for k = 1:mm
        Btemp(:,:,k) = B3(:,:,k)*Btemp(:,:,k);
      end
      % looks like we are making a whole crystal from substrate to surface and weighing that with sqrt(BB(b)).
      % The other way is to only make the film, create the weighted version. then attached the weighted version to the substrate.
       F_film(r,idx) = F_film(r,idx) + reshape(Btemp(2,1,:)./Btemp(1,1,:),[1,mm]);
       F_film(r,idx) = F_film(r,idx)*sqrt(BB(n));
             tally = tally + sqrt(BB(n));

       r = r + 1;
    end
    %refl(idx) = abs(reshape(B2(2,1,:)./B2(1,1,:),[1,mm])).^2;

   %refl(idx) = sum(abs(F_film(:,idx)).^2,1); % Obs! Intensity, not amplitude squared!
    refl(idx) = abs(sum(F_film(:,idx),1)/tally).^2;  % Sum over roughness then square
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

    r = 1;
    tally = 0;
    for n = nvector
      rho_rest  = rho_e_rough{n}(substrate_end:end);
      rho_restr = rho_rest(end:-1:1);

      [A1,A2,A3,A4] = make_A_matrix_vec(currentQ,lambda,rho_restr,dz,0,pol);

      B3 = do_matrix_propagation_optimzed_vec_winograd(A1,A2,A3,A4);
      Btemp = B2;
      for k = 1:mm
        Btemp(:,:,k) = B3(:,:,k)*Btemp(:,:,k);
      end

      F_film(r,idx) = F_film(r,idx) + reshape(B2(2,1,:)./B2(1,1,:),[1,mm]);
      F_film(r,idx) = F_film(r,idx)*sqrt(BB(n));
      tally = tally + sqrt(BB(n));
      r = r + 1;
    end
    %refl(idx) = sum(abs(F_film).^2,1); % Obs! Intensity, not amplitude squared!
    refl(idx) = abs(sum(F_film(:,idx),1)/tally).^2;  % Sum over roughness then square

  end


end

function [A1,A2,A3,A4] = make_A_matrix(Q,lambda,rho_e,dz)
    % Constants
    re = 2.814042735053330e-05;
    k0 = (2*pi/lambda);
    factor      = 8 * k0^2 * lambda^2 / (2 * pi);
    
    % Precompute layer properties
    sld   = rho_e * re / (2 * pi ); 
    delta = factor * real(sld);
    beta  = 1i*factor * imag(sld);
    %A = zeros(2,2,length(rho_e),length(Q));
    n = length(rho_e);
    m = length(Q);

    A1 = zeros(n,m);
    A2 = A1;
    A3 = A1;
    A4 = A1;
    % Compute kz and kz of adjacent layers
    for i = 1:length(Q)
      kz_all     = sqrt(Q(i).^2 - delta + beta);  % Vectorized kz
      kz_shifted = circshift(kz_all, 1);     % Shifted kz for previous layer


      denom = 1./(kz_shifted + kz_all);
      r     = (kz_shifted-kz_all).*denom;
      t     = 2*kz_shifted.*denom;

      exp_pos  = exp(1i * kz_all * dz / 2);

      invtp     = exp_pos./t;
      invtm     = 1./(t.*exp_pos);

      A1(:,i) = invtp;
      A2(:,i) = invtm.*r;
      A3(:,i) = invtp.*r;
      A4(:,i) = invtm;

    end
end

% function [A1,A2,A3,A4] = make_A_matrix_vec(Q,lambda,rho_e,dz,sigma)
%     % Constants
%     re = 2.814042735053330e-05;
%     k0 = (2*pi/lambda);
%     factor      = 8 * k0^2 * lambda^2 / (2 * pi);
% 
%     % Precompute layer properties
%     sld   = rho_e' * re / (2 * pi ); 
%     delta = factor * real(sld);
%     beta  = 1i*factor * imag(sld);
% 
%     n = length(rho_e);
%     m = length(Q);
% 
%     Qsqr = (Q(:)').^2;
% 
%     Qsqr_rep   = repmat(Qsqr,[n,1]);
%     delta_rep  = repmat(delta,[1,m]);
%     beta_rep   = repmat(beta,[1,m]);
%     kz_all     = (sqrt(Qsqr_rep - delta_rep + beta_rep)').';
% 
%     kz_shifted = circshift(kz_all,1,1);
%     denom      = 1./(kz_shifted + kz_all);
%     r          = (kz_shifted - kz_all).*denom;
%     t          = 2*kz_shifted.*denom;
%     exp_pos    = exp(kz_all * (1i*dz / 2));
% 
%     invtp      = exp_pos./t;
%     invtm      = 1./(t.*exp_pos);
% 
%     ehp = exp( ( kz_all + kz_shifted).^2*(-sigma^2/2));
%     ehm = exp( ( kz_all - kz_shifted).^2*(-sigma^2/2));
% 
%     % A1 = invtp;
%     % A2 = invtm .* r; 
%     % A3 = invtp .* r; 
%     % A4 = invtm;
% 
%     A1 = ehm.*invtp;
%     A2 = r.*ehp.*invtm;
%     A3 = r.*ehp.*invtp;
%     A4 = ehm.*invtm;
% end

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

function B = do_matrix_propagation_optimzed_vec(A) %,z)

  [m,n,o,l] = size(A);
  
   B         = zeros(2,2,l);
   B(1,1,:)  = 1;
   B(2,2,:)  = 1;
   A_sum_d   = reshape(A(1,1,:,:) + A(2,2,:,:),[o,l]);   % (a + d)
   C_sum_d   = reshape(A(2,1,:,:) + A(2,2,:,:),[o,l]);   % (c + d)
   A_sum_b   = reshape(A(1,1,:,:) + A(1,2,:,:),[o,l]);   % (a + b)
   B_minus_d = reshape(A(1,2,:,:) - A(2,2,:,:),[o,l]);   % (b - d)
   C_minus_a = reshape(A(2,1,:,:) - A(1,1,:,:),[o,l]);   % (c - a)
   

   for j = o:-1:2
    B_sum_h   = reshape(B(1,1,:) + B(2,2,:),[1,l]);   % (e + h)
    F_minus_h = reshape(B(1,2,:) - B(2,2,:),[1,l]); % (f - h)
    G_minus_e = reshape(B(2,1,:) - B(1,1,:),[1,l]); % (g - e)
    B_sum_f   = reshape(B(1,1,:) + B(1,2,:),[1,l]);   % (e + f)
    B_sum_g   = reshape(B(2,1,:) + B(2,2,:),[1,l]);   % (g + h)
    
    M1 = A_sum_d(j,:) .* B_sum_h;      % M1 = (a + d)(e + h)
    M2 = C_sum_d(j,:) .* reshape(B(1,1,:),[1,l]);       % M2 = (c + d)e
    M3 = reshape(A(1,1,j,:),[1,l]) .* F_minus_h;     % M3 = a(f - h)
    M4 = reshape(A(2,2,j,:),[1,l]) .* G_minus_e;     % M4 = d(g - e)
    M5 = reshape(A_sum_b(j,:),[1,l]) .* reshape(B(2,2,:),[1,l]);       % M5 = (a + b)h
    M6 = C_minus_a(j,:) .* B_sum_f;    % M6 = (c - a)(e + f)
    M7 = B_minus_d(j,:) .* B_sum_g;    % M7 = (b - d)(g + h)
    
    B(1,1,:) = M1 + M4 - M5 + M7;
    B(1,2,:) = M3 + M5;
    B(2,1,:) = M2 + M4;
    B(2,2,:) = M1 - M2 + M3 + M6;
   end

end

function B = do_matrix_propagation_optimzed_vec2(A1,A2,A3,A4) %,z)

  [o,l] = size(A1);
  
   
   A_sum_d   = A1 + A4;   % (a + d)
   C_sum_d   = A3 + A4;   % (c + d)
   A_sum_b   = A1 + A2;   % (a + b)
   B_minus_d = A2 - A4;   % (b - d)
   C_minus_a = A3 - A1;   % (c - a)

   B1 = ones(1,l);
   B2 = zeros(1,l);
   B3 = B2;
   B4 = B1;

   for j = o:-1:2
    B_sum_h   = B1 + B4; % (e + h)
    F_minus_h = B2 - B4; % (f - h)
    G_minus_e = B3 - B1; % (g - e)
    B_sum_f   = B1 + B2; % (e + f)
    B_sum_g   = B3 + B4; % (g + h)
    
    M1 = A_sum_d(j,:) .* B_sum_h;   
    M2 = C_sum_d(j,:) .* B1;  
    M3 = A1(j,:)      .* F_minus_h;     
    M4 = A4(j,:)      .* G_minus_e;   
    M5 = A_sum_b(j,:) .* B4;     
    M6 = C_minus_a(j,:) .* B_sum_f;   
    M7 = B_minus_d(j,:) .* B_sum_g;    
    
    B1 = M1 + M4 - M5 + M7;
    B2 = M3 + M5;
    B3 = M2 + M4;
    B4 = M1 - M2 + M3 + M6;
   end

   B(1,1,:) = B1;
   B(1,2,:) = B2;
   B(2,1,:) = B3;
   B(2,2,:) = B4;

end

function B = do_matrix_propagation_optimzed_vec_winograd(A1,A2,A3,A4) %,z)

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

   B(1,1,:) = B1;
   B(1,2,:) = B2;
   B(2,1,:) = B3;
   B(2,2,:) = B4;

end

function B = do_matrix_propagation_simple(A1,A2,A3,A4) %,z)

  [o,l] = size(A1);

  % [ A1 A2   * [ B1 B2
  %   A3 A4 ]     B3 B4 ]

  B = zeros(2,2,l);
  B(1,1,:) = 1;
  B(2,2,:) = 1;

  A(1,1,:,:) = A1;
  A(1,2,:,:) = A2;
  A(2,1,:,:) = A3;
  A(2,2,:,:) = A4;

  for k = 1:l
    for j = o:-1:2
      B(:,:,k) = A(:,:,j,k)*B(:,:,k);
    end
  end
end

function B = karstadt_schwartz_2x2(A1,A2,A3,A4)
  [o,l] = size(A1);

  % [ A1 A2   * [ B1 B2
  %   A3 A4 ]     B3 B4 ]

  B1 = ones(1,l);
  B2 = zeros(1,l);
  B3 = B2;
  B4 = B1;

  a1 = A1 + A4;
  a2 = A3 + A4;
  a3 = A1 + A2;
  a4 = A3 - A1;
  a5 = A2 - A4; 

 
  for j = o:-1:2
    % Define intermediate products using Karstadt-Schwartz method
    P1 = a1(j,:) .* (B1 + B4); 
    P2 = a2(j,:) .* B1;         
    P3 = A1(j,:) .* (B2 - B4);          
    P4 = A4(j,:) .* (B3 - B1);          
    P5 = a3(j,:) .* B4;           
    P6 = a4(j,:) .* (B1 + B2);
    P7 = a5(j,:) .* (B3 + B4);

    % Combine intermediate products to get the resulting matrix
    B1 = P1 + P4 - P5 + P7; % e11 = P1 + P4 - P5 + P7
    B2 = P3 + P5;           % e12 = P3 + P5
    B3 = P2 + P4;           % e21 = P2 + P4
    B4 = P1 + P3 - P2 + P6; % e22 = P1 + P3 - P2 + P6
  end

  B(1,1,:) = B1;
  B(1,2,:) = B2;
  B(2,1,:) = B3;
  B(2,2,:) = B4;


end
