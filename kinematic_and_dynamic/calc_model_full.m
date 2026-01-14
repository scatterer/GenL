function [I,F_cost,S_struct] = calc_model_full(FVr_temp, S_struct)

  [instrument,stack,control]  = translate_fitting_vector(FVr_temp,S_struct.instrument,S_struct.stack,S_struct.control);


  % PL     = (1+cosd(2*2)^2.*cosd(theta*2).^2)./(1+cosd(2*2).^2); %.* 1./(sind(twotheta));
  % A_film = (1-exp(-2*mu.*tot_thick./sind(theta)));

  % add background
  if instrument.background == 0 % with linear background a*Q+b
    bkg = instrument.a*instrument.Q+instrument.b;
  elseif instrument.background == 1  % with polynomial background c*(a+Q).^2+b
    bkg = instrument.c*(instrument.a+instrument.Q).^2+instrument.b;
  end

  if strcmp(S_struct.control.model,'density') 
    output = parratt_matrix_repeated_rhobuiltin(instrument.Q,instrument.wavelength,stack,control,instrument);
  elseif strcmp(S_struct.control.model,'kinematic')
    [output,S_struct.stack] = calc_kinematic_full(instrument.Q,instrument.wavelength,stack,control,instrument);
  elseif strcmp(S_struct.control.model,'darwin')
    % TODO
  end
  if instrument.resolution == 0
    output.refl = instrument.I0.*output.refl + bkg;
  else
    output.refl = instrument.I0.*GaussConv(instrument.Q,output.refl,instrument.resolution) + bkg;
  end

  I = output.refl;
 
  
  F_cost = sum( abs( log10(output.refl(:)) - log10(S_struct.I) ))./(S_struct.N-1);
  %F_cost = 0;
end