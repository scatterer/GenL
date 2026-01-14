function [instrument,stack,control]  = translate_fitting_vector(GKP_bestmem,instrument,stack,control)
 

  k = 1;
  for i = 1:length(stack)
    stack{i}.N          = GKP_bestmem(k); k = k + 1;
    stack{i}.dinterface = GKP_bestmem(k); k = k + 1;
    stack{i}.scale      = GKP_bestmem(k); k = k + 1;
    stack{i}.area_scale = GKP_bestmem(k); k = k + 1;
    stack{i}.sigma      = GKP_bestmem(k); k = k + 1;
    stack{i}.bottom_strain_amplitude = GKP_bestmem(k); k = k + 1;
    stack{i}.bottom_strain_end       = GKP_bestmem(k); k = k + 1;
    stack{i}.top_strain_amplitude = GKP_bestmem(k); k = k + 1;
    stack{i}.top_strain_end       = GKP_bestmem(k); k = k + 1;
  end


  instrument.I0 = GKP_bestmem(k); k = k + 1;
  instrument.resolution = GKP_bestmem(k); k = k + 1;
  instrument.a = GKP_bestmem(k); k = k + 1;
  instrument.b = GKP_bestmem(k); k = k + 1;
  instrument.c = GKP_bestmem(k); k = k + 1;

end