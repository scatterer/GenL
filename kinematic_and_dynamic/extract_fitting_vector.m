function [GKP_bestmem,FVR_minbound,FVR_maxbound]  = extract_fitting_vector(stack,instrument,control)
 

GKP_bestmem  = [];
FVR_minbound = [];
FVR_maxbound = []; 

  for i = 1:length(stack)
    if length(stack{i}.N) == 1
      GKP_bestmem(end+1) = stack{i}.N;
      FVR_minbound(end+1) = stack{i}.N;
      FVR_maxbound(end+1) = stack{i}.N;
    else
      GKP_bestmem(end+1) = stack{i}.N(1);
      FVR_minbound(end+1) = stack{i}.N(2);
      FVR_maxbound(end+1) = stack{i}.N(3);
    end

    if length(stack{i}.dinterface) == 1
      GKP_bestmem(end+1) = stack{i}.dinterface;
      FVR_minbound(end+1) = stack{i}.dinterface;
      FVR_maxbound(end+1) = stack{i}.dinterface;
    else
      GKP_bestmem(end+1) = stack{i}.dinterface(1);
      FVR_minbound(end+1) = stack{i}.dinterface(2);
      FVR_maxbound(end+1) = stack{i}.dinterface(3);
    end

    if length(stack{i}.scale) == 1
      GKP_bestmem(end+1) = stack{i}.scale;
      FVR_minbound(end+1) = stack{i}.scale;
      FVR_maxbound(end+1) = stack{i}.scale;
    else
      GKP_bestmem(end+1) = stack{i}.scale(1);
      FVR_minbound(end+1) = stack{i}.scale(2);
      FVR_maxbound(end+1) = stack{i}.scale(3);
    end

    if length(stack{i}.area_scale) == 1
      GKP_bestmem(end+1) = stack{i}.area_scale;
      FVR_minbound(end+1) = stack{i}.area_scale;
      FVR_maxbound(end+1) = stack{i}.area_scale;
    else
      GKP_bestmem(end+1) = stack{i}.area_scale(1);
      FVR_minbound(end+1) = stack{i}.area_scale(2);
      FVR_maxbound(end+1) = stack{i}.area_scale(3);
    end


     if length(stack{i}.sigma) == 1
      GKP_bestmem(end+1) = stack{i}.sigma;
      FVR_minbound(end+1) = stack{i}.sigma;
      FVR_maxbound(end+1) = stack{i}.sigma;
    else
      GKP_bestmem(end+1) = stack{i}.sigma(1);
      FVR_minbound(end+1) = stack{i}.sigma(2);
      FVR_maxbound(end+1) = stack{i}.sigma(3);
     end

  if length(stack{i}.bottom_strain_amplitude) == 1
    GKP_bestmem(end+1) = stack{i}.bottom_strain_amplitude;
    FVR_minbound(end+1) = stack{i}.bottom_strain_amplitude;
    FVR_maxbound(end+1) = stack{i}.bottom_strain_amplitude;
  else
    GKP_bestmem(end+1) = stack{i}.bottom_strain_amplitude(1);
    FVR_minbound(end+1) = stack{i}.bottom_strain_amplitude(2);
    FVR_maxbound(end+1) = stack{i}.bottom_strain_amplitude(3);
  end
  
  if length(stack{i}.bottom_strain_end) == 1
    GKP_bestmem(end+1) = stack{i}.bottom_strain_end;
    FVR_minbound(end+1) = stack{i}.bottom_strain_end;
    FVR_maxbound(end+1) = stack{i}.bottom_strain_end;
  else
    GKP_bestmem(end+1) = stack{i}.bottom_strain_end(1);
    FVR_minbound(end+1) = stack{i}.bottom_strain_end(2);
    FVR_maxbound(end+1) = stack{i}.bottom_strain_end(3);
  end

  if length(stack{i}.top_strain_amplitude) == 1
    GKP_bestmem(end+1) = stack{i}.top_strain_amplitude;
    FVR_minbound(end+1) = stack{i}.top_strain_amplitude;
    FVR_maxbound(end+1) = stack{i}.top_strain_amplitude;
  else
    GKP_bestmem(end+1) = stack{i}.top_strain_amplitude(1);
    FVR_minbound(end+1) = stack{i}.top_strain_amplitude(2);
    FVR_maxbound(end+1) = stack{i}.top_strain_amplitude(3);
  end

  if length(stack{i}.top_strain_end) == 1
    GKP_bestmem(end+1) = stack{i}.top_strain_end;
    FVR_minbound(end+1) = stack{i}.top_strain_end;
    FVR_maxbound(end+1) = stack{i}.top_strain_end;
  else
    GKP_bestmem(end+1) = stack{i}.top_strain_end(1);
    FVR_minbound(end+1) = stack{i}.top_strain_end(2);
    FVR_maxbound(end+1) = stack{i}.top_strain_end(3);
  end


  end

  if length(instrument.I0) == 1
    GKP_bestmem(end+1) = instrument.I0;
    FVR_minbound(end+1) = instrument.I0;
    FVR_maxbound(end+1) = instrument.I0;
  else
    GKP_bestmem(end+1) = instrument.I0(1);
    FVR_minbound(end+1) = instrument.I0(2);
    FVR_maxbound(end+1) = instrument.I0(3);
  end

  if length(instrument.resolution) == 1
    GKP_bestmem(end+1) = instrument.resolution;
    FVR_minbound(end+1) = instrument.resolution;
    FVR_maxbound(end+1) = instrument.resolution;
  else
    GKP_bestmem(end+1) = instrument.resolution(1);
    FVR_minbound(end+1) = instrument.resolution(2);
    FVR_maxbound(end+1) = instrument.resolution(3);
  end

  if length(instrument.a) == 1
    GKP_bestmem(end+1) = instrument.a;
    FVR_minbound(end+1) = instrument.a;
    FVR_maxbound(end+1) = instrument.a;
  else
    GKP_bestmem(end+1) = instrument.a(1);
    FVR_minbound(end+1) = instrument.a(2);
    FVR_maxbound(end+1) = instrument.a(3);
  end

  if length(instrument.b) == 1
    GKP_bestmem(end+1) = instrument.b;
    FVR_minbound(end+1) = instrument.b;
    FVR_maxbound(end+1) = instrument.b;
  else
    GKP_bestmem(end+1) = instrument.b(1);
    FVR_minbound(end+1) = instrument.b(2);
    FVR_maxbound(end+1) = instrument.b(3);
  end

  if length(instrument.c) == 1
    GKP_bestmem(end+1) = instrument.c;
    FVR_minbound(end+1) = instrument.c;
    FVR_maxbound(end+1) = instrument.c;
  else
    GKP_bestmem(end+1) = instrument.c(1);
    FVR_minbound(end+1) = instrument.c(2);
    FVR_maxbound(end+1) = instrument.c(3);
  end


end