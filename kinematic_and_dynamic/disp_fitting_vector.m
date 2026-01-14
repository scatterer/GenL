function disp_fitting_vector(GKP_bestmem,stack,instrument)
 

  k = 1;
  for i = 1:length(stack)
    disp([' stack ',num2str(i),' N                        : ',num2str(GKP_bestmem(k))]); k = k + 1;
    disp([' stack ',num2str(i),' dinterface               : ',num2str(GKP_bestmem(k))]); k = k + 1;
    disp([' stack ',num2str(i),' scale                    : ',num2str(GKP_bestmem(k))]); k = k + 1;
    disp([' stack ',num2str(i),' area_scale               : ',num2str(GKP_bestmem(k))]); k = k + 1;
    disp([' stack ',num2str(i),' sigma                    : ',num2str(GKP_bestmem(k))]); k = k + 1;
    disp([' stack ',num2str(i),' bottom_strain_amplitude  : ',num2str(GKP_bestmem(k))]); k = k + 1;
    disp([' stack ',num2str(i),' bottom_strain_end        : ',num2str(GKP_bestmem(k))]); k = k + 1;
    disp([' stack ',num2str(i),' top_strain_amplitude     : ',num2str(GKP_bestmem(k))]); k = k + 1;
    disp([' stack ',num2str(i),' top_strain_end           : ',num2str(GKP_bestmem(k))]); k = k + 1;
  end

  disp([' instrument  I0           : ',num2str(GKP_bestmem(k))]); k = k + 1;
  disp([' instrument  resolution   : ',num2str(GKP_bestmem(k))]); k = k + 1;
  disp([' instrument  background a : ',num2str(GKP_bestmem(k))]); k = k + 1;
  disp([' instrument  background b : ',num2str(GKP_bestmem(k))]); k = k + 1;
  disp([' instrument  background c : ',num2str(GKP_bestmem(k))]); k = k + 1;

  % instrument.I0 = GKP_bestmem(k); k = k + 1;
  % instrument.resolution = GKP_bestmem(k); k = k + 1;
  % instrument.a = GKP_bestmem(k); k = k + 1;
  % instrument.b = GKP_bestmem(k); k = k + 1;
  % instrument.c = GKP_bestmem(k); k = k + 1;

end