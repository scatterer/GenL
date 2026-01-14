close all;
clear all;

% ----- DATA -----

% state which samples are supposed to be fitted
sample = 'Example_data_10nmFe'; % sample name

% load data data
I_unit  = 1; % 1: cps, 2: counts
time    = 25; % measuring times per point

fileExt = '.txt'; % extension of data files
path    = './examples/'; % path where data is stored
file    = strcat(path, sample, fileExt);
data    = load(file);

% Skip some data points
% nothing to skip for Fe
skip    = 1; % datapoints
endskip = 1; % datapoints

m = max(data(:,2));
X = data(skip:end-endskip,1);
Y = data(skip:end-endskip,2);

% error bars on data
if I_unit == 0
    Y = Y/time; % Y1 should be in counts per second
end

uY = sqrt(Y*time)/time;

Y0 = Y;

% Add a footprint correction to the data
doFootPrint = 0;
if doFootPrint
    F = 0.5./(10*sind(X/2));
    F(sind(X/2)>(0.5/10)) = 1;
    Y  = Y.*F;
    uY = uY.*F;
    
    m = max(Y);
end

Y  = Y/m;
uY = uY/m;

% ----- SAMPLE -----

% Layer structure starting from bottom
% if x: set, if [x y z]: fit x between y and z

% layer 1
layer.direction     = 1;
layer.N             = [1e6, 1e6, 3e8];
layer.filename      = 'MgO_001_fractional.vasp';
layer.dinterface    = 0; % 0: bottom
layer.scale         = [1, 0.9, 1.1];  % A
layer.area_scale    = [1, 1, 1];      % A^2 
layer.roughness     = false;          % 
layer.sigma         = 0;              % 
layer.pre_calc_f    = [];
layer.bottom_strain_amplitude = 0.00;
layer.bottom_strain_end       = 0;
layer.top_strain_amplitude = 0;
layer.top_strain_end       = 0; % Stop straining after bottom_strain_end number of atoms (not u.c.)

stack{1} = layer;

% layer 2
layer.direction     = 1;
layer.N             = [28.5, 10, 40];
layer.filename      = 'Fe_fractional.vasp';
layer.dinterface    = [1.4, 1.3, 1.43];
layer.scale         = [1.04, 0.9, 1.06];
layer.area_scale    = [1.1927, 0.5, 1.5];
layer.roughness     = false;
layer.sigma         = 0;        % RMS roughness value
layer.bottom_strain_amplitude = [0.2, 0.004, 0.4];
layer.bottom_strain_end       = [5, 2, 20]; % Stop straining after bottom_strain_end number of atoms (not u.c.)
layer.top_strain_amplitude    = 0;
layer.top_strain_end          = 0; % Stop straining after top_strain_end number of atoms (not u.c.)

layer.pre_calc_f    = [];

stack{2} = layer;

% stack properties
control.vacuum_thick    = 5;         % Thickness of air on either side of the stack [A]
control.slices          = 50;       % Constant slice thickness throughout the stack is dz = (stack{1}.lat_par*stack{1}.scaling)/slices
control.maxQ0           = 30;        % Max Q value for sampling of the form factor for density conversion
control.stepQ0          = 0.1;       % Step size of the Q value for sampling of the form factor for density conversion
control.pol             = 2;         % 0 - sigma polarization, 1 - pi polarization, 2 - unpolarized
control.model           = 'density'; % kinematic, darwin, density 
control.dotransmission  = false;     % Calculate the transmitted intensity through the stack (beta)
control.plot_density    = false;

% strain from substrate: strained(1) = 1: tensile out-of-plane, -1: compressive out-of-plane, 0: no strain

% ----- INSTRUMENT -----

instrument.wavelength  = 1.54056; % instrument wavelength in Ångström
instrument.theta       = X(:)/2;
instrument.Q           = 4*pi/instrument.wavelength*sind(instrument.theta);
instrument.I0          = [5000.8437, 50, 1000000];    % Scale factor for the model intensity
instrument.resolution  = [0.0054205, 0.0001, 0.006]; % FWHM of a Gaussian function 
instrument.background  = 0;    % background fitting 0: linear, 1: exponential
instrument.theta_mPath = 1;    % optics on 1: incidence side, -1: detector side
instrument.theta_m     = 2;    % Monochromator angle 

% linear: a*Q+b, polynomial: c*(a+Q).^2+b

instrument.a           = [5.1469e-4, -1e-6, 1e0];    % 1st background parameter
instrument.b           = [1.2366e-7, -1e-6, 1e-1];         % 2nd background parameter 4e-7
instrument.c           = 2;                       % 3rd background parameter, irrelevant if background == 0

% ----- FITTING -----
iterations = 100; % fitting iterations

%-----------------------------End: TO DO-----------------------------

S_struct.I            = Y(:);
S_struct.uI           = uY(:);
S_struct.F_weight     = 0.6;
S_struct.F_CR         = 0.7;
S_struct.I_bnd_constr = 1;
S_struct.I_itermax    = iterations;
S_struct.F_VTR        = 0;
S_struct.I_strategy   = 1;
S_struct.I_refresh    = 1;
S_struct.I_plotting   = 1;

[GKP_bestmem,FVr_minbound,FVr_maxbound] = extract_fitting_vector(stack,instrument,control);

S_struct.I_D          = length(GKP_bestmem);
% S_struct.I_NP         = 1*S_struct.I_D; % 23
S_struct.I_NP         = 10;
S_struct.GKP_bestmem  = GKP_bestmem;
S_struct.FVr_minbound = FVr_minbound;
S_struct.FVr_maxbound = FVr_maxbound;
S_struct.N            = length(GKP_bestmem);

% ----- INSTRUMENT -----
S_struct.instrument   = instrument;

% ----- SAMPLE -----
S_struct.stack        = stack;
S_struct.control      = control;

% saving the plot?
S_struct.saving       = 1; % 1: yes, 0: no
S_struct.path         = path;
S_struct.sample       = sample;

% plot 0: 2theta, 1:Q
S_struct.plot_unit    = 0;

% access through GUI?
S_struct.GUI          = 0; % 1: yes, 0: no

PlotIt_laue_full(GKP_bestmem,0,S_struct);
[p,S_y,I_nf] = deopt_full('objfun_full','PlotIt_laue_full',S_struct);
