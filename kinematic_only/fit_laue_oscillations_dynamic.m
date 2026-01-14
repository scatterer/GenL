%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:           Anna L. Ravensburg, Gunnar K. Pálsson
% Description:      Read in your x-ray diffraction data, put in which
%                   materials the sample consists of.
%                   Run the program to fit the Laue oscillations. It will automatically fit
%                   and plot your data including the fit.
% To do as user:    The lables "TO DO" indicate, where input is required.
%                   Please fill in all your file names and parameters.
%
% Note:
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;

%% DEFINE DATA AND SAMPLE

%-----------------------------Start: TO DO-----------------------------

% state which samples are supposed to be fitted
sample = 'Example_data_10nmFe'; % sample name

% load data data
I_unit = 1; % 1: cps, 2: counts
fileExt = '.txt'; % extension of data files
path = 'C:\Users\Anna\Box\Anna\Matlab\LaueMatlab\'; % path where data is stored
file = strcat(path, sample, fileExt);
data = load(file);

% state how many iterations per fit you would like to have
itters = 10; % fitting iterations

% material
Z = [26]; % array of Z numbers of elements present [Z1, Z2, ...]
density = 2/(2.866^3); % unit cell density in atoms/Å^3
composition = [100]; % array of atomic percent of elements [x1, x2, ...]
debye_waller_coeff = -0.3328; % coefficient for debye-waller-factor in ...?

% instrument
theta_m = 2; % polarisation, 0: unpolarized, x: incident angle of e.g. Göbel mirror, ...
theta_mPath = 1; % optics on 1: incidence side, -1: detector side

% sample diffraction parameters
thickness = 100; % layer thickness in Ångström
peak = 64.92; % estimated peak position in 2Theta in degrees
range = 6; % 2Theta range in degrees around the peak to be included in the fit

% create min and max bounds
Nmax_min_per = 0.80; % minimum percentage [0,1] of sample contributing to coherent scattering
peak_min = peak*0.9965; % minimum peak position in 2Theta in degrees
peak_max = peak*1.0035; % maximum peak position in 2Theta in degrees

% does the layer grow strained?
strained = [0,0]; % [0,0]: no
% strain from substrate: strained(1) = 1: tensile out-of-plane, -1: compressive out-of-plane, 0: no strain
% strain from capping: strained2) = 1: tensile out-of-plane, -1: compressive out-of-plane, 0: no strain

% does the substrate need to be included?
substrate = 0; % 1: yes, 0: no

% does roughness need to be included?
roughness = 0; % 1: yes, 0: no

% measurement parameters
lambda   = 1.5406; % instrument wavelength in Ångströmö
time = 25; % measuring times per point

resolution = 0.005; % diffractometer resolution
amplitude = 0.0069; % scale factor of intensity
background = 0; % background fitting 0: linear, 1: exponential
% linear: a*Q+b, polynomial: c*(a+Q).^2+b
a = 0; % 1st background parameter
b = 0.1; % 2nd background parameter
c = 0; % 3rd background parameter, irrelevant if background == 0

% create min and max bounds
resolution_min = resolution/20;
resolution_max = resolution*5;
amplitude_min = 1e-3;
amplitude_max = 0.05;
a_min = 0;
a_max = 1.0;
b_min = 0;
b_max = 3;
c_min = 0; % irrelevant if background == 0
c_max = 0; % irrelevant if background == 0

% if strain, scattering from the substrate/cap, or layer roughness should be 
% included, please edit the next lines. Otherwise they are automatically put
% to 0.

if strained(1) ~= 0
    % including strain in the fitting from the substrate
    alpha_1 = 0.03; % 1st strain parameter: substrate

    % create min and max bounds
    alpha_1_min = -0.2;
    alpha_1_max = 0.2;

end

if strained(2) ~= 0
    % including strain in the fitting from the cap
    alpha_2 = 0.03; % 2nd strain parameter: cap

    % create min and max bounds
    alpha_2_min = -0.2;
    alpha_2_max = 0.2;
end

if substrate == 1
    % the substrate peak is modeled with a lorentzian profile
    I_lor    = 51; % intensity substrate peak
    w_lor    = 0.004; % width substrate peak
    peak_lor = 29.955; % substrate peak position in 2Theta in Ångström
    x0_lor   = lambda/(2*sin((peak_lor/2)/180*pi)); % d spacing for hkl peak in Ångström

    % create min and max bounds
    I_lor_min  = I_lor*0.9;
    I_lor_max  = I_lor*1.1;
    w_lor_min  = w_lor*0.1;
    w_lor_max  = w_lor*1.05;
    x0_lor_min = x0_lor * 0.999;
    x0_lor_max = x0_lor * 1.001;
end

if roughness == 1
    % including roughness in the fitting
    sigma = 3; % layer roughness in Ångström

    % create min and max bound
    sigma_min = 0.1;
    sigma_max = thickness/10;
end

%-----------------------------End: TO DO-------------------------------
% no strain?
if strained(1) == 0
    % excluding strain in the fitting
    alpha_1 = 0;
    alpha_1_min = 0;
    alpha_1_max = 0;
end

if strained(2) == 0
    % excluding strain in the fitting
    alpha_2 = 0;
    alpha_2_min = 0;
    alpha_2_max = 0;
end

% no substrate visible?
if substrate == 0
    % excluding the substrate
    I_lor  = 0;
    w_lor  = 0;
    x0_lor = 0;
    I_lor_min  = 0;
    I_lor_max  = 0;
    w_lor_min  = 0;
    w_lor_max  = 0;
    x0_lor_min = 0;
    x0_lor_max = 0;
end

% no substrate roughness included?
if roughness == 0
    % excluding roguhness
    sigma = 0;
    sigma_min = 0;
    sigma_max = 0;
end

% linear background?
if background == 0
    c = 0;
    c_min = c;
    c_max = c;
end

% calculate d and min and max bounds
d = lambda/(2*sin((peak/2)/180*pi));
d_min = lambda/(2*sin((peak_max/2)/180*pi));
d_max = lambda/(2*sin((peak_min/2)/180*pi));

% calculate average coherent scattering thickness
Nmax = thickness/d; % maximum number of atomic layers contributing
Nmax_min = Nmax * Nmax_min_per;

%% SPLITTING DATA IN X AND Y
% data range
X1 = data(data(:,1)>=(peak-range),1);
X1 = X1(X1(:,1)<=(peak+range)); % data in 2Theta in Ångström
Y1 = data(data(:,1)>=(peak-range),2);
Y1 = Y1(1:length(X1)); % intensity data

% error bars on data
if I_unit == 0
    Y1 = Y1/time; % Y1 should be in counts per second
end
uY1 = sqrt(Y1*time)/time;

%% FITTING PARAMETERS, MAX AND MIN BOUND
% [d, alpha_1, bkgLow, bkgHigh, N0, resolution, scale factor, alpha_2, bkgHighc, sigma, Ilorentz, wlorentz, xlorentz]

% start values
GKP_bestmem  = [d,      alpha_1,     a,       b,       Nmax,      resolution,        amplitude,      alpha_2,      c,       sigma,      I_lor,      w_lor,      x0_lor];

% min/max bound
FVr_minbound = [d_min, alpha_1_min,  a_min,   b_min,   Nmax_min,  resolution_min,    amplitude_min,  alpha_2_min,  c_min,   sigma_min,  I_lor_min,   w_lor_min,  x0_lor_min];
FVr_maxbound = [d_max, alpha_1_max,  a_max,   b_max,   Nmax,      resolution_max,    amplitude_max,  alpha_2_max,  c_max,   sigma_max,  I_lor_max,   w_lor_max,  x0_lor_max];   

%% CALCULATING THE REST OF THE PARAMETERS
% twotheta vector
twotheta = X1;
% Q vector
Q = 4*pi/lambda*sind(twotheta/2); % in 1/Ångtröm

% form factor
FF = Read_form_factor_coefficients(Z,lambda);
% Q dependent form factor
[f, f_sqrd_real, f_av_sqrd_real] = Form_factors(Q, FF, composition./100);

% calculate absorption coefficient
re = 2.81794*1e-15*1e10; % in Å
mu = 2*re*density*lambda*imag(f); % in 1/Å
mu = mu * 1e10; % in 1/m
    
%% FITTING SETTINGS
Eout = Q;
uI   = uY1;
Iout = Y1;
stop = 0;
N = length(FVr_minbound);

F_VTR        = stop;
I_D          = length(GKP_bestmem);
I_bnd_constr = 1;      
I_NP         = 15*I_D;
I_itermax    = itters;
F_weight     = 0.7;
F_CR         = 0.8;
I_strategy   = 1;
I_refresh    = 1;
I_plotting   = 1;

n     = length(Iout);
m     = length(GKP_bestmem);
ratio = Iout./uI;

S_struct.lambda       = lambda;
S_struct.Q            = Q(:);
S_struct.I            = Iout(:);
S_struct.L            = length(Q);
S_struct.N            = length(Q);
S_struct.f            = f;
S_struct.uI           = uI(:);
S_struct.I_NP         = I_NP;
S_struct.F_weight     = F_weight;
S_struct.F_CR         = F_CR;
S_struct.I_D          = I_D;
S_struct.FVr_minbound = FVr_minbound;
S_struct.FVr_maxbound = FVr_maxbound;
S_struct.I_bnd_constr = I_bnd_constr;
S_struct.I_itermax    = I_itermax;
S_struct.F_VTR        = F_VTR;
S_struct.I_strategy   = I_strategy;
S_struct.I_refresh    = I_refresh;
S_struct.I_plotting   = I_plotting;
S_struct.GKP_bestmem  = GKP_bestmem;

% sample parameter
S_struct.sample       = sample;
S_struct.path         = path;
S_struct.strained     = strained;
S_struct.substrate    = substrate;
S_struct.roughness    = roughness;
S_struct.background   = background;
S_struct.dwf          = debye_waller_coeff;
S_struct.Nmax         = Nmax;
S_struct.theta_m      = theta_m;
S_struct.theta_mPath  = theta_mPath;
S_struct.mu           = mu;

% saving the plot?
S_struct.saving       = 1; % 1: yes, 0: no

% plot 0: 2theta, 1:Q
S_struct.plot_unit    = 0;

% access through GUI?
S_struct.GUI          = 0; % 1: yes, 0: no

%% FITTING
[p,S_y,I_nf] = deopt('objfun_laue','PlotIt_laue',S_struct);
    
%% OUTPUT

% p vector contains final parameters
% [d, alpha_1, a, b, N0, resolution, scale factor, alpha_2, c, sigma, Ilor, wlor, xlor]

% sample
d_fit     = p(1); % d spacing in Ångström
N0_fit    = p(5); % number of coherently scattering planes

% strain
alpha_1_fit = p(2); % 1st strain parameter
alpha_2_fit = p(8); % 2nd strain parameter

% background
a_fit     = p(3); % 1st background parameter
b_fit     = p(4); % 2nd background parameter
c_fit     = p(9); % 3rd background parameter

% roughness
sigma_fit = p(10); % layer roughness

% substrate
I_lor_fit        = p(11); % intensity substrate peak
w_lor_fit        = p(12); % width substrate peak
d_sub_fit  = p(13); % d spacing substrate peak in Ångström

% machine
resolution_fit = p(6);
amplitude_fit  = p(7); 

% create stack of atoms
x = ones(1,round(N0_fit))*d_fit;
pos = cumsum(x);
pos = [0 pos];

% same array upside down
pos_ud = flip(pos,2);

% calculate substrate strain
if S_struct.strained(1) ~= 0
    strain_bottom = exp(-alpha_1*pos);
    if S_struct.strained(1) == 1
        pos = pos - strain_bottom;
    elseif S_struct.strained(1) == -1
         pos = pos + strain_bottom;
    end
end
% calculate cap strain
if S_struct.strained(2) ~= 0
    strain_top = exp(-alpha_2*pos_ud);
    if S_struct.strained(2) == 1
        pos = pos + strain_top;
    elseif S_struct.strained(2) == -1
         pos = pos - strain_top;
    end
end


dspacings = diff(pos); % array of d spacings for each atom
dav       = mean(diff(pos)); % average of d spacing over layer thickness

% quantitative estimate of sample crystal quality
coherent = N0_fit/Nmax; % how much scatters coherently?

%% SAVING
% save fit parameters
fileID = fopen(strcat(path, sample,'_LaueFit_parameters','.txt'),'w');
fprintf(fileID,strcat('Date: ', string(datetime("now")), '\n'));
fprintf(fileID,'d,    N0,    alpha_1,    alpha_2,    a,    b,    c,    sigma,     resolution,     scale factor,    I,    w,     x,     dav,     N0/N\n');
fprintf(fileID,'%.7f    %.7f    %.7f    %.7f    %.7f    %.7f    %.7f    %.7f    %.7f    %.7f    %.7f    %.7f    %.7f    %.7f    %.7f', d_fit, N0_fit, alpha_1_fit, alpha_2_fit, a_fit, b_fit, c_fit, sigma_fit, resolution_fit, amplitude_fit, I_lor_fit, w_lor_fit, d_sub_fit, dav, coherent);
fprintf(fileID,'\n');
fclose(fileID);

% save d spacing over layer thickness
fileID = fopen(strcat(path, sample,'_LaueFit_d','.txt'),'w');
for z = 1:length(dspacings)
    n = z+1;
    fprintf(fileID,'%.7f    %.7f', n, dspacings(z));
    fprintf(fileID,'\n');
end
fclose(fileID);

%% PLOTTING
% plot final plot
PlotIt_laue(p,1,S_struct);

% plot d spacing over layer thickness
fig_d = figure(2);
set(fig_d,'Color','w');

n = [2:length(dspacings)+1]; % number of layers
plot(n(2:end),dspacings(2:end),'-r','LineWidth',1.5); % plot evolution of d
set(gca,'YScale','linear');
xlabel('n'); % x axis lable
ylabel('d [Å]'); % y axis lable
title(strcat(sample)); % title
box on;
set(gca,'FontSize',14);
drawnow;

% save the graph
saveas(fig_d, strcat(path, sample,'_graph_d'), 'png');

% plot fitted parameter within its limits
fig_para = figure(3);
set(fig_para,'Color','w');

% plot all fitted parameters within their upper and lower limit
for l = 1:length(FVr_minbound)
    scatter(l,(p(l)-FVr_minbound(l))/(FVr_maxbound(l)-FVr_minbound(l)),150,'filled', 'r');
    hold on;
end
set(gca,'YScale','linear');
xlabel('Fitting parameters'); % x axis lable
xticks(linspace(1,13,13));
xticklabels({'d', 'alpha1', 'a', 'b', 'N0', 'res', 'amp', 'alpha2', 'c', 'sigma', 'I', 'w', 'x'});
ylabel('Parameter bound'); % y axis lable
title(strcat(sample)); % title
xlim([0 14]); % x axis limits
ylim([0 1]); % y axis limits
grid on;
set(gca,'FontSize',14);
drawnow;

% save the graph
saveas(fig_para, strcat(path, ...
    sample,'_graph_parameter'), 'png');
