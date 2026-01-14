% GaAs substrate 
lambda = 1.54056;
theta  = (0:0.01:90)';
Q = 4*pi/lambda*sind(theta);

layer.direction     = 3;
layer.N             = 1e8;
layer.filename      = 'GaAs_alt_fractional.vasp';
layer.dinterface    = 0;
layer.scale         = 1.001;
layer.area_scale    = 1.001;
layer.roughness     = false;
layer.sigma         = 0.0;
layer.pre_calc_f    = [];

control.vacuum_thick    = 20;        % Thickness of air on either side of the stack [A]
control.slices          = 400;       % Constant slice thickness throughout the stack is dz = (stack{1}.lat_par*stack{1}.scaling)/slices
control.maxQ0           = 75;        %75 is max % Max Q value for sampling of the form factor for density conversion
control.stepQ0          = 0.01;      % Step size of the Q value for sampling of the form factor for density conversion
control.pol             = 0;         % 0 - sigma polarization, 1 - pi polarization, 2 - unpolarized
control.model           = 'density'; % kinematic, darwin, density 
control.plot_density    = false;

instrument.theta_mPath = 1;    % optics on 1: incidence side, -1: detector side
instrument.theta_m     = 2;    % Monochromator angle 

stack{1} = layer;

output0 = parratt_matrix_repeated_rhobuiltin(Q,lambda,stack,control,instrument);
control.pol = 1;
output1 = parratt_matrix_repeated_rhobuiltin(Q,lambda,stack,control,instrument);
control.pol = 2;
output2 = parratt_matrix_repeated_rhobuiltin(Q,lambda,stack,control,instrument);

fig = figure('Color','w','Position',[1 1 800 600]);
hold on;
plot(theta,output0.refl,'-k','LineWidth',2);
plot(theta,output1.refl,'--k','LineWidth',2);
plot(theta,output2.refl,'-.k','LineWidth',2);
set(gca,'YScale','log');
legend('density \sigma','density \pi','density partial pol')
box on;
set(gca,'FontSize',16)
xlabel('\alpha_i');
ylabel('I/I_i');
set(gca,'LineWidth',1.5)