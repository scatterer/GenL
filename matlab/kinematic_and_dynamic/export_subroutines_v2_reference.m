clear;
script_dir = fileparts(mfilename('fullpath'));
repo_dir = fullfile(script_dir, '..', '..');
addpath(fullfile(script_dir, 'subroutines_v2'));
output_dir = fullfile(repo_dir, 'validation', 'matlab_v2');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

lambda = 1.54056;
instrument.theta_mPath = 1;
instrument.theta_m = 2;
control.vacuum_thick = 5;
control.slices = 100;
control.maxQ0 = 30;
control.stepQ0 = 0.1;
control.pol = 2;
control.model = 'density';
control.plot_density = false;
control.auto_convergence_check = false;

cases = {};
cases{end+1} = make_case('gaas', 20:0.02:35, ...
    {make_layer(3, 1e6, 'GaAs_alt_fractional.vasp', 0, 1.001, 1.001)});
cases{end+1} = make_case('fe_film', 58.92:0.02:68, ...
    {make_layer(1, 1e6, 'MgO_001_fractional.vasp', 0, 1, 1), ...
     make_layer(1, 28.5, 'Fe_fractional.vasp', 1.4, 1.04, 1.1927)});
cases{end+1} = make_case('w_film', 81:0.02:93, ...
    {make_layer(1, 1e6, 'Al2O3_11-20_fractional.vasp', 0, 1, 1), ...
     make_layer(3, 45, 'W_110_fractional.vasp', 1.4, 1, 1)});

fe_v = {make_layer(1, 1e6, 'MgO_001_fractional.vasp', 0.279226, 1, 1)};
for repeat = 1:11
    fe_v{end+1} = make_layer(1, 13, 'V_fractional.vasp', 1.033, 1.02552, 1); %#ok<SAGROW>
    fe_v{end+1} = make_layer(1, 2, 'Fe_fractional.vasp', 1.64691, 0.97458, 1); %#ok<SAGROW>
end
fe_v{end+1} = make_layer(1, 14, 'V_fractional.vasp', 0.997391, 1.02715, 1);
cases{end+1} = make_case('fe_v_superlattice', 50:0.02:75, fe_v);

for k = 1:numel(cases)
    current = cases{k};
    twotheta = current.twotheta(:);
    Q = 4*pi/lambda*sind(twotheta/2);
    result = parratt_matrix_repeated_rhoanalytic( ...
        Q, lambda, current.stack, control, instrument);
    writematrix([twotheta Q(:) result.refl(:)], ...
        fullfile(output_dir, [current.name '_reflectivity.csv']));
    writematrix([result.z(:) real(result.rho_e(:)) imag(result.rho_e(:))], ...
        fullfile(output_dir, [current.name '_density.csv']));
end

function item = make_case(name, twotheta, stack)
item.name = name;
item.twotheta = twotheta;
item.stack = stack;
end

function layer = make_layer(direction, N, filename, dinterface, scale, area_scale)
layer.direction = direction;
layer.N = N;
layer.filename = filename;
layer.dinterface = dinterface;
layer.scale = scale;
layer.area_scale = area_scale;
layer.roughness = false;
layer.sigma = 0;
layer.pre_calc_f = [];
layer.bottom_strain_amplitude = 0;
layer.bottom_strain_end = 0;
layer.top_strain_amplitude = 0;
layer.top_strain_end = 0;
end
