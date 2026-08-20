%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build:            2026-08-18 12:08 CEST (UTC+02:00)
% Build ID:         20260818T120800+0200
% Author:           refactor based on Gunnar Palsson / Johan Bylin code
% Description:      Repeated-unit-cell Parratt calculation using a smooth,
%                   band-limited atomic scattering-density model.
%
% User-facing defaults are intentionally simple:
%   * The published atomic form factor is used unchanged over every physically
%     accessible Q for the chosen wavelength.
%   * Above that range the ENTIRE form factor (Gaussian terms, c, f', and f'')
%     is smoothly brought to zero at control.maxQ0 (default 75/A). This avoids
%     both an infinite-Q extrapolation and the sinc ringing of a hard Q cutoff.
%   * Atomic density is never cut at one lattice parameter. Contributions are
%     evaluated only until their remaining real-space tail is negligible.
%   * Vacuum is enlarged automatically when needed to contain the retained
%     surface density tail.
%   * The continuous density transform is cached independently of the slice
%     thickness, so fitting calls only translate/add atomic kernels and average
%     them into the current slices.
%   * On first use, a 2x-slice calculation checks numerical convergence and
%     prints one short message. Subsequent fitting calls stay silent.
%
% Optional expert controls (normally leave unset):
%   density_tail_tol          default 1e-10
%   density_qpass             default 4*pi/lambda
%   density_qstop             default control.maxQ0 (or 75/A)
%   density_auto_vacuum       default true
%   auto_convergence_check    default true
%   convergence_rel_tol       default 1e-3
%   convergence_log_tol       default 1e-2 decades
%   density_internal_dq       default 0.01/A (internal cached quadrature)
%
% Existing external dependencies are kept:
%   read_poscar
%   Read_form_factor_coefficients
%   DB_prefactor
%   generate_strain
%   propagate_vectorized_chunks_opt_dual_20260818_1208
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = parratt_matrix_repeated_rhoanalytic(Q,lambda,stack,control,instrument)

build_stamp = '2026-08-18 12:08 CEST (UTC+02:00)';
build_id = '20260818T120800+0200';
output = run_core(Q,lambda,stack,control,instrument);
output.diagnostics.code_build = build_stamp;
output.diagnostics.code_build_id = build_id;
output.diagnostics.source_file = mfilename('fullpath');
run_first_use_check(Q,lambda,stack,control,instrument,output);

end

% =========================================================================
% Core calculation. This is called once normally and once more at 2x slices
% by the first-use convergence checker.
% =========================================================================
function output = run_core(Q,lambda,stack,control,instrument)

output = struct('trans',[],'absorption',[],'I',[],'refl',[],'diagnostics',struct());

if isempty(stack)
    error('parratt_matrix_repeated_rhoanalytic:EmptyStack','stack must contain at least one layer.');
end
if ~isfield(control,'slices') || control.slices < 2
    error('parratt_matrix_repeated_rhoanalytic:BadSlices','control.slices must be >= 2.');
end

n_layers = numel(stack);
validate_roughness_layout(stack);

% 1. Geometry and cached continuous atomic kernels.
layers = prepare_layers(stack,lambda);
slices = round(control.slices);
dz = layers(1).lat_par/slices;
opts = make_density_options(Q,lambda,control);

max_support = 0;
for o = 1:n_layers
    for s = 1:numel(layers(o).atoms)
        model = layers(o).atoms(s);
        model.kernel = get_atomic_kernel(model,opts);
        model.support = model.kernel.support;
        layers(o).atoms(s) = model;
        max_support = max(max_support,model.support);
    end
end

% 2. Atomic positions / strain profiles.
placements = prepare_placements(stack,layers);
totthick = placements(end).max_z;

% 3. Choose a grid. The requested vacuum is a minimum; by default it is
% enlarged automatically so that a retained surface atom is not clipped.
user_vacuum = getfield_default(control,'vacuum_thick',0);
required_vacuum = max_support + 2*dz;
if opts.density_auto_vacuum
    vacuum_target = max(user_vacuum,required_vacuum);
else
    vacuum_target = user_vacuum;
end
vacuum_slices = max(0,ceil(vacuum_target/dz));
vacuum_thick = dz*vacuum_slices;
z = -vacuum_thick:dz:(totthick + vacuum_thick);

% 4. Build the three substrate pieces used by the repeated-cell propagator.
% The lower boundary is terminated, the middle cell is exactly periodic, and
% the upper piece is terminated at the substrate/film interface. Neighboring
% substrate atoms are included automatically until their density support ends.
start_idx = slices;
idx0_end = vacuum_slices + start_idx;
idx1_beg = vacuum_slices + start_idx;
idx1_end = vacuum_slices + 2*start_idx;
substrate_end = vacuum_slices + 2*start_idx;
idx2_end = min(numel(z),vacuum_slices + 3*start_idx);

z_0 = z(1:idx0_end);
z_1 = z(idx1_beg:idx1_end);
z_rest = z(substrate_end:end);

rho_0 = substrate_density_piece(z_0,dz,layers(1),'bottom',0);
rho_1 = substrate_density_piece(z_1,dz,layers(1),'periodic',1);
rho_rest_sub = substrate_density_piece(z_rest,dz,layers(1),'top',2);

% Assemble a full working vector for films/roughness. rho_0 and rho_1 are
% already frozen above, so film tails cannot contaminate the repeated bulk cell.
rho_e = complex(zeros(size(z)));
rho_e(1:idx0_end) = rho_0;
rho_e(idx1_beg:idx1_end) = rho_1;
rho_e(substrate_end:end) = rho_rest_sub;
rho_2 = rho_e(substrate_end:idx2_end);
z_2 = z(substrate_end:idx2_end);

% 5. Add films. Every atom is evaluated to its tolerance-based support only;
% there is no lattice-parameter cutoff.
rho_e_rough = [];
rough_info = struct('nvector',[],'weights_legacy',[]);

for o = 2:n_layers
    if placements(o).is_rough
        nvector = placements(o).nvector;
        rough_weights = placements(o).weights;
        rho_e_rough = cell(max(nvector),1);
        weights_legacy = zeros(1,max(nvector));

        for ii = 1:numel(nvector)
            m = nvector(ii);
            rho_variant = rho_e;
            rho_variant = add_layer_density(rho_variant,z,dz,layers(o), ...
                placements(o).positions{ii},placements(o).sorted{ii});
            rho_e_rough{m} = rho_variant;
            weights_legacy(m) = rough_weights(ii);
        end
        rough_info.nvector = nvector;
        rough_info.weights_legacy = weights_legacy;
    else
        rho_e = add_layer_density(rho_e,z,dz,layers(o), ...
            placements(o).positions,placements(o).sorted);
    end
end

% 6. Propagation -- API-compatible with the existing propagator.
rho_0r = rho_0(end:-1:1);
rho_1r = rho_1(end:-1:1);

if isempty(rho_e_rough)
    rho_rest = rho_e(substrate_end:end);
else
    rho_rest = rho_e_rough{rough_info.nvector(end)}(substrate_end:end);
end

if getfield_default(control,'plot_density',false)
    figure(2); clf;
    plot(z_0,real(rho_0),'-r'); hold on;
    plot(z_1,real(rho_1),'-k');
    plot(z_2,real(rho_2),'-b');
    plot(z_rest,real(rho_rest),'-m');
    xlabel('z (A)'); ylabel('Re scattering density');
    drawnow;
end

final_is_rough = placements(end).is_rough;
if final_is_rough
    nvector = rough_info.nvector;
    rough_weights_legacy = rough_info.weights_legacy;

    if control.pol == 0 || control.pol == 1
        output.refl = propagate_vectorized_chunks_opt( ...
            Q,lambda,rho_0r,rho_1r,round(stack{1}.N),rho_e_rough,dz, ...
            substrate_end,nvector,rough_weights_legacy,control.pol, ...
            control.roughness_averaging);
    elseif control.pol == 2
        [refl_p,refl_s] = propagate_vectorized_chunks_opt( ...
            Q,lambda,rho_0r,rho_1r,round(stack{1}.N),rho_e_rough,dz, ...
            substrate_end,nvector,rough_weights_legacy,2, ...
            control.roughness_averaging);
        output.refl = combine_polarizations(refl_p,refl_s,instrument);
    else
        error('parratt_matrix_repeated_rhoanalytic:BadPolarization','Unknown control.pol value.');
    end
else
    rho_restr = rho_rest(end:-1:1);
    final_sigma = getfield_default(stack{end},'sigma',0);

    if control.pol == 0 || control.pol == 1
        output.refl = propagate_vectorized_chunks_opt( ...
            Q,lambda,rho_0r,rho_1r,round(stack{1}.N),rho_restr,dz, ...
            final_sigma,control.pol);
    elseif control.pol == 2
        [refl_p,refl_s] = propagate_vectorized_chunks_opt( ...
            Q,lambda,rho_0r,rho_1r,round(stack{1}.N),rho_restr,dz, ...
            final_sigma,2);
        output.refl = combine_polarizations(refl_p,refl_s,instrument);
    else
        error('parratt_matrix_repeated_rhoanalytic:BadPolarization','Unknown control.pol value.');
    end
end

output.diagnostics.dz = dz;
output.diagnostics.slices = slices;
output.diagnostics.q_data_max = opts.q_data_max;
output.diagnostics.q_physical_max = opts.q_physical_max;
output.diagnostics.q_valid_max = opts.q_valid_max;
output.diagnostics.density_qpass = opts.qpass;
output.diagnostics.density_qstop = opts.qstop;
% Backward-compatible diagnostic names used by the first-use checker.
output.diagnostics.constant_qpass = opts.qpass;
output.diagnostics.constant_qstop = opts.qstop;
output.diagnostics.max_atom_support = max_support;
output.diagnostics.density_tail_tol = opts.density_tail_tol;
output.diagnostics.density_internal_dq = opts.internal_dq;
output.diagnostics.atomic_kernel_cache = 'continuous-density primitive; independent of slices/dz';
output.diagnostics.vacuum_thick_requested = user_vacuum;
output.diagnostics.vacuum_thick_required = required_vacuum;
output.diagnostics.vacuum_thick_used = vacuum_thick;
output.diagnostics.vacuum_was_extended = vacuum_thick > user_vacuum + dz/2;
output.diagnostics.nyquist_Q = pi/dz;
output.diagnostics.substrate_lat_par = layers(1).lat_par;
output.diagnostics.min_lat_par = min([layers.lat_par]);
output.z = z;
if isempty(rho_e_rough)
    output.rho_e = rho_e;
else
    output.rho_e = rho_e_rough{rough_info.nvector(end)};
end

end

% =========================================================================
% Substrate density pieces for the repeated-cell construction.
% =========================================================================
function rho = substrate_density_piece(z,dz,layer,mode,reference_uc)

rho = complex(zeros(size(z)));
a = layer.lat_par;
zmin = min(z) - dz/2;
zmax = max(z) + dz/2;

for s = 1:numel(layer.z_uc)
    zs = layer.z_uc(s);
    atom = layer.atoms(s);
    R = atom.support;

    % Any atom whose support overlaps [zmin,zmax] must be included.
    nmin = ceil((zmin - R - zs)/a);
    nmax = floor((zmax + R - zs)/a);

    switch mode
        case 'bottom'
            % Crystal begins at unit cell 0 and extends toward +z.
            nmin = max(nmin,0);
        case 'top'
            % The explicit upper substrate cell is reference_uc; no substrate
            % atoms exist above it because the film/interface begins there.
            nmax = min(nmax,reference_uc);
        case 'periodic'
            % No termination: use every neighboring cell that overlaps.
        otherwise
            error('parratt_matrix_repeated_rhoanalytic:BadSubstrateMode', ...
                'Unknown substrate density mode.');
    end

    for n = nmin:nmax
        pos_z = zs + n*a;
        rho = add_atom_density(rho,z,dz,pos_z,layer.area,atom,R);
    end
end

end

% =========================================================================
% Automatic first-use numerical check.
% =========================================================================
function run_first_use_check(Q,lambda,stack,control,instrument,output)

if ~getfield_default(control,'auto_convergence_check',true)
    return;
end

persistent checked_keys
if isempty(checked_keys)
    checked_keys = {};
end

% Key only on numerical settings that should remain fixed throughout a fit.
% Structural fit parameters are deliberately excluded so the check is not
% repeated at every optimizer step.
D0 = output.diagnostics;
key = sprintf('build=%s|lam=%.12g|s=%d|qpass=%.12g|qstop=%.12g|tail=%.3g|pol=%d|%s', ...
    D0.code_build_id,lambda,round(control.slices),D0.constant_qpass,D0.constant_qstop, ...
    getfield_default(control,'density_tail_tol',1e-10),control.pol,stack_identity(stack));
if any(strcmp(checked_keys,key))
    return;
end

control2 = control;
control2.slices = 2*round(control.slices);
control2.plot_density = false;
control2.auto_convergence_check = false;

try
    fine = run_core(Q,lambda,stack,control2,instrument);
    checked_keys{end+1} = key;
catch ME
    warning('parratt_matrix_repeated_rhoanalytic:ConvergenceCheckFailed', ...
        'First-use convergence check could not run: %s',ME.message);
    return;
end

r0 = abs(output.refl(:));
r1 = abs(fine.refl(:));
scale = max([r0; r1; eps]);
rel_peak = max(abs(r1-r0))/scale;

floor_level = max(scale*1e-14,1e-18);
mask = max(r0,r1) > floor_level;
if any(mask)
    log_change = max(abs(log10(max(r1(mask),floor_level)) - ...
                         log10(max(r0(mask),floor_level))));
else
    log_change = 0;
end

rel_tol = getfield_default(control,'convergence_rel_tol',1e-3);
log_tol = getfield_default(control,'convergence_log_tol',1e-2);
issues = {};

if rel_peak > rel_tol || log_change > log_tol
    issues{end+1} = sprintf( ...
        ['doubling slices from %d to %d changes the curve by %.3g of peak ' ...
         'and up to %.3g decades; try control.slices = %d'], ...
        round(control.slices),control2.slices,rel_peak,log_change,control2.slices);
end

D = output.diagnostics;
transition_width = D.constant_qstop - D.constant_qpass;
if transition_width < 10
    issues{end+1} = sprintf( ...
        ['only %.2f A^-1 remains between the protected scattering range and ' ...
         'the trusted form-factor limit; this wavelength is close to the ' ...
         'high-Q limit of the atomic parameterization'],transition_width);
end

if ~getfield_default(control,'density_auto_vacuum',true) && ...
        D.vacuum_thick_requested + D.dz < D.vacuum_thick_required
    issues{end+1} = sprintf( ...
        'vacuum is too small for the retained atomic tail; use at least %.2f A', ...
        D.vacuum_thick_required);
end

if isempty(issues)
    if D.vacuum_was_extended
        vacnote = sprintf('; vacuum automatically increased to %.2f A',D.vacuum_thick_used);
    else
        vacnote = '';
    end
    fprintf(['parratt density check [%s]: OK (%d slices/cell; 2x-slice change %.2g%s).\n'], ...
        D0.code_build_id,round(control.slices),rel_peak,vacnote);
else
    msg = sprintf('parratt density check [%s]: settings need attention:\n',D0.code_build_id);
    for k = 1:numel(issues)
        msg = [msg sprintf('  - %s\n',issues{k})]; %#ok<AGROW>
    end
    warning('parratt_matrix_repeated_rhoanalytic:NumericalSettings','%s',msg);
end

end

% =========================================================================
% Layer preparation
% =========================================================================
function layers = prepare_layers(stack,lambda)

n_layers = numel(stack);
layers = repmat(struct('lat_par',[],'area',[],'z_uc',[],'atoms',[], ...
                       'type',[],'type_nr',[],'r',[]),n_layers,1);

Zt = {'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr', ...
      'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd', ...
      'In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', ...
      'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra', ...
      'Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'};

for o = 1:n_layers
    [type,type_nr,r,a1,a2,a3] = load_layer_geometry(stack{o});

    as = sqrt(stack{o}.area_scale);
    switch stack{o}.direction
        case 1
            scaling = a1*stack{o}.scale;
            area = abs(norm(cross(a2*as,a3*as)));
        case 2
            scaling = a2*stack{o}.scale;
            area = abs(norm(cross(a1*as,a3*as)));
        case 3
            scaling = a3*stack{o}.scale;
            area = abs(norm(cross(a1*as,a2*as)));
        otherwise
            error('parratt_matrix_repeated_rhoanalytic:BadDirection','stack{%d}.direction must be 1, 2, or 3.',o);
    end

    lat_par = norm(scaling);
    z_uc = r*scaling';

    n_atoms = size(r,1);
    atoms = repmat(empty_atom_model(),n_atoms,1);
    atom_counter = 0;
    species_counter = 0;

    for i = 1:numel(type)
        if isempty(type{i})
            continue;
        end
        species_counter = species_counter + 1;
        cur_Z = find(strcmp(Zt,type{i}),1,'first');
        if isempty(cur_Z)
            error('parratt_matrix_repeated_rhoanalytic:UnknownElement','Unknown element symbol "%s".',type{i});
        end

        FF = read_form_factor_coefficients_cached(cur_Z,lambda);
        dw_B = DB_prefactor(cur_Z);
        species_model = make_atom_model(FF,dw_B,type{i},cur_Z);

        n_this = type_nr(species_counter);
        for j = 1:n_this
            atom_counter = atom_counter + 1;
            atoms(atom_counter) = species_model;
        end
    end

    if atom_counter ~= n_atoms
        error('parratt_matrix_repeated_rhoanalytic:AtomCountMismatch', ...
            'POSCAR atom count (%d) and type/type_nr count (%d) disagree in layer %d.', ...
            n_atoms,atom_counter,o);
    end

    layers(o).lat_par = lat_par;
    layers(o).area = area;
    layers(o).z_uc = z_uc(:);
    layers(o).atoms = atoms;
    layers(o).type = type;
    layers(o).type_nr = type_nr;
    layers(o).r = r;
end

end

function [type,type_nr,r,a1,a2,a3] = load_layer_geometry(layer)
if isfield(layer,'pre_calc_f') && ~isempty(layer.pre_calc_f)
    pc = layer.pre_calc_f;
    a1 = pc(1).a1;
    a2 = pc(1).a2;
    a3 = pc(1).a3;
    type = pc(1).type;
    type_nr = pc(1).type_nr;
    r = pc(1).r;
    return;
end

% POSCAR geometry is static during a fit. Cache it internally so the density
% model does not depend on feeding a modified stack back through calc_model_full.
persistent geometry_keys geometry_values
if isempty(geometry_keys)
    geometry_keys = {};
    geometry_values = {};
end
key = char(layer.filename);
idx = find(strcmp(geometry_keys,key),1,'first');
if isempty(idx)
    [type,type_nr,r,a1,a2,a3] = read_poscar(layer.filename);
    G = struct('type',{type},'type_nr',type_nr,'r',r,'a1',a1,'a2',a2,'a3',a3);
    geometry_keys{end+1} = key;
    geometry_values{end+1} = G;
else
    G = geometry_values{idx};
    type = G.type;
    type_nr = G.type_nr;
    r = G.r;
    a1 = G.a1;
    a2 = G.a2;
    a3 = G.a3;
end
end


function FF = read_form_factor_coefficients_cached(Z,lambda)
% The function is called inside fitting loops. Atomic coefficients depend only
% on element and wavelength, so avoid re-reading the database every iteration.
persistent keys values
if isempty(keys)
    keys = {};
    values = {};
end
key = sprintf('%d|%.12g',Z,lambda);
idx = find(strcmp(keys,key),1,'first');
if isempty(idx)
    FF = Read_form_factor_coefficients(Z,lambda);
    keys{end+1} = key;
    values{end+1} = FF;
else
    FF = values{idx};
end
end

function model = empty_atom_model()
model = struct('symbol','','Z',[],'a',zeros(1,5),'beta',zeros(1,5), ...
               'constant',0,'dw_B',0,'kernel',struct(),'support',0);
end

function model = make_atom_model(FF,dw_B,symbol,Z)
model = empty_atom_model();
model.symbol = symbol;
model.Z = Z;
model.dw_B = dw_B;

% This mirrors Form_factors exactly for one pure element:
% f(Q) = c + f_1 - i f_2 + sum_j a_j exp[-b_j (Q/4pi)^2]
model.a = reshape(FF.a(1,1:5),1,5);
model.beta = reshape(FF.b(1,1:5),1,5) + dw_B;
model.constant = FF.c(1) + FF.f_1(1) - 1i*FF.f_2(1);

if any(model.beta <= 0)
    error('parratt_matrix_repeated_rhoanalytic:BadGaussianWidth', ...
        'Non-positive Gaussian beta encountered for %s.',symbol);
end
end

% =========================================================================
% Position / strain preparation
% =========================================================================
function placements = prepare_placements(stack,layers)

n_layers = numel(stack);
z_ss = {layers.z_uc};
lat_par = [layers.lat_par];
placements = repmat(struct('is_rough',false,'positions',[],'sorted',[], ...
                           'nvector',[],'weights',[],'max_z',[]),n_layers,1);

% Three explicit substrate unit cells are used by the repeated-UC scheme.
last_atom_z = 2*layers(1).lat_par + max(layers(1).z_uc);
placements(1).max_z = last_atom_z;

for o = 2:n_layers
    N = round(stack{o}.N);
    startz = last_atom_z + stack{o}.dinterface;
    is_rough = logical(getfield_default(stack{o},'roughness',false));
    placements(o).is_rough = is_rough;

    if is_rough
        sigma = stack{o}.sigma;
        [nvector,weights] = roughness_distribution(N,sigma);
        positions = cell(numel(nvector),1);
        sorted = cell(numel(nvector),1);
        max_z = -inf;
        for ii = 1:numel(nvector)
            m = nvector(ii);
            [positions{ii},sorted{ii}] = generate_strain( ...
                z_ss{o},lat_par,m,o,startz,stack);
            max_z = max(max_z,max(positions{ii}));
        end
        placements(o).nvector = nvector;
        placements(o).weights = weights;
        placements(o).positions = positions;
        placements(o).sorted = sorted;
        placements(o).max_z = max_z;
        last_atom_z = max_z;
    else
        [positions,sorted] = generate_strain( ...
            z_ss{o},lat_par,N,o,startz,stack);
        placements(o).positions = positions;
        placements(o).sorted = sorted;
        placements(o).max_z = max(positions);
        last_atom_z = placements(o).max_z;
    end
end

end

function [nvector,weights] = roughness_distribution(N,sigma)
nvector = max(1,round(N-3*sigma)):max(1,round(N+3*sigma));
if sigma == 0
    weights = ones(size(nvector));
else
    A = 1/(sqrt(2*pi)*sigma);
    weights = A*exp(-(nvector-N).^2/(2*sigma^2));
end
end

function validate_roughness_layout(stack)
rough = false(1,numel(stack));
for o = 1:numel(stack)
    rough(o) = isfield(stack{o},'roughness') && logical(stack{o}.roughness);
end
if rough(1)
    error('parratt_matrix_repeated_rhoanalytic:SubstrateRoughness', ...
        'The current repeated-substrate propagator does not support substrate roughness in this function.');
end
rough_idx = find(rough);
if ~isempty(rough_idx) && any(rough_idx ~= numel(stack))
    error('parratt_matrix_repeated_rhoanalytic:InternalRoughLayer', ...
        ['The existing propagate_vectorized_chunks_opt interface only carries one roughness ensemble. ' ...
         'This refactor therefore requires roughness, if present, to be on the final layer.']);
end
end

% =========================================================================
% Density options and smooth full-form-factor regularization
% =========================================================================
function opts = make_density_options(Q,lambda,control)
q_data_max = max(abs(Q(:)));
opts.density_tail_tol = getfield_default(control,'density_tail_tol',1e-10);
opts.density_auto_vacuum = getfield_default(control,'density_auto_vacuum',true);
opts.q_valid_max = getfield_default(control,'maxQ0',75);
opts.internal_dq = getfield_default(control,'density_internal_dq',0.01);

if opts.density_tail_tol <= 0 || opts.density_tail_tol >= 1
    error('parratt_matrix_repeated_rhoanalytic:BadTailTolerance', ...
        'density_tail_tol must lie between 0 and 1.');
end
if opts.internal_dq <= 0
    error('parratt_matrix_repeated_rhoanalytic:BadInternalDQ', ...
        'density_internal_dq must be positive.');
end

% Protect the complete physically accessible scattering range.  Above this
% point we do not claim knowledge that the tabulated parameterization does not
% provide: the entire form factor is tapered smoothly to zero at qstop.
q_physical_max = 4*pi/lambda;
opts.qpass = getfield_default(control,'density_qpass',q_physical_max);
opts.qpass = max(opts.qpass,q_data_max);
opts.qstop = getfield_default(control,'density_qstop',opts.q_valid_max);
opts.qstop = min(opts.qstop,opts.q_valid_max);

if opts.qpass >= opts.qstop
    error('parratt_matrix_repeated_rhoanalytic:QRegularizationRange', ...
        ['The wavelength leaves no room for the automatic high-Q taper: ' ...
         'Qpass=%.3g and Qstop=%.3g A^-1. Use a less energetic wavelength, or ' ...
         'increase the trusted form-factor limit only if justified.'], ...
        opts.qpass,opts.qstop);
end

opts.q_data_max = q_data_max;
opts.q_physical_max = q_physical_max;
end

% =========================================================================
% Density generation from cached continuous atomic kernels
% =========================================================================
function rho = add_layer_density(rho,z,dz,layer,positions,sorted_idx)

n_uc_atoms = numel(layer.z_uc);
if mod(numel(positions),n_uc_atoms) ~= 0
    error('parratt_matrix_repeated_rhoanalytic:PositionCount', ...
        'generate_strain returned a position count not divisible by atoms/unit-cell.');
end
N = numel(positions)/n_uc_atoms;

ll = 1;
for l = 1:N
    for s = 1:n_uc_atoms
        atom_index = sorted_idx(s);
        pos_z = positions(ll);
        atom = layer.atoms(atom_index);
        rho = add_atom_density(rho,z,dz,pos_z,layer.area,atom,atom.support);
        ll = ll + 1;
    end
end
end

function rho = add_atom_density(rho,z,dz,pos_z,area,model,cutoff)
if cutoff <= 0 || isempty(z)
    return;
end

% z is uniform. Compute the overlapping index range directly instead of
% scanning the complete density vector with find() for every atom.
z0 = z(1);
grid_pad = 1e-9*max(dz,1);
lower = pos_z - cutoff - dz/2 - grid_pad;
upper = pos_z + cutoff + dz/2 + grid_pad;
i0 = max(1,ceil((lower-z0)/dz) + 1);
i1 = min(numel(z),floor((upper-z0)/dz) + 1);
if i1 < i0
    return;
end
idx = i0:i1;
xc = z(idx)-pos_z;

% kernel.primitive is the primitive of the continuous band-limited atomic
% density. Therefore the exact slice average is [P(x+dz/2)-P(x-dz/2)]/dz.
% The continuous kernel is independent of dz and remains cached even if a
% lattice parameter (and hence dz) is varied by the fitting algorithm.
p_hi = primitive_interp_uniform(model.kernel,xc + dz/2);
p_lo = primitive_interp_uniform(model.kernel,xc - dz/2);
contrib = (p_hi-p_lo)/dz;
rho(idx) = rho(idx) + contrib/area;
end

function P = primitive_interp_uniform(kernel,xq)
% Fast linear interpolation on the uniformly sampled odd primitive P(-x)=-P(x).
% Outside the cached x range the primitive has already reached its asymptote,
% so clamp to its final value rather than extrapolating.
was_row = isrow(xq);
xv = xq(:);
sgn = sign(xv);
a = abs(xv);
vals = kernel.primitive(:);
n = numel(vals);
Pabs = complex(zeros(size(a)));

if n == 1 || kernel.dx <= 0
    Pabs(:) = vals(end);
else
    xmax = kernel.dx*(n-1);
    outside = a >= xmax;
    Pabs(outside) = vals(end);

    inside = ~outside;
    if any(inside)
        u = a(inside)/kernel.dx;
        k0 = floor(u);              % zero-based lower sample index
        frac = u-k0;
        k = k0 + 1;                 % MATLAB index
        k = min(k,n-1);
        Pabs(inside) = vals(k) + frac.*(vals(k+1)-vals(k));
    end
end
P = sgn.*Pabs;
P(xv == 0) = 0;
if was_row
    P = P.';
end
end

function kernel = get_atomic_kernel(model,opts)
% Continuous atomic kernels depend on species, wavelength-dependent anomalous
% terms, and the Q regularization -- but NOT on the Parratt slice thickness.
% This is important during fitting: changing a lattice parameter changes dz,
% but no expensive Q transform has to be rebuilt.
persistent keys values
if isempty(keys)
    keys = {};
    values = {};
end
astr = sprintf('%.12g,',model.a);
bstr = sprintf('%.12g,',model.beta);
key = sprintf(['Z=%d|B=%.12g|C=%.12g,%.12g|a=%s|b=%s|qp=%.12g|qs=%.12g|' ...
               'tol=%.3g|dq=%.6g'], ...
    model.Z,model.dw_B,real(model.constant),imag(model.constant),astr,bstr, ...
    opts.qpass,opts.qstop,opts.density_tail_tol,opts.internal_dq);
idx = find(strcmp(keys,key),1,'first');
if ~isempty(idx)
    kernel = values{idx};
    return;
end

kernel = build_atomic_kernel(model,opts);
keys{end+1} = key;
values{end+1} = kernel;
end

function kernel = build_atomic_kernel(model,opts)
% Transform the complete published form factor multiplied by the C-infinity
% flat-passband window. The resulting continuous density and its primitive are
% cached. Slice averages for arbitrary dz are then obtained from the primitive.

nint = ceil(opts.qstop/opts.internal_dq);
if mod(nint,2) ~= 0
    nint = nint + 1;
end
q = linspace(0,opts.qstop,nint+1);
h = opts.qstop/nint;
sw = ones(size(q));
sw(2:2:end-1) = 4;
sw(3:2:end-2) = 2;
sw = sw*(h/3);

u = (q/(4*pi)).^2;
fq = model.constant .* exp(-model.dw_B*u);
for j = 1:numel(model.a)
    if model.a(j) ~= 0
        fq = fq + model.a(j).*exp(-model.beta(j)*u);
    end
end
W = smooth_flat_window(q,opts.qpass,opts.qstop);
base = 2 * (sw .* fq .* W);

% 32 samples per shortest retained real-space period. This interpolation grid
% is a numerical cache only and does not define the Parratt slice thickness.
dx = pi/(16*opts.qstop);
rmax = 8;
for iter = 1:8
    x = 0:dx:rmax;
    [rho_x,primitive_x] = density_primitive_sum_blocks(q,base,x);
    scale = max(abs(rho_x));
    if scale == 0
        kernel.dx = dx;
        kernel.primitive = complex([0 0]);
        kernel.support = 0;
        return;
    end

    env = abs(rho_x);
    for k = numel(env)-1:-1:1
        env(k) = max(env(k),env(k+1));
    end
    idx = find(env <= opts.density_tail_tol*scale,1,'first');
    if ~isempty(idx)
        kernel.dx = dx;
        % Keep the full computed primitive out to rmax. It is cheap in memory
        % and gives a safe asymptotic value for slice edges beyond support.
        kernel.primitive = primitive_x;
        kernel.support = x(idx);
        return;
    end
    rmax = 2*rmax;
end

warning('parratt_matrix_repeated_rhoanalytic:KernelSupport', ...
    ['Atomic density tail did not reach the requested %.1e tolerance within ' ...
     '%.1f A; using that support.'],opts.density_tail_tol,rmax);
kernel.dx = dx;
kernel.primitive = primitive_x;
kernel.support = rmax;
end

function [rho_x,P_x] = density_primitive_sum_blocks(q,base,x)
% Compute, in blocks,
%   rho(x) = 2 int f(q)W(q) cos(qx) dq
%   P(x)   = int_0^x rho(t) dt
%          = 2 int f(q)W(q) sin(qx)/q dq.
% Simpson weights and the factor of two are already included in base.
rho_x = complex(zeros(size(x)));
P_x = complex(zeros(size(x)));
block = 128;
qcol = q(:);
brow = base(:).';
if numel(q) > 1
    qpos = q(2:end);
    bdiv = base(2:end)./qpos;
else
    qpos = [];
    bdiv = [];
end

for i = 1:block:numel(x)
    jj = i:min(i+block-1,numel(x));
    xb = x(jj);
    phase = qcol*xb;
    rho_x(jj) = brow*cos(phase);
    if isempty(qpos)
        P_x(jj) = base(1)*xb;
    else
        P_x(jj) = base(1)*xb + bdiv*sin(phase(2:end,:));
    end
end
end

function W = smooth_flat_window(q,qpass,qstop)
% C-infinity flat-top window. All derivatives vanish at qpass and qstop,
% giving much faster-decaying real-space tails than a rectangular cutoff.
W = ones(size(q));
W(q >= qstop) = 0;
idx = q > qpass & q < qstop;
if any(idx)
    t = (q(idx)-qpass)/(qstop-qpass);
    left = exp(-1./t);
    right = exp(-1./(1-t));
    s = left./(left+right);
    W(idx) = 1-s;
end
end

% =========================================================================
% Small helpers
% =========================================================================
function value = getfield_default(s,name,default_value)
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default_value;
end
end

function id = stack_identity(stack)
% Stable first-use-check key: include layer count and source filenames, but not
% fit parameters such as scale, N, strain or roughness that may change at every
% optimizer iteration.
parts = cell(1,numel(stack));
for k = 1:numel(stack)
    if isfield(stack{k},'filename') && ~isempty(stack{k}.filename)
        src = char(stack{k}.filename);
    elseif isfield(stack{k},'pre_calc_f') && ~isempty(stack{k}.pre_calc_f)
        pc = stack{k}.pre_calc_f;
        if isfield(pc(1),'type')
            src = strjoin(pc(1).type,',');
        else
            src = sprintf('layer%d',k);
        end
    else
        src = sprintf('layer%d',k);
    end
    parts{k} = src;
end
id = strjoin(parts,'|');
end

function refl = combine_polarizations(refl_p,refl_s,instrument)
mono = cosd(instrument.theta_m*2).^2;
refl = (refl_s + mono.*refl_p)./(1 + mono);
end
