%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:           Gunnar Palsson / dual-polarization optimization
% Build:            2026-08-18 12:08 CEST (UTC+02:00)
% Build ID:         20260818T120800+0200
% Description:      Reflection-recursion propagator with shared sigma/pi
%                   propagation for pol=2.
%
% pol = 0: sigma, one output
% pol = 1: pi,    one output
% pol = 2: calculate BOTH polarizations in one traversal. In this case:
%          [refl_p,refl_s] = propagate_vectorized_chunks_opt_dual_20260818_1208(...)
%
% For pol=2 the expensive kz=sqrt(...) and propagation exp(...) are evaluated
% only once per slice and shared by sigma and pi. The Fresnel coefficients and
% nonlinear reflection recursions remain separate, so this is exact relative to
% calling the existing sigma and pi branches independently.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [refl,refl_second] = propagate_vectorized_chunks_opt(Q,lambda,rho_0r,rho_1r,N,rho_arg,dz,varargin)

    refl_second = [];

    % Dispatch on the existing argument layouts.
    if numel(varargin) == 2
        rough_mode = false;
        rho_restr  = rho_arg;
        sigma      = varargin{1}; %#ok<NASGU>
        pol        = varargin{2};
        average_mode = '';
        coherent_mode = false;

    elseif numel(varargin) == 4 || numel(varargin) == 5
        rough_mode     = true;
        rho_e_rough    = rho_arg;
        substrate_end  = varargin{1};
        nvector        = varargin{2};
        BB             = varargin{3};
        pol            = varargin{4};

        if numel(varargin) == 5
            average_mode = varargin{5};
            if isstring(average_mode)
                average_mode = char(average_mode);
            end
            if ~ischar(average_mode)
                error('average mode must be ''coherent'' or ''incoherent''.');
            end
            average_mode = lower(strtrim(average_mode));
        else
            average_mode = 'coherent';
        end

        if ~strcmp(average_mode,'coherent') && ~strcmp(average_mode,'incoherent')
            error('average mode must be ''coherent'' or ''incoherent''.');
        end
        coherent_mode = strcmp(average_mode,'coherent');
        [rough_stacks,weights] = prepare_crystal_roughness_inputs( ...
            rho_e_rough,substrate_end,nvector,BB);

    else
        error(['Unsupported call signature. Use either the normal 9-input ' ...
               'call or the top-layer crystal-roughness 11/12-input call.']);
    end

    if pol ~= 0 && pol ~= 1 && pol ~= 2
        error('pol must be 0 (sigma), 1 (pi), or 2 (both).');
    end

    m = numel(Q);
    if m == 0
        refl = zeros(size(Q));
        if pol == 2
            refl_second = zeros(size(Q));
        end
        return;
    end

    chunk_size = choose_chunk_size(m,rough_mode);
    if pol == 2 && m > 50000
        % Dual mode carries two Mobius maps/recursions. Keep the old RAM
        % estimate conservative without changing the tested helper API.
        chunk_size = max(1,floor(0.55*chunk_size));
    end
    repeat_count = N - 2;

    if pol ~= 2
        % ---------------------------------------------------------------
        % Original single-polarization path, unchanged in arithmetic.
        % ---------------------------------------------------------------
        refl_work = zeros(1,m);
        for first = 1:chunk_size:m
            last = min(first + chunk_size - 1,m);
            idx  = first:last;
            qvec = Q(idx);
            mm   = numel(idx);

            Rbase = zeros(1,mm);
            Rbase = apply_stack_reflection_vec(qvec,lambda,rho_0r,dz,pol,Rbase);

            [h11,h12,h21,h22] = build_stack_mobius_vec( ...
                qvec,lambda,rho_1r,dz,pol);
            [p11,p12,p21,p22] = mobius_power_vec( ...
                h11,h12,h21,h22,repeat_count);
            Rbase = (p11 .* Rbase + p12) ./ (p21 .* Rbase + p22);
            clear h11 h12 h21 h22 p11 p12 p21 p22

            if ~rough_mode
                R = apply_stack_reflection_vec( ...
                    qvec,lambda,rho_restr,dz,pol,Rbase);
                refl_work(idx) = abs(R).^2;
            else
                if coherent_mode
                    average_value = complex(zeros(1,mm));
                    for kk = 1:numel(rough_stacks)
                        Rn = apply_stack_reflection_vec( ...
                            qvec,lambda,rough_stacks{kk},dz,pol,Rbase);
                        average_value = average_value + weights(kk) .* Rn;
                    end
                    refl_work(idx) = abs(average_value).^2;
                else
                    average_value = zeros(1,mm);
                    for kk = 1:numel(rough_stacks)
                        Rn = apply_stack_reflection_vec( ...
                            qvec,lambda,rough_stacks{kk},dz,pol,Rbase);
                        average_value = average_value + weights(kk) .* abs(Rn).^2;
                    end
                    refl_work(idx) = average_value;
                end
            end
        end
        refl = reshape(refl_work,size(Q));
        return;
    end

    % -------------------------------------------------------------------
    % Exact dual-polarization path. kz and propagation phase are shared.
    % Output order is [pi,sigma] to match the caller's existing naming.
    % -------------------------------------------------------------------
    refl_p_work = zeros(1,m);
    refl_s_work = zeros(1,m);

    for first = 1:chunk_size:m
        last = min(first + chunk_size - 1,m);
        idx  = first:last;
        qvec = Q(idx);
        mm   = numel(idx);

        Rbase_s = zeros(1,mm);
        Rbase_p = zeros(1,mm);
        [Rbase_s,Rbase_p] = apply_stack_reflection_dual_vec( ...
            qvec,lambda,rho_0r,dz,Rbase_s,Rbase_p);

        [s11,s12,s21,s22,p11,p12,p21,p22] = ...
            build_stack_mobius_dual_vec(qvec,lambda,rho_1r,dz);

        [sp11,sp12,sp21,sp22] = mobius_power_vec( ...
            s11,s12,s21,s22,repeat_count);
        [pp11,pp12,pp21,pp22] = mobius_power_vec( ...
            p11,p12,p21,p22,repeat_count);

        Rbase_s = (sp11 .* Rbase_s + sp12) ./ (sp21 .* Rbase_s + sp22);
        Rbase_p = (pp11 .* Rbase_p + pp12) ./ (pp21 .* Rbase_p + pp22);

        clear s11 s12 s21 s22 p11 p12 p21 p22 ...
              sp11 sp12 sp21 sp22 pp11 pp12 pp21 pp22

        if ~rough_mode
            [Rs,Rp] = apply_stack_reflection_dual_vec( ...
                qvec,lambda,rho_restr,dz,Rbase_s,Rbase_p);
            refl_s_work(idx) = abs(Rs).^2;
            refl_p_work(idx) = abs(Rp).^2;
        else
            if coherent_mode
                avg_s = complex(zeros(1,mm));
                avg_p = complex(zeros(1,mm));
                for kk = 1:numel(rough_stacks)
                    [Rns,Rnp] = apply_stack_reflection_dual_vec( ...
                        qvec,lambda,rough_stacks{kk},dz,Rbase_s,Rbase_p);
                    avg_s = avg_s + weights(kk) .* Rns;
                    avg_p = avg_p + weights(kk) .* Rnp;
                end
                refl_s_work(idx) = abs(avg_s).^2;
                refl_p_work(idx) = abs(avg_p).^2;
            else
                avg_s = zeros(1,mm);
                avg_p = zeros(1,mm);
                for kk = 1:numel(rough_stacks)
                    [Rns,Rnp] = apply_stack_reflection_dual_vec( ...
                        qvec,lambda,rough_stacks{kk},dz,Rbase_s,Rbase_p);
                    avg_s = avg_s + weights(kk) .* abs(Rns).^2;
                    avg_p = avg_p + weights(kk) .* abs(Rnp).^2;
                end
                refl_s_work(idx) = avg_s;
                refl_p_work(idx) = avg_p;
            end
        end
    end

    refl = reshape(refl_p_work,size(Q));
    refl_second = reshape(refl_s_work,size(Q));
end

function [rough_stacks,weights] = prepare_crystal_roughness_inputs( ...
    rho_e_rough,substrate_end,nvector,BB)
%PREPARE_CRYSTAL_ROUGHNESS_INPUTS Validate and normalize parent data.
%
% The current parent stores both rho_e_rough and BB using the ABSOLUTE
% integer unit-cell count as the MATLAB index.  Keep that convention here
% so no density-generation changes are required.

    if ~iscell(rho_e_rough)
        error('rho_e_rough must be a cell array.');
    end

    if ~isscalar(substrate_end) || ~isfinite(substrate_end) || ...
            substrate_end < 1 || substrate_end ~= round(substrate_end)
        error('substrate_end must be a positive integer index.');
    end

    nlist = nvector(:).';
    if isempty(nlist) || any(~isfinite(nlist)) || ...
            any(nlist < 1) || any(nlist ~= round(nlist))
        error('nvector must contain positive integer unit-cell counts.');
    end

    if max(nlist) > numel(rho_e_rough)
        error('rho_e_rough does not contain all absolute indices in nvector.');
    end
    if max(nlist) > numel(BB)
        error('BB does not contain all absolute indices in nvector.');
    end

    rough_stacks = cell(1,numel(nlist));
    for kk = 1:numel(nlist)
        n = nlist(kk);
        if isempty(rho_e_rough{n})
            error('rho_e_rough{%d} is empty.',n);
        end
        if substrate_end > numel(rho_e_rough{n})
            error('substrate_end exceeds the length of rho_e_rough{%d}.',n);
        end

        % Same operation as the old roughness propagator, but performed once
        % per realization instead of once per realization per Q chunk.
        rough_stacks{kk} = rho_e_rough{n}(end:-1:substrate_end);
    end

    % The parent constructs B(n) proportional to the Gaussian probability
    % for integer thickness n.  Since the distribution is truncated at
    % approximately +/-3 sigma and sampled only at integers, normalize the
    % discrete weights over the actually included candidates.
    weights = BB(nlist);
    weights = weights(:).';

    if ~isreal(weights) || any(~isfinite(weights)) || any(weights < 0)
        error('BB weights must be finite, real, and nonnegative.');
    end

    sw = sum(weights);
    if ~(isfinite(sw) && sw > 0)
        error('Crystal-thickness weights have zero or invalid total weight.');
    end

    weights = weights ./ sw;
end


function chunk_size = choose_chunk_size(m,rough_mode)
%CHOOSE_CHUNK_SIZE Size Q chunks from currently available physical RAM.
%
% Small and medium Q vectors are always processed in one pass; querying the
% operating system for free RAM would cost more than it could save there.
% For larger Q vectors, use a conservative peak-memory estimate and 15% of
% currently available physical RAM.

    fast_path_limit = 50000;
    if m <= fast_path_limit
        chunk_size = m;
        return;
    end

    available_bytes = get_available_ram_bytes();

    if rough_mode
        % Rbase + ensemble accumulator/Rn + Mobius build/power temporaries.
        bytes_per_q = 1152;
    else
        % Mobius build/power temporaries + final reflection stack.
        bytes_per_q = 960;
    end

    ram_fraction = 0.15;
    budget_bytes = ram_fraction * available_bytes;

    chunk_size = floor(budget_bytes / bytes_per_q);

    % Always make progress.  For the usual small/medium Q vectors, this
    % simply selects all Q points and therefore adds no chunking overhead.
    chunk_size = max(1,min(m,chunk_size));
end


function available_bytes = get_available_ram_bytes()
%GET_AVAILABLE_RAM_BYTES Best-effort cross-platform free-RAM query.
%
% Windows: MATLAB memory().
% macOS:   vm_stat, counting immediately reclaimable page classes.
% Linux:   /proc/meminfo MemAvailable.
%
% Cache the result briefly so repeated fitting calls with a very large Q do
% not launch an OS query every time.  If the OS query fails, fall back to
% 512 MiB.

    persistent cached_bytes cached_time

    tnow = now;
    if ~isempty(cached_bytes) && ~isempty(cached_time) && ...
            (tnow - cached_time)*86400 < 5
        available_bytes = cached_bytes;
        return;
    end

    available_bytes = NaN;

    if ispc
        try
            [~,systemview] = memory;
            available_bytes = systemview.PhysicalMemory.Available;
        catch
            available_bytes = NaN;
        end

    elseif ismac
        try
            [status,txt] = system('vm_stat');
            if status == 0
                tok = regexp(txt,'page size of\s+([0-9]+)\s+bytes', ...
                    'tokens','once');
                if isempty(tok)
                    page_size = 4096;
                else
                    page_size = str2double(tok{1});
                end

                pages_free        = vm_stat_pages(txt,'Pages free');
                pages_inactive    = vm_stat_pages(txt,'Pages inactive');
                pages_speculative = vm_stat_pages(txt,'Pages speculative');
                pages_purgeable   = vm_stat_pages(txt,'Pages purgeable');

                available_bytes = page_size * (pages_free + pages_inactive + ...
                    pages_speculative + pages_purgeable);
            end
        catch
            available_bytes = NaN;
        end

    elseif isunix
        try
            fid = fopen('/proc/meminfo','r');
            if fid >= 0
                txt = fread(fid,'*char')';
                fclose(fid);
                tok = regexp(txt,'MemAvailable:\s*([0-9]+)\s*kB', ...
                    'tokens','once');
                if ~isempty(tok)
                    available_bytes = str2double(tok{1}) * 1024;
                end
            end
        catch
            available_bytes = NaN;
        end
    end

    if ~(isfinite(available_bytes) && available_bytes > 0)
        available_bytes = 512 * 1024^2;
    end

    cached_bytes = available_bytes;
    cached_time  = tnow;
end


function pages = vm_stat_pages(txt,label)
%VM_STAT_PAGES Extract one page count from macOS vm_stat output.

    tok = regexp(txt,[label ':\s*([0-9]+)\.'],'tokens','once');
    if isempty(tok)
        pages = 0;
    else
        pages = str2double(tok{1});
    end
end


function R = apply_stack_reflection_vec(Q,lambda,rho_e,dz,pol,R)
%APPLY_STACK_REFLECTION_VEC Apply one finite stack directly to R.
%
% If A=[A1 A2; A3 A4] left-multiplies the existing transfer matrix, then
%
%       Rnew = (A3 + A4*R)/(A1 + A2*R).
%
% With the A matrix used in the reference code this reduces exactly to
%
%                 r + q R
%       Rnew = ---------------,       q = exp(-i*kz_cur*dz)
%                1 + r q R
%
% The layer order below is nn:-1:2, matching
% do_matrix_propagation_optimzed_vec_winograd exactly.

    re = 2.814042735053330e-05;
    k0 = 2*pi/lambda;
    factor = 8 * k0^2 * lambda^2 / (2*pi);
    phase_factor = -1i * dz;

    nn = length(rho_e);
    if nn < 2
        return;
    end

    % Preserve the conjugating transpose in the original code verbatim.
    sld   = rho_e' * re / (2*pi);
    delta = factor * real(sld);
    beta  = 1i * factor * imag(sld);

    % Preserve Q(:)' from the original as well.
    Qsqr = (Q(:)').^2;

    if pol == 1
        % Identical pi-polarization refractive-index convention to the
        % supplied make_A_matrix_vec.
        deltan = lambda^2 * real(sld) / (2*pi);
        betan  = 1i * lambda^2 * imag(sld) / (2*pi);
        nlayer = 1 - deltan - betan;
        n2 = nlayer.^2;
    end

    % IMPORTANT: the original expression
    %
    %   kz_all = (sqrt(... )').';
    %
    % conjugates sqrt(...).  conj(sqrt(...)) is its row-wise equivalent.
    kz_cur = conj(sqrt(Qsqr - delta(nn) + beta(nn)));

    for j = nn:-1:2
        kz_prev = conj(sqrt(Qsqr - delta(j-1) + beta(j-1)));

        if pol == 0
            r = (kz_prev - kz_cur) ./ (kz_prev + kz_cur);
        else
            fact1 = n2(j-1) .* kz_cur;
            fact2 = n2(j)   .* kz_prev;
            r = (fact1 - fact2) ./ (fact1 + fact2);
        end

        qphase = exp(phase_factor .* kz_cur);

        % Fused form of the reflection recurrence.  This is the same hot
        % loop as the tested optimized non-roughness version.
        x = qphase .* R;
        R = (r + x) ./ (1 + r .* x);

        kz_cur = kz_prev;
    end
end


function [h11,h12,h21,h22] = build_stack_mobius_vec(Q,lambda,rho_e,dz,pol)
%BUILD_STACK_MOBIUS_VEC Build the full reflection map of the repeated cell.

    re = 2.814042735053330e-05;
    k0 = 2*pi/lambda;
    factor = 8 * k0^2 * lambda^2 / (2*pi);
    phase_factor = -1i * dz;

    nn = length(rho_e);
    mm = numel(Q);

    h11 = ones(1,mm);
    h12 = zeros(1,mm);
    h21 = zeros(1,mm);
    h22 = ones(1,mm);

    if nn < 2
        return;
    end

    sld   = rho_e' * re / (2*pi);
    delta = factor * real(sld);
    beta  = 1i * factor * imag(sld);
    Qsqr  = (Q(:)').^2;

    if pol == 1
        deltan = lambda^2 * real(sld) / (2*pi);
        betan  = 1i * lambda^2 * imag(sld) / (2*pi);
        nlayer = 1 - deltan - betan;
        n2 = nlayer.^2;
    end

    kz_cur = conj(sqrt(Qsqr - delta(nn) + beta(nn)));

    normalize_every = 32;
    layer_counter   = 0;

    for j = nn:-1:2
        kz_prev = conj(sqrt(Qsqr - delta(j-1) + beta(j-1)));

        if pol == 0
            r = (kz_prev - kz_cur) ./ (kz_prev + kz_cur);
        else
            fact1 = n2(j-1) .* kz_cur;
            fact2 = n2(j)   .* kz_prev;
            r = (fact1 - fact2) ./ (fact1 + fact2);
        end

        qphase = exp(phase_factor .* kz_cur);
        rq     = r .* qphase;

        n11 = qphase .* h11 + r  .* h21;
        n12 = qphase .* h12 + r  .* h22;
        n21 = rq     .* h11      + h21;
        n22 = rq     .* h12      + h22;

        h11 = n11;
        h12 = n12;
        h21 = n21;
        h22 = n22;

        layer_counter = layer_counter + 1;
        if mod(layer_counter,normalize_every) == 0
            [h11,h12,h21,h22] = normalize_mobius_vec( ...
                h11,h12,h21,h22);
        end

        kz_cur = kz_prev;
    end

    [h11,h12,h21,h22] = normalize_mobius_vec( ...
        h11,h12,h21,h22);
end


function [r11,r12,r21,r22] = mobius_power_vec(a11,a12,a21,a22,k)
%MOBIUS_POWER_VEC Vectorized exponentiation by squaring over all Q values.

    mm = numel(a11);

    r11 = ones(1,mm);
    r12 = zeros(1,mm);
    r21 = zeros(1,mm);
    r22 = ones(1,mm);

    if k <= 0
        return;
    end

    [a11,a12,a21,a22] = normalize_mobius_vec( ...
        a11,a12,a21,a22);

    while k > 0
        if mod(k,2) == 1
            [r11,r12,r21,r22] = mobius_multiply_vec( ...
                r11,r12,r21,r22,a11,a12,a21,a22);
        end

        k = floor(k/2);

        if k > 0
            tr  = a11 + a22;
            off = a12 .* a21;

            b11 = a11 .* a11 + off;
            b12 = a12 .* tr;
            b21 = a21 .* tr;
            b22 = a22 .* a22 + off;

            [a11,a12,a21,a22] = normalize_mobius_vec( ...
                b11,b12,b21,b22);
        end
    end
end


function [c11,c12,c21,c22] = mobius_multiply_vec( ...
    a11,a12,a21,a22,b11,b12,b21,b22)
%MOBIUS_MULTIPLY_VEC Vectorized 2x2 multiplication over Q.

    c11 = a11 .* b11 + a12 .* b21;
    c12 = a11 .* b12 + a12 .* b22;
    c21 = a21 .* b11 + a22 .* b21;
    c22 = a21 .* b12 + a22 .* b22;

    [c11,c12,c21,c22] = normalize_mobius_vec( ...
        c11,c12,c21,c22);
end


function [a11,a12,a21,a22] = normalize_mobius_vec(a11,a12,a21,a22)
%NORMALIZE_MOBIUS_VEC Remove a common scale at each Q point.

    scale = max(max(abs(a11),abs(a12)),max(abs(a21),abs(a22)));
    scale(scale == 0) = 1;

    a11 = a11 ./ scale;
    a12 = a12 ./ scale;
    a21 = a21 ./ scale;
    a22 = a22 ./ scale;
end

function [Rs,Rp] = apply_stack_reflection_dual_vec(Q,lambda,rho_e,dz,Rs,Rp)
%APPLY_STACK_REFLECTION_DUAL_VEC Propagate sigma and pi together.
%
% Optimized for repeated use in fitting:
%   - sqrt() evaluated once over all Q/layer combinations
%   - exp()  evaluated once over all Q/layer combinations
%   - original reflection recurrence retained because profiling showed
%     that the algebraically fused version is slower in MATLAB
%   - Q dimension is stored along rows of each kz column so that
%     kz(:,j) is contiguous in memory
%
% Inputs:
%   Q       - momentum transfer values
%   lambda  - scalar wavelength
%   rho_e   - layer electron density profile
%   dz      - scalar layer spacing
%   Rs, Rp  - input reflection amplitudes
%
% Outputs:
%   Rs, Rp  - propagated reflection amplitudes

    re = 2.81404273505333e-05;

    nn = numel(rho_e);

    % Original loop is empty for nn < 3.
    if nn < 3
        return
    end

    % Preserve caller's Rs/Rp shapes.
    sizeRs = size(Rs);
    sizeRp = size(Rp);

    % Work with column vectors internally.
    Q  = Q(:);
    Rs = Rs(:);
    Rp = Rp(:);

    % ---------------------------------------------------------------
    % Layer-dependent quantities
    % ---------------------------------------------------------------

    % Preserve the original rho_e' conjugation behavior.
    sld = rho_e' .* (re / (2*pi));

    % Keep the original expression rather than changing numerical
    % evaluation order unnecessarily.
    factor = 8*(2*pi/lambda)^2 * lambda^2/(2*pi);

    delta = factor .* real(sld);
    beta  = 1i .* factor .* imag(sld);

    deltan = lambda^2 .* real(sld) ./ (2*pi);
    betan  = 1i .* lambda^2 .* imag(sld) ./ (2*pi);

    n2 = (1 - deltan - betan).^2;

    % ---------------------------------------------------------------
    % Precompute kz for every Q and layer.
    %
    % Qsqr:    [nQ x 1]
    % kz_shift:[1 x nn]
    % kz:      [nQ x nn]
    %
    % For Q=4225 and nn=690, kz is about 47 MB complex double.
    % ---------------------------------------------------------------

    Qsqr = Q.^2;

    kz_shift = (-delta + beta).';

    kz = conj(sqrt(Qsqr + kz_shift));

    % ---------------------------------------------------------------
    % Precompute phase factors.
    %
    % This was the largest speedup in profiling:
    % exp() went from ~107 s to ~15 s.
    % ---------------------------------------------------------------

    phase_factor = -1i .* dz;

    qphase_all = exp(phase_factor .* kz);

    % ---------------------------------------------------------------
    % Reflection recurrence
    %
    % IMPORTANT:
    % Preserve the indexing of the original function exactly.
    %
    % Initial current layer:
    %       nn
    %
    % First previous layer:
    %       nn-1
    % ---------------------------------------------------------------

    kz_cur = kz(:,nn);
    qphase = qphase_all(:,nn);

    for j = nn:-1:2

        kz_prev = kz(:,j-1);

        % -----------------------------------------------------------
        % Sigma / s polarization
        % -----------------------------------------------------------

        rs = (kz_prev - kz_cur) ./ ...
             (kz_prev + kz_cur);

        xs = qphase .* Rs;

        Rs = (rs + xs) ./ ...
             (1 + rs .* xs);

        % -----------------------------------------------------------
        % Pi / p polarization
        % -----------------------------------------------------------

        fact1 = n2(j-1) .* kz_cur;
        fact2 = n2(j)   .* kz_prev;

        rp = (fact1 - fact2) ./ ...
             (fact1 + fact2);

        xp = qphase .* Rp;

        Rp = (rp + xp) ./ ...
             (1 + rp .* xp);

        % -----------------------------------------------------------
        % Move to next interface.
        % -----------------------------------------------------------

        kz_cur = kz_prev;

        % No need to extract the phase after the final iteration.
        if j > 2
            qphase = qphase_all(:,j-1);
        end
    end

    % Restore original input shapes.
    Rs = reshape(Rs,sizeRs);
    Rp = reshape(Rp,sizeRp);
end


function [s11,s12,s21,s22,p11,p12,p21,p22] = ...
    build_stack_mobius_dual_vec(Q,lambda,rho_e,dz)
%BUILD_STACK_MOBIUS_DUAL_VEC Build sigma and pi maps in one slice traversal.

    re = 2.814042735053330e-05;
    factor = 16*pi;
    phase_factor = -1i * dz;

    nn = length(rho_e);
    mm = numel(Q);

    s11 = ones(1,mm); s12 = zeros(1,mm); s21 = zeros(1,mm); s22 = ones(1,mm);
    p11 = ones(1,mm); p12 = zeros(1,mm); p21 = zeros(1,mm); p22 = ones(1,mm);

    if nn < 2
        return;
    end

    sld   = rho_e' * re / (2*pi);
    delta = factor * real(sld);
    beta  = 1i * factor * imag(sld);
    Qsqr  = (Q(:)').^2;

    deltan = lambda^2 * real(sld) / (2*pi);
    betan  = 1i * lambda^2 * imag(sld) / (2*pi);
    n2 = (1 - deltan - betan).^2;

    kz_cur = conj(sqrt(Qsqr - delta(nn) + beta(nn)));

    normalize_every = 32;
    layer_counter = 0;

    for j = nn:-1:2
        kz_prev = conj(sqrt(Qsqr - delta(j-1) + beta(j-1)));

        rs = (kz_prev - kz_cur) ./ (kz_prev + kz_cur);
        fact1 = n2(j-1) .* kz_cur;
        fact2 = n2(j)   .* kz_prev;
        rp = (fact1 - fact2) ./ (fact1 + fact2);

        qphase = exp(phase_factor .* kz_cur);
        rqs = rs .* qphase;
        rqp = rp .* qphase;

        ns11 = qphase .* s11 + rs .* s21;
        ns12 = qphase .* s12 + rs .* s22;
        ns21 = rqs    .* s11 + s21;
        ns22 = rqs    .* s12 + s22;

        np11 = qphase .* p11 + rp .* p21;
        np12 = qphase .* p12 + rp .* p22;
        np21 = rqp    .* p11 + p21;
        np22 = rqp    .* p12 + p22;

        s11 = ns11; s12 = ns12; s21 = ns21; s22 = ns22;
        p11 = np11; p12 = np12; p21 = np21; p22 = np22;

        layer_counter = layer_counter + 1;
        if mod(layer_counter,normalize_every) == 0
            [s11,s12,s21,s22] = normalize_mobius_vec(s11,s12,s21,s22);
            [p11,p12,p21,p22] = normalize_mobius_vec(p11,p12,p21,p22);
        end

        kz_cur = kz_prev;
    end

    [s11,s12,s21,s22] = normalize_mobius_vec(s11,s12,s21,s22);
    [p11,p12,p21,p22] = normalize_mobius_vec(p11,p12,p21,p22);
end
