%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:           Gunnar Palsson / reflection-recursion rewrite
% Timestamp:        2026-08-14 14:04 +02:00
% Description:      One merged optimized propagator for both:
%                   (1) the normal single-stack calculation, and
%                   (2) top-layer crystal-thickness averaging.
%
% NORMAL CALL -- unchanged parent signature:
%   refl = propagate_vectorized_chunks_opt(Q,lambda,rho_0r,rho_1r,N, ...
%                                          rho_restr,dz,sigma,pol)
%
% TOP-LAYER CRYSTAL-ROUGHNESS CALL:
%   refl = propagate_vectorized_chunks_opt(Q,lambda,rho_0r,rho_1r,N, ...
%                                          rho_e_rough,dz,substrate_end, ...
%                                          nvector,BB,pol)
%
% Optional final argument for top-layer averaging, intended to be passed
% from control.roughness_averaging:
%   ...,pol,control.roughness_averaging
% where the value is 'coherent' or 'incoherent'.
%
% The normal branch below is the already-tested optimized reflection
% recursion.  The roughness branch computes the same common substrate +
% repeated-cell reflection once, applies each candidate top-layer stack to
% that complex reflection amplitude, and only then performs the ensemble
% average.
%
% CRYSTAL-THICKNESS WEIGHTING
% ---------------------------
% parratt_matrix_repeated_rhobuiltin constructs BB(n) as the Gaussian
% weight for the integer unit-cell count n.  Therefore the statistically
% correct weights are BB(n), normalized over nvector.  Do NOT take
% sqrt(BB) before averaging.
%
% coherent:
%   refl(Q) = abs(sum_n w_n * R_n(Q)).^2
%
% incoherent:
%   refl(Q) = sum_n w_n * abs(R_n(Q)).^2
%
% where sum_n w_n = 1.
%
% This also fixes the old roughness routine's B2/Btemp remainder bug: every
% candidate uses the reflection after applying its B3/top-layer stack.
%
% IMPORTANT NUMERICAL-COMPATIBILITY DETAILS
% -----------------------------------------
% Preserve both conjugations present in the reference transfer code:
%
%   sld = rho_e' * re/(2*pi)
%   kz_all = (sqrt(... )').';
%
% The old microscopic/interface sigma roughness is deliberately inactive;
% sigma remains accepted only for compatibility with the normal call.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refl = propagate_vectorized_chunks_opt(Q,lambda,rho_0r,rho_1r,N,rho_arg,dz,varargin)

    % ---------------------------------------------------------------
    % Dispatch only on the existing argument layouts.
    % ---------------------------------------------------------------
    if numel(varargin) == 2
        % Existing optimized call:
        %   ...,rho_restr,dz,sigma,pol
        rough_mode = false;
        rho_restr  = rho_arg;
        sigma      = varargin{1}; %#ok<NASGU>
        pol        = varargin{2};
        average_mode = '';
        coherent_mode = false;

    elseif numel(varargin) == 4 || numel(varargin) == 5
        % Existing roughness data, now handled by this same function:
        %   ...,rho_e_rough,dz,substrate_end,nvector,BB,pol[,mode]
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

        % Resolve the averaging mode ONCE.  There is no string comparison
        % inside the roughness-realization loop below.
        coherent_mode = strcmp(average_mode,'coherent');

        % Pre-extract/reverse every top-stack realization once.  This avoids
        % repeating rho_full(end:-1:substrate_end) for every Q chunk.
        [rough_stacks,weights] = prepare_crystal_roughness_inputs( ...
            rho_e_rough,substrate_end,nvector,BB);

    else
        error(['Unsupported call signature. Use either the normal 9-input ' ...
               'call or the top-layer crystal-roughness 11/12-input call.']);
    end

    if pol ~= 0 && pol ~= 1
        error('pol must be 0 (sigma) or 1 (pi).');
    end

    m = numel(Q);
    if m == 0
        refl = zeros(size(Q));
        return;
    end

    % Choose the largest Q chunk that fits comfortably in currently
    % available RAM.  For ordinary Q vectors this normally returns m, so
    % there is no chunking overhead.  Chunking only activates for very large
    % Q grids or low-memory situations.
    chunk_size   = choose_chunk_size(m,rough_mode);
    refl_work    = zeros(1,m);
    repeat_count = N - 2;

    for first = 1:chunk_size:m
        last = min(first + chunk_size - 1,m);
        idx  = first:last;
        qvec = Q(idx);
        mm   = numel(idx);

        % ---------------------------------------------------------------
        % Common base reflection -- identical arithmetic to the tested
        % optimized non-roughness implementation:
        %
        %       Bbase = B1^(N-2) * B0
        % ---------------------------------------------------------------
        Rbase = zeros(1,mm);
        Rbase = apply_stack_reflection_vec(qvec,lambda,rho_0r,dz,pol,Rbase);

        [h11,h12,h21,h22] = build_stack_mobius_vec( ...
            qvec,lambda,rho_1r,dz,pol);

        [p11,p12,p21,p22] = mobius_power_vec( ...
            h11,h12,h21,h22,repeat_count);

        Rbase = (p11 .* Rbase + p12) ./ (p21 .* Rbase + p22);

        % These eight full-Q arrays are no longer needed.  Releasing them
        % here lowers the peak memory of the top-layer ensemble branch.
        clear h11 h12 h21 h22 p11 p12 p21 p22

        if ~rough_mode
            % -----------------------------------------------------------
            % NORMAL BRANCH
            % Exactly the already-tested optimized final-stack operation.
            % -----------------------------------------------------------
            R = apply_stack_reflection_vec( ...
                qvec,lambda,rho_restr,dz,pol,Rbase);

            refl_work(idx) = abs(R).^2;

        else
            % -----------------------------------------------------------
            % TOP-LAYER CRYSTAL-THICKNESS ENSEMBLE
            %
            % For every candidate thickness n:
            %
            %   old transfer notation: Btemp = B3_n * Bbase
            %   reflection notation:   Rn = apply(B3_n,Rbase)
            %
            % This is where the old remainder branch accidentally used B2
            % instead of Btemp.  There is no separate B2/Btemp object here,
            % so all Q points necessarily include B3_n.
            % -----------------------------------------------------------
            % Dispatch ONCE per Q chunk.  Keeping two tight loops is faster
            % than testing a string or mode flag for every realization.
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
