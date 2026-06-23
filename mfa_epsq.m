function mfa = mfa_epsq(run, var, varargin)
% MFA_EPSQ  Extract bifractal parameters (H1, C1) from structure functions.
%
%   mfa = mfa_epsq(RUN, VAR) loads PSD scale breaks and structure function
%   output for the specified RUN and variable VAR ('c' or 'r'), fits
%   scaling exponents K(q) within the diagnosed scaling regime, and
%   computes the bifractal parameters:
%
%     H1 = K(q=1)     (smoothness / Hurst exponent)
%     C1 = dK/dq|q=1  (intermittency, via finite difference)
%
%   The Pierrehumbert (1996) hyperbolic model zeta(q) = a*q/(1 + a*q/z_inf)
%   is also fit to obtain (a, z_inf), from which alternative estimates
%   h1b and c1b are derived via Eqs. 4-5 of Witte et al. (2022, JAS).
%
%   mfa = mfa_epsq(RUN, VAR, fd) reads from directory FD.
%
%   mfa = mfa_epsq(RUN, VAR, fd, sf_mode) specifies structure function
%   mode: '1d' (default) loads sf_q<VAR>_<RUN>.mat, '2d' loads
%   sf2d_q<VAR>_<RUN>.mat. Both formats are compatible with the fitting
%   stage.
%
%   Output structure fields (appended to the loaded SF data):
%     mfa.lev.Kq  - scaling exponents [nq x nz x nt]
%     mfa.lev.h1  - H1 = K(q=1)                    [nz x nt]
%     mfa.lev.c1  - C1 = (K(1.1) - K(0.9)) / 0.2   [nz x nt]
%     mfa.lev.a   - Pierrehumbert parameter a        [nz x nt]
%     mfa.lev.zi  - Pierrehumbert parameter z_inf    [nz x nt]
%     mfa.lev.h1b - H1 from Pierrehumbert fit        [nz x nt]
%     mfa.lev.c1b - C1 from Pierrehumbert fit        [nz x nt]

%% Parse input
fd      = './';
sf_mode = '1d';
if nargin >= 3, fd      = varargin{1}; end
if nargin >= 4, sf_mode = varargin{2}; end
if nargin >= 5, outputd = varargin{3}; end


%% Build filenames and set grid spacing
switch sf_mode
    case '1d', sf_prefix = 'sf_q';
    case '2d', sf_prefix = 'sf2d_q';
    otherwise, error('sf_mode must be ''1d'' or ''2d''')
end

switch var
    case {'c', 'l', 'r'}
        sb_str = sprintf('psd_q%s_%s.mat', var, run);
        sf_str = sprintf('%s%s_%s.mat', sf_prefix, var, run);
    otherwise
        error('Variable %s not configured', var)
end

% Grid spacing: read from PSD file if available, otherwise set manually
psd_data = load(fullfile(outputd, sb_str), 'sb', 'r');
if isfield(psd_data, 'r')
    % Infer dx from the PSD length-scale array and normalised wavenumber
    % For now, use a fixed value consistent with the CM1 configuration.
    % TODO: store dx in the PSD file and read it here.
    dx = 0.05;  % km (50 m grid spacing)
else
    dx = 0.05;
end

sb = psd_data.sb;

%% Load structure function data
if ~exist(fullfile(outputd, sf_str), 'file')
    error('Structure function file %s not found', sf_str)
end
mfa = load(fullfile(outputd, sf_str));

%% Scale breaks in log2 units relative to grid spacing
l2sb = log2(sb / dx);

%% Structure function orders (must match sf_epsq_les_1step)
q  = [0.9, 1.0, 1.1];
nq = length(q);

%% Pierrehumbert (1996) model: zeta(q) = a*q / (1 + a*q/z_inf)
zfun = @(b, x) b(1) .* x ./ (1 + b(1) .* x ./ b(2));
b0   = [0.5, 2];   % initial guess for [a, z_inf]

%% Get array dimensions
szi = size(mfa.lev.sf);
nr  = szi(1);
nz  = szi(3);
if length(szi) == 4
    nt = szi(4);
elseif length(szi) == 3
    nt = 1;
else
    error('Unexpected SF array dimensions: %s', mat2str(szi))
end

%% Allocate output
Kq   = NaN(nq, nz, nt);
zinf = NaN(nz, nt);
a    = NaN(nz, nt);

%% Main analysis loop
for j = 1:nt
    for k = 1:nz

        % Extract SF and singularity measures at this level/time
        if nt == 1
            sfjk  = squeeze(mfa.lev.sf(:, :, k));
            epsqk = squeeze(mfa.lev.epsq(:, :, k));
            sb_k  = l2sb(k, :);
        else
            sfjk  = squeeze(mfa.lev.sf(:, :, k, j));
            epsqk = squeeze(mfa.lev.epsq(:, :, k, j));
            sb_k  = l2sb(k, :, j);
        end

        % Skip if no scale breaks or no valid S_1 data
        if any(isnan(sb_k)) || all(isnan(sfjk(:, q == 1)))
            continue
        end

        % Fit K(q) for each order
        for m = 1:nq
            if q(m) == 1
                % Integer order: fit log2(S_q) vs log2(r) directly
                r  = log2(1:nr)';
                ia = find(abs(r - sb_k(2)) == min(abs(r - sb_k(2))), 1, 'first');
                iz = find(abs(r - sb_k(1)) == min(abs(r - sb_k(1))), 1, 'first');
                if isempty(ia) || isempty(iz) || ia >= iz
                    continue
                end
                fit_K = polyfit(r(ia:iz), log2(sfjk(ia:iz, m)), 1);
                Kq(m, k, j) = fit_K(1);
            else
                % Non-integer order: fit -log2(eps_q) vs dyadic scale
                r  = (0 : size(epsqk, 1) - 1)';
                ia = find(abs(r - sb_k(2)) == min(abs(r - sb_k(2))), 1, 'first');
                iz = find(abs(r - sb_k(1)) == min(abs(r - sb_k(1))), 1, 'first');
                if isempty(ia) || isempty(iz) || ia >= iz
                    continue
                end
                fit_K = polyfit(r(ia:iz), -log2(epsqk(ia:iz, m)), 1);
                Kq(m, k, j) = fit_K(1);
            end
        end

        % Pierrehumbert nonlinear fit (requires all 3 K(q) values)
        Kq_vec = squeeze(Kq(:, k, j));
        if any(isnan(Kq_vec))
            continue
        end
        try
            warning('off', 'all')
            fit_z = fitnlm(q(:), Kq_vec(:), zfun, b0);
            warning('on', 'all')
            a(k, j)    = fit_z.Coefficients.Estimate(1);
            zinf(k, j) = fit_z.Coefficients.Estimate(2);
        catch
            warning('on', 'all')
            % fitnlm can fail if K(q) is non-monotonic or poorly conditioned
        end
    end
end

%% Compute bifractal parameters
% Method 1: finite-difference estimate of C1
%   H1 = K(q=1)
%   C1 = dK/dq |_{q=1} ≈ (K(1.1) - K(0.9)) / (1.1 - 0.9)
%
%   *** CRITICAL: the parenthesisation below is intentional.
%   *** The previous version had an operator-precedence bug:
%   ***   K(1.1) - K(0.9)/(1.1-0.9)  [WRONG: divides only the 2nd term]
%   *** Correct form:
%   ***   (K(1.1) - K(0.9)) / (1.1 - 0.9)

mfa.lev.Kq = Kq;
mfa.lev.h1 = squeeze(Kq(q == 1, :, :));
mfa.lev.c1 = squeeze((Kq(q == 1.1, :, :) - Kq(q == 0.9, :, :)) / (1.1 - 0.9));

% Mask out levels where K(q=2) would be negative (if available) or
% where the fit is clearly unphysical
% (With only 3 q-values we don't have K(q=2), so skip this check.)

% Method 2: Pierrehumbert fit -> (H1, C1) via Eqs. 4-5 of Witte et al. (2022)
%   H1 = a * z_inf / (a + z_inf)
%   C1 = z_inf / (z_inf/a + 1)^2
%
%   NOTE: with only 3 data points and 2 free parameters, the Pierrehumbert
%   fit has 1 degree of freedom and may be poorly constrained. The
%   finite-difference estimates (h1, c1) are the primary output; the
%   Pierrehumbert estimates (h1b, c1b) serve as a cross-check.
mfa.lev.a   = a;
mfa.lev.zi  = zinf;
mfa.lev.h1b = a .* zinf ./ (a + zinf);
mfa.lev.c1b = zinf ./ (zinf ./ a + 1).^2;

fprintf('MFA complete: %i levels, %i timesteps\n', nz, nt)

end
