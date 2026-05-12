function psd1 = psd_slope_fit_les(psd1, ql, qcflag)
% PSD_SLOPE_FIT_LES  Diagnose spectral slope and scaling regime bounds.
%
%   psd1 = psd_slope_fit_les(psd1, QL, QCFLAG) takes the PSD struct from
%   get_psd_cm1, a domain-mean cloud water profile QL, and a flag QCFLAG
%   (1 = cloud variable, 0 = other) and fits log-log spectral slopes at
%   each vertical level. The scaling regime is identified by optimising a
%   cost function that balances R^2, compensated spectral slope, and RMSE.
%
%   The lower bound of the scaling regime is fixed at r_min = 6*dx,
%   following Bryan et al. (2003) for the effective resolution of 3rd-order
%   advection schemes.
%
%   Adds fields to psd1.avgz:
%       slope  - spectral exponent beta (positive, i.e. -slope of loglog PSD)
%       r2     - coefficient of determination for the fit
%       sb     - scale break bounds [upper, lower] in physical units
%       sbn    - scale break bounds as frequency-vector indices

%% Configuration
ct = find(ql > 1e-5, 1, 'last');  % highest level with cloud
if isempty(ct), ct = 0; end

xn = pi ./ psd1.w;   % normalised length scale (wavelengths / dx)
nf = length(xn);

% Minimum span of scaling regime (in decades of wavelength)
if nf < 256
    nmin  = 5;
    dxmin = 0.1;
else
    nmin  = 5;
    dxmin = 0.5;
end

%% Only process the avgz field
fld = 'avgz';
if ~isfield(psd1, fld), return; end

nz = size(psd1.(fld).E, 2);

% Allocate output arrays
psd1.(fld).slope  = NaN(nz, 1);
psd1.(fld).r2     = NaN(nz, 1);
psd1.(fld).sb     = NaN(nz, 2);
psd1.(fld).sbn    = NaN(nz, 2);
psd1.(fld).slope2 = NaN(nz, 1);

xmj = psd1.r;     % physical length scale
fmj = psd1.w;     % angular frequency

%% Loop over vertical levels
for j = 1:ct
    % Skip levels already fitted, or below-cloud levels for cloud variables
    if ~isnan(psd1.(fld).slope(j)) || (qcflag && ql(j) < 1e-5)
        continue
    end

    Emj = psd1.(fld).E(:, j);
    nj  = length(Emj);

    % Lower bound of scaling regime: r/dx = 6 (Bryan et al. 2003)
    kl = find(abs(xn - 6) == min(abs(xn - 6)), 1);

    % Upper bound search: from index 2 downward to kl+nmin
    kf = 2;

    slopes  = NaN(kl - kf + 1, 1);
    r2s     = slopes;
    dx_span = slopes;
    slopes2 = slopes;
    rmse    = slopes;

    for ki = kf : kl + nmin
        test = log10(xmj(ki) / xmj(kl));
        if test < dxmin || kl - ki < nmin
            continue
        end
        [slope_i, ~, ~, r2_i, ~] = logfit(fmj(ki:kl), Emj(ki:kl), 'loglog');
        slopes(ki)  = slope_i;
        r2s(ki)     = r2_i;
        dx_span(ki) = test;

        [slopes2(ki), ~, MSE, ~, ~] = logfit(log10(fmj(ki:kl)), ...
            Emj(ki:kl) .* fmj(ki:kl).^(-slope_i), 'linear');
        rmse(ki) = sqrt(MSE);
    end

    % Cost function: balance R2, compensated slope flatness, and RMSE
    dist1 = sqrt((1 - r2s / max(r2s(:))).^2 + ...
                 (1 ./ log10(abs(slopes2))).^2 + ...
                 (1 ./ log10(rmse)).^2);
    mind1 = find(dist1 == min(dist1(:)));

    first = mind1;
    last  = kl;

    % Quality check: reject if no unique minimum or R2 < 0.90
    if length(first) ~= 1 || length(last) ~= 1 || ...
       r2s(dist1 == min(dist1(:))) < 0.90
        continue
    end

    psd1.(fld).sb(j, 1)  = xmj(first);
    psd1.(fld).sb(j, 2)  = xmj(last);
    psd1.(fld).sbn(j, 1) = first;
    psd1.(fld).sbn(j, 2) = last;
    psd1.(fld).slope(j)  = -slopes(first);
    psd1.(fld).r2(j)     = r2s(first);
    psd1.(fld).slope2(j) = -slopes2(first);
end

end
