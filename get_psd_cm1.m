function psd1 = get_psd_cm1(in, varargin)
% GET_PSD_CM1  Compute power spectral density of a 3D field at each level.
%
%   psd1 = get_psd_cm1(IN) computes the Welch PSD of the 3D array IN
%   (dimensions: x, y, z for a single timestep) along the x-dimension,
%   averaged over all y-transects. Returns a struct with fields:
%       avgz.E  - mean PSD at each vertical level [nf x nz]
%       w       - normalised angular frequency from pwelch
%       xn      - normalised length scale (pi./w)
%       r       - physical length scale (xn * dx), if dx is provided
%
%   psd1 = get_psd_cm1(IN, qlprof) restricts processing to levels where
%   QLPROF (domain-mean cloud water profile) exceeds 1e-5 kg/kg.
%
%   psd1 = get_psd_cm1(IN, qlprof, dx) also computes physical length
%   scales from grid spacing DX.

%% Parse dimensions and optional arguments
[nx, ny, nz] = size(in);

ql = [];
dx = [];
if nargin >= 2,  ql = varargin{1};  end
if nargin >= 3,  dx = varargin{2};  end

%% Determine highest level to analyse
if ~isempty(ql)
    kf = find(ql > 1e-5, 1, 'last');
    if isempty(kf), kf = 0; end
else
    kf = nz;
end

%% Allocate PSD output
%  pwelch with a Hann window of length nx and default nfft = max(256, 2^nextpow2(nx))
%  returns floor(nfft/2)+1 frequency points.
nfft = max(256, 2^nextpow2(nx));
nf   = floor(nfft/2) + 1;

psd1 = struct();
psd1.avgz.E = NaN(nf, nz);

%% Compute PSD level-by-level, averaged over y-transects
for k = 1:kf
    in_k = squeeze(in(:, :, k));
    if k == 1
        [psd1_x, psd1.w] = pwelch(in_k, hann(nx));
    else
        psd1_x = pwelch(in_k, hann(nx));
    end
    if any(isnan(psd1_x(:)))
        warning('NaN detected in PSD at level %i', k)
    end
    psd1.avgz.E(:, k) = mean(psd1_x, 2, 'omitnan');
end

%% Compute length scales
psd1.xn = pi ./ psd1.w;
if ~isempty(dx)
    psd1.r = psd1.xn * dx;
end

end
