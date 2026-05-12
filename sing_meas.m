function epsq = sing_meas(var, q, varargin)
% SING_MEAS  Coarse-grained singularity measure for multifractal analysis.
%
%   epsq = sing_meas(var, q) computes the qth-order singularity measure
%   of the 1D signal VAR at dyadic scales using the method of Davis et al.
%   (1994, JGR). The signal is truncated to the largest power-of-two length.
%
%   epsq = sing_meas(var, q, minx, l) restricts the finest grain to
%   separations >= MINX, where L is the grid spacing (so that the first
%   subdivision index is floor(log2(minx/l))).
%
%   The output EPSQ has length nextpow2(length(var))-1, one entry per
%   dyadic scale from finest to coarsest.
%
%   References:
%     Davis, A., A. Marshak, W. Wiscombe, and R. Cahalan, 1994:
%       Multifractal characterizations of nonstationarity and intermittency
%       in geophysical fields. J. Geophys. Res., 99, 8055-8072.

%% Parse optional arguments for minimum scale
if nargin > 2
    minx = varargin{1};
    l    = varargin{2};
    r1   = max(0, floor(log2(minx / l)));
else
    r1 = 0;
end

%% Truncate signal to largest power-of-two length
N2  = nextpow2(length(var)) - 1;   % exponent: 2^N2 points used
Npt = 2^N2;
v   = var(1:Npt);

%% First-difference (proper increments, no zero-boundary)
%  Produces Npt-1 elements; pad with NaN to keep length Npt for
%  compatibility with the non-overlapping partition below.
dv       = NaN(Npt, 1);
dv(1:end-1) = v(2:end) - v(1:end-1);

%% Normalised absolute increments: epsilon(1; x) in Eq. 7a of Davis (1994)
e1x = abs(dv) ./ mean(abs(dv), 'omitnan');

%% Coarse-grain at dyadic scales using non-overlapping partitions
%  At scale 2^n, divide the signal into Npt/2^n bins of width 2^n and
%  average e1x within each bin, then compute the qth moment.
scales = r1 : (N2 - 1);
epsq   = NaN(N2, 1);

for idx = 1:length(scales)
    n     = scales(idx);
    binsz = 2^n;
    nbins = floor(Npt / binsz);
    if nbins < 2
        continue
    end
    binmeans = NaN(nbins, 1);
    for b = 1:nbins
        j1 = (b-1)*binsz + 1;
        j2 = b*binsz;
        binmeans(b) = mean(e1x(j1:j2), 'omitnan');
    end
    epsq(idx + r1) = mean(binmeans.^q, 'omitnan');
end

end
