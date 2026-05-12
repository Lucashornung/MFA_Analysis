function epsq = sing_meas_2d(field, q)
% SING_MEAS_2D  Two-dimensional singularity measure for multifractal analysis.
%
%   epsq = sing_meas_2d(FIELD, Q) computes the qth-order singularity
%   measure of the 2D field FIELD using dyadic box-counting, following the
%   2D extension of the Davis et al. (1994) approach.
%
%   The procedure:
%     1. Compute a local measure field from the mean of absolute
%        differences in x and y (using circshift for periodic BCs):
%            e1(i,j) = ( |f(i+1,j) - f(i,j)| + |f(i,j+1) - f(i,j)| ) / 2
%        normalised by the domain mean of e1.
%
%     2. At each dyadic scale 2^n, partition the field into non-overlapping
%        boxes of size 2^n x 2^n and average e1 within each box.
%
%     3. Compute eps_q(n) = < box_mean^q > (the qth moment of box averages).
%
%   The field is truncated to the largest square with side length 2^N.
%   Output has length N = nextpow2(min(nx,ny)) - 1, one entry per scale.
%
%   Inputs:
%     FIELD  - 2D array [nx x ny]
%     Q      - moment order (scalar, typically 0.9 or 1.1)
%
%   Output:
%     epsq   - [N x 1] vector of eps_q at dyadic scales 2^0, 2^1, ..., 2^(N-1)
%
%   Reference:
%     Davis, A., A. Marshak, W. Wiscombe, and R. Cahalan, 1994:
%       Multifractal characterizations of nonstationarity and intermittency
%       in geophysical fields. J. Geophys. Res., 99, 8055-8072.

[nx, ny] = size(field);

%% Truncate to largest square power-of-two domain
N2   = min(nextpow2(nx), nextpow2(ny)) - 1;
Npt  = 2^N2;
f    = field(1:Npt, 1:Npt);

%% Compute local 2D measure field
%  Mean of absolute differences in x and y directions, using circshift
%  for periodic boundary conditions.
dx_diff = abs(circshift(f, [-1, 0]) - f);
dy_diff = abs(circshift(f, [0, -1]) - f);
e1 = (dx_diff + dy_diff) / 2;

%% Normalise to unit mean
e1_mean = mean(e1(:), 'omitnan');
if e1_mean == 0 || isnan(e1_mean)
    epsq = NaN(N2, 1);
    return
end
e1 = e1 / e1_mean;

%% Coarse-grain at dyadic scales using non-overlapping 2D boxes
%  At scale 2^n, the domain is partitioned into (Npt/2^n)^2 boxes.
%  The box average of e1 is computed, then the qth moment is taken
%  over all boxes.
epsq = NaN(N2, 1);

for n = 0:(N2 - 1)
    binsz = 2^n;
    nbins = Npt / binsz;

    if nbins < 2
        continue
    end

    % Vectorised box averaging via reshape:
    %   reshape(e1, binsz, nbins, binsz, nbins) maps element (i,j,k,l) to
    %   e1(i + (j-1)*binsz, k + (l-1)*binsz), where i,k index within a box
    %   and j,l index the box. Averaging over dims 1 and 3 gives box means.
    temp      = reshape(e1, binsz, nbins, binsz, nbins);
    box_means = squeeze(mean(mean(temp, 1), 3));   % [nbins x nbins]

    epsq(n + 1) = mean(box_means(:).^q, 'omitnan');
end

end
