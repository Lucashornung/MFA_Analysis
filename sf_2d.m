function [Sq_iso, r_bins, Sq_x, Sq_y] = sf_2d(field, q, rmax)
% SF_2D  Two-dimensional isotropic structure function with periodic BCs.
%
%   [Sq_iso, r_bins] = sf_2d(FIELD, Q) computes the qth-order isotropic
%   2D structure function of the 2D array FIELD, using circshift for
%   periodic boundary conditions:
%
%       S_q(r) = < |f(x + dx, y + dy) - f(x, y)|^q >
%
%   averaged over all displacement vectors (dx, dy) with |d| in each
%   radial bin. Radial bins are integer-spaced: bin k contains all
%   displacements with k - 0.5 <= sqrt(dx^2 + dy^2) < k + 0.5.
%
%   [Sq_iso, r_bins, Sq_x, Sq_y] = sf_2d(FIELD, Q) also returns the
%   structure function computed along x-only (dy = 0) and y-only (dx = 0)
%   displacements, for comparison with 1D transect-based estimates.
%
%   sf_2d(FIELD, Q, RMAX) specifies the maximum radial separation in grid
%   units. Default: floor(min(nx, ny) / 2), which avoids wrap-around
%   ambiguity in periodic domains.
%
%   Inputs:
%     FIELD  - 2D array [nx x ny]
%     Q      - structure function order (scalar)
%     RMAX   - maximum radial separation in grid spacings (optional)
%
%   Outputs:
%     Sq_iso - isotropic (radially averaged) S_q(r), [rmax x 1]
%     r_bins - radial bin centres in grid spacings, [rmax x 1]
%     Sq_x   - S_q along x-axis only (dy = 0), [rmax x 1]
%     Sq_y   - S_q along y-axis only (dx = 0), [rmax x 1]

[nx, ny] = size(field);

if nargin < 3
    rmax = floor(min(nx, ny) / 2);
end

%% Preallocate radial accumulation arrays
Sq_sum = zeros(rmax, 1);
Sq_cnt = zeros(rmax, 1);

%% Iterate over displacement vectors in the upper half-plane
%  Exploit symmetry: S_q(dx, dy) = S_q(-dx, -dy), so we only compute
%  displacements in the upper half-plane (dy > 0, or dy == 0 with dx > 0)
%  and count each contribution once.

for dy = 0:rmax
    for dx = -rmax:rmax
        % Skip origin
        if dx == 0 && dy == 0
            continue
        end
        % Upper half-plane only (skip lower half to avoid double-counting)
        if dy == 0 && dx < 0
            continue
        end

        r = sqrt(dx^2 + dy^2);
        rbin = round(r);
        if rbin < 1 || rbin > rmax
            continue
        end

        % Structure function for this displacement (circshift = periodic BC)
        shifted = circshift(field, [-dx, -dy]);
        sq_val  = mean(abs(field(:) - shifted(:)).^q, 'omitnan');

        Sq_sum(rbin) = Sq_sum(rbin) + sq_val;
        Sq_cnt(rbin) = Sq_cnt(rbin) + 1;
    end
end

%% Radially averaged (isotropic) structure function
Sq_iso = Sq_sum ./ max(Sq_cnt, 1);
Sq_iso(Sq_cnt == 0) = NaN;
r_bins = (1:rmax)';

%% Axis-aligned structure functions for 1D comparison
if nargout > 2
    nr_x = min(rmax, floor(nx/2));
    nr_y = min(rmax, floor(ny/2));
    Sq_x = NaN(rmax, 1);
    Sq_y = NaN(rmax, 1);

    for dx = 1:nr_x
        shifted = circshift(field, [-dx, 0]);
        Sq_x(dx) = mean(abs(field(:) - shifted(:)).^q, 'omitnan');
    end

    for dy = 1:nr_y
        shifted = circshift(field, [0, -dy]);
        Sq_y(dy) = mean(abs(field(:) - shifted(:)).^q, 'omitnan');
    end
end

end
