function sf_epsq_2d_1step(run, var, ti, varargin)
% SF_EPSQ_2D_1STEP  2D structure functions and singularity measures for one timestep.
%
%   sf_epsq_2d_1step(RUN, VAR, TI) computes 2D isotropic structure functions
%   and singularity measures for variable VAR at timestep index TI.
%
%   sf_epsq_2d_1step(RUN, VAR, TI, fd) reads from directory FD.
%
%   sf_epsq_2d_1step(RUN, VAR, TI, fd, outputd) reads stats and psd files from a seperate output directory.
%
%   Dimension ordering is determined dynamically via nc_dims / nc_read.

tStart = tic;

%% Parse input
if nargin > 3
    fd = varargin{1};
else
    fd = './';
end

if nargin > 4
    outd = varargin{2};
end

%% Structure function orders
q  = [0.9, 1.0, 1.1];
nq = length(q);

%% Locate NetCDF file and query dimensions
D = dir(fullfile(fd, ['*', run, '*.nc']));
if isempty(D)
    error('No NetCDF files matching *%s*.nc found in %s', run, fd)
end
fn = D(1).name;
fp = fullfile(fd, fn);

% NetCDF variable name
switch var
    case {'c', 'l', 'r'},  ncvar = ['q', var];
    otherwise,              ncvar = var;
end

% Dynamic dimension ordering
[dmap, dsz, dnames] = nc_dims(fp, ncvar);
nx = dsz.x;
ny = dsz.y;
nz = dsz.z;

fprintf('Grid: nx=%i, ny=%i, nz=%i (dims: %s, %s, %s, %s)\n', ...
    nx, ny, nz, dnames.x, dnames.y, dnames.z, dnames.t)

%% Derived parameters
rmax   = floor(min(nx, ny) / 2);
N_dyad = min(nextpow2(nx), nextpow2(ny)) - 1;

%% Load statistics and PSD scale breaks
switch var
    case {'l', 'c'}
        loadvars = {'cfrac', 'qcprof'};
        outnm    = sprintf('sf2d_q%s_%s', var, run);
    case 'r'
        loadvars = {'cfrac', 'rfrac', 'qcprof', 'qrprof'};
        outnm    = sprintf('sf2d_q%s_%s', var, run);
    otherwise
        error('Variable %s not configured for 2D multifractal analysis', var)
end

stats = load(fullfile(outd, sprintf('cm1_%s_stats.mat', run)), loadvars{:});
psd   = load(fullfile(outd, sprintf('psd_q%s_%s.mat', var, run)), 'slope');

cfrac = stats.cfrac;
if isfield(stats, 'rfrac'), rfrac = stats.rfrac; end
slope = psd.slope;

%% Set area fraction threshold
switch var
    case {'l', 'c'}
        thresh = 0.2;
        frac   = cfrac(:, ti);
    case 'r'
        thresh = 0.015;
        frac   = rfrac(:, ti);
end
slopet = slope(:, ti);

%% Set output filename
tstr   = sprintf('_%05i', ti);
outnmt = fullfile(outd, [outnm, tstr, '.mat']);

%% Short-circuit if no levels meet threshold
top = find(frac > thresh, 1, 'last');
if isempty(top) || max(frac) < thresh
    fprintf('No levels exceed area fraction threshold at timestep %i\n', ti)
    out.q           = q;
    out.lev.sf      = NaN(rmax, nq, nz);
    out.lev.epsq    = NaN(N_dyad, nq, nz);
    out.lev.sf_x    = NaN(rmax, nz);
    out.lev.sf_y    = NaN(rmax, nz);
    out.i           = nz;
    out.k           = nq;
    save(outnmt, '-struct', 'out')
    return
end

%% Initialise output
out.q           = q;
out.lev.sf      = NaN(rmax, nq, nz);
out.lev.epsq    = NaN(N_dyad, nq, nz);
out.lev.sf_x    = NaN(rmax, nz);
out.lev.sf_y    = NaN(rmax, nz);

%% Loop over vertical levels
fprintf('Computing 2D structure functions for timestep %i (levels 1..%i)\n', ti, top)

for ilev = 1:top

    if slopet(ilev) < 1 || slopet(ilev) > 3 || frac(ilev) < thresh
        continue
    end

    % Read full 2D field: result is always (nx, ny) regardless of file layout
    field_2d = nc_read(fp, ncvar, dmap, dsz, 't', ti, 'z', ilev);

    % q = 1: 2D isotropic structure function
    iq1 = find(q == 1);
    [Sq_iso, ~, Sq_x, Sq_y] = sf_2d(field_2d, 1, rmax);
    out.lev.sf(:, iq1, ilev) = Sq_iso;
    out.lev.sf_x(:, ilev)    = Sq_x;
    out.lev.sf_y(:, ilev)    = Sq_y;

    % q = 0.9 and 1.1: 2D singularity measures
    for iq = 1:nq
        if q(iq) == 1
            continue
        end
        epsq_out = sing_meas_2d(field_2d, q(iq));
        out.lev.epsq(:, iq, ilev) = epsq_out;
    end

    elapsed = toc(tStart);
    fprintf('  Level %3i/%i  (%.1f min elapsed)\n', ilev, top, elapsed/60)

    out.i = ilev;
    save(outnmt, '-struct', 'out')
end

fprintf('2D structure functions complete for timestep %i (%.1f min)\n', ti, toc(tStart)/60)

end
