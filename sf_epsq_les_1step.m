function sf_epsq_les_1step(run, var, ti, varargin)
% SF_EPSQ_LES_1STEP  1D structure functions for one timestep of CM1 output.
%
%   sf_epsq_les_1step(RUN, VAR, TI) computes 1D structure functions and
%   singularity measures for variable VAR at timestep index TI.
%
%   sf_epsq_les_1step(RUN, VAR, TI, fd) reads from directory FD.
%
%   Dimension ordering is determined dynamically via nc_dims / nc_read.

tStart = tic;

%% Parse input
if nargin > 3
    fd = varargin{1};
else
    fd = './';
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
nxg = dsz.x;
nz  = dsz.z;

fprintf('Grid: nx=%i, nz=%i (dims: %s=%i, %s=%i, %s=%i, %s=%i)\n', ...
    nxg, nz, dnames.x, dsz.x, dnames.y, dsz.y, dnames.z, dsz.z, dnames.t, dsz.t)

%% Load statistics and PSD results
switch var
    case {'l', 'c'}
        loadvars = {'cfrac', 'qcprof'};
        outnm    = sprintf('sf_q%s_%s', var, run);
    case 'r'
        loadvars = {'cfrac', 'rfrac', 'qcprof', 'qrprof'};
        outnm    = sprintf('sf_q%s_%s', var, run);
    otherwise
        error('Variable %s not configured for multifractal analysis', var)
end

stats = load(fullfile(fd, sprintf('cm1_%s_stats.mat', run)), loadvars{:});
psd   = load(fullfile(fd, sprintf('psd_q%s_%s.mat', var, run)), 'slope');

cfrac = stats.cfrac;
if isfield(stats, 'rfrac'), rfrac = stats.rfrac; end
slope = psd.slope;

%% Set area fraction threshold
switch var
    case {'qc', 'l', 'c'}
        thresh = 0.015;
        frac   = cfrac(:, ti);
    case {'qr', 'r'}
        thresh = 0.015;
        frac   = rfrac(:, ti);
    otherwise
        error('Set up multifractal analysis for variable: %s', var)
end
slopet = slope(:, ti);

%% Set output filename
tstr   = sprintf('_%05i', ti);
outnmt = fullfile(varargin{2}, strcat(outnm, tstr, '.mat'));

%% Short-circuit if no levels meet threshold
top = find(frac > thresh, 1, 'last');
if isempty(top) || max(frac) < thresh
    fprintf('No levels exceed area fraction threshold at timestep %i\n', ti)
    out.q        = q;
    out.lev.sf   = NaN(nxg-2, nq, nz);
    out.lev.epsq = NaN(nextpow2(nxg)-1, nq, nz);
    out.i        = nz;
    out.k        = nq;
    save(outnmt, '-struct', 'out')
    return
end

%% Load or initialise output
if exist(outnmt, 'file')
    out = load(outnmt);
    fprintf('Loaded existing output: %s\n', outnmt)
else
    out.q        = q;
    out.lev.sf   = NaN(nxg-2, nq, nz);
    out.lev.epsq = NaN(nextpow2(nxg)-1, nq, nz);
    out.i        = 0;
    out.k        = 0;
end

%% Compute structure functions
fprintf('Computing 1D structure functions for timestep %i (levels 1..%i)\n', ti, top)

for ilev = 1:top
    out.i = ilev;
    skip_level = slopet(ilev) < 1 || slopet(ilev) > 3 || frac(ilev) < thresh;

    for iq = 1:nq
        out.k = iq;
        if skip_level || all(~isnan(out.lev.sf(:, iq, ilev)))
            continue
        end

        if q(iq) == 1
            % Integer order: direct structure function
            sf_out = agg_sf_les(fp, ncvar, dmap, dsz, ilev, ti, q(iq));
            if isempty(sf_out)
                fprintf('  agg_sf_les returned empty at level %i\n', ilev)
                continue
            end
            out.lev.sf(:, iq, ilev) = sf_out;
        else
            % Non-integer order: singularity measure
            epsq_out = agg_epsq_les(fp, ncvar, dmap, dsz, ilev, ti, q(iq));
            if isempty(epsq_out)
                continue
            end
            out.lev.epsq(:, iq, ilev) = epsq_out;
        end

        save(outnmt, '-struct', 'out')
    end

    if mod(ilev, 10) == 0
        fprintf('  Level %i/%i (%.1f min elapsed)\n', ilev, top, toc(tStart)/60)
    end
end

fprintf('1D structure functions complete for timestep %i (%.1f min)\n', ti, toc(tStart)/60)

end
