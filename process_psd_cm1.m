function process_psd_cm1(run, var, varargin)
% PROCESS_PSD_CM1  Compute and fit PSDs for all timesteps of a CM1 run.
%
%   process_psd_cm1(RUN, VAR) reads the NetCDF file <RUN>.nc from
%   the current directory and computes PSDs of variable VAR ('c' for cloud
%   water, 'r' for rain water) at each vertical level and timestep.
%
%   process_psd_cm1(RUN, VAR, fd) reads from directory FD.
%
%   Dimension ordering is determined dynamically via nc_dims / nc_read.

tic

%% Parse input
if nargin > 2
    fd = varargin{1};
else
    fd = './';
end

%% Build file path and NetCDF variable name
fp = fullfile(fd, sprintf('%s.nc', run));

switch var
    case {'l', 'c', 'r'},  ncvar = ['q', var];
    otherwise,              ncvar = var;
end

%% Query dimension ordering
[dmap, dsz, dnames] = nc_dims(fp, ncvar);
nx = dsz.x;  ny = dsz.y;  nz = dsz.z;  nt = dsz.t;

%% Handle run name variants
% if contains(run, 'cm1out')
%     run = strrep(run, 'cm1out', 'cm1');
% end

%% Load domain-mean stats (for cloud-top diagnosis)
statnm   = sprintf('cm1_%s_stats.mat', run);
statpath = fullfile(fd, statnm);
if exist(statpath, 'file')
    stats = load(statpath);
    sflag = true;
else
    sflag = false;
    warning('Stats file %s not found; will compute ql profiles on the fly.', statnm)
end

%% Set output filename
switch var
    case {'l', 'r', 'c'},        outnm = sprintf('psd_q%s_%s.mat', var, run);
    case {'w', 'q', 'qc', 'qr'}, outnm = sprintf('psd_%s_%s.mat', var, run);
    otherwise, error('Unrecognised variable: %s', var)
end
outpath = fullfile(varargin{2}, outnm);

%% Check for partial output from a previous interrupted run
i_start = 1;
if exist(outpath, 'file')
    prev = load(outpath);
    if isfield(prev, 'i')
        i_start = prev.i + 1;
        fprintf('Resuming from timestep %i/%i\n', i_start, nt)
    end
end

%% Cloud-variable flag for slope fitting
is_cloud_var = any(strcmp(var, {'c', 'l'}));

%% Read grid spacing and time (using actual dimension names from file)
x  = ncread(fp, dnames.x);
dx = x(2) - x(1);

%% Initialise output structure
time = ncread(fp, dnames.t);
psdo.Time = time;

nfft = max(256, 2^nextpow2(nx));
nf   = floor(nfft/2) + 1;

first_call = true;
for i = i_start:nt
    if first_call
        first_call = false;
        if sflag, qcprof_all = stats.qcprof; end
        psdo.z = ncread(fp, dnames.z);
        if i == 1
            psdo.E     = NaN(nf, nz, nt);
            psdo.slope = NaN(nz, nt);
            psdo.sb    = NaN(nz, 2, nt);
            psdo.sbn   = NaN(nz, 2, nt);
        end
    end

    % Read 3D field: result is always (x, y, z)
    vari = nc_read(fp, ncvar, dmap, dsz, 't', i);

    % Cloud water profile for vertical bounds
    if sflag
        qlprof = qcprof_all(:, i);
    else
        if any(strcmp(var, {'c', 'l'}))
            qlprof = squeeze(mean(mean(vari, 1, 'omitnan'), 2, 'omitnan'));
        else
            ql_3d  = nc_read(fp, 'qc', dmap, dsz, 't', i);
            qlprof = squeeze(mean(mean(ql_3d, 1, 'omitnan'), 2, 'omitnan'));
        end
    end

    if sum(qlprof > 1e-5) == 0
        fprintf('  t=%i: no cloud, skipping\n', i)
        continue
    end

    % Compute PSD and fit spectral slope
    psdi = get_psd_cm1(vari, qlprof, dx);
    psdi = psd_slope_fit_les(psdi, qlprof, is_cloud_var);

    % Store results
    psdo.E(:, :, i)     = psdi.avgz.E;
    psdo.slope(:, i)    = psdi.avgz.slope;
    psdo.sb(:, :, i)    = psdi.avgz.sb;
    psdo.sbn(:, :, i)   = psdi.avgz.sbn;
    if i == i_start
        psdo.w = psdi.w;
        psdo.r = psdi.r;
        psdo.k = psdi.xn;
    end
    psdo.i = i;

    elapsed = toc;
    if elapsed > 18000 || i == nt || i == i_start
        save(fullfile(varargin{2}, outnm), '-struct', 'psdo');
        fprintf('  Checkpoint at t=%i (%.1f min elapsed)\n', i, elapsed/60)
    end
end

fprintf('PSD processing complete for %s/%s (%i timesteps)\n', run, var, nt)

end
