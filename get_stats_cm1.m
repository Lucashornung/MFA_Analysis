function get_stats_cm1(run, varargin)
% GET_STATS_CM1  Compute domain-mean profiles and area fractions from CM1 output.
%
%   get_stats_cm1(RUN) reads the NetCDF file <RUN>.nc from the
%   current directory and computes vertical profiles of domain-mean cloud
%   water (qcprof), rain water (qrprof), and their horizontal area
%   fractions (cfrac, rfrac). Output is saved to cm1_<RUN>_stats.mat.
%
%   get_stats_cm1(RUN, fd) reads from directory FD instead of './'.
%
%    
%    get_stats_cm1(RUN, fd, outputd) writes the output file to the outputd directory 
%    NOTE!!!! THIS FUNCTION WILL FAIL IN THE CASE THAT THERE IS AN OUTPUTD BUT NO FD PARAMETER
%
%   Dimension ordering is determined dynamically from the NetCDF metadata
%   via nc_dims / nc_read, so the code is agnostic to whether the file
%   uses (t,x,y,z), (x,y,z,t), or any other permutation.

%% Parse input
if nargin > 1
    fd = varargin{1};
else
    fd = './';
end
fp = fullfile(fd, sprintf('%s.nc', run));

%% Select variables to process
if contains(run, 'sdm')
    vars = {'c', 'r', 'out3', 'out4', 'out5', 'out6'};
else
    vars = {'c', 'r'};
end

%% Initialise output structure (z coordinate set after first nc_dims call)
out = struct();

%% Loop over variables
for j = 1:length(vars)

    % Map short variable name to profile/fraction field names
    switch vars{j}
        case {'l', 'c'},       prf_nme = 'qcprof'; frc_nme = 'cfrac';
        case 'r',              prf_nme = 'qrprof'; frc_nme = 'rfrac';
        case 'q',              prf_nme = 'qvprof'; frc_nme = '';
        case 'r1',             prf_nme = 'M1prof'; frc_nme = '';
        case 'r2',             prf_nme = 'M2prof'; frc_nme = '';
        case {'out3', 'r5'},   prf_nme = 'M5prof'; frc_nme = '';
        case {'out4', 'r6'},   prf_nme = 'M6prof'; frc_nme = 'M6frac';
        case 'out5',           prf_nme = 'Ntprof'; frc_nme = '';
        case 'ncd',            prf_nme = 'Ncprof'; frc_nme = '';
        case {'out6', 'ncr'},  prf_nme = 'Nrprof'; frc_nme = '';
        otherwise
            prf_nme = sprintf('%sprof', vars{j});
            frc_nme = '';
    end

    % NetCDF variable name: 'c' -> 'qc', 'r' -> 'qr', etc.
    switch vars{j}
        case {'l', 'c', 'r'},  loadvar = ['q', vars{j}];
        otherwise,              loadvar = vars{j};
    end

    % Query dimension ordering for this variable
    [dmap, dsz, dnames] = nc_dims(fp, loadvar);
    nx = dsz.x;  ny = dsz.y;  nz = dsz.z;  nt = dsz.t;

    % Store time and z coordinates on first variable (use actual names from file)
    if ~isfield(out, 't')
        out.t = ncread(fp, dnames.t);
        out.z = ncread(fp, dnames.z);
    end

    % Allocate output arrays if not yet present
    if ~isfield(out, prf_nme)
        out.(prf_nme) = NaN(nz, nt);
        if ~isempty(frc_nme)
            out.(frc_nme) = NaN(nz, nt);
        end
    end

    % Loop over time steps
    for i = 1:nt
        % Read 3D field: result is always (x, y, z) regardless of file layout
        in = nc_read(fp, loadvar, dmap, dsz, 't', i);

        switch prf_nme
            case 'M6prof'
                in(in == 0) = NaN;
                in = 10 * log10(in);
        end

        % Domain-mean profile: mean over x and y -> (nz,)
        out.(prf_nme)(:, i) = squeeze(mean(mean(in, 1, 'omitnan'), 2, 'omitnan'));

        % Area fraction profile
        if ~isempty(frc_nme)
            for k = 1:nz
                qk = in(:, :, k);
                switch prf_nme
                    case 'M6prof'
                        out.(frc_nme)(k, i) = sum(~isnan(qk(:))) / numel(qk);
                    otherwise
                        out.(frc_nme)(k, i) = sum(qk(:) > 1e-5) / numel(qk);
                end
            end
        end
    end
end

%% Save
save(fullfile(varargin{2}, sprintf('cm1_%s_stats.mat', run)), '-struct', 'out')
fprintf('Saved cm1_%s_stats.mat (%i timesteps, %i levels)\n', run, nt, nz)

end
