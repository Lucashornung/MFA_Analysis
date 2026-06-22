% RUN_MFA_PIPELINE  Driver script for multifractal analysis of CM1 LES output.
%
% Author: Mikael Witte (with annotation assistance from Claude)
% Last updated: 17 April 2026
%
%   This script runs the four-stage pipeline from CM1 NetCDF output to
%   bifractal parameters (H1, C1) for cloud and/or rain water mixing ratio.
%
%   Pipeline stages:
%     1. get_stats_cm1      - domain-mean profiles and area fractions
%     2. process_psd_cm1    - power spectral density and scaling regime
%     3. sf_epsq_*_1step    - structure functions (per timestep)
%                             1D: sf_epsq_les_1step (transect-based, as in Witte et al. 2022)
%                             2D: sf_epsq_2d_1step  (isotropic radial, level-by-level)
%     4. mfa_epsq           - fit scaling exponents -> (H1, C1)
%
%   Expected input:
%     A single NetCDF file named <run>.nc with dimensions (t,x,y,z)
%     and variables qc (cloud water) and qr (rain water) in kg/kg.
%
%   References:
%     Witte, M. K., H. Morrison, A. B. Davis, and J. Teixeira, 2022:
%       Limitations of bin and bulk microphysics in reproducing the observed
%       spatial structure of light precipitation. J. Atmos. Sci., 79, 161-178.

%% ========================================================================
%  USER CONFIGURATION
%  ========================================================================

run = 'cm1out_000002';       % run identifier (matches filename <run>.nc)
fd  = '~/Documents/Work/sct_les'; % directory containing the NetCDF file
outputd = ''                 % directory for output files
var = 'c';                  % variable to analyse: 'c' (cloud) or 'r' (rain)

% Structure function mode:
%   '1d' - 1D transects averaged over y (as in Witte et al. 2022 JAS)
%   '2d' - 2D isotropic with periodic-BC circshift (preferred for LES)
sf_mode = '1d';

% Timesteps to analyse for structure functions (Stage 3).
% These should correspond to quasi-steady-state output (skip spin-up).
% Set to [] to process all timesteps found in the stats file.
ti_list = [];

%% ========================================================================
%  STAGE 1: Compute domain-mean statistics
%  ========================================================================
fprintf('\n=== STAGE 1: get_stats_cm1 ===\n')

statsfile = fullfile(outputd, sprintf('cm1_%s_stats.mat', run));
if exist(statsfile, 'file')
    fprintf('Stats file already exists: %s\n', statsfile)
else
    get_stats_cm1(run, fd, outputd);
end

%% ========================================================================
%  STAGE 2: Compute PSDs and diagnose scaling regime
%  ========================================================================
fprintf('\n=== STAGE 2: process_psd_cm1 ===\n')

psdfile = fullfile(outputd, sprintf('psd_q%s_%s.mat', var, run));
if exist(psdfile, 'file')
    fprintf('PSD file already exists: %s\n', psdfile)
else
    process_psd_cm1(run, var, fd, outputd);
end

%% ========================================================================
%  STAGE 3: Compute structure functions (one call per timestep)
%  ========================================================================
fprintf('\n=== STAGE 3: sf_epsq_les_1step ===\n')

% Determine timestep list if not specified
if isempty(ti_list)
    stats   = load(statsfile, 'qcprof');
    [~, nt] = size(stats.qcprof);
    ti_list = 1:nt;
end

% Dispatch to 1D or 2D structure function driver
switch sf_mode
    case '1d'
        sf_prefix = 'sf_q';
        for ti = ti_list
            fprintf('\n--- Timestep %i/%i (1D) ---\n', ti, ti_list(end))
            sf_epsq_les_1step(run, var, ti, fd, outputd);
        end
    case '2d'
        sf_prefix = 'sf2d_q';
        for ti = ti_list
            fprintf('\n--- Timestep %i/%i (2D) ---\n', ti, ti_list(end))
            sf_epsq_2d_1step(run, var, ti, fd, outputd);
        end
    otherwise
        error('sf_mode must be ''1d'' or ''2d''')
end

% Merge per-timestep SF files into a single file for mfa_epsq
fprintf('\nMerging per-timestep structure function files...\n')
merge_sf_files(run, var, ti_list, outputd, sf_prefix);

%% ========================================================================
%  STAGE 4: Extract bifractal parameters
%  ========================================================================
fprintf('\n=== STAGE 4: mfa_epsq ===\n')

mfa = mfa_epsq(run, var, fd, sf_mode);

% Save final output
outfile = fullfile(outputd, sprintf('mfa_%s%s_%s.mat', sf_prefix, var, run));
save(outfile, '-struct', 'mfa')
fprintf('\nFinal output saved to %s\n', outfile)

% Print summary statistics
h1_mean = mean(mfa.lev.h1(:), 'omitnan');
c1_mean = mean(mfa.lev.c1(:), 'omitnan');
fprintf('\n=== RESULTS SUMMARY ===\n')
fprintf('  Mode:     %s\n', upper(sf_mode))
fprintf('  Variable: q%s\n', var)
fprintf('  Domain-mean H1 = %.3f  (smoothness)\n', h1_mean)
fprintf('  Domain-mean C1 = %.3f  (intermittency)\n', c1_mean)
fprintf('  For reference:  passive scalar in turbulence -> (H1, C1) ~ (0.33, 0.05)\n')

%% ========================================================================
%  STAGE 5 (optional): Compare 1D vs 2D results
%  ========================================================================
if strcmp(sf_mode, '2d')
    fprintf('\n=== 1D vs 2D comparison available ===\n')
    fprintf('  The 2D output files (sf2d_*) contain axis-aligned S_1(r)\n')
    fprintf('  in fields lev.sf_x and lev.sf_y for direct comparison\n')
    fprintf('  with the 1D transect-based results.\n')
    fprintf('  Run again with sf_mode = ''1d'' to generate 1D output,\n')
    fprintf('  then use compare_1d_2d() below.\n')
end

%% ========================================================================
%  HELPER: Merge per-timestep SF output into one file
%  ========================================================================
function merge_sf_files(run, var, ti_list, fd, prefix, outputd)
    outnm_base  = sprintf('%s%s_%s', prefix, var, run);
    merged_file = fullfile(outputd, [outnm_base, '.mat']);

    first = true;
    for ti = ti_list
        tstr    = sprintf('_%05i', ti);
        step_fn = fullfile(outputd, [outnm_base, tstr, '.mat']);
        if ~exist(step_fn, 'file')
            warning('Missing SF file for timestep %i: %s', ti, step_fn)
            continue
        end
        step = load(step_fn);

        if first
            first = false;
            nt = length(ti_list);
            [nr, nq, nz] = size(step.lev.sf);
            merged.q    = step.q;
            merged.t    = NaN(nt, 1);
            merged.lev.sf   = NaN(nr, nq, nz, nt);
            merged.lev.epsq = NaN(size(step.lev.epsq, 1), nq, nz, nt);
            % Carry 2D directional fields if present
            if isfield(step.lev, 'sf_x')
                merged.lev.sf_x = NaN(nr, nz, nt);
                merged.lev.sf_y = NaN(nr, nz, nt);
            end
        end

        idx = find(ti_list == ti);
        if isfield(step, 'time')
            merged.t(idx) = step.time;
        end
        merged.lev.sf(:, :, :, idx)   = step.lev.sf;
        if isfield(step.lev, 'epsq')
            merged.lev.epsq(:, :, :, idx) = step.lev.epsq;
        end
        if isfield(step.lev, 'sf_x') && isfield(merged.lev, 'sf_x')
            merged.lev.sf_x(:, :, idx) = step.lev.sf_x;
            merged.lev.sf_y(:, :, idx) = step.lev.sf_y;
        end
    end

    if ~first
        save(merged_file, '-struct', 'merged')
        fprintf('Merged %i timesteps into %s\n', length(ti_list), merged_file)
    else
        warning('No SF files found to merge')
    end
end
