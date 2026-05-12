function out = agg_epsq_les(fp, ncvar, dmap, dsz, klev, ti, q)
% AGG_EPSQ_LES  Aggregate singularity measures over one vertical level.
%
%   out = agg_epsq_les(FP, NCVAR, DMAP, DSZ, KLEV, TI, Q)
%
%   Reads horizontal transects of variable NCVAR at vertical level KLEV and
%   timestep TI from NetCDF file FP, computes the qth-order singularity
%   measure via sing_meas for each y-transect, and returns the y-average.
%
%   Dimension ordering is handled via DMAP/DSZ from nc_dims.
%
%   Inputs:
%     FP    - full path to NetCDF file
%     NCVAR - NetCDF variable name (e.g., 'qc', 'qr')
%     DMAP  - dimension map struct from nc_dims
%     DSZ   - dimension size struct from nc_dims
%     KLEV  - vertical level index
%     TI    - time index
%     Q     - singularity measure order (typically 0.9 or 1.1)
%
%   Output:
%     out   - [nextpow2(nx)-1 x 1] vector, or [] if all NaN.

nx = dsz.x;
ny = dsz.y;

epsq = [];

%% Loop over y-transects
for iy = 1:ny
    % Read one x-transect: result is (nx,) regardless of file layout
    transect = nc_read(fp, ncvar, dmap, dsz, 't', ti, 'y', iy, 'z', klev);

    if iy == 1
        nr   = nextpow2(length(transect)) - 1;
        epsq = NaN(nr, ny);
    end
    epsq(:, iy) = sing_meas(transect, q);
end

%% Average over all y-transects
if isempty(epsq) || all(isnan(epsq(:)))
    out = [];
else
    out = mean(epsq, 2, 'omitnan');
end

end
