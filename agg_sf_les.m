function out = agg_sf_les(fp, ncvar, dmap, dsz, klev, ti, q)
% AGG_SF_LES  Aggregate qth-order structure functions over one vertical level.
%
%   out = agg_sf_les(FP, NCVAR, DMAP, DSZ, KLEV, TI, Q)
%
%   Reads horizontal transects of variable NCVAR at vertical level KLEV and
%   timestep TI from NetCDF file FP, computes the qth-order structure
%   function S_q(r) for each y-transect, and returns the y-averaged result.
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
%     Q     - structure function order (typically 1 or 2)
%
%   Output:
%     out   - [nx-2 x 1] vector of S_q(r), or [] if all NaN.

nx = dsz.x;
ny = dsz.y;
nr = nx - 2;
sf = NaN(nr, ny);

%% Loop over y-transects
for iy = 1:ny
    % Read one x-transect: result is (nx,) regardless of file layout
    transect = nc_read(fp, ncvar, dmap, dsz, 't', ti, 'y', iy, 'z', klev);

    sfr = NaN(nr, 1);
    for r = 1:nr
        sfr(r) = mean(abs(transect((r+1):nx) - transect(1:nx-r)).^q, 'omitnan');
    end
    sf(:, iy) = sfr;
end

%% Average over all y-transects
if all(isnan(sf(:)))
    out = [];
else
    out = mean(sf, 2, 'omitnan');
end

end
