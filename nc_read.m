function data = nc_read(fp, varname, dmap, dim_sizes, varargin)
% NC_READ  Read a NetCDF variable and permute to standard (x,y,z) order.
%
%   data = nc_read(FP, VARNAME, DMAP, DIM_SIZES, 'dim', idx, ...)
%
%   Reads variable VARNAME from file FP using the dimension map DMAP
%   (from nc_dims). Name-value pairs specify which indices to select
%   (same syntax as nc_slice). The result is permuted so that the spatial
%   dimensions are always in (x, y, z) order regardless of file layout,
%   and singleton dimensions from slicing are squeezed out.
%
%   Examples:
%     % Read a 3D volume at timestep ti -> result is (x, y, z)
%     vol = nc_read(fp, 'qc', dmap, dsz, 't', ti);
%
%     % Read a 2D horizontal field at timestep ti, level kz -> result is (x, y)
%     field = nc_read(fp, 'qc', dmap, dsz, 't', ti, 'z', kz);
%
%     % Read a 1D transect at timestep ti, y-index iy, level kz -> result is (x,)
%     row = nc_read(fp, 'qc', dmap, dsz, 't', ti, 'y', iy, 'z', kz);

%% Build start/count from nc_slice
[start, count] = nc_slice(dmap, dim_sizes, varargin{:});

%% Read raw data (dimension order matches file)
raw = ncread(fp, varname, start, count);

%% Build permutation to reorder to (x, y, z, t)
target_order = [dmap.x, dmap.y, dmap.z, dmap.t];
[~, perm] = sort(target_order);

%% Permute and squeeze
data = squeeze(permute(raw, perm));

end
