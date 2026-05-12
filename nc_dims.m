function [dmap, dim_sizes, dnames] = nc_dims(fp, varname)
% NC_DIMS  Query dimension ordering, sizes, and names of a NetCDF variable.
%
%   [dmap, dsz, dnames] = nc_dims(FP, VARNAME) returns three structs
%   describing the dimensions of variable VARNAME in the NetCDF file FP.
%   Each struct has fields x, y, z, t corresponding to the logical role
%   of each dimension.
%
%   dmap   — 1-based position of each dimension in the variable's
%            dimension list. For a file with dimensions (time, ni, nj, nk):
%                dmap.t = 1, dmap.x = 2, dmap.y = 3, dmap.z = 4
%
%   dsz    — Size (length) of each dimension:
%                dsz.t = 33, dsz.x = 128, dsz.y = 128, dsz.z = 100
%
%   dnames — Actual dimension names as they appear in the file:
%                dnames.t = 'time', dnames.x = 'ni', dnames.y = 'nj',
%                dnames.z = 'nk'
%            Useful for diagnostics, coordinate variable reads, and
%            passing to other NetCDF tools.
%
%   Recognised dimension name aliases (case-insensitive):
%       x: 'x', 'ni', 'xh', 'west_east', 'lon', 'rlon'
%       y: 'y', 'nj', 'yh', 'south_north', 'lat', 'rlat'
%       z: 'z', 'nk', 'zh', 'bottom_top', 'lev', 'level', 'height'
%       t: 't', 'time', 'nt', 'times'
%
%   Example:
%       [dmap, dsz, dnames] = nc_dims('cm1out.nc', 'qc');
%       fprintf('x-dimension is called "%s" with %i points\n', dnames.x, dsz.x);
%       % Read a 2D field at timestep 5, level 20:
%       field = nc_read('cm1out.nc', 'qc', dmap, dsz, 't', 5, 'z', 20);
%
%   See also: nc_slice, nc_read

%% Alias tables (case-insensitive matching)
x_aliases = {'x', 'ni', 'xh', 'west_east', 'lon', 'rlon'};
y_aliases = {'y', 'nj', 'yh', 'south_north', 'lat', 'rlat'};
z_aliases = {'z', 'nk', 'zh', 'bottom_top', 'lev', 'level', 'height'};
t_aliases = {'t', 'time', 'nt', 'times'};

%% Query variable info from NetCDF metadata
vinfo = ncinfo(fp, varname);
ndims_var = length(vinfo.Dimensions);

dmap      = struct('x', [], 'y', [], 'z', [], 't', []);
dim_sizes = struct('x', [], 'y', [], 'z', [], 't', []);
dnames    = struct('x', '', 'y', '', 'z', '', 't', '');

for i = 1:ndims_var
    raw_name = vinfo.Dimensions(i).Name;
    dname    = lower(raw_name);
    dlen     = vinfo.Dimensions(i).Length;

    if any(strcmp(dname, x_aliases))
        dmap.x      = i;
        dim_sizes.x = dlen;
        dnames.x    = raw_name;
    elseif any(strcmp(dname, y_aliases))
        dmap.y      = i;
        dim_sizes.y = dlen;
        dnames.y    = raw_name;
    elseif any(strcmp(dname, z_aliases))
        dmap.z      = i;
        dim_sizes.z = dlen;
        dnames.z    = raw_name;
    elseif any(strcmp(dname, t_aliases))
        dmap.t      = i;
        dim_sizes.t = dlen;
        dnames.t    = raw_name;
    else
        warning('nc_dims: unrecognised dimension "%s" at position %i in variable "%s"', ...
                raw_name, i, varname)
    end
end

%% Validate: all four dimensions must be identified for a 4D variable
missing = {};
if isempty(dmap.x), missing{end+1} = 'x'; end
if isempty(dmap.y), missing{end+1} = 'y'; end
if isempty(dmap.z), missing{end+1} = 'z'; end
if isempty(dmap.t), missing{end+1} = 't'; end

if ~isempty(missing) && ndims_var == 4
    error(['nc_dims: could not identify dimensions {%s} in variable "%s".\n' ...
           '  Found: %s\n' ...
           '  Add missing aliases to the tables in nc_dims.m.'], ...
        strjoin(missing, ', '), varname, ...
        strjoin({vinfo.Dimensions.Name}, ', '))
end

end
