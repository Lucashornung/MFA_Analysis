function [start, count] = nc_slice(dmap, dim_sizes, varargin)
% NC_SLICE  Build start/count vectors for ncread from a dimension map.
%
%   [start, count] = nc_slice(dmap, dim_sizes, 'dim1', idx1, 'dim2', idx2, ...)
%
%   By default, each dimension selects ALL elements (start=1, count=size).
%   Name-value pairs override specific dimensions to select a single index
%   (count=1) or a range.
%
%   Examples:
%     % Read a full 2D horizontal field at one timestep and level:
%     [s, c] = nc_slice(dmap, dsz, 't', ti, 'z', kz);
%     data = squeeze(ncread(fp, varname, s, c));
%
%     % Read a single x-transect at one y, z, and t:
%     [s, c] = nc_slice(dmap, dsz, 't', ti, 'y', iy, 'z', kz);
%     data = squeeze(ncread(fp, varname, s, c));
%
%     % Read a 3D (x,y,z) volume at one timestep:
%     [s, c] = nc_slice(dmap, dsz, 't', ti);
%     data = squeeze(ncread(fp, varname, s, c));

ndims_var = 0;
fields = {'x', 'y', 'z', 't'};
for i = 1:length(fields)
    if ~isempty(dmap.(fields{i}))
        ndims_var = max(ndims_var, dmap.(fields{i}));
    end
end

%% Defaults: read everything
start = ones(1, ndims_var);
count = ones(1, ndims_var);

for i = 1:length(fields)
    f = fields{i};
    if ~isempty(dmap.(f)) && ~isempty(dim_sizes.(f))
        count(dmap.(f)) = dim_sizes.(f);
    end
end

%% Apply overrides from name-value pairs
for i = 1:2:length(varargin)
    dim_name = varargin{i};
    dim_val  = varargin{i+1};

    if ~isfield(dmap, dim_name)
        error('nc_slice: unknown dimension name "%s"', dim_name)
    end
    pos = dmap.(dim_name);
    if isempty(pos)
        error('nc_slice: dimension "%s" not present in this variable', dim_name)
    end

    if isscalar(dim_val)
        start(pos) = dim_val;
        count(pos) = 1;
    elseif length(dim_val) == 2
        % [start_idx, count] range
        start(pos) = dim_val(1);
        count(pos) = dim_val(2);
    end
end

end
