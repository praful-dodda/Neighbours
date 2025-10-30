function [psub, zsub, dsub, nsub, index] = neighbours_stug_index(p0, grid_data, nmax, dmax)
% neighbours_stug_index - Index-based neighbor search for uniformly gridded data
%
% Optimized for uniform grids by working entirely in index space during search.
% Only computes actual distances at the end for output.
%
% Key advantages:
% - No distance calculations during search (pure integer arithmetic)
% - Handles anisotropic grids (different dx, dy)
% - Ellipsoidal space-time expansion
% - Efficient edge case handling
%
% SYNTAX:
%
% [psub, zsub, dsub, nsub, index] = neighbours_stug_index(p0, grid_data, nmax, dmax);
%
% INPUT:
%
% p0         1 by 3 (or 6)  [lon, lat, time] or [lon, lat, time, idx_x, idx_y, idx_t]
% grid_data  struct or char Structure with grid data or path to NetCDF file
%            Must contain: .x, .y, .time, .Lon, .Lat, .Z
% nmax       scalar         maximum number of neighbors to return
% dmax       1 by 3         [max_spatial_dist, max_temporal_dist, spacetime_metric]
%
% OUTPUT:
%
% psub       m by 6         [lon, lat, time, idx_x, idx_y, idx_t] for m neighbors
% zsub       m by 1         data values at selected neighbors
% dsub       m by 2         [spatial_distance, temporal_distance] from p0
% nsub       scalar         number of neighbors returned (m ≤ nmax)
% index      m by 1         linear indices into Z array

%% Input handling
if ischar(grid_data) || isstring(grid_data)
    % NetCDF file path provided
    data.x = ncread(grid_data, 'x');
    data.y = ncread(grid_data, 'y');
    data.time = ncread(grid_data, 'time');
    data.Lon = ncread(grid_data, 'Lon');
    data.Lat = ncread(grid_data, 'Lat');
    data.Z = ncread(grid_data, 'BENZ');
else
    data = grid_data;
end

% Extract target coordinates
lon0 = p0(1);
lat0 = p0(2);
t0 = p0(3);

% Get grid dimensions
nx = length(data.x);
ny = length(data.y);
nt = length(data.time);

% Compute or extract grid indices for p0
if length(p0) >= 6
    idx_x = round(p0(4));
    idx_y = round(p0(5));
    idx_t = round(p0(6));
else
    % Find nearest grid point to p0
    [~, idx_x] = min(abs(data.x - lon0));
    [~, idx_y] = min(abs(data.y - lat0));
    [~, idx_t] = min(abs(data.time - t0));
end

% Ensure indices are within bounds
idx_x = max(1, min(nx, idx_x));
idx_y = max(1, min(ny, idx_y));
idx_t = max(1, min(nt, idx_t));

% Compute grid spacing (assume uniform)
dx = abs(median(diff(data.x)));
dy = abs(median(diff(data.y)));
dt = abs(median(diff(data.time)));

% Handle edge case: empty data
if isempty(data.Z) || all(isnan(data.Z(:)))
    psub = [];
    zsub = [];
    dsub = [];
    nsub = 0;
    index = [];
    return
end

%% Step 1: Convert dmax to index space (grid cells)
% This allows us to work purely with integer indices
rx_max = ceil(dmax(1) / dx);  % Maximum cells in x direction
ry_max = ceil(dmax(1) / dy);  % Maximum cells in y direction
rt_max = ceil(dmax(2) / dt);  % Maximum cells in time direction

%% Step 2: Determine anisotropy and expansion strategy
% If dx and dy are very different, use elliptical expansion
anisotropy_ratio = dy / dx;
is_anisotropic = (anisotropy_ratio < 0.5 || anisotropy_ratio > 2.0);

% Estimate NaN ratio to determine search strategy
sample_size = min(1000, numel(data.Z));
sample_indices = randperm(numel(data.Z), sample_size);
estimated_nanratio = sum(isnan(data.Z(sample_indices))) / sample_size;

%% Step 3: Calculate required number of shells
% Estimate how many shells needed to find nmax neighbors
% For a square shell at distance r, approximate number of points = 8*r (perimeter)
% For circle, area = π*r², so r_needed ≈ sqrt(nmax/π)

if is_anisotropic
    % Ellipse area = π * a * b where a, b are semi-axes
    % Adjust for anisotropy
    r_spatial_estimate = sqrt(nmax / (pi * (1 - estimated_nanratio)));
    rx_shells_needed = ceil(r_spatial_estimate);
    ry_shells_needed = ceil(r_spatial_estimate * anisotropy_ratio);
else
    % Isotropic case
    r_shells_needed = ceil(sqrt(nmax / (pi * (1 - estimated_nanratio))));
    rx_shells_needed = r_shells_needed;
    ry_shells_needed = r_shells_needed;
end

% Ensure we don't exceed dmax bounds
rx_shells_needed = min(rx_shells_needed, rx_max);
ry_shells_needed = min(ry_shells_needed, ry_max);

%% Step 4: Spatial expansion in index space (NO distance calculations yet)
candidates_ix = [];
candidates_iy = [];
candidates_it = [];

% Expand spatially first, then temporally for each spatial location
% This creates an ellipsoid in space-time

for r_shell = 0:max(rx_shells_needed, ry_shells_needed)
    % Generate spatial shell points
    % For each shell, we check points at Chebyshev distance r_shell

    for di = -r_shell:r_shell
        for dj = -r_shell:r_shell
            % Only process points on shell boundary
            if max(abs(di), abs(dj)) ~= r_shell
                continue
            end

            ix = idx_x + di;
            iy = idx_y + dj;

            % Check spatial bounds
            if ix < 1 || ix > nx || iy < 1 || iy > ny
                continue
            end

            % Check if within elliptical spatial bounds (for anisotropic case)
            if is_anisotropic
                % Ellipse equation: (x/a)² + (y/b)² ≤ 1
                ellipse_dist = (di / rx_shells_needed)^2 + (dj / ry_shells_needed)^2;
                if ellipse_dist > 1.1  % Small tolerance
                    continue
                end
            end

            % Now expand temporally for this spatial location
            % Use ellipsoid: spatial_dist²/dmax(1)² + temporal_dist²/dmax(2)² ≤ 1
            % In index space: (di*dx)²/dmax(1)² + (dj*dy)²/dmax(1)² + (dk*dt)²/dmax(2)² ≤ 1

            spatial_dist_idx = sqrt((di * dx)^2 + (dj * dy)^2);

            % Calculate maximum temporal range for this spatial location
            if spatial_dist_idx <= dmax(1)
                % How much "budget" left for temporal distance?
                % ellipsoid: s²/a² + t²/b² ≤ 1  =>  t² ≤ b²(1 - s²/a²)
                temporal_budget_squared = dmax(2)^2 * (1 - (spatial_dist_idx / dmax(1))^2);
                temporal_budget = sqrt(max(0, temporal_budget_squared));
                dt_max_cells = floor(temporal_budget / dt);
            else
                % Outside spatial range, skip
                continue
            end

            % Expand temporally
            for dk = -dt_max_cells:dt_max_cells
                it = idx_t + dk;

                % Check temporal bounds
                if it < 1 || it > nt
                    continue
                end

                % Verify ellipsoidal constraint (space-time metric)
                spacetime_dist_idx = sqrt((di * dx)^2 + (dj * dy)^2 + (dmax(3) * dk * dt)^2);
                combined_constraint = (spatial_dist_idx / dmax(1))^2 + ...
                                      (abs(dk * dt) / dmax(2))^2;

                if combined_constraint > 1
                    continue
                end

                % Check if data exists (not NaN)
                if ~isnan(data.Z(ix, iy, it))
                    candidates_ix = [candidates_ix; ix];
                    candidates_iy = [candidates_iy; iy];
                    candidates_it = [candidates_it; it];
                end
            end
        end
    end

    % Early termination if we have enough candidates
    if length(candidates_ix) >= nmax * 1.5
        break
    end
end

% Handle case where no candidates found
if isempty(candidates_ix)
    psub = [];
    zsub = [];
    dsub = [];
    nsub = 0;
    index = [];
    return
end

%% Step 5: ONLY NOW compute actual distances (for output and final sorting)
n_candidates = length(candidates_ix);

% Extract coordinates for candidates
if ndims(data.Lon) == 3
    lon_candidates = zeros(n_candidates, 1);
    lat_candidates = zeros(n_candidates, 1);
    for i = 1:n_candidates
        lon_candidates(i) = data.Lon(candidates_ix(i), candidates_iy(i), candidates_it(i));
        lat_candidates(i) = data.Lat(candidates_ix(i), candidates_iy(i), candidates_it(i));
    end
else
    % 2D Lon/Lat arrays
    lon_candidates = zeros(n_candidates, 1);
    lat_candidates = zeros(n_candidates, 1);
    for i = 1:n_candidates
        lon_candidates(i) = data.Lon(candidates_ix(i), candidates_iy(i));
        lat_candidates(i) = data.Lat(candidates_ix(i), candidates_iy(i));
    end
end

time_candidates = data.time(candidates_it);

% Compute actual distances
spatial_dist = sqrt((lon_candidates - lon0).^2 + (lat_candidates - lat0).^2);
temporal_dist = abs(time_candidates - t0);

% Apply dmax constraints (should already be satisfied, but verify)
valid_mask = (spatial_dist <= dmax(1)) & (temporal_dist <= dmax(2));
valid_idx = find(valid_mask);

if isempty(valid_idx)
    psub = [];
    zsub = [];
    dsub = [];
    nsub = 0;
    index = [];
    return
end

% Filter to valid points
candidates_ix = candidates_ix(valid_idx);
candidates_iy = candidates_iy(valid_idx);
candidates_it = candidates_it(valid_idx);
lon_candidates = lon_candidates(valid_idx);
lat_candidates = lat_candidates(valid_idx);
time_candidates = time_candidates(valid_idx);
spatial_dist = spatial_dist(valid_idx);
temporal_dist = temporal_dist(valid_idx);

% Compute combined space-time distance
spacetime_dist = spatial_dist + dmax(3) * temporal_dist;

% Sort by space-time distance with deterministic tie-breaking
sort_matrix = [spacetime_dist, candidates_ix, candidates_iy, candidates_it];
[~, sort_idx] = sortrows(sort_matrix);
n_return = min(nmax, length(sort_idx));
select_idx = sort_idx(1:n_return);

%% Prepare outputs
psub = [lon_candidates(select_idx), lat_candidates(select_idx), time_candidates(select_idx), ...
        candidates_ix(select_idx), candidates_iy(select_idx), candidates_it(select_idx)];

% Extract data values
zsub = zeros(n_return, 1);
index = zeros(n_return, 1);
for i = 1:n_return
    ix = candidates_ix(select_idx(i));
    iy = candidates_iy(select_idx(i));
    it = candidates_it(select_idx(i));
    zsub(i) = data.Z(ix, iy, it);
    index(i) = sub2ind([nx, ny, nt], ix, iy, it);
end

dsub = [spatial_dist(select_idx), temporal_dist(select_idx)];
nsub = n_return;

end
