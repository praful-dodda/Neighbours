function [psub, zsub, dsub, nsub, index] = neighbours_stug_index(p0, grid_data, nmax, dmax)
% neighbours_stug_index - Optimized index-based neighbor search for uniformly gridded data
%
% SYNTAX:
%   [psub, zsub, dsub, nsub, index] = neighbours_stug_index(p0, grid_data, nmax, dmax);
%
% INPUT:
%   p0         1 by 3 (or 6)  [lon, lat, time] or [lon, lat, time, idx_x, idx_y, idx_t]
%   grid_data  struct         Structure with grid data
%              Must contain: .x, .y, .time, .Lon, .Lat, .Z
%   nmax       scalar         maximum number of neighbors to return
%   dmax       1 by 3         [max_spatial_dist_km, max_temporal_dist, spacetime_weight_km_per_timeunit]
%                            Set to 0 for no constraint in that dimension
%
% OUTPUT:
%   psub       m by 6         [lon, lat, time, idx_x, idx_y, idx_t] for m neighbors
%   zsub       m by 1         data values at selected neighbors
%   dsub       m by 2         [spatial_distance_km, temporal_distance] from p0
%   nsub       scalar         number of neighbors returned (m ≤ nmax)
%   index      m by 1         linear indices into Z array
%
% NOTES:
%   - Uses Haversine distance for lat/lon (handles spherical Earth geometry)
%   - Handles antimeridian crossing (dateline at ±180°)
%   - Validates grid uniformity (rejects non-uniform grids)
%   - Pre-filters by distance constraints before sorting
%   - Shell-based iterative expansion for efficiency

%% ========================================================================
%% PHASE 0: INPUT VALIDATION AND GRID INFORMATION
%% ========================================================================

% Load grid data if path provided
if ischar(grid_data)
    error('NetCDF loading not implemented. Please pass struct directly.');
end

% Extract grid information
x = grid_data.x(:);      % Longitude grid vector
y = grid_data.y(:);      % Latitude grid vector
t = grid_data.time(:);   % Time grid vector
Z = grid_data.Z;         % Data array [nx, ny, nt]

% Grid dimensions
nx = length(x);
ny = length(y);
nt = length(t);
total_cells = nx * ny * nt;

% Validate that we have enough dimensions
if nx < 2 || ny < 2 || nt < 1
    error('Grid must have at least 2x2x1 dimensions.');
end

%% Validate Grid Uniformity
dx_vec = diff(x);
dy_vec = diff(y);

% Check uniformity in x (tolerance 1%)
if std(dx_vec) / mean(abs(dx_vec)) > 0.01
    error('Grid spacing in x is not uniform (std/mean > 1%%). This function requires uniform grids.');
end

% Check uniformity in y (tolerance 1%)
if std(dy_vec) / mean(abs(dy_vec)) > 0.01
    error('Grid spacing in y is not uniform (std/mean > 1%%). This function requires uniform grids.');
end

% Check uniformity in time (if more than 1 time step)
if nt > 1
    dt_vec = diff(t);
    if std(dt_vec) / mean(abs(dt_vec)) > 0.01
        error('Grid spacing in time is not uniform (std/mean > 1%%). This function requires uniform grids.');
    end
    dt = mean(abs(dt_vec));
else
    dt = 1;  % Single time step: set dummy value
end

% Compute uniform grid spacing
dx = mean(abs(dx_vec));
dy = mean(abs(dy_vec));

%% ========================================================================
%% PHASE 1: INITIALIZATION
%% ========================================================================

% Extract POI coordinates
lon_poi = p0(1);
lat_poi = p0(2);
time_poi = p0(3);

% Check if indices are provided, otherwise compute them
if length(p0) >= 6
    i_center = round(p0(4));
    j_center = round(p0(5));
    k_center = round(p0(6));
else
    % Compute center indices (MATLAB uses 1-based indexing)
    i_center = round((lon_poi - x(1)) / dx) + 1;
    j_center = round((lat_poi - y(1)) / dy) + 1;

    if nt > 1
        k_center = round((time_poi - t(1)) / dt) + 1;
    else
        k_center = 1;
    end

    % Clamp to valid range
    i_center = max(1, min(i_center, nx));
    j_center = max(1, min(j_center, ny));
    k_center = max(1, min(k_center, nt));
end

% Check if grid has enough points
if nmax > total_cells
    warning('nmax (%d) exceeds total grid cells (%d). Returning all valid cells.', nmax, total_cells);
    nmax = total_cells;
end

% Parse distance constraints and spacetime weight
max_spatial_dist = dmax(1);    % km (0 = no limit)
max_temporal_dist = dmax(2);   % time units (0 = no limit)

if length(dmax) >= 3 && dmax(3) > 0
    spacetime_weight = dmax(3);  % km per time unit
else
    spacetime_weight = 0;
end

%% ========================================================================
%% PHASE 2: ESTIMATE INITIAL SEARCH RADIUS
%% ========================================================================

SAFETY_FACTOR = 1.5;
EXPANSION_RATE = 1.3;
MAX_ITERATIONS = 10;
VERBOSE = false;  % Set to true for progress reporting

% Adjust safety factor if near edges
edge_proximity = min([i_center, nx - i_center + 1, j_center, ny - j_center + 1]);
if edge_proximity < 5
    SAFETY_FACTOR = 2.0;  % More aggressive near edges
end

% Volume-based estimation for initial radius
cells_needed = nmax * SAFETY_FACTOR;

% Spherical approximation in normalized space
r_unit = (3 * cells_needed / (4 * pi))^(1/3);

% Compute anisotropy scaling factors
% Account for spacetime weight in temporal dimension
if spacetime_weight > 0
    s_min = min([dx, dy, dt * spacetime_weight]);
    alpha_x = sqrt(dx / s_min);
    alpha_y = sqrt(dy / s_min);
    alpha_t = sqrt((dt * spacetime_weight) / s_min);
else
    s_min = min([dx, dy, dt]);
    alpha_x = sqrt(dx / s_min);
    alpha_y = sqrt(dy / s_min);
    alpha_t = sqrt(dt / s_min);
end

% Initial search radii in index space
radius_i = max(1, ceil(r_unit * alpha_x));
radius_j = max(1, ceil(r_unit * alpha_y));
radius_k = max(1, ceil(r_unit * alpha_t));

% If distance constraints provided, limit initial radius
if max_spatial_dist > 0
    % Convert spatial distance to grid indices (conservative estimate)
    % Use Haversine to estimate how many degrees correspond to max distance
    % At equator: 1 degree ≈ 111 km, but varies with latitude
    % Use conservative estimate based on latitude
    km_per_deg_lon = 111.0 * cos(deg2rad(lat_poi));
    km_per_deg_lat = 111.0;

    if km_per_deg_lon > 0
        max_radius_i = ceil(max_spatial_dist / (km_per_deg_lon * dx));
    else
        max_radius_i = inf;
    end
    max_radius_j = ceil(max_spatial_dist / (km_per_deg_lat * dy));

    radius_i = min(radius_i, max_radius_i);
    radius_j = min(radius_j, max_radius_j);
end

if max_temporal_dist > 0 && nt > 1
    max_radius_temporal = ceil(max_temporal_dist / dt);
    radius_k = min(radius_k, max_radius_temporal);
end

% Special case: 2D grid (no time dimension)
if nt == 1
    radius_k = 0;
end

%% ========================================================================
%% PHASE 3: SETUP FOR ITERATIVE SEARCH
%% ========================================================================

% Determine if ellipsoidal filtering is beneficial
anisotropy_ratios = [dx/dy, dx/dt, dy/dx, dy/dt, dt/dx, dt/dy];
max_anisotropy = max(anisotropy_ratios);
use_ellipsoid = (max_anisotropy > 5);  % Only for highly anisotropic grids

if VERBOSE
    fprintf('Grid: %dx%dx%d = %d cells\n', nx, ny, nt, total_cells);
    fprintf('Spacing: dx=%.3f°, dy=%.3f°, dt=%.3f\n', dx, dy, dt);
    fprintf('Center: (%d,%d,%d)\n', i_center, j_center, k_center);
    fprintf('Initial radius: (%d,%d,%d)\n', radius_i, radius_j, radius_k);
    fprintf('Anisotropy: %.2f, Ellipsoid filter: %d\n', max_anisotropy, use_ellipsoid);
end

% Pre-allocate for worst-case scenario
max_possible = min(total_cells, nmax * 10);  % Upper bound on candidates
candidates_i = zeros(max_possible, 1);
candidates_j = zeros(max_possible, 1);
candidates_k = zeros(max_possible, 1);
n_candidates = 0;

% Iteration tracking
iteration = 0;
prev_radius_i = 0;
prev_radius_j = 0;
prev_radius_k = 0;

%% ========================================================================
%% PHASE 4: ITERATIVE SHELL-BASED SEARCH
%% ========================================================================

while n_candidates < nmax && iteration < MAX_ITERATIONS

    % Compute search bounds with edge handling
    i_min = max(1, i_center - radius_i);
    i_max = min(nx, i_center + radius_i);
    j_min = max(1, j_center - radius_j);
    j_max = min(ny, j_center + radius_j);
    k_min = max(1, k_center - radius_k);
    k_max = min(nt, k_center + radius_k);

    % Determine search region (full box or just shell)
    search_full = (iteration == 0);

    % Pre-compute squared radius values for ellipsoid check
    if use_ellipsoid && radius_i > 0 && radius_j > 0 && radius_k > 0
        radius_i_sq = radius_i^2;
        radius_j_sq = radius_j^2;
        radius_k_sq = radius_k^2;
    end

    % Iterate through the region
    for k = k_min:k_max
        delta_k = k - k_center;
        delta_k_sq = delta_k^2;

        % Shell check for k dimension
        k_is_shell = search_full || (abs(delta_k) >= prev_radius_k);

        for j = j_min:j_max
            delta_j = j - j_center;
            delta_j_sq = delta_j^2;

            % Shell check for j dimension
            j_is_shell = k_is_shell || (abs(delta_j) >= prev_radius_j);

            for i = i_min:i_max
                delta_i = i - i_center;

                % Shell check: skip if this cell was in previous iteration
                if ~search_full && ~j_is_shell && abs(delta_i) < prev_radius_i
                    continue;
                end

                % Ellipsoidal filtering
                if use_ellipsoid && radius_i > 0 && radius_j > 0 && radius_k > 0
                    % Avoid divisions by using squared form
                    % (Δi/ri)² + (Δj/rj)² + (Δk/rk)² ≤ 1
                    % Multiply through: Δi²·rj²·rk² + Δj²·ri²·rk² + Δk²·ri²·rj² ≤ ri²·rj²·rk²
                    threshold = radius_i_sq * radius_j_sq * radius_k_sq;
                    weighted_sum = delta_i^2 * radius_j_sq * radius_k_sq + ...
                                   delta_j_sq * radius_i_sq * radius_k_sq + ...
                                   delta_k_sq * radius_i_sq * radius_j_sq;

                    if weighted_sum > threshold
                        continue;  % Outside ellipsoid
                    end
                end

                % Check if valid data exists (not NaN)
                if isnan(Z(i, j, k))
                    continue;
                end

                % Store candidate (distances computed later)
                n_candidates = n_candidates + 1;

                % Check if we need to expand arrays (safety)
                if n_candidates > length(candidates_i)
                    % Double array size
                    new_size = length(candidates_i) * 2;
                    candidates_i = [candidates_i; zeros(new_size - length(candidates_i), 1)];
                    candidates_j = [candidates_j; zeros(new_size - length(candidates_j), 1)];
                    candidates_k = [candidates_k; zeros(new_size - length(candidates_k), 1)];
                end

                candidates_i(n_candidates) = i;
                candidates_j(n_candidates) = j;
                candidates_k(n_candidates) = k;

            end
        end
    end

    % Check termination conditions
    if n_candidates >= nmax
        if VERBOSE
            fprintf('Iteration %d: Found sufficient candidates (%d >= %d)\n', ...
                    iteration, n_candidates, nmax);
        end
        break;
    end

    % Check if we've hit grid boundaries and can't expand further
    if i_min == 1 && i_max == nx && ...
       j_min == 1 && j_max == ny && ...
       k_min == 1 && k_max == nt
        if VERBOSE
            fprintf('Iteration %d: Reached full grid, stopping.\n', iteration);
        end
        break;
    end

    % Early termination: if we have enough candidates and further expansion
    % is unlikely to improve results
    if n_candidates >= nmax * 0.8 && iteration > 0
        % Estimate if next expansion will find closer neighbors
        current_max_index_dist = sqrt(radius_i^2 + radius_j^2 + radius_k^2);
        next_radius_i = ceil(radius_i * EXPANSION_RATE);
        next_radius_j = ceil(radius_j * EXPANSION_RATE);
        next_radius_k = ceil(radius_k * EXPANSION_RATE);
        next_max_index_dist = sqrt(next_radius_i^2 + next_radius_j^2 + next_radius_k^2);

        if next_max_index_dist > current_max_index_dist * 2
            if VERBOSE
                fprintf('Iteration %d: Early termination (expansion unlikely to help)\n', iteration);
            end
            break;
        end
    end

    % Progress reporting
    if VERBOSE
        fprintf('Iteration %d: Found %d/%d candidates, radius=(%d,%d,%d)\n', ...
                iteration, n_candidates, nmax, radius_i, radius_j, radius_k);
    end

    % Store current bounds for next iteration
    prev_radius_i = radius_i;
    prev_radius_j = radius_j;
    prev_radius_k = radius_k;

    % Adaptive expansion based on deficit
    deficit = nmax - n_candidates;
    deficit_ratio = deficit / max(nmax, 1);

    if deficit_ratio > 0.5
        % Large deficit: expand more aggressively
        expansion = 1.5;
    elseif deficit_ratio > 0.2
        % Moderate deficit: standard expansion
        expansion = EXPANSION_RATE;
    else
        % Small deficit: gentle expansion
        expansion = 1.15;
    end

    % Apply expansion with distance constraint checks
    radius_i = ceil(radius_i * expansion);
    radius_j = ceil(radius_j * expansion);

    if nt > 1
        radius_k = ceil(radius_k * expansion);
    end

    % Respect distance constraints
    if max_spatial_dist > 0
        km_per_deg_lon = 111.0 * cos(deg2rad(lat_poi));
        km_per_deg_lat = 111.0;

        if km_per_deg_lon > 0
            max_radius_i = ceil(max_spatial_dist / (km_per_deg_lon * dx));
            radius_i = min(radius_i, max_radius_i);
        end

        max_radius_j = ceil(max_spatial_dist / (km_per_deg_lat * dy));
        radius_j = min(radius_j, max_radius_j);
    end

    if max_temporal_dist > 0 && nt > 1
        max_radius_temporal = ceil(max_temporal_dist / dt);
        radius_k = min(radius_k, max_radius_temporal);
    end

    iteration = iteration + 1;
end

% Trim candidates to actual size
candidates_i = candidates_i(1:n_candidates);
candidates_j = candidates_j(1:n_candidates);
candidates_k = candidates_k(1:n_candidates);

% Check if any candidates found
if n_candidates == 0
    warning('No valid neighbors found in search region.');
    psub = [];
    zsub = [];
    dsub = [];
    nsub = 0;
    index = [];
    return;
end

%% ========================================================================
%% PHASE 5: COMPUTE DISTANCES AND PRE-FILTER
%% ========================================================================

% Convert indices to coordinates (vectorized)
lon_candidates = x(candidates_i);
lat_candidates = y(candidates_j);
time_candidates = t(candidates_k);

% Calculate spatial distances using Haversine formula with antimeridian handling
spatial_dist = haversine_vectorized(lon_poi, lat_poi, lon_candidates, lat_candidates);

% Calculate temporal distances
temporal_dist = abs(time_candidates - time_poi);

% PRE-FILTER by distance constraints BEFORE sorting
% This reduces the array size for sorting, improving performance
pre_filter_mask = true(n_candidates, 1);

if max_spatial_dist > 0
    pre_filter_mask = pre_filter_mask & (spatial_dist <= max_spatial_dist);
end

if max_temporal_dist > 0
    pre_filter_mask = pre_filter_mask & (temporal_dist <= max_temporal_dist);
end

% Apply pre-filter
n_filtered = sum(pre_filter_mask);

if n_filtered == 0
    warning('No neighbors satisfy distance constraints (dmax=[%.1f km, %.1f time]).',...
            max_spatial_dist, max_temporal_dist);
    psub = [];
    zsub = [];
    dsub = [];
    nsub = 0;
    index = [];
    return;
end

% Apply filter to all arrays
candidates_i = candidates_i(pre_filter_mask);
candidates_j = candidates_j(pre_filter_mask);
candidates_k = candidates_k(pre_filter_mask);
lon_candidates = lon_candidates(pre_filter_mask);
lat_candidates = lat_candidates(pre_filter_mask);
time_candidates = time_candidates(pre_filter_mask);
spatial_dist = spatial_dist(pre_filter_mask);
temporal_dist = temporal_dist(pre_filter_mask);

if VERBOSE && n_filtered < n_candidates
    fprintf('Pre-filtered: %d -> %d candidates (%.1f%% reduction)\n', ...
            n_candidates, n_filtered, 100*(1 - n_filtered/n_candidates));
end

n_candidates = n_filtered;

%% ========================================================================
%% PHASE 6: SELECT NEAREST NEIGHBORS
%% ========================================================================

% Apply spacetime metric for combined distance
if spacetime_weight > 0
    % Combined spacetime distance with proper weighting
    combined_dist = sqrt(spatial_dist.^2 + (spacetime_weight * temporal_dist).^2);
else
    % Use spatial distance only for sorting
    combined_dist = spatial_dist;
end

% Sort by combined distance and select top nmax
n_select = min(nmax, n_candidates);

[~, sort_idx] = sort(combined_dist);
selected_idx = sort_idx(1:n_select);
nsub = n_select;

%% ========================================================================
%% PHASE 7: PREPARE OUTPUT (VECTORIZED)
%% ========================================================================

% Extract final results
final_i = candidates_i(selected_idx);
final_j = candidates_j(selected_idx);
final_k = candidates_k(selected_idx);

% Build output arrays (pre-allocated and vectorized)
psub = zeros(nsub, 6);
psub(:, 1) = x(final_i);           % Longitude
psub(:, 2) = y(final_j);           % Latitude
psub(:, 3) = t(final_k);           % Time
psub(:, 4) = final_i;              % Index x
psub(:, 5) = final_j;              % Index y
psub(:, 6) = final_k;              % Index t

% Extract data values (vectorized using linear indexing)
% Fast linear index calculation (MATLAB column-major order)
linear_idx = final_i + (final_j - 1) * nx + (final_k - 1) * nx * ny;
zsub = Z(linear_idx);

% Distance matrix
dsub = zeros(nsub, 2);
dsub(:, 1) = spatial_dist(selected_idx);   % Spatial distance (km)
dsub(:, 2) = temporal_dist(selected_idx);  % Temporal distance

% Return linear indices
index = linear_idx;

if VERBOSE
    fprintf('Completed: %d neighbors found in %d iterations\n', nsub, iteration);
end

end

%% ========================================================================
%% HELPER FUNCTIONS
%% ========================================================================

function dist = haversine_vectorized(lon1, lat1, lon2, lat2)
% HAVERSINE_VECTORIZED - Compute great circle distance on Earth with antimeridian handling
%
% Vectorized implementation for arrays of points.
% Handles the antimeridian (dateline) crossing at ±180° longitude.
% Returns distance in kilometers.
%
% INPUT:
%   lon1, lat1: scalar (reference point in degrees)
%   lon2, lat2: vectors (candidate points in degrees)
%
% OUTPUT:
%   dist: vector of distances in km
%
% NOTES:
%   - Properly handles longitude wrapping at ±180°
%   - Uses WGS84 mean Earth radius (6371 km)
%   - Numerically stable for antipodal points

    R = 6371.0;  % Earth radius in kilometers (mean WGS84)

    % Convert to radians
    lon1_rad = deg2rad(lon1);
    lat1_rad = deg2rad(lat1);
    lon2_rad = deg2rad(lon2);
    lat2_rad = deg2rad(lat2);

    % Calculate differences
    dlat = lat2_rad - lat1_rad;
    dlon = lon2_rad - lon1_rad;

    % Handle antimeridian crossing
    % Wrap dlon to [-π, π] for shortest path
    % If dlon > π, we should go the other way around the globe
    dlon = mod(dlon + pi, 2*pi) - pi;

    % Haversine formula
    % Numerically stable for all distances including antipodal points
    a = sin(dlat/2).^2 + cos(lat1_rad) .* cos(lat2_rad) .* sin(dlon/2).^2;

    % Use atan2 for numerical stability (better than asin near antipodal points)
    c = 2 * atan2(sqrt(a), sqrt(1-a));

    dist = R * c;
end
