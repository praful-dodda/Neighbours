% example_stg_to_stug_usage.m
% Examples of how to use stg_to_stug conversion function

%% Example 1: Basic usage with auto resolution
fprintf('=== Example 1: Basic Usage ===\n');

% Assuming you have loaded: KS.softdata
% grid_data = stg_to_stug(KS.softdata);

% This will:
% - Auto-calculate optimal resolution based on station density
% - Use linear interpolation
% - Show progress
% - Take ~2-10 seconds for 78k stations × 19 times

%% Example 2: Specify custom resolution
fprintf('\n=== Example 2: Custom Resolution ===\n');

% High resolution (smaller grid spacing = more grid points)
% grid_data_hires = stg_to_stug(KS.softdata, 0.05);  % 0.05 degree spacing

% Low resolution (larger grid spacing = fewer grid points, faster)
% grid_data_lores = stg_to_stug(KS.softdata, 0.2);   % 0.2 degree spacing

%% Example 3: Advanced options
fprintf('\n=== Example 3: Advanced Options ===\n');

% grid_data = stg_to_stug(KS.softdata, ...
%     'resolution', 0.1, ...           % Grid spacing
%     'method', 'natural', ...         % Natural neighbor interpolation
%     'extrapolation', 'none', ...     % NaN outside convex hull
%     'verbose', true, ...             % Show progress
%     'includeVariance', true);        % Also interpolate variance

%% Example 4: Custom bounds
fprintf('\n=== Example 4: Custom Spatial Bounds ===\n');

% Specify exact region to grid
% grid_data = stg_to_stug(KS.softdata, ...
%     'resolution', 0.1, ...
%     'bounds', [-100, -90, 25, 35]);  % [lon_min lon_max lat_min lat_max]

%% Example 5: Complete workflow with neighbor search
fprintf('\n=== Example 5: Complete Workflow ===\n');

% Step 1: Convert STG to STUG (one-time cost)
fprintf('Step 1: Converting STG to STUG...\n');
% grid_data = stg_to_stug(KS.softdata, 0.1);

% Step 2: Use with neighbours_stug_optimized
fprintf('Step 2: Using with neighbour search...\n');
% p0 = [-95.0, 30.0, 2017.5];  % Point of interest [lon, lat, time]
% nmax = 50;
% dmax = [500.0, 365, 0.01];  % [spatial_km, temporal_days, weight]
%
% [psub, zsub, dsub, nsub] = neighbours_stug_optimized(p0, grid_data, nmax, dmax);

% Step 3: Or use with neighbours_stug_index (for Haversine distance)
% [psub, zsub, dsub, nsub] = neighbours_stug_index(p0, grid_data, nmax, dmax);

%% Example 6: Resolution selection guide
fprintf('\n=== Example 6: Resolution Selection Guide ===\n');

% Rule of thumb for resolution selection:
%
% 1. Check your station density:
%    mean_nn_dist = estimate_station_spacing(KS.softdata.sMS);
%
% 2. Choose resolution based on use case:
%    - Preserve detail: resolution = 0.3 × mean_nn_dist
%    - Balance speed/detail: resolution = 0.5 × mean_nn_dist (default)
%    - Faster processing: resolution = 1.0 × mean_nn_dist
%
% 3. Memory consideration:
%    grid_points = (lon_range/res) × (lat_range/res) × n_times
%    memory_MB = grid_points × 4 / 1024^2
%
% Example: 50° × 30° region, 19 times
%    res=0.05: 1000×600×19 = 11.4M points → 45 MB
%    res=0.10:  500×300×19 =  2.9M points → 11 MB
%    res=0.20:  250×150×19 =  0.7M points →  3 MB

%% Helper function: Estimate station spacing
function mean_dist = estimate_station_spacing(sMS)
    % Estimate mean nearest neighbor distance
    n = size(sMS, 1);
    if n > 1000
        % Sample for speed
        sample_idx = randperm(n, 1000);
        sample_coords = sMS(sample_idx, :);
        search_coords = sMS;
    else
        sample_coords = sMS;
        search_coords = sMS;
    end

    nn_dists = zeros(size(sample_coords, 1), 1);
    for i = 1:size(sample_coords, 1)
        dists = sqrt(sum((search_coords - repmat(sample_coords(i, :), size(search_coords, 1), 1)).^2, 2));
        dists(dists == 0) = inf;
        nn_dists(i) = min(dists);
    end

    mean_dist = mean(nn_dists);
    fprintf('Estimated mean station spacing: %.4f degrees (%.1f km at mid-latitudes)\n', ...
            mean_dist, mean_dist * 111);
end

%% Example 7: Batch processing multiple datasets
fprintf('\n=== Example 7: Batch Processing ===\n');

% Process multiple datasets with same settings
% datasets = {KS.softdata, other_data1, other_data2};
% grid_datasets = cell(size(datasets));
%
% for i = 1:length(datasets)
%     fprintf('Processing dataset %d/%d...\n', i, length(datasets));
%     grid_datasets{i} = stg_to_stug(datasets{i}, ...
%         'resolution', 0.1, ...
%         'verbose', false);  % Suppress individual progress
% end

%% Example 8: Quality checking the conversion
fprintf('\n=== Example 8: Quality Checking ===\n');

% After conversion, check interpolation quality:
% grid_data = stg_to_stug(KS.softdata, 0.1);
%
% % 1. Check fill ratio
% fprintf('Fill ratio: %.1f%%\n', grid_data.metadata.fill_ratio * 100);
%
% % 2. Visualize one time slice
% figure;
% t_idx = 10;  % Choose a time index
% imagesc(grid_data.x, grid_data.y, grid_data.Z(:,:,t_idx)');
% axis xy; colorbar;
% title(sprintf('Interpolated grid at time %d', t_idx));
% xlabel('Longitude'); ylabel('Latitude');
%
% % 3. Overlay original station locations
% hold on;
% scatter(KS.softdata.sMS(:,1), KS.softdata.sMS(:,2), 10, 'r', 'filled');
% legend('Interpolated grid', 'Original stations');

%% Example 9: Memory-efficient processing for large datasets
fprintf('\n=== Example 9: Memory-Efficient Processing ===\n');

% For very large datasets, process in chunks:
% function grid_data = convert_large_stg(stg_data, resolution)
%     % Split into spatial chunks
%     lon_range = [min(stg_data.sMS(:,1)), max(stg_data.sMS(:,1))];
%     n_chunks = 4;  % Process in 4 spatial chunks
%
%     chunk_bounds = linspace(lon_range(1), lon_range(2), n_chunks + 1);
%
%     for i = 1:n_chunks
%         % Extract chunk
%         chunk_mask = stg_data.sMS(:,1) >= chunk_bounds(i) & ...
%                      stg_data.sMS(:,1) < chunk_bounds(i+1);
%
%         chunk_data.sMS = stg_data.sMS(chunk_mask, :);
%         chunk_data.tME = stg_data.tME;
%         chunk_data.Xms = stg_data.Xms(chunk_mask, :);
%
%         % Convert chunk
%         grid_chunk = stg_to_stug(chunk_data, 'resolution', resolution);
%
%         % Merge chunks (implementation depends on use case)
%         if i == 1
%             grid_data = grid_chunk;
%         else
%             % Concatenate grids
%         end
%     end
% end

%% Comparison: STG vs STUG approaches
fprintf('\n=== STG vs STUG Comparison ===\n');

% Direct STG approach (NO conversion):
% - Use: neighbours_stg(p0, KS.softdata, nmax, dmax)
% - Time: 5-20 ms per query
% - Memory: Original data only
% - Accuracy: Exact station measurements
% - Best for: Single/few queries, need exact values

% STUG approach (WITH conversion):
% - Step 1: grid_data = stg_to_stug(KS.softdata, 0.1)  [2-10 seconds]
% - Step 2: neighbours_stug_index(p0, grid_data, nmax, dmax) [5-10 ms]
% - Time: One-time conversion + fast queries
% - Memory: 2-5x original (depends on resolution)
% - Accuracy: Interpolation introduces error
% - Best for: Many queries (100+), want gap filling

fprintf('\nRecommendation: Use STG approach unless you need >100 queries\n');
fprintf('on the same dataset or require gap-filled gridded data.\n');
