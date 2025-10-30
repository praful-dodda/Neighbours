% compare_neighbours_implementations.m
% Comprehensive comparison between neighbours_stug_index.m and neighbours.m
%
% This script:
% 1. Creates test grids with various characteristics
% 2. Calls both functions with equivalent inputs
% 3. Compares outputs for correctness
% 4. Measures and compares execution times
% 5. Validates that both functions find the same neighbors

clear; close all;

fprintf('=== Comparison: neighbours_stug_index vs neighbours.m ===\n\n');

%% Test Configuration
% Realistic dimensions for global monthly soft data:
% - High spatial resolution (nx, ny: 50-300)
% - Low temporal dimension (nt: 12-24 months)
% - Various NaN ratios representing data availability

test_configs = {
    struct('name', 'Small Grid - 12 months (10% NaN)', 'nx', 50, 'ny', 50, 'nt', 12, 'nan_ratio', 0.1, 'nmax', 15)
    struct('name', 'Medium Grid - 24 months (30% NaN)', 'nx', 100, 'ny', 100, 'nt', 24, 'nan_ratio', 0.3, 'nmax', 20)
    struct('name', 'Medium Grid - 18 months (70% NaN)', 'nx', 100, 'ny', 100, 'nt', 18, 'nan_ratio', 0.7, 'nmax', 20)
    struct('name', 'Large Grid - 24 months (30% NaN)', 'nx', 200, 'ny', 200, 'nt', 24, 'nan_ratio', 0.3, 'nmax', 25)
    struct('name', 'Very Sparse - 24 months (95% NaN)', 'nx', 150, 'ny', 150, 'nt', 24, 'nan_ratio', 0.95, 'nmax', 10)
    struct('name', 'Very Large Grid - 24 months (50% NaN)', 'nx', 300, 'ny', 300, 'nt', 24, 'nan_ratio', 0.5, 'nmax', 30)
    struct('name', 'Very Large Grid - 24 months (0% NaN)', 'nx', 300, 'ny', 300, 'nt', 24, 'nan_ratio', 0, 'nmax', 30)
};

% Storage for results
n_tests = length(test_configs);
results = struct();

%% Run all test configurations
for test_idx = 1:n_tests
    cfg = test_configs{test_idx};

    fprintf('Test %d/%d: %s\n', test_idx, n_tests, cfg.name);
    fprintf('  Grid: %dx%dx%d, NaN ratio: %.1f%%, nmax: %d\n', ...
        cfg.nx, cfg.ny, cfg.nt, cfg.nan_ratio*100, cfg.nmax);

    %% Create uniform grid
    grid_data = create_test_grid(cfg.nx, cfg.ny, cfg.nt, cfg.nan_ratio);

    %% Convert grid to vector format for neighbours.m
    [c_vector, Z_vector] = grid_to_vector(grid_data);

    fprintf('  Grid created: %d total points, %d non-NaN values\n', ...
        cfg.nx*cfg.ny*cfg.nt, length(Z_vector));

    %% Test point (center of grid)
    p0_lon = (grid_data.x(1) + grid_data.x(end)) / 2;
    p0_lat = (grid_data.y(1) + grid_data.y(end)) / 2;
    p0_time = (grid_data.time(1) + grid_data.time(end)) / 2;
    p0 = [p0_lon, p0_lat, p0_time];

    %% Parameters
    nmax = cfg.nmax;
    dmax = [5.0, 10.0, 0.1];  % [spatial_max, temporal_max, metric]

    %% Call neighbours.m (reference implementation)
    fprintf('  Running neighbours.m (reference)...\n');
    tic;
    [csub_ref, Zsub_ref, dsub_ref, nsub_ref, index_ref] = ...
        neighbours(p0, c_vector, Z_vector, nmax, dmax);
    time_ref = toc;
    fprintf('    Time: %.4f sec, Found: %d neighbors\n', time_ref, nsub_ref);

    %% Call neighbours_stug_index.m
    fprintf('  Running neighbours_stug_index.m...\n');
    tic;
    [psub_opt, zsub_opt, dsub_opt, nsub_opt, index_opt] = ...
        neighbours_stug_index(p0, grid_data, nmax, dmax);
    time_opt = toc;
    fprintf('    Time: %.4f sec, Found: %d neighbors\n', time_opt, nsub_opt);

    %% Compare results
    fprintf('  Comparing results...\n');

    % Extract just spatial coordinates from optimized version
    csub_opt = psub_opt(:, 1:3);  % [lon, lat, time]

    % Compare number of neighbors
    if nsub_ref == nsub_opt
        fprintf('    ✓ Same number of neighbors: %d\n', nsub_ref);
    else
        fprintf('    ✗ DIFFERENT number of neighbors: ref=%d, opt=%d\n', nsub_ref, nsub_opt);
    end

    % Compare neighbor sets (allowing for different ordering)
    if nsub_ref > 0 && nsub_opt > 0
        [match_quality, detailed_comparison] = compare_neighbor_sets(...
            csub_ref, Zsub_ref, dsub_ref, ...
            csub_opt, zsub_opt, dsub_opt, ...
            p0, dmax);

        fprintf('    Neighbor set match: %.1f%%\n', match_quality.percent_match);

        if match_quality.percent_match >= 95
            fprintf('    ✓ Neighbor sets match well\n');
        elseif match_quality.percent_match >= 80
            fprintf('    ⚠ Neighbor sets mostly match\n');
        else
            fprintf('    ✗ SIGNIFICANT DIFFERENCE in neighbor sets\n');
        end

        % Compare distances
        if match_quality.distances_match
            fprintf('    ✓ Distances match for common neighbors\n');
        else
            fprintf('    ⚠ Some distance discrepancies (max diff: %.6f)\n', ...
                match_quality.max_distance_diff);
        end

        % Compare values
        if match_quality.values_match
            fprintf('    ✓ Values match for common neighbors\n');
        else
            fprintf('    ⚠ Some value discrepancies\n');
        end

        % Compare distance distributions (more robust than exact neighbor matching)
        fprintf('  Distance Distribution Comparison:\n');
        [dist_comparison] = compare_distance_distributions(dsub_ref, dsub_opt, dmax);
        fprintf('    Mean distance diff:  %.6f (%.2f%%)\n', ...
            dist_comparison.mean_diff, dist_comparison.mean_percent_diff);
        fprintf('    Median distance diff: %.6f (%.2f%%)\n', ...
            dist_comparison.median_diff, dist_comparison.median_percent_diff);
        fprintf('    Distribution KS test: p-value = %.4f', dist_comparison.ks_pvalue);
        if dist_comparison.ks_pvalue > 0.05
            fprintf(' ✓ (distributions match)\n');
        else
            fprintf(' ⚠ (distributions differ)\n');
        end
    end

    %% Performance comparison
    speedup = time_ref / time_opt;
    fprintf('  Performance:\n');
    fprintf('    neighbours.m:           %.4f sec\n', time_ref);
    fprintf('    neighbours_stug_opt:    %.4f sec\n', time_opt);
    fprintf('    Speedup:                %.2fx\n', speedup);

    if speedup > 1
        fprintf('    ✓ Optimized version is faster\n');
    elseif speedup < 0.9
        fprintf('    ⚠ Optimized version is SLOWER (unexpected)\n');
    else
        fprintf('    ≈ Similar performance\n');
    end

    %% Store results
    results(test_idx).config = cfg;
    results(test_idx).nsub_ref = nsub_ref;
    results(test_idx).nsub_opt = nsub_opt;
    results(test_idx).time_ref = time_ref;
    results(test_idx).time_opt = time_opt;
    results(test_idx).speedup = speedup;
    results(test_idx).match_quality = match_quality;
    results(test_idx).csub_ref = csub_ref;  % Store for visualization
    results(test_idx).csub_opt = csub_opt;
    results(test_idx).p0 = p0;
    if nsub_ref > 0 && nsub_opt > 0
        results(test_idx).dist_comparison = dist_comparison;
    end

    fprintf('\n');
end

%% Summary Report
fprintf('=== SUMMARY ===\n\n');

fprintf('Performance Summary:\n');
fprintf('%-30s | %8s | %8s | %8s | %8s\n', ...
    'Test', 'Ref(s)', 'Opt(s)', 'Speedup', 'Match%%');
fprintf('%s\n', repmat('-', 80, 1));

for i = 1:n_tests
    fprintf('%-30s | %8.4f | %8.4f | %7.2fx | %7.1f%%\n', ...
        results(i).config.name, ...
        results(i).time_ref, ...
        results(i).time_opt, ...
        results(i).speedup, ...
        results(i).match_quality.percent_match);
end

fprintf('\n');

% Overall statistics
all_speedups = [results.speedup];
all_matches = [results.match_quality];
avg_speedup = mean(all_speedups);
median_speedup = median(all_speedups);
min_speedup = min(all_speedups);
max_speedup = max(all_speedups);
avg_match = mean([all_matches.percent_match]);

fprintf('Overall Statistics:\n');
fprintf('  Average speedup:  %.2fx\n', avg_speedup);
fprintf('  Median speedup:   %.2fx\n', median_speedup);
fprintf('  Min speedup:      %.2fx\n', min_speedup);
fprintf('  Max speedup:      %.2fx\n', max_speedup);
fprintf('  Average match:    %.1f%%\n', avg_match);

fprintf('\n');

% Correctness validation
if avg_match >= 95
    fprintf('✓ CORRECTNESS: Both functions produce equivalent results\n');
else
    fprintf('✗ CORRECTNESS: Significant differences detected\n');
end

if avg_speedup > 1
    fprintf('✓ PERFORMANCE: Optimized version is %.2fx faster on average\n', avg_speedup);
else
    fprintf('⚠ PERFORMANCE: Optimized version not faster\n');
end

fprintf('\n=== Comparison Complete ===\n');

%% Save results for visualization
fprintf('\nSaving results to compare_results.mat...\n');
save('compare_results.mat', 'results', 'test_configs');
fprintf('Results saved. Use visualize_accuracy.m to generate accuracy plots.\n');

%% Helper Functions

function grid_data = create_test_grid(nx, ny, nt, nan_ratio)
    % Create a uniform grid with specified dimensions and NaN ratio

    grid_data.x = linspace(-100, -80, nx)';
    grid_data.y = linspace(25, 35, ny)';
    grid_data.time = linspace(2400, 2500, nt)';

    [X_grid, Y_grid] = meshgrid(grid_data.x, grid_data.y);
    grid_data.Lon = X_grid';
    grid_data.Lat = Y_grid';

    % Create data with spatial and temporal patterns
    rng(42);  % Reproducible
    grid_data.Z = zeros(nx, ny, nt);
    for it = 1:nt
        for ix = 1:nx
            for iy = 1:ny
                % Some pattern plus noise
                grid_data.Z(ix, iy, it) = ...
                    sin(ix/10) * cos(iy/10) * exp(-it/50) + 0.1*randn();
            end
        end
    end

    % Add NaNs
    if nan_ratio > 0
        nan_mask = rand(nx, ny, nt) < nan_ratio;
        grid_data.Z(nan_mask) = NaN;
    end
end

function [c_vector, Z_vector] = grid_to_vector(grid_data)
    % Convert grid format to vector format for neighbours.m

    nx = length(grid_data.x);
    ny = length(grid_data.y);
    nt = length(grid_data.time);

    % Pre-allocate
    c_all = zeros(nx*ny*nt, 3);
    Z_all = zeros(nx*ny*nt, 1);

    idx = 1;
    for it = 1:nt
        for iy = 1:ny
            for ix = 1:nx
                c_all(idx, 1) = grid_data.Lon(ix, iy);
                c_all(idx, 2) = grid_data.Lat(ix, iy);
                c_all(idx, 3) = grid_data.time(it);
                Z_all(idx) = grid_data.Z(ix, iy, it);
                idx = idx + 1;
            end
        end
    end

    % Remove NaN values
    valid_idx = ~isnan(Z_all);
    c_vector = c_all(valid_idx, :);
    Z_vector = Z_all(valid_idx);
end

function [match_quality, details] = compare_neighbor_sets(...
    csub_ref, Zsub_ref, dsub_ref, ...
    csub_opt, zsub_opt, dsub_opt, ...
    p0, dmax)
    % Compare two sets of neighbors, allowing for different ordering

    match_quality = struct();
    details = struct();

    n_ref = size(csub_ref, 1);
    n_opt = size(csub_opt, 1);

    if n_ref == 0 || n_opt == 0
        match_quality.percent_match = 0;
        match_quality.distances_match = false;
        match_quality.values_match = false;
        match_quality.max_distance_diff = NaN;
        return;
    end

    % For each reference neighbor, find if it exists in optimized set
    tolerance_spatial = 1e-6;  % Tolerance for coordinate matching
    tolerance_distance = 1e-6; % Tolerance for distance matching
    tolerance_value = 1e-9;    % Tolerance for value matching

    matched_ref = false(n_ref, 1);
    matched_opt = false(n_opt, 1);
    distance_diffs = zeros(n_ref, 1);
    value_diffs = zeros(n_ref, 1);

    for i = 1:n_ref
        % Find this neighbor in optimized set
        coord_diff = sqrt(sum((csub_opt - repmat(csub_ref(i,:), n_opt, 1)).^2, 2));
        [min_diff, match_idx] = min(coord_diff);

        if min_diff < tolerance_spatial
            % Found a match
            matched_ref(i) = true;
            matched_opt(match_idx) = true;

            % Compare distances
            if size(dsub_ref, 2) == 1
                % Old format: single distance
                dist_ref = dsub_ref(i);
                dist_opt = sqrt(dsub_opt(match_idx,1)^2 + (dmax(3)*dsub_opt(match_idx,2))^2);
            else
                % New format: [spatial, temporal]
                dist_ref = sqrt(dsub_ref(i,1)^2 + (dmax(3)*dsub_ref(i,2))^2);
                dist_opt = sqrt(dsub_opt(match_idx,1)^2 + (dmax(3)*dsub_opt(match_idx,2))^2);
            end
            distance_diffs(i) = abs(dist_ref - dist_opt);

            % Compare values
            value_diffs(i) = abs(Zsub_ref(i) - zsub_opt(match_idx));
        end
    end

    % Calculate match statistics
    n_matched = sum(matched_ref);
    match_quality.percent_match = (n_matched / max(n_ref, n_opt)) * 100;
    match_quality.n_ref = n_ref;
    match_quality.n_opt = n_opt;
    match_quality.n_matched = n_matched;

    % Distance comparison
    max_distance_diff = max(distance_diffs(matched_ref));
    match_quality.max_distance_diff = max_distance_diff;
    match_quality.distances_match = max_distance_diff < tolerance_distance;

    % Value comparison
    max_value_diff = max(value_diffs(matched_ref));
    match_quality.max_value_diff = max_value_diff;
    match_quality.values_match = max_value_diff < tolerance_value;

    % Store details
    details.matched_ref = matched_ref;
    details.matched_opt = matched_opt;
    details.distance_diffs = distance_diffs;
    details.value_diffs = value_diffs;
end

function [dist_comparison] = compare_distance_distributions(dsub_ref, dsub_opt, dmax)
    % Compare statistical properties of distance distributions
    % This is more robust than exact neighbor matching for uniform grids

    dist_comparison = struct();

    % Extract spatial distances (first column)
    if size(dsub_ref, 2) >= 1
        spatial_ref = dsub_ref(:, 1);
        spatial_opt = dsub_opt(:, 1);
    else
        spatial_ref = dsub_ref;
        spatial_opt = dsub_opt;
    end

    % Compute combined space-time distance
    if size(dsub_ref, 2) >= 2
        combined_ref = dsub_ref(:,1) + dmax(3) * dsub_ref(:,2);
        combined_opt = dsub_opt(:,1) + dmax(3) * dsub_opt(:,2);
    else
        combined_ref = spatial_ref;
        combined_opt = spatial_opt;
    end

    % Compare means
    mean_ref = mean(combined_ref);
    mean_opt = mean(combined_opt);
    dist_comparison.mean_diff = abs(mean_ref - mean_opt);
    dist_comparison.mean_percent_diff = 100 * dist_comparison.mean_diff / max(mean_ref, 1e-10);

    % Compare medians
    median_ref = median(combined_ref);
    median_opt = median(combined_opt);
    dist_comparison.median_diff = abs(median_ref - median_opt);
    dist_comparison.median_percent_diff = 100 * dist_comparison.median_diff / max(median_ref, 1e-10);

    % Compare standard deviations
    std_ref = std(combined_ref);
    std_opt = std(combined_opt);
    dist_comparison.std_diff = abs(std_ref - std_opt);

    % Kolmogorov-Smirnov test (tests if distributions are from same continuous distribution)
    [~, dist_comparison.ks_pvalue] = kstest2(combined_ref, combined_opt);

    % Store raw stats
    dist_comparison.mean_ref = mean_ref;
    dist_comparison.mean_opt = mean_opt;
    dist_comparison.median_ref = median_ref;
    dist_comparison.median_opt = median_opt;
    dist_comparison.std_ref = std_ref;
    dist_comparison.std_opt = std_opt;
end
