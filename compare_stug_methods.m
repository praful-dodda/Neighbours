% compare_stug_methods.m
% Compare all three neighbour search methods for uniformly gridded data:
% 1. neighbours.m (reference implementation)
% 2. neighbours_stug_optimized.m (adaptive strategy with shell expansion)
% 3. neighbours_stug_index.m (pure index-based search)
%
% This script validates correctness and measures performance across
% various grid configurations and data sparsity levels.

clear; close all;
fprintf('=== Comparison: Three STUG Methods ===\n');
fprintf('  1. neighbours.m (reference)\n');
fprintf('  2. neighbours_stug_optimized.m (adaptive)\n');
fprintf('  3. neighbours_stug_index.m (index-based)\n\n');

%% Test scenarios
% Format: {nx, ny, nt, nan_ratio, nmax, description}
scenarios = {
    {100, 100, 12, 0.1, 20, 'Small grid, 10% NaN, 12 months'}
    {200, 200, 24, 0.3, 25, 'Medium grid, 30% NaN, 24 months'}
    {200, 200, 24, 0.7, 25, 'Medium grid, 70% NaN, 24 months'}
    {300, 300, 24, 0.3, 30, 'Large grid, 30% NaN, 24 months'}
    {200, 200, 24, 0.95, 15, 'Sparse grid, 95% NaN, 24 months'}
    {300, 300, 24, 0.0, 30, 'Large dense grid, 0% NaN, 24 months'}
};

n_scenarios = length(scenarios);
results = struct();

%% Run comparison tests
for s = 1:n_scenarios
    % Extract scenario parameters
    scenario = scenarios{s};
    nx = scenario{1};
    ny = scenario{2};
    nt = scenario{3};
    nan_ratio = scenario{4};
    nmax = scenario{5};
    description = scenario{6};

    fprintf('\n=== Test %d/%d: %s ===\n', s, n_scenarios, description);
    fprintf('Grid: %dx%dx%d = %d points\n', nx, ny, nt, nx*ny*nt);
    fprintf('NaN ratio: %.1f%%, nmax: %d\n', nan_ratio*100, nmax);

    % Create test grid
    grid_data.x = linspace(-180, 180, nx)';
    grid_data.y = linspace(-90, 90, ny)';
    grid_data.time = (1:nt)';

    [X, Y] = meshgrid(grid_data.x, grid_data.y);
    grid_data.Lon = X';
    grid_data.Lat = Y';

    % Create synthetic data
    rng(42 + s);  % For reproducibility
    grid_data.Z = randn(nx, ny, nt);

    % Add NaN values
    nan_mask = rand(nx, ny, nt) < nan_ratio;
    grid_data.Z(nan_mask) = NaN;

    n_valid = sum(~isnan(grid_data.Z(:)));
    fprintf('Valid points: %d (%.1f%%)\n', n_valid, 100*(1-nan_ratio));

    % Test point (center of grid, middle time)
    p0 = [0, 0, ceil(nt/2)];
    dmax = [5.0, 6.0, 0.5];

    fprintf('\nTest point: [%.1f, %.1f, %.1f]\n', p0(1), p0(2), p0(3));
    fprintf('dmax: [%.1f spatial, %.1f temporal, %.1f metric]\n', dmax(1), dmax(2), dmax(3));

    %% Method 1: neighbours.m (reference)
    fprintf('\nMethod 1: neighbours.m (reference)...\n');

    % Prepare data in vector format for neighbours.m
    [I, J, K] = ndgrid(1:nx, 1:ny, 1:nt);
    lon_vec = grid_data.Lon(:);
    lat_vec = grid_data.Lat(:);
    time_vec = repmat(grid_data.time, nx*ny, 1);
    z_vec = grid_data.Z(:);

    % Remove NaN values
    valid_idx = ~isnan(z_vec);
    c_vector = [lon_vec(valid_idx), lat_vec(valid_idx), time_vec(valid_idx)];
    z_vector = z_vec(valid_idx);

    tic;
    [csub_ref, zsub_ref, dsub_ref, nsub_ref, ~] = ...
        neighbours(p0, c_vector, z_vector, nmax, dmax);
    time_ref = toc;

    fprintf('  Time: %.4f sec\n', time_ref);
    fprintf('  Found: %d neighbors\n', nsub_ref);
    if nsub_ref > 0
        fprintf('  Distance range: [%.4f, %.4f] spatial, [%.4f, %.4f] temporal\n', ...
            min(dsub_ref(:,1)), max(dsub_ref(:,1)), ...
            min(dsub_ref(:,2)), max(dsub_ref(:,2)));
    end

    %% Method 2: neighbours_stug_optimized.m
    fprintf('\nMethod 2: neighbours_stug_optimized.m (adaptive)...\n');

    tic;
    [psub_opt, zsub_opt, dsub_opt, nsub_opt, ~] = ...
        neighbours_stug_optimized(p0, grid_data, nmax, dmax);
    time_opt = toc;

    fprintf('  Time: %.4f sec\n', time_opt);
    fprintf('  Found: %d neighbors\n', nsub_opt);
    if nsub_opt > 0
        csub_opt = psub_opt(:, 1:3);
        fprintf('  Distance range: [%.4f, %.4f] spatial, [%.4f, %.4f] temporal\n', ...
            min(dsub_opt(:,1)), max(dsub_opt(:,1)), ...
            min(dsub_opt(:,2)), max(dsub_opt(:,2)));
    end

    %% Method 3: neighbours_stug_index.m
    fprintf('\nMethod 3: neighbours_stug_index.m (index-based)...\n');

    tic;
    [psub_idx, zsub_idx, dsub_idx, nsub_idx, ~] = ...
        neighbours_stug_index(p0, grid_data, nmax, dmax);
    time_idx = toc;

    fprintf('  Time: %.4f sec\n', time_idx);
    fprintf('  Found: %d neighbors\n', nsub_idx);
    if nsub_idx > 0
        csub_idx = psub_idx(:, 1:3);
        fprintf('  Distance range: [%.4f, %.4f] spatial, [%.4f, %.4f] temporal\n', ...
            min(dsub_idx(:,1)), max(dsub_idx(:,1)), ...
            min(dsub_idx(:,2)), max(dsub_idx(:,2)));
    end

    %% Compare results
    fprintf('\n--- Comparison ---\n');

    % Compare optimized vs reference
    if nsub_ref > 0 && nsub_opt > 0
        [match_opt] = compare_neighbor_sets(csub_ref, csub_opt, 1e-6);
        fprintf('Optimized vs Reference:\n');
        fprintf('  Neighbor match: %.1f%%\n', match_opt.percent_match);
        fprintf('  Distance comparison: mean diff = %.6f\n', ...
            abs(mean(dsub_ref(:,1)) - mean(dsub_opt(:,1))));
    else
        match_opt.percent_match = 0;
        fprintf('Optimized vs Reference: Cannot compare (one returned 0 neighbors)\n');
    end

    % Compare index-based vs reference
    if nsub_ref > 0 && nsub_idx > 0
        [match_idx] = compare_neighbor_sets(csub_ref, csub_idx, 1e-6);
        fprintf('Index-based vs Reference:\n');
        fprintf('  Neighbor match: %.1f%%\n', match_idx.percent_match);
        fprintf('  Distance comparison: mean diff = %.6f\n', ...
            abs(mean(dsub_ref(:,1)) - mean(dsub_idx(:,1))));
    else
        match_idx.percent_match = 0;
        fprintf('Index-based vs Reference: Cannot compare (one returned 0 neighbors)\n');
    end

    % Compare index-based vs optimized
    if nsub_opt > 0 && nsub_idx > 0
        [match_opt_idx] = compare_neighbor_sets(csub_opt, csub_idx, 1e-6);
        fprintf('Index-based vs Optimized:\n');
        fprintf('  Neighbor match: %.1f%%\n', match_opt_idx.percent_match);
    else
        match_opt_idx.percent_match = 0;
        fprintf('Index-based vs Optimized: Cannot compare (one returned 0 neighbors)\n');
    end

    %% Performance comparison
    fprintf('\n--- Performance ---\n');
    fprintf('Reference:    %.4f sec (baseline)\n', time_ref);
    fprintf('Optimized:    %.4f sec (%.2fx speedup)\n', time_opt, time_ref/time_opt);
    fprintf('Index-based:  %.4f sec (%.2fx speedup)\n', time_idx, time_ref/time_idx);

    if time_idx < time_opt
        fprintf('Winner: Index-based is %.2fx faster than optimized\n', time_opt/time_idx);
    else
        fprintf('Winner: Optimized is %.2fx faster than index-based\n', time_idx/time_opt);
    end

    %% Store results
    results(s).config.nx = nx;
    results(s).config.ny = ny;
    results(s).config.nt = nt;
    results(s).config.nan_ratio = nan_ratio;
    results(s).config.nmax = nmax;
    results(s).config.name = description;

    results(s).nsub_ref = nsub_ref;
    results(s).nsub_opt = nsub_opt;
    results(s).nsub_idx = nsub_idx;

    results(s).time_ref = time_ref;
    results(s).time_opt = time_opt;
    results(s).time_idx = time_idx;

    results(s).speedup_opt = time_ref / time_opt;
    results(s).speedup_idx = time_ref / time_idx;
    results(s).relative_speedup = time_opt / time_idx;

    if nsub_ref > 0 && nsub_opt > 0
        results(s).match_opt_vs_ref = match_opt.percent_match;
    else
        results(s).match_opt_vs_ref = NaN;
    end

    if nsub_ref > 0 && nsub_idx > 0
        results(s).match_idx_vs_ref = match_idx.percent_match;
    else
        results(s).match_idx_vs_ref = NaN;
    end

    if nsub_opt > 0 && nsub_idx > 0
        results(s).match_idx_vs_opt = match_opt_idx.percent_match;
    else
        results(s).match_idx_vs_opt = NaN;
    end
end

%% Summary
fprintf('\n\n=== OVERALL SUMMARY ===\n\n');

% Create summary table
fprintf('%-30s | %8s | %8s | %8s | %8s | %8s | %6s | %6s\n', ...
    'Scenario', 'Ref(s)', 'Opt(s)', 'Idx(s)', 'Opt_x', 'Idx_x', 'M_opt', 'M_idx');
fprintf('%s\n', repmat('-', 120, 1));

for s = 1:n_scenarios
    fprintf('%-30s | %8.4f | %8.4f | %8.4f | %7.2fx | %7.2fx | %5.1f%% | %5.1f%%\n', ...
        results(s).config.name, ...
        results(s).time_ref, ...
        results(s).time_opt, ...
        results(s).time_idx, ...
        results(s).speedup_opt, ...
        results(s).speedup_idx, ...
        results(s).match_opt_vs_ref, ...
        results(s).match_idx_vs_ref);
end

fprintf('\n');

% Statistics
speedups_opt = [results.speedup_opt];
speedups_idx = [results.speedup_idx];
relative_speedups = [results.relative_speedup];

fprintf('Optimized Method Statistics:\n');
fprintf('  Average speedup: %.2fx\n', mean(speedups_opt));
fprintf('  Median speedup:  %.2fx\n', median(speedups_opt));
fprintf('  Range: %.2fx - %.2fx\n\n', min(speedups_opt), max(speedups_opt));

fprintf('Index-based Method Statistics:\n');
fprintf('  Average speedup: %.2fx\n', mean(speedups_idx));
fprintf('  Median speedup:  %.2fx\n', median(speedups_idx));
fprintf('  Range: %.2fx - %.2fx\n\n', min(speedups_idx), max(speedups_idx));

fprintf('Index-based vs Optimized:\n');
fprintf('  Average: %.2fx', mean(relative_speedups));
if mean(relative_speedups) < 1
    fprintf(' (index-based is faster)\n');
elseif mean(relative_speedups) > 1
    fprintf(' (optimized is faster)\n');
else
    fprintf(' (comparable)\n');
end
fprintf('  Scenarios where index-based wins: %d/%d\n', sum(relative_speedups > 1), n_scenarios);

% Accuracy
match_opt_avg = nanmean([results.match_opt_vs_ref]);
match_idx_avg = nanmean([results.match_idx_vs_ref]);

fprintf('\nAccuracy (vs Reference):\n');
fprintf('  Optimized average match: %.1f%%\n', match_opt_avg);
fprintf('  Index-based average match: %.1f%%\n', match_idx_avg);

%% Verdict
fprintf('\n--- Recommendations ---\n');

if mean(speedups_idx) > mean(speedups_opt) * 1.1
    fprintf('✓ Index-based method is significantly faster (%.1f%% improvement)\n', ...
        100 * (mean(speedups_idx) / mean(speedups_opt) - 1));
    fprintf('  Recommended for most use cases\n');
elseif mean(speedups_opt) > mean(speedups_idx) * 1.1
    fprintf('✓ Optimized method is significantly faster (%.1f%% improvement)\n', ...
        100 * (mean(speedups_opt) / mean(speedups_idx) - 1));
    fprintf('  Recommended for most use cases\n');
else
    fprintf('✓ Both methods offer comparable performance\n');
    fprintf('  Choose based on specific use case:\n');
    fprintf('  - Index-based: Better for anisotropic grids\n');
    fprintf('  - Optimized: Better for isotropic grids with adaptive strategy\n');
end

if match_idx_avg >= 95 && match_opt_avg >= 95
    fprintf('✓ Both methods are accurate (>95%% match with reference)\n');
elseif match_idx_avg < 90 || match_opt_avg < 90
    fprintf('⚠ Warning: Some accuracy concerns detected\n');
    fprintf('  Review test cases with low match percentages\n');
end

% Save results
save('compare_stug_results.mat', 'results', 'scenarios');
fprintf('\nResults saved to: compare_stug_results.mat\n');

fprintf('\n=== Comparison Complete ===\n');

%% Helper function
function [match_info] = compare_neighbor_sets(csub1, csub2, tolerance)
    % Compare two sets of neighbor coordinates
    % Returns percentage match and detailed comparison

    n1 = size(csub1, 1);
    n2 = size(csub2, 1);

    % Find matches
    n_matches = 0;
    for i = 1:n1
        dists = sqrt(sum((csub2 - repmat(csub1(i,:), n2, 1)).^2, 2));
        if any(dists < tolerance)
            n_matches = n_matches + 1;
        end
    end

    match_info.n1 = n1;
    match_info.n2 = n2;
    match_info.n_matches = n_matches;
    match_info.percent_match = 100 * n_matches / max(n1, n2);
end
