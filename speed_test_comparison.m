% speed_test_comparison.m
% Focused speed comparison between neighbours_stug_optimized.m and neighbours.m
%
% This script runs multiple iterations to get reliable timing measurements
% and shows the performance benefits of the optimized implementation.

clear; close all;

fprintf('=== Speed Test: neighbours_stug_optimized vs neighbours.m ===\n\n');

%% Test scenarios
% Realistic scenarios for global monthly soft data:
% - Global grids typically: 0.5-1 degree resolution = 360x180 to 720x360
% - Monthly data: 12-24 months typical
% - High spatial, low temporal dimensions

scenarios = {
    % [nx, ny, nt, nan_ratio, nmax, description]
    {50,  50,  12, 0.1,  15, 'Small grid - 1 year, low NaN'}
    {100, 100, 24, 0.3,  20, 'Medium grid - 2 years, 30%% NaN'}
    {150, 150, 24, 0.5,  20, 'Large grid - 2 years, 50%% NaN'}
    {200, 200, 24, 0.3,  25, 'Very large - 2 years, 30%% NaN'}
    {200, 200, 18, 0.7,  25, 'Very large - 18 months, sparse (70%% NaN)'}
    {300, 300, 24, 0.5,  30, 'Global-scale grid (300x300x24)'}
};

n_scenarios = length(scenarios);
n_iterations = 5;  % Number of runs per scenario for averaging

% Storage for results
timing_results = struct();

%% Run benchmarks
fprintf('Running %d iterations per scenario for reliable timing...\n\n', n_iterations);

for s_idx = 1:n_scenarios
    scenario = scenarios{s_idx};
    nx = scenario{1};
    ny = scenario{2};
    nt = scenario{3};
    nan_ratio = scenario{4};
    nmax = scenario{5};
    desc = scenario{6};

    fprintf('Scenario %d/%d: %s\n', s_idx, n_scenarios, desc);
    fprintf('  Grid: %dx%dx%d = %d points\n', nx, ny, nt, nx*ny*nt);
    fprintf('  NaN ratio: %.0f%%, nmax: %d\n', nan_ratio*100, nmax);

    %% Create test grid
    grid_data = struct();
    grid_data.x = linspace(-100, -80, nx)';
    grid_data.y = linspace(25, 35, ny)';
    grid_data.time = linspace(2400, 2500, nt)';

    [X_grid, Y_grid] = meshgrid(grid_data.x, grid_data.y);
    grid_data.Lon = X_grid';
    grid_data.Lat = Y_grid';

    % Create data with pattern
    rng(42);
    grid_data.Z = zeros(nx, ny, nt);
    for it = 1:nt
        for ix = 1:nx
            for iy = 1:ny
                grid_data.Z(ix, iy, it) = sin(ix/10) * cos(iy/10) + 0.1*randn();
            end
        end
    end

    % Add NaNs
    if nan_ratio > 0
        nan_mask = rand(nx, ny, nt) < nan_ratio;
        grid_data.Z(nan_mask) = NaN;
    end

    %% Convert to vector format for neighbours.m
    c_vector = [];
    Z_vector = [];
    for it = 1:nt
        for iy = 1:ny
            for ix = 1:nx
                if ~isnan(grid_data.Z(ix, iy, it))
                    c_vector = [c_vector; ...
                        grid_data.Lon(ix, iy), grid_data.Lat(ix, iy), grid_data.time(it)];
                    Z_vector = [Z_vector; grid_data.Z(ix, iy, it)];
                end
            end
        end
    end

    fprintf('  Non-NaN points: %d\n', length(Z_vector));

    %% Test point (center)
    p0 = [(grid_data.x(1)+grid_data.x(end))/2, ...
          (grid_data.y(1)+grid_data.y(end))/2, ...
          (grid_data.time(1)+grid_data.time(end))/2];

    %% Parameters
    dmax = [5.0, 15.0, 0.1];

    %% Benchmark neighbours.m
    fprintf('  Benchmarking neighbours.m...\n');
    times_ref = zeros(n_iterations, 1);
    nsub_ref = 0;

    for iter = 1:n_iterations
        tic;
        [~, ~, ~, nsub_ref, ~] = neighbours(p0, c_vector, Z_vector, nmax, dmax);
        times_ref(iter) = toc;
    end

    avg_time_ref = mean(times_ref);
    std_time_ref = std(times_ref);
    fprintf('    Average: %.4f ± %.4f sec (%d neighbors)\n', ...
        avg_time_ref, std_time_ref, nsub_ref);

    %% Benchmark neighbours_stug_optimized.m
    fprintf('  Benchmarking neighbours_stug_optimized.m...\n');
    times_opt = zeros(n_iterations, 1);
    nsub_opt = 0;

    for iter = 1:n_iterations
        tic;
        [~, ~, ~, nsub_opt, ~] = neighbours_stug_optimized(p0, grid_data, nmax, dmax);
        times_opt(iter) = toc;
    end

    avg_time_opt = mean(times_opt);
    std_time_opt = std(times_opt);
    fprintf('    Average: %.4f ± %.4f sec (%d neighbors)\n', ...
        avg_time_opt, std_time_opt, nsub_opt);

    %% Calculate speedup
    speedup = avg_time_ref / avg_time_opt;
    fprintf('  Speedup: %.2fx', speedup);
    if speedup > 10
        fprintf(' 🚀 (excellent!)');
    elseif speedup > 5
        fprintf(' ⚡ (very good)');
    elseif speedup > 2
        fprintf(' ✓ (good)');
    elseif speedup > 1
        fprintf(' + (faster)');
    else
        fprintf(' ⚠ (slower)');
    end
    fprintf('\n');

    % Time saved
    time_saved = avg_time_ref - avg_time_opt;
    fprintf('  Time saved: %.4f sec (%.1f%% reduction)\n', ...
        time_saved, (1 - avg_time_opt/avg_time_ref)*100);

    %% Store results
    timing_results(s_idx).description = desc;
    timing_results(s_idx).grid_size = [nx, ny, nt];
    timing_results(s_idx).total_points = nx*ny*nt;
    timing_results(s_idx).valid_points = length(Z_vector);
    timing_results(s_idx).nan_ratio = nan_ratio;
    timing_results(s_idx).nmax = nmax;
    timing_results(s_idx).time_ref_avg = avg_time_ref;
    timing_results(s_idx).time_ref_std = std_time_ref;
    timing_results(s_idx).time_opt_avg = avg_time_opt;
    timing_results(s_idx).time_opt_std = std_time_opt;
    timing_results(s_idx).speedup = speedup;
    timing_results(s_idx).nsub_ref = nsub_ref;
    timing_results(s_idx).nsub_opt = nsub_opt;

    fprintf('\n');
end

%% Summary Table
fprintf('=== PERFORMANCE SUMMARY ===\n\n');

fprintf('%-30s | %10s | %10s | %10s | %8s | %8s\n', ...
    'Scenario', 'Grid Pts', 'Valid Pts', 'Ref Time', 'Opt Time', 'Speedup');
fprintf('%s\n', repmat('-', 100, 1));

for i = 1:n_scenarios
    r = timing_results(i);
    fprintf('%-30s | %10d | %10d | %8.4fs | %8.4fs | %7.2fx\n', ...
        r.description, r.total_points, r.valid_points, ...
        r.time_ref_avg, r.time_opt_avg, r.speedup);
end

fprintf('\n');

%% Statistical Summary
all_speedups = [timing_results.speedup];
fprintf('Speedup Statistics:\n');
fprintf('  Mean:     %.2fx\n', mean(all_speedups));
fprintf('  Median:   %.2fx\n', median(all_speedups));
fprintf('  Min:      %.2fx\n', min(all_speedups));
fprintf('  Max:      %.2fx\n', max(all_speedups));
fprintf('  Std Dev:  %.2f\n', std(all_speedups));

fprintf('\n');

%% Scaling Analysis
fprintf('Scaling Analysis (how performance changes with grid size):\n');

grid_sizes = [timing_results.total_points];
times_ref = [timing_results.time_ref_avg];
times_opt = [timing_results.time_opt_avg];

[~, sort_idx] = sort(grid_sizes);

fprintf('  Grid Size  ->  neighbours.m  |  optimized  |  speedup\n');
fprintf('  %s\n', repmat('-', 60, 1));
for i = sort_idx
    fprintf('  %9d  ->  %8.4f sec  |  %8.4f sec  |  %6.2fx\n', ...
        grid_sizes(i), times_ref(i), times_opt(i), ...
        times_ref(i)/times_opt(i));
end

fprintf('\n');

%% Efficiency by NaN ratio
fprintf('Performance vs NaN Ratio:\n');

nan_ratios_unique = unique([timing_results.nan_ratio]);
for nr = nan_ratios_unique
    idx = find([timing_results.nan_ratio] == nr);
    avg_speedup_for_ratio = mean([timing_results(idx).speedup]);
    fprintf('  NaN ratio %.0f%%: average speedup %.2fx\n', ...
        nr*100, avg_speedup_for_ratio);
end

fprintf('\n');

%% Memory Efficiency Estimate
fprintf('Memory Efficiency:\n');
fprintf('  neighbours.m creates vector of all valid points: O(n)\n');
fprintf('  neighbours_stug_optimized.m works with candidates only: O(nmax)\n');

largest_scenario = timing_results(n_scenarios);
memory_ref_est = largest_scenario.valid_points * 4 * 8;  % 4 arrays × 8 bytes
memory_opt_est = largest_scenario.nmax * 3 * 8;  % 3 small arrays × 8 bytes

fprintf('  For largest test (%s):\n', largest_scenario.description);
fprintf('    neighbours.m estimated memory:       ~%.2f MB\n', memory_ref_est/1e6);
fprintf('    neighbours_stug_opt estimated memory: ~%.4f MB\n', memory_opt_est/1e6);
fprintf('    Memory reduction:                     ~%.1fx\n', memory_ref_est/memory_opt_est);

fprintf('\n');

%% Recommendations
fprintf('=== RECOMMENDATIONS ===\n\n');

fprintf('Use neighbours_stug_optimized.m when:\n');
fprintf('  ✓ Data is on a uniform grid (constant spacing)\n');
fprintf('  ✓ Grid size > 100×100×100 points\n');
fprintf('  ✓ NaN ratio > 30%% (sparse data)\n');
fprintf('  ✓ Need to query many points (speedup compounds)\n');
fprintf('  ✓ Memory constraints are important\n');

fprintf('\n');

fprintf('Use neighbours.m when:\n');
fprintf('  ✓ Data at arbitrary irregular locations\n');
fprintf('  ✓ Small datasets (<10,000 points)\n');
fprintf('  ✓ Need multi-variable indexing\n');
fprintf('  ✓ Grid is not uniform\n');

fprintf('\n=== Speed Test Complete ===\n');

%% Save results
save('speed_test_results.mat', 'timing_results', 'scenarios');
fprintf('\nResults saved to: speed_test_results.mat\n');
