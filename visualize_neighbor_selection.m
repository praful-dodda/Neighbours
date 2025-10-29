% visualize_neighbor_selection.m
% Create comprehensive visualizations showing:
% - Grid data distribution
% - Selected neighbors by both methods
% - Overlap and differences
% - Spatial and temporal views
%
% Prerequisites: Run compare_neighbours_implementations.m first

clear; close all;

fprintf('=== Neighbor Selection Visualization ===\n\n');

% Load comparison results
if ~exist('compare_results.mat', 'file')
    error('Please run compare_neighbours_implementations.m first to generate compare_results.mat');
end

load('compare_results.mat', 'results');
n_tests = length(results);

fprintf('Loaded comparison results for %d test scenarios\n', n_tests);

% Select a representative test case to visualize
% Pick medium grid with moderate NaN ratio
test_idx = 2;  % Adjust as needed
if test_idx > n_tests
    test_idx = min(2, n_tests);
end

cfg = results(test_idx).config;
fprintf('Visualizing test case: %s\n', cfg.name);
fprintf('  Grid: %dx%dx%d, NaN ratio: %.1f%%, nmax: %d\n\n', ...
    cfg.nx, cfg.ny, cfg.nt, cfg.nan_ratio*100, cfg.nmax);

% Recreate the grid for this test
grid_data = create_test_grid(cfg.nx, cfg.ny, cfg.nt, cfg.nan_ratio);

% Extract neighbor data
csub_ref = results(test_idx).csub_ref;
csub_opt = results(test_idx).csub_opt;
p0 = results(test_idx).p0;

fprintf('Creating visualizations...\n');

%% Figure 1: Spatial view (lon-lat) at target time
figure('Position', [100, 100, 1400, 500]);

% Find time slice closest to p0
[~, t_idx] = min(abs(grid_data.time - p0(3)));

subplot(1, 3, 1);
% Plot all valid data points at this time
Z_slice = grid_data.Z(:, :, t_idx);
valid_mask = ~isnan(Z_slice);
[valid_x, valid_y] = find(valid_mask);
scatter(grid_data.x(valid_x), grid_data.y(valid_y), 20, Z_slice(valid_mask), 'filled', 'MarkerFaceAlpha', 0.4);
hold on;
% Target point
plot(p0(1), p0(2), 'r*', 'MarkerSize', 20, 'LineWidth', 3);
% Reference neighbors at this time (if any)
ref_at_time = csub_ref(abs(csub_ref(:,3) - p0(3)) < abs(grid_data.time(2) - grid_data.time(1)), :);
if ~isempty(ref_at_time)
    plot(ref_at_time(:,1), ref_at_time(:,2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
end
colorbar;
xlabel('Longitude');
ylabel('Latitude');
title(sprintf('neighbours.m selection\n(green circles = selected at t≈%.0f)', p0(3)));
grid on;
axis tight;
hold off;

subplot(1, 3, 2);
% Same for optimized
scatter(grid_data.x(valid_x), grid_data.y(valid_y), 20, Z_slice(valid_mask), 'filled', 'MarkerFaceAlpha', 0.4);
hold on;
plot(p0(1), p0(2), 'r*', 'MarkerSize', 20, 'LineWidth', 3);
opt_at_time = csub_opt(abs(csub_opt(:,3) - p0(3)) < abs(grid_data.time(2) - grid_data.time(1)), :);
if ~isempty(opt_at_time)
    plot(opt_at_time(:,1), opt_at_time(:,2), 'mo', 'MarkerSize', 10, 'LineWidth', 2);
end
colorbar;
xlabel('Longitude');
ylabel('Latitude');
title(sprintf('neighbours\\_stug\\_optimized selection\n(magenta circles = selected at t≈%.0f)', p0(3)));
grid on;
axis tight;
hold off;

subplot(1, 3, 3);
% Overlap view
scatter(grid_data.x(valid_x), grid_data.y(valid_y), 20, Z_slice(valid_mask), 'filled', 'MarkerFaceAlpha', 0.2);
hold on;
plot(p0(1), p0(2), 'r*', 'MarkerSize', 20, 'LineWidth', 3);
% Find common neighbors
tolerance = 1e-6;
common_neighbors = [];
ref_only = [];
opt_only = [];

for i = 1:size(ref_at_time, 1)
    dists = sqrt(sum((opt_at_time - repmat(ref_at_time(i,:), size(opt_at_time,1), 1)).^2, 2));
    if any(dists < tolerance)
        common_neighbors = [common_neighbors; ref_at_time(i,:)];
    else
        ref_only = [ref_only; ref_at_time(i,:)];
    end
end

for i = 1:size(opt_at_time, 1)
    dists = sqrt(sum((ref_at_time - repmat(opt_at_time(i,:), size(ref_at_time,1), 1)).^2, 2));
    if ~any(dists < tolerance)
        opt_only = [opt_only; opt_at_time(i,:)];
    end
end

% Plot different categories
if ~isempty(common_neighbors)
    plot(common_neighbors(:,1), common_neighbors(:,2), 'go', 'MarkerSize', 12, 'LineWidth', 2.5);
end
if ~isempty(ref_only)
    plot(ref_only(:,1), ref_only(:,2), 'bs', 'MarkerSize', 10, 'LineWidth', 2);
end
if ~isempty(opt_only)
    plot(opt_only(:,1), opt_only(:,2), 'r^', 'MarkerSize', 10, 'LineWidth', 2);
end
legend({'Data', 'Target', 'Both methods', 'Reference only', 'Optimized only'}, 'Location', 'best');
colorbar;
xlabel('Longitude');
ylabel('Latitude');
title(sprintf('Overlap analysis (t≈%.0f)', p0(3)));
grid on;
axis tight;
hold off;

sgtitle(sprintf('Spatial View: %s', cfg.name), 'FontSize', 14, 'FontWeight', 'bold');

saveas(gcf, sprintf('neighbor_spatial_view_test%d.png', test_idx));
fprintf('  Saved: neighbor_spatial_view_test%d.png\n', test_idx);

%% Figure 2: 3D space-time view
figure('Position', [150, 150, 1200, 900]);

% All data points (sample for clarity)
sample_ratio = min(1, 1000 / sum(~isnan(grid_data.Z(:))));
subplot(2, 2, 1);
[ix, iy, it] = ind2sub(size(grid_data.Z), find(~isnan(grid_data.Z)));
n_sample = floor(length(ix) * sample_ratio);
sample_idx = randperm(length(ix), n_sample);
scatter3(grid_data.x(ix(sample_idx)), grid_data.y(iy(sample_idx)), ...
    grid_data.time(it(sample_idx)), 10, grid_data.Z(~isnan(grid_data.Z)), 'filled', 'MarkerFaceAlpha', 0.3);
hold on;
plot3(p0(1), p0(2), p0(3), 'r*', 'MarkerSize', 20, 'LineWidth', 3);
colorbar;
xlabel('Lon');
ylabel('Lat');
zlabel('Time');
title(sprintf('All data (%.0f%% sample)', sample_ratio*100));
grid on;
view(45, 30);
hold off;

% Reference method neighbors
subplot(2, 2, 2);
scatter3(grid_data.x(ix(sample_idx)), grid_data.y(iy(sample_idx)), ...
    grid_data.time(it(sample_idx)), 5, [0.8 0.8 0.8], 'filled', 'MarkerFaceAlpha', 0.1);
hold on;
plot3(p0(1), p0(2), p0(3), 'r*', 'MarkerSize', 20, 'LineWidth', 3);
plot3(csub_ref(:,1), csub_ref(:,2), csub_ref(:,3), 'go', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'g');
% Draw lines to target
for i = 1:size(csub_ref, 1)
    plot3([p0(1), csub_ref(i,1)], [p0(2), csub_ref(i,2)], [p0(3), csub_ref(i,3)], ...
        'g-', 'LineWidth', 0.5, 'Color', [0, 0.8, 0, 0.3]);
end
xlabel('Lon');
ylabel('Lat');
zlabel('Time');
title(sprintf('neighbours.m (%d neighbors)', size(csub_ref,1)));
grid on;
view(45, 30);
legend({'Data', 'Target', 'Selected'}, 'Location', 'best');
hold off;

% Optimized method neighbors
subplot(2, 2, 3);
scatter3(grid_data.x(ix(sample_idx)), grid_data.y(iy(sample_idx)), ...
    grid_data.time(it(sample_idx)), 5, [0.8 0.8 0.8], 'filled', 'MarkerFaceAlpha', 0.1);
hold on;
plot3(p0(1), p0(2), p0(3), 'r*', 'MarkerSize', 20, 'LineWidth', 3);
plot3(csub_opt(:,1), csub_opt(:,2), csub_opt(:,3), 'mo', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'm');
% Draw lines to target
for i = 1:size(csub_opt, 1)
    plot3([p0(1), csub_opt(i,1)], [p0(2), csub_opt(i,2)], [p0(3), csub_opt(i,3)], ...
        'm-', 'LineWidth', 0.5, 'Color', [1, 0, 1, 0.3]);
end
xlabel('Lon');
ylabel('Lat');
zlabel('Time');
title(sprintf('neighbours\\_stug\\_optimized (%d neighbors)', size(csub_opt,1)));
grid on;
view(45, 30);
legend({'Data', 'Target', 'Selected'}, 'Location', 'best');
hold off;

% Overlap view
subplot(2, 2, 4);
scatter3(grid_data.x(ix(sample_idx)), grid_data.y(iy(sample_idx)), ...
    grid_data.time(it(sample_idx)), 5, [0.8 0.8 0.8], 'filled', 'MarkerFaceAlpha', 0.1);
hold on;
plot3(p0(1), p0(2), p0(3), 'r*', 'MarkerSize', 20, 'LineWidth', 3);

% Find common neighbors in 3D
common_3d = [];
ref_only_3d = [];
opt_only_3d = [];

for i = 1:size(csub_ref, 1)
    dists = sqrt(sum((csub_opt - repmat(csub_ref(i,:), size(csub_opt,1), 1)).^2, 2));
    if any(dists < tolerance)
        common_3d = [common_3d; csub_ref(i,:)];
    else
        ref_only_3d = [ref_only_3d; csub_ref(i,:)];
    end
end

for i = 1:size(csub_opt, 1)
    dists = sqrt(sum((csub_ref - repmat(csub_opt(i,:), size(csub_ref,1), 1)).^2, 2));
    if ~any(dists < tolerance)
        opt_only_3d = [opt_only_3d; csub_opt(i,:)];
    end
end

if ~isempty(common_3d)
    plot3(common_3d(:,1), common_3d(:,2), common_3d(:,3), ...
        'go', 'MarkerSize', 12, 'LineWidth', 2.5, 'MarkerFaceColor', 'g');
end
if ~isempty(ref_only_3d)
    plot3(ref_only_3d(:,1), ref_only_3d(:,2), ref_only_3d(:,3), ...
        'bs', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'b');
end
if ~isempty(opt_only_3d)
    plot3(opt_only_3d(:,1), opt_only_3d(:,2), opt_only_3d(:,3), ...
        'r^', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'r');
end
xlabel('Lon');
ylabel('Lat');
zlabel('Time');
title(sprintf('Overlap: %d common, %d ref-only, %d opt-only', ...
    size(common_3d,1), size(ref_only_3d,1), size(opt_only_3d,1)));
grid on;
view(45, 30);
legend({'Data', 'Target', 'Both', 'Ref only', 'Opt only'}, 'Location', 'best');
hold off;

sgtitle(sprintf('3D Space-Time View: %s', cfg.name), 'FontSize', 14, 'FontWeight', 'bold');

saveas(gcf, sprintf('neighbor_3d_view_test%d.png', test_idx));
fprintf('  Saved: neighbor_3d_view_test%d.png\n', test_idx);

%% Figure 3: Distance distributions
figure('Position', [200, 200, 1200, 800]);

if isfield(results(test_idx), 'dist_comparison')
    dc = results(test_idx).dist_comparison;

    % Spatial distances
    subplot(2, 3, 1);
    dsub_ref = results(test_idx).csub_ref;  % Need to recompute
    dsub_opt = results(test_idx).csub_opt;
    % Compute distances
    spatial_dist_ref = sqrt((csub_ref(:,1) - p0(1)).^2 + (csub_ref(:,2) - p0(2)).^2);
    spatial_dist_opt = sqrt((csub_opt(:,1) - p0(1)).^2 + (csub_opt(:,2) - p0(2)).^2);

    histogram(spatial_dist_ref, 15, 'FaceColor', 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    hold on;
    histogram(spatial_dist_opt, 15, 'FaceColor', 'm', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    xlabel('Spatial Distance');
    ylabel('Count');
    title('Spatial Distance Distribution');
    legend({'neighbours.m', 'optimized'});
    grid on;
    hold off;

    % Temporal distances
    subplot(2, 3, 2);
    temporal_dist_ref = abs(csub_ref(:,3) - p0(3));
    temporal_dist_opt = abs(csub_opt(:,3) - p0(3));

    histogram(temporal_dist_ref, 15, 'FaceColor', 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    hold on;
    histogram(temporal_dist_opt, 15, 'FaceColor', 'm', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    xlabel('Temporal Distance');
    ylabel('Count');
    title('Temporal Distance Distribution');
    legend({'neighbours.m', 'optimized'});
    grid on;
    hold off;

    % Combined space-time distance
    subplot(2, 3, 3);
    dmax = [5.0, 15.0, 0.1];  % Default from test
    combined_dist_ref = spatial_dist_ref + dmax(3) * temporal_dist_ref;
    combined_dist_opt = spatial_dist_opt + dmax(3) * temporal_dist_opt;

    histogram(combined_dist_ref, 15, 'FaceColor', 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    hold on;
    histogram(combined_dist_opt, 15, 'FaceColor', 'm', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    xlabel('Combined Space-Time Distance');
    ylabel('Count');
    title('Combined Distance Distribution');
    legend({'neighbours.m', 'optimized'});
    grid on;
    hold off;

    % CDF comparison
    subplot(2, 3, 4);
    [f_ref, x_ref] = ecdf(combined_dist_ref);
    [f_opt, x_opt] = ecdf(combined_dist_opt);
    plot(x_ref, f_ref, 'g-', 'LineWidth', 2);
    hold on;
    plot(x_opt, f_opt, 'm--', 'LineWidth', 2);
    xlabel('Distance');
    ylabel('Cumulative Probability');
    title('Cumulative Distribution');
    legend({'neighbours.m', 'optimized'});
    grid on;
    hold off;

    % Q-Q plot
    subplot(2, 3, 5);
    qqplot(combined_dist_ref, combined_dist_opt);
    title('Q-Q Plot');
    xlabel('neighbours.m quantiles');
    ylabel('optimized quantiles');
    grid on;

    % Statistics comparison
    subplot(2, 3, 6);
    axis off;
    stats_text = {
        '\bf Distribution Statistics'
        ''
        sprintf('\\bf Mean distances:')
        sprintf('  neighbours.m: %.4f', dc.mean_ref)
        sprintf('  optimized:    %.4f', dc.mean_opt)
        sprintf('  Difference:   %.4f (%.2f%%)', dc.mean_diff, dc.mean_percent_diff)
        ''
        sprintf('\\bf Median distances:')
        sprintf('  neighbours.m: %.4f', dc.median_ref)
        sprintf('  optimized:    %.4f', dc.median_opt)
        sprintf('  Difference:   %.4f (%.2f%%)', dc.median_diff, dc.median_percent_diff)
        ''
        sprintf('\\bf Std deviation:')
        sprintf('  neighbours.m: %.4f', dc.std_ref)
        sprintf('  optimized:    %.4f', dc.std_opt)
        sprintf('  Difference:   %.4f', dc.std_diff)
        ''
        sprintf('\\bf KS Test:')
        sprintf('  p-value: %.4f', dc.ks_pvalue)
        if dc.ks_pvalue > 0.05
            '  → Distributions match ✓'
        else
            '  → Distributions differ'
        end
    };
    text(0.1, 0.9, stats_text, 'VerticalAlignment', 'top', 'FontSize', 10, 'Interpreter', 'tex');
end

sgtitle(sprintf('Distance Distribution Analysis: %s', cfg.name), 'FontSize', 14, 'FontWeight', 'bold');

saveas(gcf, sprintf('neighbor_distributions_test%d.png', test_idx));
fprintf('  Saved: neighbor_distributions_test%d.png\n', test_idx);

fprintf('\n=== Visualization Complete ===\n');
fprintf('Generated 3 figures for test case %d\n', test_idx);
fprintf('  - Spatial view (2D slice at target time)\n');
fprintf('  - 3D space-time view\n');
fprintf('  - Distance distribution analysis\n');

%% Helper function
function grid_data = create_test_grid(nx, ny, nt, nan_ratio)
    % Recreate test grid (same as in compare script)
    grid_data.x = linspace(-100, -80, nx)';
    grid_data.y = linspace(25, 35, ny)';
    grid_data.time = linspace(2400, 2500, nt)';

    [X_grid, Y_grid] = meshgrid(grid_data.x, grid_data.y);
    grid_data.Lon = X_grid';
    grid_data.Lat = Y_grid';

    rng(42);
    grid_data.Z = zeros(nx, ny, nt);
    for it = 1:nt
        for ix = 1:nx
            for iy = 1:ny
                grid_data.Z(ix, iy, it) = ...
                    sin(ix/10) * cos(iy/10) * exp(-it/50) + 0.1*randn();
            end
        end
    end

    if nan_ratio > 0
        nan_mask = rand(nx, ny, nt) < nan_ratio;
        grid_data.Z(nan_mask) = NaN;
    end
end
