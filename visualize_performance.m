% visualize_performance.m
% Create visualizations comparing neighbours.m and neighbours_stug_optimized.m
%
% Run this after speed_test_comparison.m to visualize the results

clear; close all;

% Load results from speed test
if ~exist('speed_test_results_index.mat', 'file')
    error('Please run speed_test_comparison.m first to generate results');
end

load('speed_test_results.mat', 'timing_results', 'scenarios');

fprintf('Creating performance visualizations...\n');

%% Extract data
n_tests = length(timing_results);
grid_sizes = [timing_results.total_points];
valid_points = [timing_results.valid_points];
nan_ratios = [timing_results.nan_ratio];
times_ref = [timing_results.time_ref_avg];
times_opt = [timing_results.time_opt_avg];
speedups = [timing_results.speedup];
descriptions = {timing_results.description};

%% Figure 1: Execution Time Comparison
figure('Position', [100, 100, 1200, 400]);

subplot(1, 3, 1);
x = 1:n_tests;
bar(x, [times_ref', times_opt']);
set(gca, 'XTickLabel', descriptions, 'XTickLabelRotation', 45);
legend({'neighbours.m', 'neighbours\_stug\_opt'}, 'Location', 'northwest');
ylabel('Time (seconds)');
title('Execution Time Comparison');
grid on;

subplot(1, 3, 2);
bar(x, speedups);
set(gca, 'XTickLabel', descriptions, 'XTickLabelRotation', 45);
ylabel('Speedup Factor');
title('Speedup (higher is better)');
grid on;
hold on;
plot([0, n_tests+1], [1, 1], 'r--', 'LineWidth', 1.5);  % Reference line at 1x
hold off;

subplot(1, 3, 3);
semilogy(grid_sizes, times_ref, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(grid_sizes, times_opt, 's-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Total Grid Points');
ylabel('Time (seconds, log scale)');
title('Scalability with Grid Size');
legend({'neighbours.m', 'neighbours\_stug\_opt'}, 'Location', 'northwest');
grid on;
hold off;

sgtitle('Performance Comparison: neighbours.m vs neighbours\_stug\_optimized.m', ...
    'FontSize', 14, 'FontWeight', 'bold');

% Save figure
saveas(gcf, 'performance_comparison.png');
fprintf('  Saved: performance_comparison.png\n');

%% Figure 2: Speedup Analysis
figure('Position', [150, 150, 1000, 600]);

subplot(2, 2, 1);
scatter(grid_sizes, speedups, 100, nan_ratios*100, 'filled');
xlabel('Grid Size (total points)');
ylabel('Speedup Factor');
title('Speedup vs Grid Size');
colorbar;
caxis([0 100]);
ylabel(colorbar, 'NaN Ratio (%)');
grid on;

subplot(2, 2, 2);
scatter(nan_ratios*100, speedups, 100, grid_sizes, 'filled');
xlabel('NaN Ratio (%)');
ylabel('Speedup Factor');
title('Speedup vs NaN Ratio');
colorbar;
ylabel(colorbar, 'Grid Size');
grid on;

subplot(2, 2, 3);
time_saved = times_ref - times_opt;
bar(x, time_saved);
set(gca, 'XTickLabel', descriptions, 'XTickLabelRotation', 45);
ylabel('Time Saved (seconds)');
title('Absolute Time Savings');
grid on;

subplot(2, 2, 4);
percent_reduction = (1 - times_opt./times_ref) * 100;
bar(x, percent_reduction);
set(gca, 'XTickLabel', descriptions, 'XTickLabelRotation', 45);
ylabel('% Time Reduction');
title('Relative Time Savings');
grid on;

sgtitle('Speedup Analysis', 'FontSize', 14, 'FontWeight', 'bold');

% Save figure
saveas(gcf, 'speedup_analysis.png');
fprintf('  Saved: speedup_analysis.png\n');

%% Figure 3: Detailed Performance Metrics
figure('Position', [200, 200, 1200, 500]);

subplot(1, 2, 1);
% Stacked bar showing time breakdown
time_data = [times_opt', (times_ref - times_opt)'];
bar(x, time_data, 'stacked');
set(gca, 'XTickLabel', descriptions, 'XTickLabelRotation', 45);
ylabel('Time (seconds)');
legend({'Optimized Time', 'Time Saved'}, 'Location', 'northwest');
title('Time Breakdown: Optimized vs Saved');
grid on;

subplot(1, 2, 2);
% Efficiency: points processed per second
points_per_sec_ref = valid_points ./ times_ref;
points_per_sec_opt = valid_points ./ times_opt;
bar(x, [points_per_sec_ref', points_per_sec_opt']);
set(gca, 'XTickLabel', descriptions, 'XTickLabelRotation', 45);
ylabel('Valid Points Processed per Second');
legend({'neighbours.m', 'neighbours\_stug\_opt'}, 'Location', 'northwest');
title('Processing Throughput');
grid on;

sgtitle('Performance Metrics', 'FontSize', 14, 'FontWeight', 'bold');

% Save figure
saveas(gcf, 'performance_metrics.png');
fprintf('  Saved: performance_metrics.png\n');

%% Summary Statistics Display
figure('Position', [250, 250, 600, 400]);
axis off;

% Find the best scenario index before building the summary text
idx_best = find(speedups == max(speedups), 1);

summary_text = {
    '\bf Performance Summary Statistics'
    ''
    sprintf('\\bf Total Test Scenarios: %d', n_tests)
    ''
    '\bf Speedup:'
    sprintf('  Mean:     %.2fx', mean(speedups))
    sprintf('  Median:   %.2fx', median(speedups))
    sprintf('  Range:    %.2fx - %.2fx', min(speedups), max(speedups))
    ''
    '\bf Time Savings:'
    sprintf('  Average time saved: %.4f sec', mean(times_ref - times_opt))
    sprintf('  Average reduction:  %.1f%%', mean(percent_reduction))
    ''
    '\bf Best Performance:'
    sprintf('  Scenario: %s', descriptions{idx_best})
    sprintf('  Speedup:  %.2fx', speedups(idx_best))
    sprintf('  From %.4fs to %.4fs', times_ref(idx_best), times_opt(idx_best))
};

text(0.1, 0.9, summary_text, 'VerticalAlignment', 'top', ...
    'FontSize', 11, 'Interpreter', 'tex');

% Save figure
saveas(gcf, 'performance_summary.png');
fprintf('  Saved: performance_summary.png\n');

fprintf('\nVisualization complete! Generated 4 figures:\n');
fprintf('  1. performance_comparison.png - Overall comparison\n');
fprintf('  2. speedup_analysis.png - Detailed speedup metrics\n');
fprintf('  3. performance_metrics.png - Throughput analysis\n');
fprintf('  4. performance_summary.png - Summary statistics\n');
