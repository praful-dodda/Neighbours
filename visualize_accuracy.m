% visualize_accuracy.m
% Create visualizations comparing accuracy of neighbours_stug_optimized.m
% vs neighbours.m reference implementation
%
% This script visualizes:
% - Neighbor match quality
% - Distance accuracy
% - Value accuracy
% - Spatial distribution of selected neighbors
%
% Prerequisites: Run compare_neighbours_implementations.m first

clear; close all;

fprintf('=== Accuracy Visualization ===\n\n');

% Check if comparison results exist
if ~exist('compare_results.mat', 'file')
    fprintf('No saved comparison results found.\n');
    fprintf('Running compare_neighbours_implementations.m first...\n\n');

    % Run comparison and save results
    compare_neighbours_implementations;

    % The comparison script should save results
    % If not, we'll create a simple version
    if ~exist('compare_results.mat', 'file')
        error('Please run compare_neighbours_implementations.m and save results to compare_results.mat');
    end
end

% Load comparison results
load('compare_results.mat', 'results');

fprintf('Loaded comparison results for %d test scenarios\n', length(results));
fprintf('Creating accuracy visualizations...\n\n');

n_tests = length(results);

%% Figure 1: Match Quality Overview
figure('Position', [100, 100, 1200, 800]);

% Extract match statistics
match_percents = zeros(n_tests, 1);
max_dist_diffs = zeros(n_tests, 1);
max_value_diffs = zeros(n_tests, 1);
test_names = cell(n_tests, 1);

for i = 1:n_tests
    match_percents(i) = results(i).match_quality.percent_match;
    max_dist_diffs(i) = results(i).match_quality.max_distance_diff;
    max_value_diffs(i) = results(i).match_quality.max_value_diff;
    test_names{i} = results(i).config.name;
end

% Subplot 1: Match percentage
subplot(2, 3, 1);
bar(1:n_tests, match_percents);
set(gca, 'XTickLabel', test_names, 'XTickLabelRotation', 45);
ylabel('Match Percentage (%)');
title('Neighbor Set Agreement');
ylim([0, 105]);
hold on;
plot([0, n_tests+1], [100, 100], 'g--', 'LineWidth', 2);
plot([0, n_tests+1], [95, 95], 'r--', 'LineWidth', 1);
hold off;
grid on;
legend({'Match %', 'Perfect', '95% Threshold'}, 'Location', 'southeast');

% Subplot 2: Distance accuracy
subplot(2, 3, 2);
semilogy(1:n_tests, max_dist_diffs, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
set(gca, 'XTickLabel', test_names, 'XTickLabelRotation', 45);
ylabel('Max Distance Difference (log scale)');
title('Distance Accuracy');
grid on;

% Subplot 3: Value accuracy
subplot(2, 3, 3);
semilogy(1:n_tests, max_value_diffs, 's-', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.8, 0.4, 0]);
set(gca, 'XTickLabel', test_names, 'XTickLabelRotation', 45);
ylabel('Max Value Difference (log scale)');
title('Value Accuracy');
grid on;

% Subplot 4: Number of neighbors comparison
subplot(2, 3, 4);
nsub_ref = [results.nsub_ref];
nsub_opt = [results.nsub_opt];
bar(1:n_tests, [nsub_ref', nsub_opt']);
set(gca, 'XTickLabel', test_names, 'XTickLabelRotation', 45);
ylabel('Number of Neighbors');
title('Neighbors Found');
legend({'neighbours.m', 'optimized'}, 'Location', 'northwest');
grid on;

% Subplot 5: Correlation between number of neighbors
subplot(2, 3, 5);
plot(nsub_ref, nsub_opt, 'o', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
max_n = max([nsub_ref, nsub_opt]);
plot([0, max_n], [0, max_n], 'r--', 'LineWidth', 1.5);
hold off;
xlabel('neighbours.m (# neighbors)');
ylabel('optimized (# neighbors)');
title('Neighbor Count Correlation');
grid on;
axis equal tight;

% Subplot 6: Overall accuracy score
subplot(2, 3, 6);
% Compute overall accuracy score (0-100)
accuracy_scores = match_percents;
distance_ok = (max_dist_diffs < 1e-6);
value_ok = (max_value_diffs < 1e-9);
accuracy_scores(~distance_ok) = accuracy_scores(~distance_ok) * 0.95;
accuracy_scores(~value_ok) = accuracy_scores(~value_ok) * 0.95;

bar(1:n_tests, accuracy_scores);
set(gca, 'XTickLabel', test_names, 'XTickLabelRotation', 45);
ylabel('Overall Accuracy Score');
title('Combined Accuracy Metric');
ylim([0, 105]);
hold on;
plot([0, n_tests+1], [100, 100], 'g--', 'LineWidth', 2);
hold off;
grid on;

sgtitle('Accuracy Comparison: neighbours\_stug\_optimized vs neighbours.m', ...
    'FontSize', 14, 'FontWeight', 'bold');

% Save figure
saveas(gcf, 'accuracy_comparison.png');
fprintf('  Saved: accuracy_comparison.png\n');

%% Figure 2: Detailed Match Quality Analysis
figure('Position', [150, 150, 1000, 600]);

% Prepare data for heatmap-style visualization
subplot(2, 2, 1);
perfect_match = (match_percents >= 99.9);
good_match = (match_percents >= 95) & (match_percents < 99.9);
partial_match = (match_percents < 95);

match_categories = zeros(n_tests, 3);
match_categories(:, 1) = perfect_match;
match_categories(:, 2) = good_match;
match_categories(:, 3) = partial_match;

bar(1:n_tests, match_categories, 'stacked');
set(gca, 'XTickLabel', test_names, 'XTickLabelRotation', 45);
ylabel('Match Quality');
title('Match Quality Categories');
legend({'Perfect (≥99.9%)', 'Good (95-99.9%)', 'Partial (<95%)'}, ...
    'Location', 'northwest');
ylim([0, 1.2]);
grid on;

% Subplot 2: Distance error distribution
subplot(2, 2, 2);
categories = categorical({'Excellent\n(<10^{-10})', 'Very Good\n(<10^{-8})', ...
    'Good\n(<10^{-6})', 'Acceptable\n(<10^{-4})', 'Poor\n(≥10^{-4})'});
categories = reordercats(categories, {'Excellent\n(<10^{-10})', 'Very Good\n(<10^{-8})', ...
    'Good\n(<10^{-6})', 'Acceptable\n(<10^{-4})', 'Poor\n(≥10^{-4})'});

dist_excellent = sum(max_dist_diffs < 1e-10);
dist_very_good = sum(max_dist_diffs >= 1e-10 & max_dist_diffs < 1e-8);
dist_good = sum(max_dist_diffs >= 1e-8 & max_dist_diffs < 1e-6);
dist_acceptable = sum(max_dist_diffs >= 1e-6 & max_dist_diffs < 1e-4);
dist_poor = sum(max_dist_diffs >= 1e-4);

dist_counts = [dist_excellent, dist_very_good, dist_good, dist_acceptable, dist_poor];
bar(categories, dist_counts);
ylabel('Number of Test Cases');
title('Distance Accuracy Distribution');
grid on;

% Subplot 3: Value error distribution
subplot(2, 2, 3);
val_excellent = sum(max_value_diffs < 1e-12);
val_very_good = sum(max_value_diffs >= 1e-12 & max_value_diffs < 1e-10);
val_good = sum(max_value_diffs >= 1e-10 & max_value_diffs < 1e-9);
val_acceptable = sum(max_value_diffs >= 1e-9 & max_value_diffs < 1e-6);
val_poor = sum(max_value_diffs >= 1e-6);

val_counts = [val_excellent, val_very_good, val_good, val_acceptable, val_poor];
bar(categories, val_counts);
ylabel('Number of Test Cases');
title('Value Accuracy Distribution');
grid on;

% Subplot 4: Summary statistics text
subplot(2, 2, 4);
axis off;

summary_stats = {
    '\bf Accuracy Summary'
    ''
    sprintf('\\bf Total Test Cases: %d', n_tests)
    ''
    '\bf Match Quality:'
    sprintf('  Perfect matches:  %d (%.1f%%)', sum(perfect_match), 100*sum(perfect_match)/n_tests)
    sprintf('  Good matches:     %d (%.1f%%)', sum(good_match), 100*sum(good_match)/n_tests)
    sprintf('  Partial matches:  %d (%.1f%%)', sum(partial_match), 100*sum(partial_match)/n_tests)
    sprintf('  Average match:    %.2f%%', mean(match_percents))
    ''
    '\bf Distance Accuracy:'
    sprintf('  Excellent (<%s): %d', '10^{-10}', dist_excellent)
    sprintf('  Max error:         %.2e', max(max_dist_diffs))
    ''
    '\bf Value Accuracy:'
    sprintf('  Excellent (<%s): %d', '10^{-12}', val_excellent)
    sprintf('  Max error:         %.2e', max(max_value_diffs))
};

text(0.1, 0.9, summary_stats, 'VerticalAlignment', 'top', ...
    'FontSize', 10, 'Interpreter', 'tex');

sgtitle('Detailed Accuracy Analysis', 'FontSize', 14, 'FontWeight', 'bold');

% Save figure
saveas(gcf, 'accuracy_detailed.png');
fprintf('  Saved: accuracy_detailed.png\n');

%% Figure 3: Accuracy vs Performance Trade-off
figure('Position', [200, 200, 800, 600]);

speedups = [results.speedup];

% Main plot: Accuracy vs Speedup
scatter(speedups, match_percents, 100, 1:n_tests, 'filled');
xlabel('Speedup Factor');
ylabel('Match Percentage (%)');
title('Accuracy vs Performance Trade-off');
colorbar('Ticks', 1:n_tests, 'TickLabels', test_names);
grid on;
hold on;

% Add reference lines
plot([min(speedups), max(speedups)], [100, 100], 'g--', 'LineWidth', 1.5);
plot([min(speedups), max(speedups)], [95, 95], 'r--', 'LineWidth', 1);

% Annotate points
for i = 1:n_tests
    if match_percents(i) < 99 || speedups(i) > 35
        text(speedups(i), match_percents(i), sprintf('  %d', i), ...
            'FontSize', 8, 'FontWeight', 'bold');
    end
end

hold off;
ylim([90, 105]);

% Save figure
saveas(gcf, 'accuracy_vs_performance.png');
fprintf('  Saved: accuracy_vs_performance.png\n');

%% Figure 4: Spatial neighbor comparison (for one test case)
% Pick a representative test case (medium grid, moderate NaN ratio)
test_idx = min(2, n_tests);  % Second test or first if only one

fprintf('\nCreating spatial neighbor comparison for test case: %s\n', ...
    results(test_idx).config.name);

% Note: This requires storing neighbor coordinates in results
% For now, create a placeholder
figure('Position', [250, 250, 1000, 400]);

subplot(1, 2, 1);
text(0.5, 0.5, {'\bf Spatial Neighbor Comparison', '', ...
    'To enable this visualization, modify', ...
    'compare_neighbours_implementations.m to save:', ...
    '- results(i).csub_ref', ...
    '- results(i).csub_opt', ...
    '- results(i).p0', ...
    '', 'Then this plot will show:', ...
    '- Target point location', ...
    '- Neighbors selected by each method', ...
    '- Differences highlighted'}, ...
    'HorizontalAlignment', 'center', 'FontSize', 11);
axis off;

subplot(1, 2, 2);
text(0.5, 0.5, {'\bf Instructions:', '', ...
    '1. Update compare script to save coordinates', ...
    '2. Re-run comparison', ...
    '3. Re-run this visualization', ...
    '', 'The spatial plot will show neighbor', ...
    'positions in lon-lat-time space'}, ...
    'HorizontalAlignment', 'center', 'FontSize', 10);
axis off;

sgtitle('Spatial Neighbor Comparison (Placeholder)', 'FontSize', 14);

% Save figure
saveas(gcf, 'spatial_comparison.png');
fprintf('  Saved: spatial_comparison.png\n');

%% Overall Summary
fprintf('\n=== Visualization Summary ===\n\n');

fprintf('Generated 4 accuracy visualization figures:\n');
fprintf('  1. accuracy_comparison.png - Match quality, distance/value accuracy\n');
fprintf('  2. accuracy_detailed.png - Detailed accuracy distributions\n');
fprintf('  3. accuracy_vs_performance.png - Accuracy-speedup trade-off\n');
fprintf('  4. spatial_comparison.png - Spatial neighbor comparison\n');
fprintf('\n');

fprintf('Key Findings:\n');
fprintf('  Average match quality: %.2f%%\n', mean(match_percents));
fprintf('  Tests with >99%% match: %d/%d\n', sum(match_percents >= 99), n_tests);
fprintf('  Tests with >95%% match: %d/%d\n', sum(match_percents >= 95), n_tests);
fprintf('  Max distance error: %.2e\n', max(max_dist_diffs));
fprintf('  Max value error: %.2e\n', max(max_value_diffs));
fprintf('\n');

% Overall verdict
if mean(match_percents) >= 99.5 && max(max_dist_diffs) < 1e-6
    fprintf('✓ EXCELLENT: Optimized version matches reference implementation\n');
    fprintf('  Both functions find essentially identical neighbors\n');
elseif mean(match_percents) >= 95 && max(max_dist_diffs) < 1e-4
    fprintf('✓ GOOD: Optimized version is accurate\n');
    fprintf('  Minor differences within acceptable tolerances\n');
else
    fprintf('⚠ WARNING: Some discrepancies detected\n');
    fprintf('  Review test cases with <95%% match\n');
end

fprintf('\n=== Accuracy Visualization Complete ===\n');
