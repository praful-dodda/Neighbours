% compare_neighbours_index.m
% Compare neighbours_stug_index against neighbours.m (reference)
% and neighbours_stug_optimized. Saves detailed per-scenario results.

clear; close all;

fprintf('=== Comparison: neighbours_stug_index vs reference/optimized ===\n');

% Configuration
% Set RUN_HEAVY=1 in the environment to include very large scenarios.
run_heavy = strcmp(getenv('RUN_HEAVY'), '1');
repeats = 3;       % timing repeats (median)

% Test scenarios: {nx, ny, nt, nan_ratio, nmax, description}
scenarios = {
    {50, 50, 12, 0, 15, 'Small dense, 0% NaN'}
    {50, 50, 12, 0, 15, 'Small, 30% NaN'}
    {100, 100, 24, 0, 20, 'Medium, 10% NaN'}
    {100, 100, 24, 0, 20, 'Medium, 70% NaN'}
};

if run_heavy
    scenarios = [scenarios, {{150,150,50,0.0,25,'Large dense, 0% NaN'}, {200,200,24,0.3,25,'Very Large, 30% NaN (heavy)'}}];
end

n_scenarios = length(scenarios);
results = struct();

% Tolerances for matching
tolerance_spatial = 1e-6;

for s = 1:n_scenarios
    sc = scenarios{s};
    nx = sc{1}; ny = sc{2}; nt = sc{3}; nan_ratio = sc{4}; nmax = sc{5}; desc = sc{6};

    fprintf('\n--- Scenario %d/%d: %s ---\n', s, n_scenarios, desc);
    fprintf('Grid: %dx%dx%d, NaN ratio: %.1f%%, nmax: %d\n', nx, ny, nt, nan_ratio*100, nmax);

    % Create grid
    grid_data = create_test_grid(nx, ny, nt, nan_ratio);

    % Test point: center
    p0 = [(grid_data.x(1)+grid_data.x(end))/2, (grid_data.y(1)+grid_data.y(end))/2, (grid_data.time(1)+grid_data.time(end))/2];

    dmax = [5.0, 10.0, 0.1];

    % Prepare vector format for neighbours.m
    [c_vector, Z_vector] = grid_to_vector(grid_data);

    % Reference: neighbours.m
    times_ref = zeros(repeats,1);
    for r = 1:repeats
        tic; [csub_ref, zsub_ref, dsub_ref, nsub_ref, index_ref] = neighbours(p0, c_vector, Z_vector, nmax, dmax); times_ref(r) = toc;
    end
    time_ref = median(times_ref);
    fprintf('  neighbours.m:   %.4fs, n=%d\n', time_ref, nsub_ref);

    % Optimized
    times_opt = zeros(repeats,1);
    for r = 1:repeats
        tic; [psub_opt, zsub_opt, dsub_opt, nsub_opt, index_opt] = neighbours_stug_optimized(p0, grid_data, nmax, dmax); times_opt(r) = toc;
    end
    time_opt = median(times_opt);
    fprintf('  optimized:      %.4fs, n=%d\n', time_opt, nsub_opt);

    % Index-based
    times_idx = zeros(repeats,1);
    for r = 1:repeats
        tic; [psub_idx, zsub_idx, dsub_idx, nsub_idx, index_idx] = neighbours_stug_index(p0, grid_data, nmax, dmax); times_idx(r) = toc;
    end
    time_idx = median(times_idx);
    fprintf('  index-based:    %.4fs, n=%d\n', time_idx, nsub_idx);

    % Compare neighbor sets (reference vs index)
    percent_match_idx = 0; max_dist_diff = NaN; max_value_diff = NaN;
    if nsub_ref>0 && nsub_idx>0
        [match_info_idx, details_idx] = compare_neighbor_sets(csub_ref, zsub_ref, dsub_ref, psub_idx(:,1:3), zsub_idx, dsub_idx, dmax);
        percent_match_idx = match_info_idx.percent_match;
        max_dist_diff = match_info_idx.max_distance_diff;
        max_value_diff = match_info_idx.max_value_diff;
    end

    % Store results
    results(s).config = struct('nx',nx,'ny',ny,'nt',nt,'nan_ratio',nan_ratio,'nmax',nmax,'desc',desc);
    results(s).time_ref = time_ref;
    results(s).time_opt = time_opt;
    results(s).time_idx = time_idx;
    results(s).nsub_ref = nsub_ref;
    results(s).nsub_opt = nsub_opt;
    results(s).nsub_idx = nsub_idx;
    results(s).percent_match_idx = percent_match_idx;
    results(s).max_dist_diff = max_dist_diff;
    results(s).max_value_diff = max_value_diff;

    % Print quick verdict
    fprintf('  Match (index vs ref): %.1f%%, max dist diff: %g, max value diff: %g\n', percent_match_idx, max_dist_diff, max_value_diff);
    fprintf('  Speedups: opt=%.2fx, idx=%.2fx\n', time_ref/time_opt, time_ref/time_idx);
end

% Save results
save('compare_neighbours_index_results.mat', 'results', 'scenarios');
fprintf('\nSaved results to compare_neighbours_index_results.mat\n');

%% Helper functions
function grid_data = create_test_grid(nx, ny, nt, nan_ratio)
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
                grid_data.Z(ix, iy, it) = sin(ix/10) * cos(iy/10) * exp(-it/50) + 0.1*randn();
            end
        end
    end

    if nan_ratio > 0
        nan_mask = rand(nx, ny, nt) < nan_ratio;
        grid_data.Z(nan_mask) = NaN;
    end
end

function [c_vector, Z_vector] = grid_to_vector(grid_data)
    nx = length(grid_data.x);
    ny = length(grid_data.y);
    nt = length(grid_data.time);

    c_all = zeros(nx*ny*nt, 3);
    Z_all = zeros(nx*ny*nt,1);
    idx = 1;
    for it = 1:nt
        for iy = 1:ny
            for ix = 1:nx
                c_all(idx,1) = grid_data.Lon(ix,iy);
                c_all(idx,2) = grid_data.Lat(ix,iy);
                c_all(idx,3) = grid_data.time(it);
                Z_all(idx) = grid_data.Z(ix,iy,it);
                idx = idx + 1;
            end
        end
    end
    valid_idx = ~isnan(Z_all);
    c_vector = c_all(valid_idx,:);
    Z_vector = Z_all(valid_idx);
end

function [match_quality, details] = compare_neighbor_sets(...
    csub_ref, Zsub_ref, dsub_ref, csub_cmp, Zsub_cmp, dsub_cmp, dmax)
    match_quality = struct(); details = struct();

    n_ref = size(csub_ref,1);
    n_cmp = size(csub_cmp,1);
    if n_ref==0 || n_cmp==0
        match_quality.percent_match = 0; match_quality.distances_match = false; match_quality.values_match = false; match_quality.max_distance_diff = NaN; match_quality.max_value_diff = NaN; return;
    end

    tolerance_spatial = 1e-6;
    matched_ref = false(n_ref,1);
    distance_diffs = zeros(n_ref,1);
    value_diffs = zeros(n_ref,1);

    for i=1:n_ref
        coord_diff = sqrt(sum((csub_cmp - repmat(csub_ref(i,:), n_cmp,1)).^2,2));
        [min_diff, idx] = min(coord_diff);
        if min_diff < tolerance_spatial
            matched_ref(i) = true;
            % compute combined distances
            if size(dsub_ref,2) == 1
                dist_ref = dsub_ref(i);
                dist_cmp = sqrt(dsub_cmp(idx,1)^2 + (dmax(3)*dsub_cmp(idx,2))^2);
            else
                dist_ref = sqrt(dsub_ref(i,1)^2 + (dmax(3)*dsub_ref(i,2))^2);
                dist_cmp = sqrt(dsub_cmp(idx,1)^2 + (dmax(3)*dsub_cmp(idx,2))^2);
            end
            distance_diffs(i) = abs(dist_ref - dist_cmp);
            value_diffs(i) = abs(Zsub_ref(i) - Zsub_cmp(idx));
        end
    end

    n_matched = sum(matched_ref);
    match_quality.percent_match = (n_matched / max(n_ref, n_cmp)) * 100;
    match_quality.n_ref = n_ref; match_quality.n_cmp = n_cmp; match_quality.n_matched = n_matched;
    match_quality.max_distance_diff = max(distance_diffs(matched_ref));
    match_quality.max_value_diff = max(value_diffs(matched_ref));
    match_quality.distances_match = match_quality.max_distance_diff < 1e-6;
    match_quality.values_match = match_quality.max_value_diff < 1e-9;
    details.matched_ref = matched_ref; details.distance_diffs = distance_diffs; details.value_diffs = value_diffs;
end
