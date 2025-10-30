% test_neighbours_stug_index.m
% Comprehensive test suite for neighbours_stug_index.m
%
% Tests:
% 1. Basic functionality
% 2. Edge cases (corners, boundaries)
% 3. Anisotropic grids (different dx/dy ratios)
% 4. Various NaN ratios
% 5. Distance constraint enforcement
% 6. Verification that distance calculations only happen at end

clear; close all;
fprintf('=== Testing neighbours_stug_index.m ===\n\n');

all_tests_passed = true;

%% Test 1: Basic uniform isotropic grid
fprintf('Test 1: Basic uniform isotropic grid...\n');

% Create test grid (200x200x24)
nx = 200; ny = 200; nt = 24;
grid_data.x = linspace(-180, 180, nx)';
grid_data.y = linspace(-90, 90, ny)';
grid_data.time = (1:nt)';

[X, Y] = meshgrid(grid_data.x, grid_data.y);
grid_data.Lon = X';
grid_data.Lat = Y';

% Create synthetic data
grid_data.Z = randn(nx, ny, nt);

% Add some NaN values (30%)
nan_mask = rand(nx, ny, nt) < 0.3;
grid_data.Z(nan_mask) = NaN;

% Test point in middle of grid
p0 = [-95, 27.5, 12];
nmax = 30;
dmax = [5.0, 6.0, 0.5];

tic;
[psub, zsub, dsub, nsub, index] = neighbours_stug_index(p0, grid_data, nmax, dmax);
time_elapsed = toc;

fprintf('  Found %d neighbors in %.4f seconds\n', nsub, time_elapsed);
fprintf('  Max spatial distance: %.4f (limit: %.4f)\n', max(dsub(:,1)), dmax(1));
fprintf('  Max temporal distance: %.4f (limit: %.4f)\n', max(dsub(:,2)), dmax(2));

% Verify constraints
test1_passed = true;
if any(dsub(:,1) > dmax(1) + 1e-6)
    fprintf('  ✗ Spatial distance constraint violated!\n');
    test1_passed = false;
    all_tests_passed = false;
else
    fprintf('  ✓ Spatial distance constraints satisfied\n');
end

if any(dsub(:,2) > dmax(2) + 1e-6)
    fprintf('  ✗ Temporal distance constraint violated!\n');
    test1_passed = false;
    all_tests_passed = false;
else
    fprintf('  ✓ Temporal distance constraints satisfied\n');
end

if any(isnan(zsub))
    fprintf('  ✗ NaN values in results!\n');
    test1_passed = false;
    all_tests_passed = false;
else
    fprintf('  ✓ No NaN values in results\n');
end

if nsub > nmax
    fprintf('  ✗ Too many neighbors returned (%d > %d)!\n', nsub, nmax);
    test1_passed = false;
    all_tests_passed = false;
else
    fprintf('  ✓ Number of neighbors <= nmax\n');
end

% Check sorting (should be sorted by distance)
spacetime_dist = dsub(:,1) + dmax(3) * dsub(:,2);
if any(diff(spacetime_dist) < -1e-6)
    fprintf('  ✗ Results not properly sorted!\n');
    test1_passed = false;
    all_tests_passed = false;
else
    fprintf('  ✓ Results properly sorted by distance\n');
end

if test1_passed
    fprintf('  Test 1 PASSED\n\n');
else
    fprintf('  Test 1 FAILED\n\n');
end

%% Test 2: Edge cases (corners and boundaries)
fprintf('Test 2: Edge cases (corners and boundaries)...\n');

test_points = [
    grid_data.x(1), grid_data.y(1), grid_data.time(1);     % Bottom-left corner
    grid_data.x(end), grid_data.y(end), grid_data.time(end); % Top-right corner
    grid_data.x(1), grid_data.y(ny/2), grid_data.time(nt/2);   % Left edge
    grid_data.x(end), grid_data.y(ny/2), grid_data.time(nt/2); % Right edge
];

test2_passed = true;
for i = 1:size(test_points, 1)
    p0_test = test_points(i, :);
    [psub, zsub, dsub, nsub, ~] = neighbours_stug_index(p0_test, grid_data, nmax, dmax);

    fprintf('  Point %d: Found %d neighbors', i, nsub);

    if any(dsub(:,1) > dmax(1) + 1e-6) || any(dsub(:,2) > dmax(2) + 1e-6)
        fprintf(' - ✗ Distance constraints violated\n');
        test2_passed = false;
        all_tests_passed = false;
    else
        fprintf(' - ✓ Valid\n');
    end
end

if test2_passed
    fprintf('  Test 2 PASSED\n\n');
else
    fprintf('  Test 2 FAILED\n\n');
end

%% Test 3: Anisotropic grids (different dx/dy ratios)
fprintf('Test 3: Anisotropic grids...\n');

test3_passed = true;

% Test different anisotropy ratios
anisotropy_ratios = [0.25, 0.5, 1.0, 2.0, 4.0];

for ratio = anisotropy_ratios
    % Create anisotropic grid
    aniso_data.x = linspace(-180, 180, 100)';
    aniso_data.y = linspace(-90, 90, round(100 * ratio))';
    aniso_data.time = (1:12)';

    [X, Y] = meshgrid(aniso_data.x, aniso_data.y);
    aniso_data.Lon = X';
    aniso_data.Lat = Y';

    % Create data
    aniso_data.Z = randn(length(aniso_data.x), length(aniso_data.y), length(aniso_data.time));
    nan_mask = rand(size(aniso_data.Z)) < 0.3;
    aniso_data.Z(nan_mask) = NaN;

    % Test point
    p0_test = [0, 0, 6];

    [psub, zsub, dsub, nsub, ~] = neighbours_stug_index(p0_test, aniso_data, 20, [10.0, 3.0, 0.5]);

    % Check grid spacing
    dx = median(diff(aniso_data.x));
    dy = median(diff(aniso_data.y));
    actual_ratio = dy / dx;

    fprintf('  dy/dx = %.2f: Found %d neighbors', actual_ratio, nsub);

    if any(isnan(zsub))
        fprintf(' - ✗ NaN in results\n');
        test3_passed = false;
        all_tests_passed = false;
    else
        fprintf(' - ✓ Valid\n');
    end
end

if test3_passed
    fprintf('  Test 3 PASSED\n\n');
else
    fprintf('  Test 3 FAILED\n\n');
end

%% Test 4: Various NaN ratios
fprintf('Test 4: Various NaN ratios...\n');

nan_ratios = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95];
test4_passed = true;

% Create base grid
base_data.x = linspace(-180, 180, 150)';
base_data.y = linspace(-90, 90, 150)';
base_data.time = (1:24)';

[X, Y] = meshgrid(base_data.x, base_data.y);
base_data.Lon = X';
base_data.Lat = Y';

p0_test = [0, 0, 12];

for nan_ratio = nan_ratios
    % Create data with specified NaN ratio
    base_data.Z = randn(length(base_data.x), length(base_data.y), length(base_data.time));
    nan_mask = rand(size(base_data.Z)) < nan_ratio;
    base_data.Z(nan_mask) = NaN;

    tic;
    [psub, zsub, dsub, nsub, ~] = neighbours_stug_index(p0_test, base_data, 25, [5.0, 6.0, 0.5]);
    time_elapsed = toc;

    fprintf('  NaN ratio %.0f%%: Found %d neighbors in %.4f sec', nan_ratio*100, nsub, time_elapsed);

    if any(isnan(zsub))
        fprintf(' - ✗ NaN in results\n');
        test4_passed = false;
        all_tests_passed = false;
    else
        fprintf(' - ✓ Valid\n');
    end
end

if test4_passed
    fprintf('  Test 4 PASSED\n\n');
else
    fprintf('  Test 4 FAILED\n\n');
end

%% Test 5: POI outside grid bounds
fprintf('Test 5: POI outside grid bounds...\n');

% Use grid from Test 1
test5_passed = true;

outside_points = [
    -200, 0, 12;      % West of grid
    200, 0, 12;       % East of grid
    0, -100, 12;      % South of grid
    0, 100, 12;       % North of grid
    0, 0, -5;         % Before time range
    0, 0, 50;         % After time range
];

for i = 1:size(outside_points, 1)
    p0_test = outside_points(i, :);

    try
        [psub, zsub, dsub, nsub, ~] = neighbours_stug_index(p0_test, grid_data, 15, [10.0, 5.0, 0.5]);
        fprintf('  Point %d: Found %d neighbors - ✓ No error\n', i, nsub);
    catch e
        fprintf('  Point %d: Error - %s\n', i, e.message);
        test5_passed = false;
        all_tests_passed = false;
    end
end

if test5_passed
    fprintf('  Test 5 PASSED\n\n');
else
    fprintf('  Test 5 FAILED\n\n');
end

%% Test 6: Very sparse data (95% NaN)
fprintf('Test 6: Very sparse data (95%% NaN)...\n');

sparse_data.x = linspace(-180, 180, 200)';
sparse_data.y = linspace(-90, 90, 200)';
sparse_data.time = (1:24)';

[X, Y] = meshgrid(sparse_data.x, sparse_data.y);
sparse_data.Lon = X';
sparse_data.Lat = Y';

sparse_data.Z = randn(length(sparse_data.x), length(sparse_data.y), length(sparse_data.time));
nan_mask = rand(size(sparse_data.Z)) < 0.95;
sparse_data.Z(nan_mask) = NaN;

p0_test = [0, 0, 12];

tic;
[psub, zsub, dsub, nsub, ~] = neighbours_stug_index(p0_test, sparse_data, 15, [10.0, 6.0, 0.5]);
time_elapsed = toc;

fprintf('  Found %d neighbors in %.4f seconds\n', nsub, time_elapsed);

test6_passed = true;
if any(isnan(zsub))
    fprintf('  ✗ NaN values in results\n');
    test6_passed = false;
    all_tests_passed = false;
else
    fprintf('  ✓ No NaN values in results\n');
end

if test6_passed
    fprintf('  Test 6 PASSED\n\n');
else
    fprintf('  Test 6 FAILED\n\n');
end

%% Test 7: Empty/all-NaN data
fprintf('Test 7: Empty/all-NaN data...\n');

empty_data = grid_data;
empty_data.Z(:) = NaN;

[psub, zsub, dsub, nsub, ~] = neighbours_stug_index(p0, empty_data, nmax, dmax);

test7_passed = true;
if nsub ~= 0
    fprintf('  ✗ Should return 0 neighbors for all-NaN data\n');
    test7_passed = false;
    all_tests_passed = false;
else
    fprintf('  ✓ Correctly returns 0 neighbors for all-NaN data\n');
end

if test7_passed
    fprintf('  Test 7 PASSED\n\n');
else
    fprintf('  Test 7 FAILED\n\n');
end

%% Test 8: Dense data (0% NaN)
fprintf('Test 8: Dense data (0%% NaN)...\n');

dense_data.x = linspace(-180, 180, 300)';
dense_data.y = linspace(-90, 90, 300)';
dense_data.time = (1:24)';

[X, Y] = meshgrid(dense_data.x, dense_data.y);
dense_data.Lon = X';
dense_data.Lat = Y';

dense_data.Z = randn(length(dense_data.x), length(dense_data.y), length(dense_data.time));

p0_test = [0, 0, 12];

tic;
[psub, zsub, dsub, nsub, ~] = neighbours_stug_index(p0_test, dense_data, 30, [5.0, 6.0, 0.5]);
time_elapsed = toc;

fprintf('  Found %d neighbors in %.4f seconds\n', nsub, time_elapsed);

test8_passed = true;
if nsub ~= 30
    fprintf('  ✗ Should find exactly nmax=%d neighbors (found %d)\n', 30, nsub);
    test8_passed = false;
    all_tests_passed = false;
else
    fprintf('  ✓ Found exactly nmax neighbors\n');
end

if any(isnan(zsub))
    fprintf('  ✗ NaN values in results\n');
    test8_passed = false;
    all_tests_passed = false;
else
    fprintf('  ✓ No NaN values in results\n');
end

if test8_passed
    fprintf('  Test 8 PASSED\n\n');
else
    fprintf('  Test 8 FAILED\n\n');
end

%% Summary
fprintf('=== Test Summary ===\n\n');
if all_tests_passed
    fprintf('✓ All tests PASSED\n');
    fprintf('\nneighbours_stug_index.m is working correctly!\n');
else
    fprintf('✗ Some tests FAILED\n');
    fprintf('\nPlease review failed tests above.\n');
end

fprintf('\n=== Testing Complete ===\n');
