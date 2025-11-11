# Data Reformatting Guide: STG to STUG

## Overview

This guide explains how to convert your data to STUG format for use with `neighbours_stug_optimized.m` or `neighbours_stug_index.m`.

**Important:** There are TWO different functions depending on your data type!

---

## Quick Decision: Which Function Should I Use?

### ✅ Use `reformat_stg_to_stug.m` if:
- Your data is **already on a uniform grid**
- Grid points have **regular spacing** (constant dx, dy)
- You just need to **reshape** the data structure
- Example: Satellite data, model output stored in vectorized format

### ✅ Use `stg_to_stug.m` if:
- Your data is from **irregular monitoring stations**
- Station locations are **not uniformly spaced**
- You need to **interpolate** onto a uniform grid
- Example: Weather station network, environmental sensors

---

## Option 1: Reformat Uniform Grid (NO Interpolation)

### For Data ALREADY on Uniform Grid

**Your case:** `KS.softdata` with uniform grid stored as vectors

```matlab
% Quick usage
grid_data = reformat_stg_to_stug(KS.softdata);

% Use with neighbor search
[psub, zsub, dsub, nsub] = neighbours_stug_index(p0, grid_data, nmax, dmax);
```

### What This Does:

1. **Validates** grid is uniformly spaced
2. **Extracts** unique x and y values
3. **Reshapes** data from [N×T] to [nx×ny×nt]
4. **Creates** proper STUG structure

**Speed:** < 1 second
**Accuracy:** Perfect (no data loss)
**Memory:** Minimal overhead

### Example Output:

```
=== Reformatting STG to STUG Structure ===
Input: 78048 grid points × 19 time steps
Detected grid dimensions: 272 (lon) × 287 (lat)
Longitude spacing: 0.100000° (uniform ✓)
Latitude spacing: 0.100000° (uniform ✓)
Reshaping data from [78048×19] to [272×287×19]...

=== Reformatting Complete ===
Output grid: 272 × 287 × 19
Grid spacing: 0.100000° × 0.100000°
Valid data points: 1482912 / 1482912 (100.0%)
Structure is ready for neighbours_stug_* functions
```

### Error Cases:

```matlab
% If grid is NOT uniform, you'll get a clear error:
Error: Grid is NOT uniformly spaced in longitude!
Mean spacing: 0.100000 degrees
Std deviation: 0.005000 degrees
Relative variation: 5.00% (tolerance: 0.0001%)
This function requires uniform grid spacing.
Use stg_to_stug.m for interpolation if you have irregular stations.
```

### Input Requirements:

```matlab
% Your data must have:
stg_data.sMS   % [N×2] Grid coordinates (lon, lat) - MUST be uniform!
stg_data.tME   % [1×T] Time values
stg_data.Xms   % [N×T] Data values
stg_data.Xvs   % [N×T] Variances (optional)

% Where N = nx × ny (complete rectangular grid)
```

### Output Structure:

```matlab
grid_data.x      % [nx×1] Longitude vector
grid_data.y      % [ny×1] Latitude vector
grid_data.time   % [nt×1] Time vector
grid_data.Lon    % [nx×ny] Longitude mesh
grid_data.Lat    % [nx×ny] Latitude mesh
grid_data.Z      % [nx×ny×nt] Data array
grid_data.Zvar   % [nx×ny×nt] Variance (if provided)
grid_data.metadata  % Reformatting information
```

---

## Option 2: Interpolate Irregular Stations (WITH Interpolation)

### For Irregular Station Data

**Use case:** Monitoring stations at random locations

```matlab
% Interpolate onto uniform grid
grid_data = stg_to_stug(station_data, 0.1);  % 0.1° resolution

% Use with neighbor search
[psub, zsub, dsub, nsub] = neighbours_stug_index(p0, grid_data, nmax, dmax);
```

### What This Does:

1. **Determines** spatial extent of stations
2. **Creates** uniform grid at specified resolution
3. **Interpolates** station data onto grid
4. **Fills** gaps between stations

**Speed:** 2-10 seconds (depends on resolution and station count)
**Accuracy:** Interpolation introduces error
**Memory:** 2-5x original data size

See `STG_TO_STUG_GUIDE.md` for detailed documentation.

---

## Complete Workflow for Uniform Grid

### Step-by-Step

```matlab
%% Step 1: Load your data
% KS.softdata with:
%   sMS: [78048×2] uniform grid coordinates
%   tME: [1×19] time values
%   Xms: [78048×19] data

%% Step 2: Reformat to STUG structure
fprintf('Reformatting data...\n');
tic;
grid_data = reformat_stg_to_stug(KS.softdata);
fprintf('Reformatting took %.3f seconds\n\n', toc);

%% Step 3: Inspect the result
fprintf('Grid information:\n');
fprintf('  Dimensions: %d × %d × %d\n', ...
        length(grid_data.x), length(grid_data.y), length(grid_data.time));
fprintf('  Spacing: %.6f° × %.6f°\n', ...
        grid_data.metadata.spacing(1), grid_data.metadata.spacing(2));
fprintf('  Total points: %d\n', numel(grid_data.Z));

%% Step 4: Use with neighbor search

p0 = [-95.0, 30.0, 2017.5];  % Point of interest
nmax = 50;
dmax = [500.0, 365, 0.01];   % [spatial_km, temporal, weight]

% Option A: Use neighbours_stug_optimized
[psub, zsub, dsub, nsub] = neighbours_stug_optimized(p0, grid_data, nmax, dmax);

% Option B: Use neighbours_stug_index (Haversine distance)
[psub, zsub, dsub, nsub] = neighbours_stug_index(p0, grid_data, nmax, dmax);

fprintf('\nFound %d neighbors\n', nsub);
fprintf('Distance range: %.1f - %.1f km\n', min(dsub(:,1)), max(dsub(:,1)));

%% Step 5: (Optional) Visualize
figure;
subplot(1,2,1);
% Show grid structure
imagesc(grid_data.x, grid_data.y, grid_data.Z(:,:,10)');
axis xy; colorbar;
title('Reformatted grid (time index 10)');
xlabel('Longitude'); ylabel('Latitude');

subplot(1,2,2);
% Show neighbors
scatter3(psub(:,1), psub(:,2), psub(:,3), 50, zsub, 'filled');
hold on;
plot3(p0(1), p0(2), p0(3), 'r*', 'MarkerSize', 15, 'LineWidth', 2);
colorbar;
title('Selected neighbors');
xlabel('Longitude'); ylabel('Latitude'); zlabel('Time');
grid on;
```

---

## Function Comparison

| Feature | reformat_stg_to_stug | stg_to_stug |
|---------|---------------------|-------------|
| **Input data** | Uniform grid | Irregular stations |
| **Operation** | Reshape only | Interpolation |
| **Speed** | < 1 second | 2-10 seconds |
| **Accuracy** | Perfect | Interpolation error |
| **Memory** | Minimal | 2-5x original |
| **Data loss** | None | None, but adds interpolated values |
| **Use when** | Grid already uniform | Need to create grid from stations |

---

## Validation: Is My Grid Uniform?

### Quick Check

```matlab
% Extract unique coordinates
unique_lon = unique(KS.softdata.sMS(:, 1));
unique_lat = unique(KS.softdata.sMS(:, 2));

% Check spacing
dx = diff(unique_lon);
dy = diff(unique_lat);

fprintf('Longitude spacing:\n');
fprintf('  Mean: %.6f°\n', mean(dx));
fprintf('  Std:  %.6f°\n', std(dx));
fprintf('  Variation: %.2f%%\n', 100 * std(dx) / abs(mean(dx)));

fprintf('Latitude spacing:\n');
fprintf('  Mean: %.6f°\n', mean(dy));
fprintf('  Std:  %.6f°\n', std(dy));
fprintf('  Variation: %.2f%%\n', 100 * std(dy) / abs(mean(dy)));

% Grid is uniform if variation < 0.01%
if std(dx) / abs(mean(dx)) < 1e-4 && std(dy) / abs(mean(dy)) < 1e-4
    fprintf('\n✓ Grid is UNIFORM - use reformat_stg_to_stug\n');
else
    fprintf('\n✗ Grid is IRREGULAR - use stg_to_stug\n');
end
```

### Expected Output for Uniform Grid:

```
Longitude spacing:
  Mean: 0.100000°
  Std:  0.000000°
  Variation: 0.00%

Latitude spacing:
  Mean: 0.100000°
  Std:  0.000000°
  Variation: 0.00%

✓ Grid is UNIFORM - use reformat_stg_to_stug
```

---

## Advanced Options

### Custom Tolerance

```matlab
% Default tolerance: 1e-6 (0.0001% variation allowed)
grid_data = reformat_stg_to_stug(KS.softdata);

% Relaxed tolerance (for grids with minor numerical noise)
grid_data = reformat_stg_to_stug(KS.softdata, 'tolerance', 1e-4);

% Strict tolerance
grid_data = reformat_stg_to_stug(KS.softdata, 'tolerance', 1e-8);
```

### Suppress Output

```matlab
% Quiet mode (no console output)
grid_data = reformat_stg_to_stug(KS.softdata, 'verbose', false);
```

---

## Common Issues

### Issue 1: "Grid validation failed - incomplete grid"

**Problem:** Number of points doesn't match grid dimensions

```
Error: Grid validation failed!
Unique lon points: 272, Unique lat points: 287
Expected grid size: 272 × 287 = 78064 points
Actual points: 78048
Grid appears to be incomplete or irregular.
```

**Solution:**
- Check if some grid points are missing
- Verify grid is truly rectangular
- Look for duplicate coordinates

```matlab
% Find missing points
unique_lon = unique(KS.softdata.sMS(:, 1));
unique_lat = unique(KS.softdata.sMS(:, 2));
expected_points = length(unique_lon) * length(unique_lat);
actual_points = size(KS.softdata.sMS, 1);

fprintf('Missing %d points\n', expected_points - actual_points);
```

### Issue 2: "Grid is NOT uniformly spaced"

**Problem:** Grid spacing varies too much

**Solution:**
- If variation is small (<1%), increase tolerance
- If variation is large, use `stg_to_stug.m` instead

```matlab
% Option 1: Increase tolerance (if variation is truly minor)
grid_data = reformat_stg_to_stug(KS.softdata, 'tolerance', 1e-3);

% Option 2: Use interpolation (if grid is irregular)
grid_data = stg_to_stug(KS.softdata, 0.1);
```

### Issue 3: "Dimensions don't match"

**Problem:** sMS and Xms have different number of points

**Solution:**
- Verify data integrity
- Check if data was partially loaded

```matlab
fprintf('sMS points: %d\n', size(KS.softdata.sMS, 1));
fprintf('Xms points: %d\n', size(KS.softdata.Xms, 1));
% These must match!
```

---

## Performance Tips

### For Large Datasets

```matlab
% If you have very large grids (millions of points):

% 1. Suppress verbose output
grid_data = reformat_stg_to_stug(KS.softdata, 'verbose', false);

% 2. Clear large variables after reformatting
clear KS;  % If no longer needed

% 3. Save reformatted data for reuse
save('grid_data.mat', 'grid_data', '-v7.3');
```

### Memory Estimation

```matlab
% Before reformatting, estimate memory:
nx = length(unique(KS.softdata.sMS(:, 1)));
ny = length(unique(KS.softdata.sMS(:, 2)));
nt = length(KS.softdata.tME);

memory_mb = (nx * ny * nt * 4) / 1024^2;  % single precision = 4 bytes
fprintf('Expected memory: %.2f MB\n', memory_mb);
```

---

## Summary

### For Uniform Grid Data (Your Case!)

```matlab
% Simple reformatting
grid_data = reformat_stg_to_stug(KS.softdata);

% Then use with STUG functions
[psub, zsub, dsub, nsub] = neighbours_stug_index(p0, grid_data, nmax, dmax);
```

**Benefits:**
- ✅ Fast (< 1 second)
- ✅ No data loss
- ✅ No interpolation error
- ✅ Minimal memory overhead

### For Irregular Stations

```matlab
% Interpolation required
grid_data = stg_to_stug(station_data, 0.1);

% Then use with STUG functions
[psub, zsub, dsub, nsub] = neighbours_stug_index(p0, grid_data, nmax, dmax);
```

See `STG_TO_STUG_GUIDE.md` for interpolation details.
