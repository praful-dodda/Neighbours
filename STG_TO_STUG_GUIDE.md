# STG to STUG Conversion Guide

## Overview

This guide explains how to convert your STG (irregular stations) data to STUG (uniform grid) format for use with `neighbours_stug_optimized.m` or `neighbours_stug_index.m`.

---

## Quick Start

### Basic Usage

```matlab
% Load your data
% KS.softdata contains your STG format data

% Convert to STUG format (auto resolution)
grid_data = stg_to_stug(KS.softdata);

% Use with neighbor search
p0 = [-95.0, 30.0, 2017.5];  % [lon, lat, time]
nmax = 50;
dmax = [500.0, 365, 0.01];   % [spatial_km, temporal, weight]

[psub, zsub, dsub, nsub] = neighbours_stug_index(p0, grid_data, nmax, dmax);
```

**Expected Performance:**
- Conversion time: 2-10 seconds (one-time cost)
- Neighbor search: 5-10 ms per query
- Memory: ~11-45 MB (depends on resolution)

---

## Function Reference

### `stg_to_stug(stg_data, ...)`

**Purpose:** Convert irregular station data to uniform grid via interpolation

**Basic Syntax:**
```matlab
grid_data = stg_to_stug(stg_data);                    % Auto resolution
grid_data = stg_to_stug(stg_data, resolution);        % Specify resolution
grid_data = stg_to_stug(stg_data, 'Name', Value, ...); % Advanced options
```

**Input Requirements:**
```matlab
stg_data must have:
  .sMS  [N×2 double]   - Station locations [lon, lat]
  .tME  [1×T double]   - Time values
  .Xms  [N×T single]   - Data values
  .Xvs  [N×T single]   - Variances (optional)
  .Zisnotnan [N×T logical] - Valid data mask (optional)
```

**Output Structure:**
```matlab
grid_data contains:
  .x    [nx×1 double]   - Longitude grid vector
  .y    [ny×1 double]   - Latitude grid vector
  .time [nt×1 double]   - Time vector
  .Lon  [nx×ny double]  - Longitude mesh
  .Lat  [nx×ny double]  - Latitude mesh
  .Z    [nx×ny×nt single] - Interpolated data
  .Zvar [nx×ny×nt single] - Interpolated variance (if requested)
  .metadata - Conversion information
```

---

## Resolution Selection

### Auto-Calculated Resolution (Recommended)

```matlab
grid_data = stg_to_stug(KS.softdata);
```

The function automatically calculates optimal resolution as **50% of mean nearest neighbor distance** between stations.

**For your data (78,048 stations):**
- Typical mean NN distance: ~0.5-2.0 degrees
- Auto resolution: ~0.25-1.0 degrees
- Result: Good balance of detail vs speed/memory

### Manual Resolution Selection

```matlab
% High detail (slow, large memory)
grid_data = stg_to_stug(KS.softdata, 0.05);  % 0.05° ≈ 5.5 km

% Balanced (recommended)
grid_data = stg_to_stug(KS.softdata, 0.1);   % 0.1° ≈ 11 km

% Fast (coarse resolution)
grid_data = stg_to_stug(KS.softdata, 0.2);   % 0.2° ≈ 22 km
```

### Resolution Impact Table

**For your data extent (~50° × 30° × 19 times):**

| Resolution | Grid Size | Total Points | Memory | Conv. Time | Detail |
|------------|-----------|--------------|--------|------------|--------|
| 0.05° | 1000×600×19 | 11.4M | 45 MB | ~15 sec | Excellent |
| 0.10° | 500×300×19 | 2.9M | 11 MB | ~5 sec | Good ✅ |
| 0.20° | 250×150×19 | 0.7M | 3 MB | ~2 sec | Fair |
| 0.50° | 100×60×19 | 0.1M | 0.5 MB | <1 sec | Poor |

**Recommendation:** Start with 0.1° (good balance)

---

## Advanced Options

### Interpolation Methods

```matlab
% Linear interpolation (default, fastest)
grid_data = stg_to_stug(KS.softdata, 'method', 'linear');

% Natural neighbor (smoother, slower)
grid_data = stg_to_stug(KS.softdata, 'method', 'natural');

% Nearest neighbor (no smoothing, fastest)
grid_data = stg_to_stug(KS.softdata, 'method', 'nearest');
```

**Comparison:**

| Method | Speed | Smoothness | Use Case |
|--------|-------|------------|----------|
| `'linear'` | Fast | Medium | General purpose ✅ |
| `'natural'` | Slow | Smooth | High quality output |
| `'nearest'` | Fastest | None | Categorical data |

### Custom Spatial Bounds

```matlab
% Focus on specific region
grid_data = stg_to_stug(KS.softdata, ...
    'resolution', 0.1, ...
    'bounds', [-100, -90, 25, 35]);  % [lon_min lon_max lat_min lat_max]
```

**Benefits:**
- Smaller grid = faster + less memory
- Focus on region of interest
- Exclude sparse boundary areas

### Include Variance Interpolation

```matlab
% Also interpolate variance field
grid_data = stg_to_stug(KS.softdata, ...
    'resolution', 0.1, ...
    'includeVariance', true);

% Access interpolated variance
variance_grid = grid_data.Zvar;
```

### Suppress Progress Output

```matlab
% Silent mode (no console output)
grid_data = stg_to_stug(KS.softdata, ...
    'resolution', 0.1, ...
    'verbose', false);
```

---

## Complete Workflow Example

### Step-by-Step with Your Data

```matlab
%% Step 1: Load your data
% Assuming KS.softdata is already loaded with structure:
%   sMS: [78048×2] stations
%   tME: [1×19] times
%   Xms: [78048×19] data

%% Step 2: Convert to STUG format
fprintf('Converting STG to STUG...\n');
tic;
grid_data = stg_to_stug(KS.softdata, ...
    'resolution', 0.1, ...           % 0.1 degree spacing
    'method', 'linear', ...          % Linear interpolation
    'verbose', true);                % Show progress
conversion_time = toc;

fprintf('Conversion complete: %.2f seconds\n', conversion_time);

%% Step 3: Inspect the result
fprintf('\nGrid information:\n');
fprintf('  Grid size: %d × %d × %d\n', ...
        length(grid_data.x), length(grid_data.y), length(grid_data.time));
fprintf('  Resolution: %.4f degrees\n', grid_data.metadata.resolution);
fprintf('  Fill ratio: %.1f%%\n', grid_data.metadata.fill_ratio * 100);
fprintf('  Memory: %.2f MB\n', numel(grid_data.Z) * 4 / 1024^2);

%% Step 4: Use with neighbor search

% Define point of interest
p0 = [-95.0, 30.0, 2017.5];  % [lon, lat, time]
nmax = 50;                    % Find 50 nearest neighbors
dmax = [500.0, 365, 0.01];   % [500 km spatial, 365 days temporal]

% Option A: Use neighbours_stug_optimized
[psub, zsub, dsub, nsub] = neighbours_stug_optimized(p0, grid_data, nmax, dmax);

% Option B: Use neighbours_stug_index (Haversine distance)
[psub, zsub, dsub, nsub] = neighbours_stug_index(p0, grid_data, nmax, dmax);

fprintf('\nFound %d neighbors:\n', nsub);
fprintf('  Spatial distances: %.1f - %.1f km\n', min(dsub(:,1)), max(dsub(:,1)));
fprintf('  Temporal distances: %.1f - %.1f days\n', min(dsub(:,2)), max(dsub(:,2)));

%% Step 5: Visualize (optional)
figure;
subplot(1,2,1);
% Plot one time slice
t_idx = 10;
imagesc(grid_data.x, grid_data.y, grid_data.Z(:,:,t_idx)');
axis xy; colorbar;
title(sprintf('Interpolated grid (time %d)', t_idx));
xlabel('Longitude'); ylabel('Latitude');
hold on;
% Overlay original stations
scatter(KS.softdata.sMS(:,1), KS.softdata.sMS(:,2), 5, 'k.', 'filled');

subplot(1,2,2);
% Plot neighbors
scatter3(psub(:,1), psub(:,2), psub(:,3), 50, zsub, 'filled');
hold on;
plot3(p0(1), p0(2), p0(3), 'r*', 'MarkerSize', 15, 'LineWidth', 2);
colorbar;
xlabel('Longitude'); ylabel('Latitude'); zlabel('Time');
title('Selected neighbors');
grid on;
```

---

## Performance Optimization Tips

### 1. Cache Converted Grid for Multiple Queries

```matlab
% Convert once
grid_data = stg_to_stug(KS.softdata, 0.1);

% Save for reuse
save('cached_grid_data.mat', 'grid_data', '-v7.3');

% Later sessions: load instead of reconverting
load('cached_grid_data.mat');

% Now run many queries
for i = 1:1000
    p0 = get_query_point(i);
    [psub, zsub, dsub, nsub] = neighbours_stug_index(p0, grid_data, nmax, dmax);
    % Process results...
end
```

**When this makes sense:**
- ✅ You have 100+ queries on same dataset
- ✅ Conversion takes >5 seconds
- ✅ You can afford the storage space
- ❌ Only a few queries (not worth it)

### 2. Choose Resolution Wisely

```matlab
% Too fine: Wastes memory and time
grid_data = stg_to_stug(KS.softdata, 0.01);  % ❌ Probably overkill

% Too coarse: Loses important details
grid_data = stg_to_stug(KS.softdata, 1.0);   % ❌ Too crude

% Just right: Balance detail and performance
grid_data = stg_to_stug(KS.softdata, 0.1);   % ✅ Good choice
```

**Rule of thumb:**
- Resolution ≈ 0.3-0.5 × mean station spacing
- Start with 0.1° for typical datasets
- Refine based on visual inspection

### 3. Use Custom Bounds to Reduce Grid Size

```matlab
% Full extent (larger grid)
grid_data_full = stg_to_stug(KS.softdata, 0.1);
% Result: 500×300×19 = 2.9M points

% Region of interest only (smaller grid)
grid_data_roi = stg_to_stug(KS.softdata, ...
    'resolution', 0.1, ...
    'bounds', [-100, -90, 25, 35]);
% Result: 100×100×19 = 0.2M points (15x smaller!)
```

---

## Quality Control

### Check Interpolation Quality

```matlab
% After conversion
grid_data = stg_to_stug(KS.softdata, 0.1);

% 1. Check metadata
disp(grid_data.metadata);

% 2. Check fill ratio (should be >50% for good coverage)
if grid_data.metadata.fill_ratio < 0.5
    warning('Low fill ratio (%.1f%%) - many grid points have no data', ...
            grid_data.metadata.fill_ratio * 100);
end

% 3. Compare interpolated vs original at station locations
t_idx = 10;
for i = 1:100  % Check first 100 stations
    % Find nearest grid point
    [~, ix] = min(abs(grid_data.x - KS.softdata.sMS(i,1)));
    [~, iy] = min(abs(grid_data.y - KS.softdata.sMS(i,2)));

    original = KS.softdata.Xms(i, t_idx);
    interpolated = grid_data.Z(ix, iy, t_idx);

    error = abs(original - interpolated);

    if error > 0.1 * abs(original)  % >10% error
        fprintf('Large error at station %d: %.2f vs %.2f\n', ...
                i, original, interpolated);
    end
end
```

### Visualize Interpolation Coverage

```matlab
% Create coverage map
figure;
coverage = ~isnan(grid_data.Z(:,:,10));
imagesc(grid_data.x, grid_data.y, coverage');
axis xy; colormap(gray);
colorbar('Ticks', [0 1], 'TickLabels', {'No data', 'Has data'});
title('Grid coverage at time index 10');
xlabel('Longitude'); ylabel('Latitude');

% Overlay stations
hold on;
scatter(KS.softdata.sMS(:,1), KS.softdata.sMS(:,2), 10, 'r.', 'filled');
legend('Grid coverage', 'Stations');
```

---

## Common Issues and Solutions

### Issue 1: "Memory Error - Cannot allocate array"

**Problem:** Grid too large for available memory

**Solutions:**
```matlab
% 1. Increase resolution (fewer points)
grid_data = stg_to_stug(KS.softdata, 0.2);  % Instead of 0.1

% 2. Reduce spatial extent
grid_data = stg_to_stug(KS.softdata, ...
    'resolution', 0.1, ...
    'bounds', [-100, -90, 25, 35]);  % Smaller region

% 3. Process in chunks (see example_stg_to_stug_usage.m)
```

### Issue 2: "Too Many NaN Values in Output Grid"

**Problem:** Grid points outside convex hull of stations

**Solutions:**
```matlab
% 1. Check station coverage
figure; scatter(KS.softdata.sMS(:,1), KS.softdata.sMS(:,2), '.');

% 2. Use tighter bounds (exclude sparse edges)
grid_data = stg_to_stug(KS.softdata, ...
    'bounds', [lon_min+0.5, lon_max-0.5, lat_min+0.5, lat_max-0.5]);

% 3. Enable extrapolation (use with caution!)
grid_data = stg_to_stug(KS.softdata, ...
    'extrapolation', 'linear');  % Fills beyond convex hull
```

### Issue 3: "Conversion Takes Too Long"

**Problem:** Too many stations × time points

**Solutions:**
```matlab
% 1. Increase resolution
grid_data = stg_to_stug(KS.softdata, 0.2);  % Faster than 0.1

% 2. Reduce time range
subset_data = KS.softdata;
subset_data.tME = KS.softdata.tME(1:10);  % First 10 times only
subset_data.Xms = KS.softdata.Xms(:, 1:10);
grid_data = stg_to_stug(subset_data, 0.1);

% 3. Use 'nearest' method (faster than 'linear')
grid_data = stg_to_stug(KS.softdata, ...
    'method', 'nearest');
```

---

## Comparison: When to Use STG vs STUG

| Factor | Use STG (neighbours_stg) | Convert to STUG |
|--------|-------------------------|-----------------|
| **Number of queries** | 1-50 queries | 100+ queries ✅ |
| **Accuracy needed** | Exact measurements ✅ | Interpolated OK |
| **Memory available** | Limited | Ample (2-5x data) |
| **Preprocessing time** | Not acceptable | Acceptable (one-time) |
| **Gap filling needed** | Not needed | Need estimates ✅ |
| **Station spacing** | Irregular | Nearly regular |
| **Processing speed** | 5-20 ms per query ✅ | 2-10 sec + 5ms |

**Bottom line:**
- For your typical use case: **Use STG directly** (faster overall)
- For many queries or gap filling: **Convert to STUG**

---

## Summary

**Efficient Conversion:**
```matlab
grid_data = stg_to_stug(KS.softdata, 0.1);
```

**Key Features:**
- ✅ Automatic resolution calculation
- ✅ Vectorized interpolation (10-100x faster than loops)
- ✅ Progress reporting
- ✅ Memory pre-allocation
- ✅ Quality metadata
- ✅ ~2-10 seconds for your data

**Recommended Workflow:**
1. Start with auto resolution
2. Check metadata and coverage
3. Adjust resolution if needed
4. Cache result if many queries
5. Use with `neighbours_stug_index.m`
