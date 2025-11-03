# Data Structure Analysis & Reformatting Strategy

## Current Data Structure Analysis

```matlab
Your data structure:
  sMS: [78048×2 double]           % Spatial monitoring sites (lat/lon pairs)
  tME: [1×19 double]              % Time values (19 time points)
  Xms: [78048×19 single]          % Data mean values (stations × time)
  Xvs: [78048×19 single]          % Data variance (stations × time)
  p: [1482912×3 double]           % Full space-time points [lon, lat, time]
  z: [1482912×1 single]           % Data values in vector format
  vs: [1482912×1 single]          % Variance in vector format
  Zisnotnan: [78048×19 logical]   % Valid data mask
  nanratio: 0                     % No NaN values
  index_stg_to_stv: [1482912×1 double] % Grid-to-vector mapping

Total points: 78048 stations × 19 time points = 1,482,912 space-time observations
```

---

## ⚠️ CRITICAL MISMATCH IDENTIFIED

### Your Data Format: **STG (Space-Time Grid)**
- **Spatial**: 78,048 **irregular** monitoring stations (not uniformly spaced)
- **Temporal**: 19 time points (regular)
- **Structure**: Stations × Time grid

### Function You're Considering: **neighbours_stug_optimized.m**
- **Expected format**: STUG (Space-Time **Uniformly Gridded**)
- **Requirements**:
  - Uniform spatial grid (constant dx, dy)
  - Regular grid vectors: x, y, time
  - 3D array: Z[nx, ny, nt]

### The Problem:
**Your spatial coordinates (`sMS`) are irregular station locations, NOT a uniform grid!**

Example of what this means:
```
Uniform Grid (STUG):          Irregular Stations (STG - YOUR DATA):
  x: [-180, -179, -178, ...]    sMS: [(-95.3, 29.7),
  y: [-90, -89, -88, ...]             (-94.2, 31.5),
  Regular spacing: dx=1°, dy=1°        (-96.8, 28.3),
                                       ...]
                                  Irregular locations!
```

---

## 🎯 RECOMMENDED SOLUTION: Use the Correct Function

### Option 1: **Use `neighbours_stg.m` (BEST CHOICE)** ✅

**Why this is correct:**
- ✅ Designed specifically for your data format (STG)
- ✅ Handles irregular spatial stations
- ✅ Already optimized for station × time structure
- ✅ **NO reformatting needed** - use data as-is!

**Your data already has the right format:**
```matlab
% Your current structure maps directly to neighbours_stg inputs:
data_stg.sMS = ans.sMS;          % [78048×2] spatial coordinates
data_stg.tME = ans.tME;          % [1×19] time values
data_stg.Xms = ans.Xms;          % [78048×19] data values
data_stg.Xvs = ans.Xvs;          % [78048×19] variances (optional)

% Call directly:
[psub, zsub, dsub, nsub] = neighbours_stg(p0, data_stg, nmax, dmax);
```

**Performance:**
- Optimized for this exact use case
- No data conversion overhead
- Handles irregular stations correctly
- Fast for moderate number of stations (78k is fine)

**Verdict:** ✅ **Use this - it's the right tool for your data**

---

## Option 2: Force into STUG Format (NOT RECOMMENDED) ❌

If you absolutely must use `neighbours_stug_optimized.m`, you'd need to:

### Approach A: Interpolate onto Uniform Grid

**Process:**
1. Determine bounding box of all stations
2. Create uniform grid (e.g., 0.1° × 0.1° resolution)
3. Interpolate irregular station data onto grid
4. Reformat into STUG structure

**Code sketch:**
```matlab
% 1. Determine spatial extent
lon_min = min(ans.sMS(:,1));
lon_max = max(ans.sMS(:,1));
lat_min = min(ans.sMS(:,2));
lat_max = max(ans.sMS(:,2));

% 2. Create uniform grid (choose resolution)
resolution = 0.1;  % degrees
grid_data.x = (lon_min:resolution:lon_max)';
grid_data.y = (lat_min:resolution:lat_max)';
grid_data.time = ans.tME';

nx = length(grid_data.x);
ny = length(grid_data.y);
nt = length(grid_data.time);

% 3. Create meshgrid
[X, Y] = meshgrid(grid_data.x, grid_data.y);
grid_data.Lon = X';
grid_data.Lat = Y';

% 4. Interpolate for each time step
grid_data.Z = NaN(nx, ny, nt);
for t = 1:nt
    % Get valid stations at this time
    valid_idx = ans.Zisnotnan(:, t);

    if sum(valid_idx) > 0
        % Interpolate using scatteredInterpolant
        F = scatteredInterpolant(ans.sMS(valid_idx, 1), ...
                                 ans.sMS(valid_idx, 2), ...
                                 double(ans.Xms(valid_idx, t)), ...
                                 'linear', 'none');

        % Evaluate on uniform grid
        for i = 1:nx
            for j = 1:ny
                grid_data.Z(i, j, t) = F(grid_data.x(i), grid_data.y(j));
            end
        end
    end
end

% Now can use neighbours_stug_optimized
[psub, zsub, dsub, nsub] = neighbours_stug_optimized(p0, grid_data, nmax, dmax);
```

**Problems with this approach:**
- ❌ **Information loss** - interpolation introduces error
- ❌ **Computational overhead** - interpolation is expensive
- ❌ **Memory overhead** - grid may have many more points than stations
- ❌ **Artificial data** - creates values where no measurements exist
- ❌ **Grid size issues** - too coarse loses detail, too fine wastes memory

**Example of the problem:**
```
Original: 78,048 stations (irregular locations)
Grid at 0.1° resolution over USA (~25°×50°):
  nx = 250, ny = 500, nt = 19
  Total grid points = 250 × 500 × 19 = 2,375,000

You'd create 2.4M grid points from 1.5M station observations!
And most grid cells may have no nearby stations.
```

---

## Option 3: Wrapper Function (MIDDLE GROUND)

Create a thin wrapper that reformats on-the-fly:

```matlab
function [psub, zsub, dsub, nsub] = neighbours_stg_wrapper(p0, stg_data, nmax, dmax)
    % Wrapper that accepts STG format and calls neighbours_stg
    %
    % This keeps your calling code simple while using the right function

    % Validate input format
    if ~isfield(stg_data, 'sMS') || ~isfield(stg_data, 'tME') || ~isfield(stg_data, 'Xms')
        error('Input must be in STG format with fields: sMS, tME, Xms');
    end

    % Call the appropriate function
    [psub, zsub, dsub, nsub] = neighbours_stg(p0, stg_data, nmax, dmax);
end
```

**Advantages:**
- ✅ Keeps your calling code unchanged
- ✅ Uses correct function internally
- ✅ No data reformatting overhead
- ✅ Clear interface

---

## 📊 Performance Comparison

| Approach | Preprocessing Time | Memory Usage | Accuracy | Recommended |
|----------|-------------------|--------------|----------|-------------|
| **Use neighbours_stg** | 0 ms (none) | Original data only | Perfect | ✅ **YES** |
| **Interpolate to STUG** | 5-30 seconds | 2-10x original | Interpolation error | ❌ NO |
| **Wrapper function** | 0 ms (none) | Original data only | Perfect | ✅ YES |

---

## 🎯 FINAL RECOMMENDATION

### **DO THIS:** ✅
```matlab
% Option 1: Use neighbours_stg directly (BEST)
[psub, zsub, dsub, nsub] = neighbours_stg(p0, ans, nmax, dmax);
```

### **DON'T DO THIS:** ❌
```matlab
% Don't force STG data into STUG format
% - Waste of computation
% - Loss of accuracy
% - Unnecessary complexity
```

---

## Why This Matters

### Your Data Characteristics:
```
78,048 stations × 19 time points = 1,482,912 observations
nanratio = 0 (all data valid)
Spatial: Irregular monitoring sites (likely air quality, weather, or environmental sensors)
Temporal: Regular time series (19 snapshots)
```

### Function Match:
```
neighbours_stg.m:
  ✅ Designed for: irregular stations × regular time
  ✅ Your data: irregular stations × regular time
  ✅ PERFECT MATCH!

neighbours_stug_optimized.m:
  ❌ Designed for: uniform grid × regular time
  ❌ Your data: irregular stations × regular time
  ❌ WRONG FUNCTION!
```

---

## If You Absolutely Must Use STUG (Not Recommended)

### When to Reformat:

**Reformat BEFORE calling function if:**
- ✅ You'll make **multiple queries** on the same dataset (amortize interpolation cost)
- ✅ Your stations are **nearly regular** (interpolation error is minimal)
- ✅ You need **gap filling** (want to estimate values between stations)

**DON'T reformat if:**
- ❌ Single or few queries (overhead not worth it)
- ❌ Stations are highly irregular (large interpolation errors)
- ❌ You have the correct function available (`neighbours_stg`)

### Reformatting Strategy (if you insist):

```matlab
function grid_data = stg_to_stug(stg_data, resolution)
    % Convert STG to STUG format via interpolation
    %
    % INPUT:
    %   stg_data: struct with sMS, tME, Xms
    %   resolution: grid spacing in degrees (e.g., 0.1)
    %
    % OUTPUT:
    %   grid_data: struct with x, y, time, Lon, Lat, Z (STUG format)

    % 1. Determine spatial extent
    lon_range = [min(stg_data.sMS(:,1)), max(stg_data.sMS(:,1))];
    lat_range = [min(stg_data.sMS(:,2)), max(stg_data.sMS(:,2))];

    % 2. Create uniform grid
    grid_data.x = (lon_range(1):resolution:lon_range(2))';
    grid_data.y = (lat_range(1):resolution:lat_range(2))';
    grid_data.time = stg_data.tME';

    nx = length(grid_data.x);
    ny = length(grid_data.y);
    nt = length(grid_data.time);

    % 3. Create coordinate arrays
    [X, Y] = meshgrid(grid_data.x, grid_data.y);
    grid_data.Lon = X';
    grid_data.Lat = Y';

    % 4. Interpolate data
    grid_data.Z = NaN(nx, ny, nt, 'single');

    fprintf('Interpolating %d time steps onto %dx%d grid...\n', nt, nx, ny);

    for t = 1:nt
        % Get valid data at this time
        if isfield(stg_data, 'Zisnotnan')
            valid_idx = stg_data.Zisnotnan(:, t);
        else
            valid_idx = ~isnan(stg_data.Xms(:, t));
        end

        n_valid = sum(valid_idx);
        if n_valid > 3  % Need at least 3 points for interpolation
            % Create interpolant
            F = scatteredInterpolant(...
                stg_data.sMS(valid_idx, 1), ...
                stg_data.sMS(valid_idx, 2), ...
                double(stg_data.Xms(valid_idx, t)), ...
                'linear', 'none');  % 'none' = NaN outside convex hull

            % Interpolate
            grid_data.Z(:, :, t) = reshape(...
                F(grid_data.Lon(:), grid_data.Lat(:)), ...
                nx, ny);
        end

        if mod(t, 5) == 0
            fprintf('  Completed %d/%d time steps\n', t, nt);
        end
    end

    fprintf('Done. Grid size: %dx%dx%d = %d cells\n', ...
            nx, ny, nt, nx*ny*nt);
end

% Usage:
grid_data = stg_to_stug(ans, 0.1);  % 0.1 degree resolution
[psub, zsub, dsub, nsub] = neighbours_stug_optimized(p0, grid_data, nmax, dmax);
```

**Estimated performance:**
```
For your data (78k stations, 19 times):
  Interpolation time: 10-30 seconds (one-time cost)
  Memory: ~2-5x your current data size
  Accuracy: Depends on station density and spatial patterns

neighbor_stg (correct function):
  Call time: 5-20 milliseconds per query
  Memory: Uses original data (no overhead)
  Accuracy: Exact (uses actual station measurements)
```

---

## Summary & Decision Matrix

| Scenario | Recommendation | Rationale |
|----------|---------------|-----------|
| **You have irregular stations** | ✅ Use `neighbours_stg.m` | Designed for this format |
| **You need exact station values** | ✅ Use `neighbours_stg.m` | No interpolation error |
| **You need speed** | ✅ Use `neighbours_stg.m` | No preprocessing |
| **Single/few queries** | ✅ Use `neighbours_stg.m` | Not worth reformatting |
| **You want to fill gaps** | ⚠️ Interpolate to STUG | But understand tradeoffs |
| **Many queries (100+)** | ⚠️ Maybe interpolate | If stations nearly regular |
| **You insist on using STUG** | ❌ Interpolate BEFORE | Amortize one-time cost |

---

## Conclusion

### **ANSWER: Neither - Use the Right Function!**

**Don't reformat at all. Use `neighbours_stg.m` which is designed for your data format.**

Your data is already in the optimal format (STG). The only reformatting needed is minimal field mapping:

```matlab
% Minimal reformatting (just structure field names if needed):
data_stg.sMS = ans.sMS;
data_stg.tME = ans.tME;
data_stg.Xms = ans.Xms;

% Call the correct function:
[psub, zsub, dsub, nsub] = neighbours_stg(p0, data_stg, nmax, dmax);
```

**If you absolutely must use STUG:**
- ⚠️ Reformat **BEFORE** calling the function (one-time interpolation cost)
- But understand this introduces error and computational overhead
- Only worth it if you have 100+ queries on the same dataset

**Bottom line:** The question assumes you should use `neighbours_stug_optimized.m`, but the real answer is you shouldn't - use `neighbours_stg.m` instead!
