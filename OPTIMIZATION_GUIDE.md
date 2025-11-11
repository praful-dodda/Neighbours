# Optimization Guide for neighbours_stug.m

## Overview

This document describes the optimization strategy for neighbor search in uniformly gridded space-time data, particularly for soft data which can have very large grids with many missing (NaN) values.

## Problem Statement

### Original Implementation Issues (`neighbours_stug.m`)

The original implementation had several inefficiencies:

1. **Arbitrary search radius**: Used `ceil(sqrt(nmax))` - no mathematical justification
2. **Full ndgrid creation**: Created complete coordinate arrays for entire search region
3. **Brute-force distance computation**: Computed distances for ALL points before filtering
4. **Memory intensive**: For large grids, created massive temporary arrays
5. **Ignored uniform spacing**: Didn't exploit the fact that grid spacing is constant

### Example of Inefficiency

For a 100×100×1000 grid with 30% NaN:
- Original: Creates arrays of size ~1000+ points, computes all distances
- Memory: ~8-16 MB temporary arrays per query
- Wasted computation: >90% of distance calculations discarded

## Optimization Strategy

### Core Insights

1. **Uniform Grid Properties**:
   - Grid spacing (dx, dy, dt) is constant
   - Can convert distances to grid indices directly
   - Integer arithmetic much faster than floating-point distance calculations

2. **Spatial Locality**:
   - Neighbors are typically close in grid-index space
   - Can use grid indices as distance proxy
   - Only compute actual distances for final candidates

3. **Early Termination**:
   - Stop searching when enough neighbors found
   - Use shells of increasing radius
   - Skip regions that can't contain valid neighbors

### Implementation: Hybrid Approach

#### Strategy A: Shell-Based Expansion (for large grids)

**When to use**: Candidate volume > 1000 points OR nmax << candidate volume

**Algorithm**:
```
1. Start at target grid point (idx_x, idx_y, idx_t)
2. Expand outward in "shells" of increasing Chebyshev distance
3. For each shell:
   a. Generate grid points at that distance
   b. Check if within bounds and dmax constraints
   c. Check if data exists (not NaN)
   d. Add valid points to candidates
4. Stop when:
   - Have >= 2×nmax candidates (for sorting), OR
   - Exceeded maximum possible search radius
5. Compute actual distances for candidates
6. Sort and return top nmax
```

**Benefits**:
- Memory: O(nmax) instead of O(search_volume)
- Early termination possible
- Cache-friendly (processes nearby points together)

**Complexity**:
- Best case: O(nmax log nmax) when data is dense
- Worst case: O(R³) where R is search radius in grid units

#### Strategy B: Index-Based Search (for small grids)

**When to use**: Candidate volume ≤ 1000 points

**Algorithm**:
```
1. Convert dmax to grid units:
   rx_max = ceil(dmax(1) / dx)
   ry_max = ceil(dmax(1) / dy)
   rt_max = ceil(dmax(2) / dt)

2. For each time slice in [idx_t - rt_max, idx_t + rt_max]:
   a. Compute temporal distance
   b. Calculate remaining spatial budget
   c. Determine effective spatial search radius
   d. For each (ix, iy) in spatial range:
      - Quick grid-distance check
      - If valid and not NaN, add to candidates

3. Compute actual distances for candidates
4. Sort and return top nmax
```

**Benefits**:
- Simple, readable code
- No complex shell generation
- Good for small search regions

**Complexity**: O(rx × ry × rt) where all are small

#### Decision Logic

```matlab
candidate_volume = (2*rx_max + 1) * (2*ry_max + 1) * (2*rt_max + 1)
use_shell_expansion = (candidate_volume > 1000) || (nmax < 0.1 * candidate_volume)
```

## Performance Comparison

### Theoretical Analysis

| Aspect | Original | Optimized |
|--------|----------|-----------|
| Memory | O(rx × ry × rt) | O(nmax) |
| Distance calcs | All points in box | Only candidates |
| Grid indices | Created via ndgrid | Direct calculation |
| Early exit | No | Yes (shell method) |

### Expected Speedup

**Scenario 1: Large sparse grid (100×100×1000, 70% NaN, nmax=20)**
- Original: ~0.5-1.0 sec
- Optimized: ~0.01-0.05 sec
- **Speedup: 10-50×**

**Scenario 2: Medium grid (50×50×100, 30% NaN, nmax=30)**
- Original: ~0.1-0.2 sec
- Optimized: ~0.005-0.02 sec
- **Speedup: 5-20×**

**Scenario 3: Small dense grid (30×30×30, 10% NaN, nmax=50)**
- Original: ~0.01-0.05 sec
- Optimized: ~0.005-0.01 sec
- **Speedup: 2-5×**

## Key Differences from neighbours_stg.m

| Feature | neighbours_stg.m | neighbours_stug_optimized.m |
|---------|------------------|------------------------------|
| **Grid Type** | Irregular spatial stations + regular time | Fully uniform 3D grid |
| **Search Strategy** | Filter MS→ME→combine | Coupled 3D grid search |
| **Distance Calculation** | All MS, all ME | Only final candidates |
| **Memory Model** | O(nMS + nME + logical grid) | O(nmax) |
| **Best Use Case** | Monitoring station networks | Regular satellite/model grids |
| **Handles NaN** | Pre-computed Zisnotnan | On-the-fly checks |
| **Spatial Structure** | Arbitrary station locations | Uniform spacing exploited |

### When to Use Each

- **neighbours.m**: Arbitrary point clouds, multi-variable indexing, irregular spacing
- **neighbours_stg.m**: Monitoring stations (irregular space) + regular time, moderate NaN
- **neighbours_stug_optimized.m**: Fully uniform grids (satellite, model data), potentially huge with high NaN ratio

## Implementation Details

### Grid Spacing Calculation

```matlab
dx = median(diff(data.x))  % Robust to slight irregularities
dy = median(diff(data.y))
dt = median(diff(data.time))
```

Uses `median` instead of `mean` to be robust to edge irregularities.

### Space-Time Distance Metric

Following the convention from other functions:
```
spacetime_distance = spatial_distance + dmax(3) × temporal_distance
```

This allows different weighting of spatial vs temporal proximity.

### Coordinate System Handling

The function handles two coordinate systems:
1. **Physical coordinates**: (Lon, Lat, time) - what user provides
2. **Grid indices**: (ix, iy, it) - used for efficient search

Conversion happens at boundaries:
- Input: Physical coords → find grid indices
- Search: Works in grid index space
- Output: Returns both physical coords AND grid indices

### NaN Handling

Direct check during search:
```matlab
if ~isnan(Z(ix, iy, it))
    % Valid data point
end
```

No pre-computation needed (unlike neighbours_stg.m which uses Zisnotnan matrix).

## Input/Output Specification

### Input

```matlab
p0: [lon, lat, time] or [lon, lat, time, idx_x, idx_y, idx_t]
    If indices provided, uses them directly (faster)
    Otherwise computes nearest grid point

grid_data: struct with fields:
    .x      nx×1    x-coordinates (can be longitude or grid units)
    .y      ny×1    y-coordinates (can be latitude or grid units)
    .time   nt×1    time coordinates
    .Lon    nx×ny   longitude at each grid point (or nx×ny×nt)
    .Lat    nx×ny   latitude at each grid point (or nx×ny×nt)
    .Z      nx×ny×nt   data values (NaN for missing)

    OR: path to NetCDF file with these variables

nmax: scalar - maximum neighbors to return

dmax: [spatial_max, temporal_max, metric]
    spatial_max: maximum spatial distance
    temporal_max: maximum temporal distance
    metric: space-time coupling factor
```

### Output

```matlab
psub:  m×6  [lon, lat, time, idx_x, idx_y, idx_t]
zsub:  m×1  data values
dsub:  m×2  [spatial_distance, temporal_distance]
nsub:  scalar  number of neighbors (m ≤ nmax)
index: m×1  linear indices into Z array
```

All outputs sorted by increasing space-time distance.

## Testing Strategy

### Correctness Tests

1. **Basic functionality**: Find neighbors, verify distances
2. **Edge cases**: Corners, boundaries, grid edges
3. **NaN handling**: Various ratios (0%, 30%, 70%, 95%)
4. **Distance constraints**: Verify dmax enforcement
5. **Sorting**: Check space-time distance ordering
6. **Empty data**: All NaN, very sparse grids

### Performance Tests

1. **Scalability**: Grid sizes 30³ to 150³
2. **NaN ratio impact**: 0% to 95% missing
3. **nmax variation**: 5 to 200 neighbors
4. **Search radius**: Small to large dmax values

### Validation Against Reference

Compare with `neighbours.m` (general implementation):
- Same inputs should give same neighbors (possibly different order)
- Distances should match
- All constraints satisfied

## Usage Examples

### Example 1: Basic Usage

```matlab
% Setup grid
grid_data.x = linspace(-100, -90, 100)';
grid_data.y = linspace(25, 30, 100)';
grid_data.time = (1:200)';

[X, Y] = meshgrid(grid_data.x, grid_data.y);
grid_data.Lon = X';
grid_data.Lat = Y';
grid_data.Z = randn(100, 100, 200);
grid_data.Z(rand(100,100,200) < 0.3) = NaN;  % 30% missing

% Query
p0 = [-95, 27.5, 100];  % [lon, lat, time]
nmax = 20;
dmax = [2.0, 10.0, 0.1];

[psub, zsub, dsub, nsub, index] = ...
    neighbours_stug_optimized(p0, grid_data, nmax, dmax);

fprintf('Found %d neighbors\n', nsub);
fprintf('Closest neighbor at distance: %.3f spatial, %.3f temporal\n', ...
    dsub(1,1), dsub(1,2));
```

### Example 2: With NetCDF File

```matlab
p0 = [-90.1, 27.1, 2404];
nmax = 25;
dmax = [3.0, 15.0, 0.1];

[psub, zsub, dsub, nsub, index] = ...
    neighbours_stug_optimized(p0, 'data.nc', nmax, dmax);
```

### Example 3: Pre-computed Indices

```matlab
% If you know the grid indices, provide them for speed
idx_x = 50;
idx_y = 45;
idx_t = 100;
p0 = [-95, 27.5, 100, idx_x, idx_y, idx_t];

[psub, zsub, dsub, nsub, index] = ...
    neighbours_stug_optimized(p0, grid_data, nmax, dmax);
```

## Maintenance and Extensions

### Potential Enhancements

1. **Adaptive shell expansion**: Vary shell thickness based on NaN density
2. **Spatial indexing**: Pre-build KD-tree for very large grids with repeated queries
3. **Parallel processing**: Parallelize shell generation for very large grids
4. **Anisotropic grids**: Handle dx ≠ dy with elliptical search regions
5. **Variable metrics**: Support different distance metrics (haversine for lat/lon)

### Code Quality

- **Modularity**: Separate helper functions for each strategy
- **Documentation**: Extensive comments explaining algorithm
- **Error handling**: Validates inputs, handles edge cases
- **Flexibility**: Supports multiple input formats

## References

### Related Functions

- `neighbours.m`: General neighbor search for arbitrary point clouds
- `neighbours_stg.m`: Optimized for station×time grids
- `coord2dist.m`: Euclidean distance matrix computation

### Key Concepts

- **Space-time kriging**: Why neighbor selection matters
- **Soft data**: Interval or probabilistic data, often dense in space-time
- **BME (Bayesian Maximum Entropy)**: Framework using these functions
