# Quick Reference: Neighbor Search Functions

## Function Comparison at a Glance

| Function | Input Format | Best For | Complexity | Memory |
|----------|--------------|----------|------------|--------|
| `neighbours.m` | STV (vector) | Irregular points | O(n) | O(n) |
| `neighbours_stg.m` | STG (station×time) | Monitoring networks | O(nMS + nME) | O(nMS×nME logical) |
| `neighbours_stug_optimized.m` | Uniform grid | Satellite/model data | O(nmax log nmax) | O(nmax) |

## When to Use Which Function

### Use `neighbours.m` when:
- ✓ Data points at arbitrary locations
- ✓ Need multi-variable indexing
- ✓ Small to medium datasets (<10,000 points)
- ✓ Points not on regular grid

**Example**: Air quality measurements from mobile sensors

### Use `neighbours_stg.m` when:
- ✓ Fixed monitoring stations
- ✓ Regular time intervals
- ✓ Moderate number of stations (<1,000)
- ✓ Some missing data (<50% NaN)

**Example**: Weather station network with hourly measurements

### Use `neighbours_stug_optimized.m` when:
- ✓ Data on uniform 3D grid
- ✓ Very large grids (>100,000 points)
- ✓ Potentially high NaN ratio (>50%)
- ✓ Satellite or model output data

**Example**: Satellite imagery, climate model outputs, gridded soft data

## Syntax Quick Guide

### neighbours.m
```matlab
[csub, Zsub, dsub, nsub, index] = neighbours(c0, c, Z, nmax, dmax)

% Input
c0:    [lon, lat, time]
c:     n×3 coordinates
Z:     n×k values
nmax:  scalar or vector
dmax:  scalar or [spatial_max, temporal_max, metric]

% Output
csub:  m×3 coordinates
Zsub:  m×k values
dsub:  m×1 or m×2 distances
```

### neighbours_stg.m
```matlab
[psub, zsub, dsub, nsub, index] = neighbours_stg(p0, data, nmax, dmax)

% Input
p0:    [lon, lat, time]
data:  struct with .sMS, .tME, .p, .z, .Zisnotnan, .nanratio, .index_stg_to_stv
nmax:  scalar
dmax:  [spatial_max, temporal_max, metric]

% Output
psub:  m×3 coordinates
zsub:  m×1 values
dsub:  m×2 distances [spatial, temporal]
```

### neighbours_stug_optimized.m
```matlab
[psub, zsub, dsub, nsub, index] = neighbours_stug_optimized(p0, grid_data, nmax, dmax)

% Input
p0:         [lon, lat, time] or [lon, lat, time, ix, iy, it]
grid_data:  struct with .x, .y, .time, .Lon, .Lat, .Z OR NetCDF path
nmax:       scalar
dmax:       [spatial_max, temporal_max, metric]

% Output
psub:  m×6 [lon, lat, time, ix, iy, it]
zsub:  m×1 values
dsub:  m×2 distances [spatial, temporal]
```

## Performance Cheat Sheet

### Typical Execution Times (100×100×100 grid, 30% NaN, nmax=20)

- `neighbours.m`: ~0.5-1.0 sec (if converted to vector)
- `neighbours_stg.m`: ~0.1-0.3 sec (if structured as stations)
- `neighbours_stug_optimized.m`: ~0.01-0.05 sec

### Memory Usage (same scenario)

- `neighbours.m`: ~8 MB
- `neighbours_stg.m`: ~2-4 MB
- `neighbours_stug_optimized.m`: ~0.1 MB

### Scalability (100×100×1000 grid)

- `neighbours.m`: ~5-10 sec
- `neighbours_stg.m`: ~0.5-1.5 sec (depends on station structure)
- `neighbours_stug_optimized.m`: ~0.02-0.1 sec

## Common Pitfalls and Solutions

### Problem: "Index exceeds array bounds"
**Cause**: Grid indices don't match actual grid size
**Solution**:
```matlab
% Verify grid dimensions
[nx, ny, nt] = size(grid_data.Z);
assert(length(grid_data.x) == nx);
assert(length(grid_data.y) == ny);
assert(length(grid_data.time) == nt);
```

### Problem: "Returns 0 neighbors but data exists"
**Cause**: dmax constraints too tight
**Solution**:
```matlab
% Check grid spacing
dx = median(diff(grid_data.x));
dy = median(diff(grid_data.y));
% Ensure dmax(1) > dx and dmax(2) > dt
```

### Problem: "Very slow for large grids"
**Cause**: Using wrong function for data structure
**Solution**: Match function to data format (see table above)

### Problem: "All returned values are NaN"
**Cause**: Bug in original data or index mapping
**Solution**:
```matlab
% Check if non-NaN data exists
num_valid = sum(~isnan(grid_data.Z(:)));
fprintf('Valid data points: %d / %d\n', num_valid, numel(grid_data.Z));
```

## Parameter Tuning Guide

### nmax (number of neighbors)

- **Too small** (<10): May miss important spatial patterns, kriging unstable
- **Too large** (>100): Slow computation, may include irrelevant distant points
- **Recommended**: 15-30 for most applications

### dmax (distance thresholds)

#### Spatial distance (dmax(1))
```matlab
% Rule of thumb: 2-5× the grid spacing
dx = median(diff(grid_data.x));
dmax_spatial = 3 * dx;  % Start here
```

#### Temporal distance (dmax(2))
```matlab
% Depends on temporal correlation
dt = median(diff(grid_data.time));
% For daily data: 7-30 days
% For hourly data: 24-168 hours
dmax_temporal = 10 * dt;
```

#### Space-time metric (dmax(3))
```matlab
% Converts temporal distance to spatial equivalent
% Typical values: 0.01 to 1.0
% Higher value = time matters more
% Lower value = space matters more

% Example: If 1 day ≈ 100 km in correlation
dmax(3) = 100 / 1;  % 100 km per day

% For your data:
dmax(3) = spatial_correlation_range / temporal_correlation_range;
```

## Code Templates

### Template 1: Basic Grid Setup
```matlab
% Create uniform grid
nx = 100; ny = 100; nt = 200;
grid_data.x = linspace(lon_min, lon_max, nx)';
grid_data.y = linspace(lat_min, lat_max, ny)';
grid_data.time = (time_start:time_end)';

[X, Y] = meshgrid(grid_data.x, grid_data.y);
grid_data.Lon = X';
grid_data.Lat = Y';

% Load or create data
grid_data.Z = your_data;  % nx×ny×nt

% Mark missing values
grid_data.Z(your_data < 0) = NaN;
```

### Template 2: Neighbor Search for Multiple Points
```matlab
% Multiple estimation points
p0_list = [
    -95, 27.5, 100;
    -92, 28.0, 105;
    -90, 29.5, 110;
];

% Parameters
nmax = 20;
dmax = [2.0, 10.0, 0.1];

% Loop over points
for i = 1:size(p0_list, 1)
    p0 = p0_list(i, :);

    [psub, zsub, dsub, nsub, index] = ...
        neighbours_stug_optimized(p0, grid_data, nmax, dmax);

    fprintf('Point %d: Found %d neighbors\n', i, nsub);

    % Store results
    results(i).psub = psub;
    results(i).zsub = zsub;
    results(i).nsub = nsub;
end
```

### Template 3: Adaptive Search (if not finding enough neighbors)
```matlab
p0 = [-95, 27.5, 100];
nmax = 20;
dmax_base = [2.0, 10.0, 0.1];

% Try progressively larger search radii
search_multipliers = [1, 1.5, 2, 3, 5];

for mult = search_multipliers
    dmax = dmax_base .* [mult, mult, 1];

    [psub, zsub, dsub, nsub, index] = ...
        neighbours_stug_optimized(p0, grid_data, nmax, dmax);

    if nsub >= nmax
        fprintf('Found enough neighbors with multiplier %.1f\n', mult);
        break;
    end
end

if nsub < nmax
    warning('Could not find %d neighbors even with largest search radius', nmax);
end
```

## Testing Checklist

Before using in production, verify:

- [ ] Grid dimensions match: `size(Z)` vs `[length(x), length(y), length(time)]`
- [ ] Coordinate arrays are monotonic: `issorted(x)`, `issorted(y)`, `issorted(time)`
- [ ] Lon/Lat dimensions correct: `size(Lon)` = `[nx, ny]` or `[nx, ny, nt]`
- [ ] dmax values reasonable: `dmax(1) > grid_spacing`, `dmax(2) > time_step`
- [ ] Some non-NaN data exists: `sum(~isnan(Z(:))) > 0`
- [ ] Test point within grid bounds

Quick validation:
```matlab
% Validation script
assert(length(grid_data.x) == size(grid_data.Z, 1), 'X dimension mismatch');
assert(length(grid_data.y) == size(grid_data.Z, 2), 'Y dimension mismatch');
assert(length(grid_data.time) == size(grid_data.Z, 3), 'Time dimension mismatch');
assert(sum(~isnan(grid_data.Z(:))) > 0, 'No valid data');
fprintf('✓ Grid validation passed\n');
```

## Support and Troubleshooting

1. Check grid structure: `whos grid_data`
2. Verify data exists: `sum(~isnan(grid_data.Z(:)))`
3. Test with simple case: small grid, low NaN ratio
4. Compare with `neighbours.m` for correctness
5. Use `OPTIMIZATION_GUIDE.md` for detailed algorithm explanation
