# Testing Guide

This guide explains how to run and interpret the various test scripts for the neighbor search implementations.

## Quick Start

To validate correctness and measure performance improvements:

```matlab
% 1. Test correctness
test_neighbours_stug_optimized

% 2. Compare with reference implementation
compare_neighbours_implementations

% 3. Detailed speed test
speed_test_comparison

% 4. (Optional) Generate visualizations
visualize_performance
```

## Test Scripts Overview

### 1. test_neighbours_stug_optimized.m

**Purpose:** Validate the optimized implementation in isolation

**What it tests:**
- Basic neighbor finding functionality
- Edge cases (corners, boundaries, grid edges)
- Various NaN ratios (0%, 10%, 30%, 50%, 70%, 90%, 95%)
- Distance constraint enforcement (dmax)
- Different nmax values
- Proper sorting by space-time distance
- Empty/sparse data handling
- Input format flexibility (3 or 6 element p0)

**Expected output:**
```
=== Testing neighbours_stug_optimized ===

Test 1: Basic uniform grid with synthetic data...
  Optimized version found 20 neighbors in 0.0234 seconds
  Max spatial distance: 1.9876 (limit: 2.0000)
  Max temporal distance: 4.8912 (limit: 5.0000)
  ✓ Distance constraints satisfied
  ✓ No NaN values in results
  ✓ Number of neighbors <= nmax
  Test 1 PASSED

Test 2: Edge cases (corners and boundaries)...
  Point 1: Found 8 neighbors
  Point 2: Found 7 neighbors
  Point 3: Found 12 neighbors
  Point 4: Found 11 neighbors
  Test 2 PASSED

...

=== All tests completed successfully! ===
```

**When to run:** After any code changes, before committing

---

### 2. compare_neighbours_implementations.m

**Purpose:** Validate that optimized version produces the same results as neighbours.m

**What it compares:**
- Number of neighbors found
- Neighbor coordinates (spatial and temporal)
- Distance values
- Data values at neighbors
- Execution time

**Test scenarios:**
1. Small Dense Grid (30³, 10% NaN)
2. Medium Grid 30% NaN (50³)
3. Medium Grid 70% NaN (50³)
4. Large Grid (100³, 30% NaN)
5. Very Sparse Grid (80³, 95% NaN)

**Expected output:**
```
=== Comparison: neighbours_stug_optimized vs neighbours.m ===

Test 1/5: Small Dense Grid
  Grid: 30x30x30, NaN ratio: 10.0%, nmax: 15
  Grid created: 27000 total points, 24300 non-NaN values
  Running neighbours.m (reference)...
    Time: 0.1245 sec, Found: 15 neighbors
  Running neighbours_stug_optimized.m...
    Time: 0.0087 sec, Found: 15 neighbors
  Comparing results...
    ✓ Same number of neighbors: 15
    Neighbor set match: 100.0%
    ✓ Neighbor sets match well
    ✓ Distances match for common neighbors
    ✓ Values match for common neighbors
  Performance:
    neighbours.m:           0.1245 sec
    neighbours_stug_opt:    0.0087 sec
    Speedup:                14.31x
    ✓ Optimized version is faster

...

=== SUMMARY ===

Performance Summary:
Test                           |  Ref(s) |  Opt(s) |  Speedup |  Match%
--------------------------------------------------------------------------------
Small Dense Grid               |   0.1245 |   0.0087 |   14.31x |   100.0%
Medium Grid 30% NaN            |   0.3421 |   0.0156 |   21.93x |   100.0%
Medium Grid 70% NaN            |   0.2987 |   0.0098 |   30.48x |   100.0%
Large Grid                     |   1.2341 |   0.0423 |   29.17x |   100.0%
Very Sparse Grid               |   0.1876 |   0.0045 |   41.69x |   100.0%

Overall Statistics:
  Average speedup:  27.52x
  Median speedup:   29.17x
  Min speedup:      14.31x
  Max speedup:      41.69x
  Average match:    100.0%

✓ CORRECTNESS: Both functions produce equivalent results
✓ PERFORMANCE: Optimized version is 27.52x faster on average
```

**Interpreting results:**
- **Match % ≥ 95%**: Excellent, both find essentially the same neighbors
- **Match % 80-95%**: Good, minor differences (possibly due to ties in distance)
- **Match % < 80%**: Investigate - may indicate a bug

**When to run:** After implementing algorithm changes, before releasing

---

### 3. speed_test_comparison.m

**Purpose:** Detailed performance benchmarking with reliable statistics

**What it measures:**
- Average execution time over multiple iterations
- Standard deviation of timing
- Speedup factor
- Time saved
- Processing throughput (points/second)

**Test scenarios:**
1. Small grid, low NaN (40³, 10%)
2. Medium grid, 30% NaN (50³)
3. Medium grid, 50% NaN (70³)
4. Large grid, 30% NaN (100³)
5. Large grid, 70% NaN (100³)
6. Very large grid (150×150×100)

**Each scenario runs 5 iterations** for statistical reliability

**Expected output:**
```
=== Speed Test: neighbours_stug_optimized vs neighbours.m ===

Running 5 iterations per scenario for reliable timing...

Scenario 1/6: Small grid, low NaN
  Grid: 40x40x40 = 64000 points
  NaN ratio: 10%, nmax: 15
  Non-NaN points: 57600
  Benchmarking neighbours.m...
    Average: 0.0876 ± 0.0034 sec (15 neighbors)
  Benchmarking neighbours_stug_optimized.m...
    Average: 0.0045 ± 0.0002 sec (15 neighbors)
  Speedup: 19.47x ⚡ (very good)
  Time saved: 0.0831 sec (94.9% reduction)

...

=== PERFORMANCE SUMMARY ===

Scenario                       |  Grid Pts |  Valid Pts |   Ref Time |   Opt Time |  Speedup
------------------------------------------------------------------------------------------------
Small grid, low NaN            |     64000 |      57600 |  0.0876s   |  0.0045s   |   19.47x
Medium grid, 30% NaN           |    125000 |      87500 |  0.1567s   |  0.0089s   |   17.60x
Medium grid, 50% NaN           |    343000 |     171500 |  0.3214s   |  0.0134s   |   23.99x
Large grid, 30% NaN            |   1000000 |     700000 |  1.2341s   |  0.0456s   |   27.06x
Large grid, 70% NaN            |   1000000 |     300000 |  0.9876s   |  0.0234s   |   42.21x
Very large grid                |   2250000 |    1125000 |  3.4567s   |  0.0876s   |   39.46x

Speedup Statistics:
  Mean:     28.30x
  Median:   25.53x
  Min:      17.60x
  Max:      42.21x
  Std Dev:  9.67

...

Results saved to: speed_test_results.mat
```

**Interpreting results:**
- **Speedup > 10x**: Excellent optimization
- **Speedup 5-10x**: Very good
- **Speedup 2-5x**: Good
- **Speedup 1-2x**: Modest improvement
- **Speedup < 1x**: Slower (investigate!)

**When to run:**
- Benchmark new hardware
- After significant algorithm changes
- For performance documentation
- Before/after MATLAB version upgrades

---

### 4. visualize_performance.m

**Purpose:** Create publication-quality performance plots

**Prerequisites:** Must run `speed_test_comparison.m` first (generates `speed_test_results.mat`)

**Generates 4 figures:**

1. **performance_comparison.png**
   - Bar chart of execution times
   - Speedup comparison
   - Scalability with grid size (log scale)

2. **speedup_analysis.png**
   - Speedup vs grid size scatter plot
   - Speedup vs NaN ratio scatter plot
   - Absolute time savings
   - Relative time savings (%)

3. **performance_metrics.png**
   - Time breakdown (optimized vs saved)
   - Processing throughput (points/sec)

4. **performance_summary.png**
   - Text summary with key statistics

**Expected output:**
```
Creating performance visualizations...
  Saved: performance_comparison.png
  Saved: speedup_analysis.png
  Saved: performance_metrics.png
  Saved: performance_summary.png

Visualization complete! Generated 4 figures:
  1. performance_comparison.png - Overall comparison
  2. speedup_analysis.png - Detailed speedup metrics
  3. performance_metrics.png - Throughput analysis
  4. performance_summary.png - Summary statistics
```

**When to run:**
- For documentation
- For presentations
- For publication figures
- After completing speed tests

---

### 5. benchmark_neighbours_stug.m

**Purpose:** Quick performance check (optimized version only, no comparison)

**What it tests:**
- Medium grid baseline (100³)
- Scalability (30³ to 150³)
- Effect of NaN ratio (0% to 95%)
- Effect of nmax (5 to 200)
- Effect of dmax (search radius)

**When to run:**
- Quick sanity check
- Parameter tuning
- When neighbours.m not available

---

## Troubleshooting

### Problem: "File not found" errors

**Solution:**
```matlab
% Ensure you're in the Neighbours directory
cd /path/to/Neighbours

% Or add to path
addpath('/path/to/Neighbours')
```

### Problem: Tests fail with "Index exceeds array bounds"

**Possible causes:**
- Grid dimensions mismatch
- Invalid p0 coordinates

**Solution:** Check grid data structure:
```matlab
whos grid_data
assert(length(grid_data.x) == size(grid_data.Z, 1))
assert(length(grid_data.y) == size(grid_data.Z, 2))
assert(length(grid_data.time) == size(grid_data.Z, 3))
```

### Problem: Match percentage < 100% in comparison

**Possible causes:**
- Numerical precision differences
- Different handling of ties (same distance)
- Small coordinate rounding errors

**If match > 95%:** Usually acceptable, minor differences
**If match < 95%:** Investigate specific cases

**Debug:**
```matlab
% Add detailed output in compare_neighbours_implementations.m
% Check the detailed_comparison structure
```

### Problem: Speedup < 1 (slower!)

**Possible causes:**
- Very small grids (overhead dominates)
- Dense data (no NaN, no early termination benefit)
- Old MATLAB version

**Debug:**
```matlab
% Profile the code
profile on
speed_test_comparison
profile viewer
```

### Problem: Out of memory

**Solutions:**
1. Reduce grid size in tests
2. Reduce number of iterations
3. Close other applications
4. Use smaller test scenarios

---

## Recommended Testing Workflow

### Before Committing Code Changes

```matlab
% 1. Quick correctness check
test_neighbours_stug_optimized

% 2. Compare with reference
compare_neighbours_implementations

% 3. If all pass, commit
```

### For Performance Documentation

```matlab
% 1. Full speed test
speed_test_comparison

% 2. Generate plots
visualize_performance

% 3. Include PNG files in documentation
```

### For Debugging

```matlab
% 1. Enable verbose output
dbstop if error

% 2. Run failing test
test_neighbours_stug_optimized

% 3. Inspect variables at breakpoint
```

---

## Expected Timing

On typical modern hardware (2020+):

| Test Script | Expected Runtime |
|-------------|------------------|
| test_neighbours_stug_optimized | ~5-10 seconds |
| compare_neighbours_implementations | ~30-60 seconds |
| speed_test_comparison | ~2-5 minutes |
| visualize_performance | ~10-20 seconds |
| benchmark_neighbours_stug | ~1-3 minutes |

**Total testing time:** ~5-10 minutes for complete validation

---

## Interpreting Performance Results

### What counts as "good" performance?

**For uniform grids (neighbours_stug_optimized):**
- Small grids (< 50³): 5-15× speedup
- Medium grids (50-100³): 15-30× speedup
- Large grids (> 100³): 25-50× speedup
- Very sparse (> 70% NaN): 30-100× speedup

### When is speedup lower?

- **Small grids:** Fixed overhead is larger fraction
- **Dense data (low NaN ratio):** Less opportunity for early termination
- **Large nmax:** Need to search more thoroughly
- **Very large dmax:** Search region is larger

### When is speedup higher?

- **Large grids:** Optimization benefits compound
- **High NaN ratio:** Shell expansion terminates early
- **Small nmax:** Find enough neighbors quickly
- **Moderate dmax:** Efficient search region

---

## Continuous Integration

For automated testing:

```matlab
% run_all_tests.m (create this for CI)

success = true;

try
    test_neighbours_stug_optimized;
catch e
    fprintf('Correctness test FAILED: %s\n', e.message);
    success = false;
end

try
    compare_neighbours_implementations;
catch e
    fprintf('Comparison test FAILED: %s\n', e.message);
    success = false;
end

if success
    fprintf('\nAll tests PASSED ✓\n');
    exit(0);
else
    fprintf('\nSome tests FAILED ✗\n');
    exit(1);
end
```

---

## Support

If tests fail or results are unexpected:

1. Check MATLAB version (tested on R2019b+)
2. Verify grid data structure format
3. Check available memory
4. Review test output for specific errors
5. Enable debugging: `dbstop if error`
6. Consult OPTIMIZATION_GUIDE.md for algorithm details
