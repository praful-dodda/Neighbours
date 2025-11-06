# Comprehensive Repository Analysis: Neighbours

## Executive Summary

This repository contains neighbor search algorithms optimized for different space-time data formats, along with kriging interpolation functions. The codebase has evolved through multiple iterations, resulting in some redundancy.

**Total Files:** 45 (35 MATLAB scripts, 10 Markdown docs)
**Redundant Files Identified:** 10+ legacy/duplicate files

---

## 📁 CORE NEIGHBOR SEARCH FUNCTIONS (Production)

### 1. **`neighbours.m`** - General Purpose (Reference Implementation)
**Status:** ✅ Production, Reference Standard

**Purpose:** Original neighbor search for arbitrary point clouds (space-time vector format)

**Data Format:**
- Input: Irregular points in space-time (any locations, any times)
- Structure: `c` [n×d] coordinates, `Z` [n×k] values

**Key Features:**
- General-purpose algorithm (works for ANY data arrangement)
- Supports both spatial-only and space-time cases
- Handles variable indexing (different variables at different locations)
- Euclidean distance for spatial, combined space-time metric

**Algorithm:**
1. Compute distances from point of interest to all data points
2. Filter by spatial distance (dmax(1)) and temporal distance (dmax(2))
3. Combine using space-time metric: spatial + dmax(3) × temporal
4. Sort by combined distance, select nmax closest

**Complexity:** O(n) where n = total data points
**Best For:** Irregular point clouds, mobile sensors, arbitrary sampling

**Drawbacks:**
- Doesn't exploit grid structure
- Computes all distances even for gridded data

---

### 2. **`neighbours_stg.m`** - Station × Time Grid Format
**Status:** ✅ Production, Optimized for Station Networks

**Purpose:** Optimized for monitoring station data with irregular spatial locations but regular time series

**Data Format:**
- Input: STG (Space-Time Grid) = Stations × Time
- Structure: `sMS` [nMS×2] station locations, `tME` [1×nME] times, `Z` [nMS×nME] data
- Example: Weather stations, air quality monitors

**Key Features:**
- Exploits regular time structure (no need to search time for each station)
- Pre-computed index mappings (stg_to_stv conversion)
- Handles missing data (NaN) efficiently
- Returns both stg and stv index formats

**Algorithm:**
1. For each monitoring station:
   - Check if within spatial distance dmax(1)
   - For that station, check which time steps within dmax(2)
   - Add valid (station, time) pairs to candidates
2. Compute combined space-time distances
3. Sort and select nmax closest

**Complexity:** O(nMS × log(nME)) - much better than O(nMS × nME) for general approach
**Best For:** Monitoring networks, fixed stations with time series

**Performance:** ~5-20x faster than `neighbours.m` for station data

---

### 3. **`neighbours_stug_optimized.m`** - Uniform Grid with Adaptive Strategy
**Status:** ✅ Production, Highly Optimized

**Purpose:** Optimized for **uniformly gridded** space-time data (satellite, model output)

**Data Format:**
- Input: STUG (Space-Time Uniformly Gridded)
- Structure: `.x` [nx×1], `.y` [ny×1], `.time` [nt×1] vectors, `.Z` [nx×ny×nt] array
- Must have constant spacing: dx, dy, dt

**Key Features:**
- **Adaptive strategy selection:**
  - Shell-based expansion for large search volumes
  - Direct index search for small volumes or sparse data
- NaN ratio estimation (samples 1000 random points)
- Avoids shell expansion if >90% NaN (switches to direct search)
- Deterministic tie-breaking for consistent results
- Pre-filtering by distance constraints

**Algorithm (Shell-based mode):**
1. Find nearest grid point to POI
2. Search in expanding shells (only NEW cells each iteration)
3. Check NaN, compute distances for valid candidates
4. Apply ellipsoid filter if anisotropic
5. Sort and select nmax closest

**Complexity:** O(k × V_shell + m log m) where k=iterations (2-3), V_shell=shell volume
**Best For:** Large uniform grids (satellite data, climate models)

**Performance:** 10-100x faster than `neighbours.m` for uniform grids
**Speedup factors by scenario:**
- Small grid (50³), 30% NaN: 15x
- Large grid (150³), 50% NaN: 50x
- Sparse (95% NaN): 0.9x (switches to direct mode)

---

### 4. **`neighbours_stug_index.m`** - Pure Index-Based with Haversine
**Status:** ✅ Production, Geographic Correctness

**Purpose:** Index-based search for uniform grids with **geographic distance** (Haversine)

**Data Format:**
- Same as `neighbours_stug_optimized.m` (STUG format)
- But uses **lat/lon coordinates** → great circle distance

**Key Features:**
- **Haversine distance formula** (spherical Earth geometry)
- **Antimeridian handling** (dateline at ±180°)
- Latitude-aware spatial constraints (adjusts for cos(latitude))
- Numerical stability for antipodal points (uses atan2)
- Grid uniformity validation (1% tolerance)
- Pre-filtering by distance constraints before sorting
- Ellipsoidal filtering for anisotropic grids

**Critical Geographic Correctness:**
```
Antimeridian crossing:
  POI: 179.5°E, Point: -179.5°W
  Without fix: dlon = -359° → 39,890 km ❌
  With fix: dlon = 1° → 111 km ✅

Latitude adjustment:
  At equator: 1° = 111 km
  At 60°N: 1° = 55 km (factor of 2!)
  At pole: All longitudes equidistant
```

**Algorithm:**
1. Validate grid uniformity (errors if not uniform)
2. Convert dmax from km to grid cells (latitude-aware)
3. Estimate required search radius (spherical approximation)
4. Shell-based expansion in index space
5. For each spatial location, calculate temporal budget (ellipsoid)
6. **Distance calculation ONLY at end** (for selected candidates)
7. Apply Haversine formula with antimeridian wrapping

**Complexity:** O(k × V_shell + nmax) - minimal distance calculations!
**Best For:** Global ocean/atmosphere data, any lat/lon gridded data

**Key Innovation:** Only ~20-30 Haversine distance calculations (vs 1000s in other methods)

---

## 🔧 UTILITY & CONVERSION FUNCTIONS

### 5. **`coord2dist.m`** - Distance Calculation
**Status:** ✅ Utility Function

**Purpose:** Compute Euclidean distance matrices using Kronecker products

**Key Feature:** Vectorized distance calculation (faster than loops)
**Used By:** `neighbours.m`, `neighbours_stg.m`

---

### 6. **`reformat_stg_to_stug.m`** - Uniform Grid Reformatter
**Status:** ✅ NEW - Production Ready

**Purpose:** Reformat data that is **already on uniform grid** from STG to STUG structure

**Important Distinction:**
- Does **NOT interpolate** (no data loss)
- Only **validates uniformity** and **reshapes** arrays
- Errors if grid is not uniform

**Algorithm:**
1. Extract unique x and y coordinates from vectorized format
2. Check spacing uniformity (default tolerance: 1e-6 = 0.0001%)
3. Validate grid completeness (nx × ny = N total points)
4. Reshape from [N×T] to [nx×ny×nt]
5. Create grid vectors and coordinate meshes

**Input:**
```
stg_data.sMS: [N×2] grid coordinates (must be uniform!)
stg_data.tME: [1×T] time values
stg_data.Xms: [N×T] data values
```

**Output:**
```
grid_data.x, .y, .time: grid vectors
grid_data.Lon, .Lat: coordinate meshes
grid_data.Z: [nx×ny×nt] data array
```

**Performance:** < 1 second (just array manipulation)
**Use Case:** When you have uniform grid stored in vectorized format (user's case!)

---

### 7. **`stg_to_stug.m`** - Irregular Station Interpolator
**Status:** ✅ Production

**Purpose:** **Interpolate** irregular monitoring stations onto uniform grid

**Important Distinction:**
- Does **interpolate** (introduces error but fills gaps)
- Creates uniform grid from non-uniform data
- Uses `scatteredInterpolant` for spatial interpolation

**Algorithm:**
1. Determine spatial extent of stations
2. Create uniform grid at specified resolution
3. For each time step:
   - Build `scatteredInterpolant` from valid stations
   - Evaluate on entire grid (vectorized)
   - Store in 3D array

**Input:** Irregular station locations (STG with non-uniform sMS)
**Output:** Uniform grid (STUG format)

**Performance:** 2-10 seconds (depends on resolution and station count)
**Memory:** 2-5x original data
**Accuracy:** Interpolation error (values between stations are estimated)

**Use Case:** When you have irregular monitoring stations and need uniform grid

---

## 📊 KRIGING FUNCTIONS

### 8. **`krigingME_stg.m`** - Kriging for STG Format
**Status:** ✅ Production

**Purpose:** Bayesian Maximum Entropy kriging for station × time data

**Uses:** `neighbours_stg.m` for neighbor selection
**Key Feature:** Handles missing data in space-time grid format

---

### 9. **`krigingME_stug.m`** - Kriging for STUG Format
**Status:** ✅ NEW - Production

**Purpose:** Bayesian Maximum Entropy kriging for uniform grids

**Uses:** `neighbours_stug_index.m` for neighbor selection
**Key Feature:** Optimized for uniformly gridded soft data

---

## 🧪 TEST SCRIPTS

### 10. **`test_neighbours_stg.m`**
**Purpose:** Test suite for `neighbours_stg.m`
**Tests:** Basic functionality, edge cases, NaN handling

---

### 11. **`test_neighbours_stug_optimized.m`**
**Purpose:** Standalone tests for `neighbours_stug_optimized.m`
**Tests:**
- Basic uniform grids
- Edge cases (corners, boundaries)
- Various NaN ratios (0% to 99%)
- Distance constraint enforcement
- Sorting verification

---

### 12. **`test_neighbours_stug_index.m`**
**Purpose:** Comprehensive tests for `neighbours_stug_index.m`
**Tests:**
- Basic functionality
- Anisotropic grids (dx/dy ratios: 0.25, 0.5, 1.0, 2.0, 4.0)
- Various NaN ratios
- Edge cases (corners, boundaries, outside grid)
- Dense (0% NaN) and sparse (95% NaN) data
- Empty/all-NaN data

---

### 13. **`test_krigingME_stg.m`**
**Purpose:** Test suite for `krigingME_stg.m`

---

## 📈 COMPARISON & BENCHMARKING SCRIPTS

### 14. **`compare_neighbours_implementations.m`**
**Purpose:** Direct comparison between `neighbours_stug_optimized.m` and `neighbours.m`

**Features:**
- Validates both find same neighbors
- Compares distances and values
- Measures speedup for each scenario
- Tests 5+ different grid configurations
- Distribution comparison (KS test, mean/median differences)
- Stores results for visualization

**Test Scenarios:**
1. Small grid (100×100×12)
2. Medium grid (200×200×24, 30% NaN)
3. Medium grid (200×200×24, 70% NaN)
4. Large grid (300×300×24)
5. Sparse grid (200×200×24, 95% NaN)
6. Large dense grid (300×300×24, 0% NaN)

---

### 15. **`compare_stug_methods.m`**
**Purpose:** Three-way comparison of all STUG methods

**Compares:**
1. `neighbours.m` (reference)
2. `neighbours_stug_optimized.m` (adaptive)
3. `neighbours_stug_index.m` (index-based)

**Output:**
- Performance comparison (which is fastest)
- Accuracy validation (match percentages)
- Recommendations based on use case

---

### 16. **`compare_neighbours_index.m`**
**Status:** ⚠️ Appears to be duplicate/alternative comparison

---

### 17. **`speed_test_comparison.m`**
**Purpose:** Detailed speed benchmarks with multiple iterations

**Features:**
- Multiple iterations (5+) per scenario for reliable timing
- Tests 6 different grid sizes and NaN ratios
- Statistical analysis (mean, median, std dev)
- Scaling analysis
- Memory estimates
- Saves results to `speed_test_results.mat`

---

### 18. **`speed_test_comparison_index.m`**
**Status:** ⚠️ Appears to be variant for index-based method

---

### 19. **`speed_test_krigingME.m`**
**Purpose:** Performance benchmarks for kriging functions

---

### 20. **`benchmark_neighbours_stug.m`**
**Purpose:** Quick performance check for optimized version only

**Tests:**
- Scalability with grid size
- Performance vs NaN ratio
- Effect of nmax and dmax parameters

---

## 📊 VISUALIZATION SCRIPTS

### 21. **`visualize_performance.m`**
**Purpose:** Generate performance visualization plots

**Requires:** Results from `speed_test_comparison.m`

**Generates:**
- Execution time comparison charts
- Speedup analysis graphs
- Throughput metrics
- Saves 4 PNG files

---

### 22. **`visualize_accuracy.m`**
**Purpose:** Generate accuracy comparison plots

**Requires:** Results from `compare_neighbours_implementations.m`

**Generates:**
- Match quality analysis (% neighbors that match)
- Distance accuracy (max error)
- Value accuracy (max error)
- Accuracy vs performance trade-off
- Saves 4 PNG files

---

### 23. **`visualize_neighbor_selection.m`**
**Purpose:** Visualize spatial distribution of selected neighbors

**Requires:** Results from `compare_neighbours_implementations.m`

**Creates 3 figures:**
1. 2D spatial view (lon-lat slice) with neighbor overlays
2. 3D space-time view with connection lines
3. Distance distribution analysis (histograms, CDFs, Q-Q plots)

**Shows:**
- Common neighbors (both methods selected)
- Reference-only neighbors
- Optimized-only neighbors
- Overlap statistics

---

### 24. **`verify_no_distance_calculations.m`**
**Purpose:** Verify `neighbours_stug_index.m` has no distance calcs during search

**Analysis:**
1. Code structure analysis (parses source to find distance calculations)
2. Algorithmic proof (confirms only index arithmetic during search)
3. Visual comparison (operation counts by phase)
4. Creates verification plots

**Confirms:** Only nmax distance calculations (at the end), not during search

---

## 📘 EXAMPLE/USAGE SCRIPTS

### 25. **`example_reformat_usage.m`**
**Purpose:** 11 examples of using `reformat_stg_to_stug.m`

**Examples:**
- Basic reformatting
- Custom tolerance
- Pre-validation techniques
- Visualization
- Error handling
- Batch processing
- Data integrity checks

---

### 26. **`example_stg_to_stug_usage.m`**
**Purpose:** 9 examples of using `stg_to_stug.m`

**Examples:**
- Basic usage with auto resolution
- Custom resolution selection
- Advanced options (method, bounds, variance)
- Complete workflow
- Resolution selection guide
- Batch processing
- Quality checking
- Memory-efficient processing
- STG vs STUG comparison

---

### 27. **`sample_krigingME_case.m`**
**Purpose:** Example workflow for kriging with ME

---

## 🧪 LEGACY TEST SCRIPTS

### 28. **`testNeighbours.m`**
**Status:** 🟡 Legacy - Original test script

**Purpose:** Original test for `neighbours.m`
**Note:** Still functional but less comprehensive than newer tests

---

### 29. **`testNeighbours_nonan.m`**
**Status:** 🟡 Legacy - Variant without NaN data

---

### 30. **`testKrigingME.m`**
**Status:** 🟡 Legacy

---

### 31. **`testKrigingME_stug.m`**
**Status:** 🟡 Legacy/newer variant

---

## 🔧 UTILITY SCRIPTS

### 32. **`createNetcdf.m`**
**Purpose:** Helper to create NetCDF files for testing/data exchange

---

## 🗑️ REDUNDANT/DEPRECATED FILES

### Legacy STG Versions (Deprecated)

#### 33. **`neighbours_stg_v1.m`**
**Status:** ❌ REDUNDANT - Bug fix version

**History:** Fixed edge case indexing issues
**Replaced by:** `neighbours_stg.m` (current production version)
**Action:** Should be moved to `/archive` or deleted

---

#### 34. **`neighbours_stg_v2.m`**
**Status:** ❌ REDUNDANT - Bug fix version

**History:** Fixed empty matrix handling
**Replaced by:** `neighbours_stg.m`
**Action:** Should be moved to `/archive` or deleted

---

#### 35. **`neighbours_stg_v3.m`**
**Status:** ❌ REDUNDANT - Additional edge case fixes

**History:** Additional edge case handling
**Replaced by:** `neighbours_stg.m`
**Action:** Should be moved to `/archive` or deleted

---

### Original STUG Implementation

#### 36. **`neighbours_stug.m`**
**Status:** ❌ REDUNDANT - Original baseline implementation

**Purpose:** Original uniform grid neighbor search (inefficient)

**Issues:**
- Uses arbitrary `ceil(sqrt(nmax))` for search radius
- Creates full ndgrid for entire search region
- Computes distances for ALL points before filtering
- Memory intensive O(search_volume) vs O(nmax)

**Replaced by:**
- `neighbours_stug_optimized.m` (10-100x faster)
- `neighbours_stug_index.m` (geographic correctness)

**Historical Value:** Baseline for performance comparisons
**Action:** Move to `/archive` - no longer recommended for production

---

## 📚 DOCUMENTATION FILES

### Production Guides

#### 37. **`README.md`**
**Status:** ✅ Main documentation

**Contents:**
- Overview of all functions
- Quick start examples
- Performance comparisons
- Testing guide
- When to use which function

---

#### 38. **`QUICK_REFERENCE.md`**
**Status:** ✅ Quick lookup guide

**Contents:**
- Function comparison table
- Usage examples
- Parameter tuning guide

---

#### 39. **`OPTIMIZATION_GUIDE.md`**
**Status:** ✅ Deep dive into optimizations

**Contents:**
- Detailed strategy explanation
- Algorithm details
- Complexity analysis
- Performance benchmarks

---

#### 40. **`TESTING_GUIDE.md`**
**Status:** ✅ How to run tests

**Contents:**
- Test suite overview
- How to interpret results
- Adding new tests

---

#### 41. **`STUG_INDEX_ALGORITHM.md`**
**Status:** ✅ Index-based algorithm documentation

**Contents:**
- Algorithm overview
- Step-by-step breakdown
- Anisotropy handling
- Complexity analysis
- Verification of design goals

---

#### 42. **`CRITIQUE_STUG_INDEX.md`**
**Status:** ✅ Code review document

**Contents:**
- Comprehensive code analysis (9.9/10 grade)
- Line-by-line review
- Performance analysis
- Best practices observed
- Minor issues identified
- Deployment recommendations

---

### Decision & Strategy Guides

#### 43. **`DATA_REFORMATTING_ANALYSIS.md`**
**Status:** ✅ Analysis of user's data format

**Contents:**
- Data structure breakdown
- STG vs STUG format comparison
- Why data doesn't need interpolation
- Decision matrix

---

#### 44. **`QUICK_DECISION_GUIDE.md`**
**Status:** ✅ Decision tree

**Contents:**
- Quick decision tree
- Visual format comparisons
- Common mistakes to avoid
- Ready-to-use code

---

#### 45. **`REFORMATTING_GUIDE.md`**
**Status:** ✅ Complete reformatting guide

**Contents:**
- When to use `reformat_stg_to_stug.m` vs `stg_to_stug.m`
- Complete workflow examples
- Quality control procedures
- Common issues and solutions

---

#### 46. **`STG_TO_STUG_GUIDE.md`**
**Status:** ✅ Interpolation guide

**Contents:**
- How to use `stg_to_stug.m`
- Resolution selection strategies
- Performance optimization
- Quality control

---

## 🚨 REDUNDANCY ANALYSIS

### Definite Redundancies (Should Remove/Archive)

1. **`neighbours_stg_v1.m`** ❌ - Superseded by `neighbours_stg.m`
2. **`neighbours_stg_v2.m`** ❌ - Superseded by `neighbours_stg.m`
3. **`neighbours_stg_v3.m`** ❌ - Superseded by `neighbours_stg.m`
4. **`neighbours_stug.m`** ❌ - Superseded by optimized versions
   - Keep only for historical/benchmark reference

### Potential Redundancies (Need Verification)

5. **`compare_neighbours_index.m`** ⚠️ - May duplicate `compare_stug_methods.m`
6. **`speed_test_comparison_index.m`** ⚠️ - May duplicate `speed_test_comparison.m`

### Legacy but Potentially Useful

7. **`testNeighbours.m`** 🟡 - Original test, less comprehensive
8. **`testNeighbours_nonan.m`** 🟡 - Specific test case
9. **`testKrigingME.m`** 🟡 - May still be useful

---

## 📊 SUMMARY STATISTICS

| Category | Count | Notes |
|----------|-------|-------|
| **Core Functions** | 4 | Production-ready neighbor search |
| **Kriging Functions** | 2 | BME kriging implementations |
| **Conversion Functions** | 3 | Format conversion & reformatting |
| **Test Scripts** | 7 | Comprehensive test suites |
| **Comparison Scripts** | 5 | Performance & accuracy validation |
| **Visualization Scripts** | 4 | Plot generation |
| **Example Scripts** | 3 | Usage demonstrations |
| **Utility Functions** | 2 | Helper functions |
| **Redundant/Legacy** | 4-6 | Should archive or remove |
| **Documentation** | 10 | Guides and references |
| **Total** | 45 | |

---

## 🎯 RECOMMENDATIONS

### Immediate Actions

1. **Move to `/archive` folder:**
   - `neighbours_stg_v1.m`
   - `neighbours_stg_v2.m`
   - `neighbours_stg_v3.m`
   - `neighbours_stug.m` (keep note it's baseline)

2. **Investigate and consolidate:**
   - `compare_neighbours_index.m` vs `compare_stug_methods.m`
   - `speed_test_comparison_index.m` vs `speed_test_comparison.m`

3. **Update README.md to clearly mark:**
   - Production functions vs legacy
   - Which files users should actually use

### Code Organization Suggestion

```
/Neighbours
  /core              # 4 production neighbor functions
  /kriging           # 2 kriging functions
  /utilities         # coord2dist, createNetcdf, converters
  /tests             # All test scripts
  /benchmarks        # All comparison & speed test scripts
  /visualization     # All visualization scripts
  /examples          # Example usage scripts
  /docs              # All markdown documentation
  /archive           # Deprecated versions (v1, v2, v3, original stug)
```

### Documentation Consolidation

The 10 markdown files could potentially be consolidated:
- Keep: README.md, QUICK_REFERENCE.md
- Merge strategy guides: DATA_REFORMATTING_ANALYSIS + QUICK_DECISION_GUIDE + REFORMATTING_GUIDE → single "DATA_FORMAT_GUIDE.md"
- Merge algorithm docs: OPTIMIZATION_GUIDE + STUG_INDEX_ALGORITHM → "ALGORITHMS.md"

---

## 🏆 PRODUCTION-READY CORE

**For Users:** These are the only functions you need to know:

1. **`neighbours.m`** - For arbitrary point clouds
2. **`neighbours_stg.m`** - For monitoring stations × time
3. **`neighbours_stug_optimized.m`** - For uniform grids (general)
4. **`neighbours_stug_index.m`** - For uniform grids (lat/lon with Haversine)
5. **`reformat_stg_to_stug.m`** - To convert uniform grid format (no interpolation)
6. **`stg_to_stug.m`** - To interpolate irregular stations to grid

Everything else is either:
- Testing/validation
- Documentation
- Examples
- Legacy/deprecated
