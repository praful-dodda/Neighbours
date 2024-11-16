# neighbours_stg.m 

This repository contains MATLAB scripts that implement an optimized algorithm to find the neighbours of a given point in space-time grid.

This repository contains different versions of a neighbours_stg.m script and test scripts.

## Files

- `neighbours_stg.m`: the main script that was tested
- `neighbours_stg_v1.m`: the first version of the script where there was a bug in the case where the point of interest is at the edge of the grid and the number of neighbours being less than nmax, the code was trying to increase the sub-domain by space and time. There the indexing was incorrect.
- `neighbours_stg_v2.m`: the next version of the script where was a bug. Where you are interested in size of a matrix (and it being empty), you are not able to do logical indexing. So, I changed such that if the size of the matrix is empty, replace it with 0 since you are interested only in the size of the matrix.

## Instructions

1. Clone the repository: `git clone`
2. Run the test script: `test_neighbours_stg.m`

