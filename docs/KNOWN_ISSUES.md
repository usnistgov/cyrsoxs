# Known Issues

No active issues are currently documented here.

## Recently Resolved

### 2D DC-component replacement for single-slice morphologies

The `Nz == 1` DC-component replacement bug in `replaceDCComponentWithAverage(...)`
has been fixed.

The current behavior is:

1. `Nz == 1` uses a 4-neighbor in-plane average.
2. `Nz > 1` retains the existing 3D six-neighbor behavior.
3. The DC voxel is no longer allowed to alias itself through the wrapped z-neighbor path.
4. The out-of-bounds `z=+1` read for single-slice inputs is eliminated.

The nearby polarization-kernel bounds checks were also tightened from
`threadID > numVoxels` to `threadID >= numVoxels`.

Validation used the NRSS pytest physics suite:

1. `tests/validation/test_2d_disk_contrast_scaling.py`
2. `tests/validation/test_analytical_2d_disk_form_factor.py`
3. `tests/validation/test_analytical_sphere_form_factor.py`
4. `tests/validation/test_sphere_contrast_scaling.py`

Result: 6 tests passed against the local pybind build.
