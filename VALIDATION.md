# Validation of the initial OctotigerII project

Host checks performed with GCC 13.3.0, C++20, on Linux x86-64. CMake configured
and built the actual new project; no old application stubs or mocked physics
were substituted. Numerical libraries were compiled with `-Wall -Wextra
-Wpedantic`.

## Completed

- Release serial CPU build: **10/10 CTest checks passed**.
- Debug serial build with AddressSanitizer and UndefinedBehaviorSanitizer:
  **10/10 checks passed**. Leak detection was disabled for this execution
  environment; leak checking is not claimed.
- Serial CPU plus Silo 4.11/HDF5 1.10.10: **11/11 checks passed**.
- HPX 1.11.0 plus Silo, two worker threads on one locality:
  **11/11 checks passed**, including all seven ordinary problems.
- Two independent HPX processes/localities using TCP and two threads each:
  Sod (64 cells), streaming (64 cells), and collapse (512 cells) completed.
  All final output fields agreed with their one-locality runs to the checked
  tolerance of `2e-13` relative plus `1e-20` absolute. This exercised remote
  component creation, snapshots, field serialization, advances, gravity
  assignment, and source-kick actions. It is a single-host network check,
  not a multi-node cluster or scaling test.

The small test suite consists of two numerical checks, the seven ordinary
example problems, an application/output check when Python is present, and a
Silo readback check when Silo is enabled. The numerical and application checks
cover:

- Hydro agreement between one block and a tiled mesh; periodic conservation
  of all five conserved fields; positive density and pressure.
- M1 agreement between one block and tiled meshes in 1D, 2D, and 3D;
  periodic conservation of all four fields; `E >= 0`, `|F| <= c E`.
- Reduced-speed streaming compared with a translated analytic Gaussian;
  physical output flux remains `F=cE` rather than `cHat E`.
- The gravity kick preserves internal energy and applies the specified
  momentum impulse.
- The new production FMM tree driver compared with an independent direct
  sum over an 8³ nonuniform density fixture. Both M2L and P2P paths execute;
  the reference uses the same cell-centered point sources and excludes self.
- FMM error improves with expansion order; potential and acceleration scale
  correctly when cell width changes.
- Unknown options, invalid mesh sizes and dimensions, unsupported gravity
  boundaries, invalid light speed, and premature maximum-step termination
  return failure.
- Actual CSV output contains the expected time, cell count, density, and
  physical gravity fields.
- Silo files are opened and read back in 1D/2D/3D. Coordinates, multimeshes,
  zone centering, data, cycle, time, and physical flux conversion are checked.
  The same output file is written twice to verify `DB_CLOBBER`.

FMM relative RMS errors (opening angle 0.5):

| p | Potential | Acceleration |
| --- | ---: | ---: |
| 3 | 2.025507e-05 | 7.413578e-04 |
| 4 | 3.701440e-06 | 1.836427e-04 |
| 5 | 1.558004e-06 | 4.481547e-05 |

Streaming energy relative L1 errors, 128 cells, t=0.2 s:

| cHat/c | Relative L1 |
| --- | ---: |
| 1 | 1.418577e-03 |
| 0.25 | 4.121192e-04 |

These finite-resolution checks establish basic correctness of this first
implementation. They do not establish astrophysical production accuracy,
mesh convergence of every example, gravitational energy conservation, or
large-scale performance. The gravity test compares the discrete point-mass
problem, not an exact continuum sphere potential.
