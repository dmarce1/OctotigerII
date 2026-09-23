# Gravity image validation

Validated with GCC 13.3, Release builds, GoogleTest, and HPX 1.11 from the available
profiling-enabled installation. All lengths, masses, potentials and accelerations
use the existing CGS interfaces.

| Build / check | Result |
| --- | --- |
| gravity-sphere, serial full CTest | 101 / 101 passed |
| gravity-sphere, HPX full CTest | 111 / 111 passed |
| Two-locality TCP gravity suite | Passed, including image fields and exchange counts |
| Three-locality TCP gravity suite | Passed, including image fields and exchange counts |
| Collapse, three periodic axes, one coupled step | Completed |
| Collapse, lower reflecting x / periodic y / open z, one coupled step | Completed |

The CTest logs include isolated regression tests, new independent Ewald kernel
checks, all seven periodic-axis masks, single and paired reflections, mixed
boundaries, root corrections, and partitioned repeated solves. Derivatives are
checked through degree 20 (the existing maximum multipole order is 10).

The sparse 8³-cell image-force convergence test compares p=3 and p=8 at theta=0.5.
For each of its eleven periodic/reflected geometries, the p=8 maximum force error
is less than 2e-5 times the maximum reference force component and less than 3%
of the p=3 error. The high-order kernel tensor check varies alpha and quadrature
resolution, including large open-coordinate separations. Independent references
also cover the 1P long direct image sum, 2P Fourier slab solution, the 3P cubic
lattice constant, and finite mirrored copies without an Ewald kernel.

The application direct reference shares the Ewald pair kernel but performs no
multipole approximation. It is supplemented by the independent checks above.
These are correctness checks on small domains and local TCP processes, not an
HPC scaling benchmark. Reproduce through the normal CMake/CTest build; the
commands and boundary semantics are in ../../gravity-images.md.
