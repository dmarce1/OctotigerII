# Source provenance

This new project was assembled on 2026-09-22 from:

- `code.tar(20260922-154024).gz`: modular dimensional mesh, unsplit finite-volume
  integrator, ideal-gas HLLC hydro, source-free M1 transport, scalar/vector
  support, and halo exchange mathematics.
- `octotiger-diagonal-fmm.tar.gz`: the separate diagonal Cartesian FMM
  `fmm.hpp` and `fmm.cpp`. The attached code archive did not contain these
  files. Only these numerical operators were imported from the earlier
  archive; its legacy grid adapter was not imported.

The root CMake project, block runtime, minimal Subgrid application API,
cell-octree gravity driver, options, examples, outputs, documentation, and
small validation suite are new. The namespace/include root is `octotigerII`.
No existing test directories or legacy application sources were copied.
The original AUTHORS and Boost Software License are retained.

The 2026-09-26 radiation-coupling work also consulted the existing OctoI
checkout at `/home/dmarce1/octotiger/src/octotiger`, branch `dominic-wip`,
commit `412c80e773bbad8c3151995e1be542b1456945e0`. Its
`octotiger/radiation/m1.hpp`, `src/radiation/rad_grid.cpp`,
`RADIATION_SO.md`, and `doc/grey-opacity-step-04.md` contain the earlier
Skinner–Ostriker gas/radiation exchange work. That implementation did not
carry into the source-free OctoII import. It supplied equations and design
context; OctoII's typed conservative source solve, forced source stages,
thick-cell flux interpolation, runtime integration, and regression tests
are new implementations. The older documentation describes first-order
source splitting and an isotropic velocity-pressure approximation; those
limitations are not silently adopted here.

The diagonal operator is the earlier separation-aligned plane-wave
factorization, with p+1 Gauss–Laguerre nodes and 2p+1 angular samples. It is not
the optimized six-direction aggregation implementation, and no performance
advantage over the old solver is asserted. See the comments in
`octotigerII/gravity/diagonal/fmm.hpp` and `src/gravity/diagonal/fmm.cpp`.

References carried by the numerical implementation:

- Skinner & Ostriker (2013), M1 radiation transport:
  https://arxiv.org/abs/1306.0010
- Greengard & Rokhlin (1997), plane-wave diagonal translation, section 7;
  the retained code implements the C_XL D C_MX factorization.
- HPX component and CMake integration:
  https://docs.hpx.dev/latest/singlehtml/index.html

The standalone `lib/helmholtz` library translates F. X. Timmes's user-supplied
`helmholtz.tbz` evaluator and the direct electron–positron routines in his
separately retrieved `timmes.tbz`. Original Fortran is preserved in C++ comments.
The separate table-generation driver and its higher-derivative completion
stencil are new OctoII work. Scientific attribution, source SHA-256 hashes,
and upstream licensing status are recorded in
[the library notice](lib/helmholtz/NOTICE.md); the original Timmes code is not
relicensed by the project's Boost license. See
[the validation record](lib/helmholtz/VALIDATION.md) for the distinction between
translation agreement and unresolved full-domain thermal accuracy.
