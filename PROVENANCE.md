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
