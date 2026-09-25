# Composition and conservative gravity work — 2026-09-24

GCC 13.3.0, Release, serial builds of all three dimensions: **504/504 tests passed**.
Full build completed without compiler warnings. The uploaded working-tree
changes, including the existing dual-energy implementation, were retained.

New coverage includes all 118 elements and case-insensitive aliases, effective
mixture A/Z, validation, INI and CLI overrides, no allocated active rho column,
tracer exclusion from density, nonuniform conservative species fluxes and AMR
split/merge, runtime reflux and uniform composition preservation, and Silo
readback of material and tracer fields.

Gravity work tests retain the original momentum kicks and verify that the
final energy correction leaves momentum exactly unchanged. Uniform and mixed
level meshes, isolated and periodic boundaries, with and without species,
are tested for three consecutive steps. Coupling error is checked within
4e-13 of the energy L1 scale, independently subtracting the measured potential
operator reciprocity defect. On these tests the uncorrected physical energy
drift stays below 4e-13 for uniform grids and 1e-9 for adaptive grids. These
are small regression fixtures, not long-term production accuracy claims.

The polytrope regression exposed nonpositive conservative E in its dilute
atmosphere. Dual energy now accepts finite E when its entropy auxiliary
supplies admissible thermodynamics; a dedicated test verifies positive pressure
and temperature and no clipping or auxiliary reset. Disabled dual energy still
requires admissible total-energy thermodynamics. All existing hydro, radiation,
AMR, gravity, storage, output, polytrope and RT regression checks pass.

Reproduce with installed Boost, Silo/HDF5 and GoogleTest dependencies:

```bash
cmake -S . -B build -DOCTOTIGERII_WITH_HPX=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
ctest --test-dir build --output-on-failure -j4
```

The attached [CTest log](ctest.log) records the complete run. HPX was not
installed in this environment, so distributed execution and the updated HPX
serialization checks were not run. The existing CMake distributed-test
registration includes the new Application test targets when HPX is enabled.

Known limits: first-order concentration transport; the existing ideal-gas EOS
uses its configured mean molecular weight; adaptive FMM reciprocity error is
not redistributed; gravitational binding-energy changes at regridding are not
compensated. See [mass fractions](../../mass-fractions.md) and
[gravity energy](../../gravity-energy.md).
