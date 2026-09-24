# GoogleTest and numerical regression testing

`OCTOTIGERII_BUILD_TESTS=ON` builds real GoogleTest suites. CTest discovers each
`TEST`, fixture, typed test and parameterized case separately. Assertions remain
active in Release builds. Numerical comparisons use analytic values, independent
direct sums, conservation laws, symmetry, convergence and decomposition checks.

## Run the current build

From the source directory after building all three dimensions:

```bash
ctest --test-dir release --output-on-failure -j 2
ctest --test-dir release --output-on-failure -L unit
ctest --test-dir release --output-on-failure -L integration
ctest --test-dir release --output-on-failure -R 'Convergence|Transport'
```

For a specific GoogleTest case or repeated/shuffled execution:

```bash
release/3d/tests/hydroChecks-3d --gtest_filter='Hydro.*'
release/3d/tests/storageChecks-3d --gtest_shuffle --gtest_repeat=10 --hpx:threads=2 --hpx:bind=none
```

HPX flags are needed only with HPX builds. Test discovery (`--gtest_list_tests`)
does not initialize HPX. Runtime cases initialize HPX once per process and run
GoogleTest on the console locality; worker localities serve the tested actions.
No change to production `main` is required. Temporary INI and Silo files use
unique directories so concurrent or repeated tests cannot overwrite one another.
The `gravityParallelChecks benchmark [level] [order]` command remains available.

To write a JUnit report:

```bash
ctest --test-dir release --output-on-failure --output-junit test-results.xml
```

## Coverage

| Layer | Checks |
| --- | --- |
| Units and math | CGS conversion/constants, compile-time dimensional type restrictions, state dimensions, vector arithmetic, aliasing, integer division, split/concatenate |
| Mesh and storage layout | Ghost/interior indexing, complete face traversal, cell geometry, dimensional child topology, invalid coordinates, time metadata, exact partition capacities, bounds and overflow |
| Options | INI/CLI precedence, aliases, generated help, all supported multipole orders, invalid and nonfinite values, deterministic sampling controls |
| Hydro | Primitive/conserved conversion, independent energy/flux values, sound speeds, identical states, stationary contact, supersonic upwinding, reflection, positivity limiter, invalid EOS/floors |
| Radiation | Isotropic and streaming M1 limits, oblique closure symmetry/trace, causal characteristic speeds, HLL upwinding, physical CGS flux, reduced light speed, realizability, vacuum and bounded roundoff repair |
| Finite-volume evolution | Uniform-state preservation, nonuniform periodic conservation, admissibility, all three limiters, face/edge/corner boundary rules, CFL scaling, invalid timesteps, boundary time tags |
| Runtime and storage | Tiled/single-block agreement, periodic conservation, stage failure/rollback/retry, task accounting, independent fields/banks, type checks, retained buffers, retirement |
| AMR | Conservative split/merge, Morton ordering, mass cap and finest level, periodic signal padding, early regrid cadence, evolved shadows, deep 2:1 face/edge/corner balance, mixed-level reflux conservation, transactional failure, adaptive FMM/direct image comparisons |
| Gravity | Independent point-mass fields, self exclusion, vacuum, superposition, geometric scaling, multipole-order convergence against independent direct sums, worker/block reproducibility, interaction counts, failure/retry |
| Verification | Exact Sod/star states, shell-integrated sphere/Gaussian reference, streaming wraparound, injected error norms, accuracy gates, spatial convergence, reproducible sampling and uncertainty |
| Output and communication | Uniform and mixed-level Silo geometry/CGS fields/readback/overwrite, leaf refinement levels, CLI/report checks, typed HPX archive roundtrip, zero-copy buffers, 2- and 3-locality execution |

## All dimensions and module configurations

One configure/build produces all three dimensions and their GoogleTest suites.
Problem-specific tests select their problem at runtime. The application tests
exercise every problem, including errors for unavailable modules or dimensions.

```bash
./build.sh release -j 12
ctest --test-dir release --output-on-failure -j 2
python3 tests/run_matrix.py --backend serial --module-matrix -j 4
python3 tests/run_matrix.py --backend hpx -j 4 --cmake-arg=-DHPX_DIR=/path/to/HPX
```

The optional module matrix uses all eight HYDRO/RADIATION/GRAVITY ON/OFF
combinations. Each configuration builds all dimensions and writes a JUnit report.
`--build-type Debug`, `--unit-only`, `--build-root PATH`, and repeatable
`--cmake-arg=-D...` options are supported.

## Dependencies and distributed runs

GoogleTest 1.12+ is found with `find_package`; otherwise the default fallback
fetches GoogleTest 1.14.0 with a SHA-256 checked archive. No package manager or
administrator rights are required. Disable fallback with
`-DOCTOTIGERII_FETCH_GOOGLETEST=OFF`, select an installation with `GTest_ROOT`,
or use `-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=/path/to/googletest` offline.
`OCTOTIGERII_BUILD_TESTS=OFF` avoids GoogleTest discovery/download entirely.

For an HPX build:

```bash
ctest --test-dir release --output-on-failure -L distributed
```

These tests launch independent HPX TCP localities on the same host, using two
threads per locality. They validate communication/correctness, not inter-node
cluster performance. `OCTOTIGERII_TEST_THREADS` changes thread count. Distributed
CTest entries run serially to avoid resource oversubscription. They require a
working HPX TCP parcelport and local socket access.
