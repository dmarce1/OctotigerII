# GoogleTest and numerical regression testing

`OCTOTIGERII_BUILD_TESTS=ON` builds real GoogleTest suites. CTest discovers each
`TEST`, fixture, typed test and parameterized case separately. Assertions remain
active in Release builds. Numerical comparisons use analytic values, independent
direct sums, conservation laws, symmetry, convergence and decomposition checks.

## Run the current build

From the source directory after building the selected problem:

```bash
ctest --test-dir release/sod/1d --output-on-failure -j 2
ctest --test-dir release/sod/1d --output-on-failure -L unit
ctest --test-dir release/sod/1d --output-on-failure -L integration
ctest --test-dir release/sod/1d --output-on-failure -R 'Convergence|Transport'
```

For a specific GoogleTest case or repeated/shuffled execution:

```bash
release/sod/1d/tests/hydroChecks --gtest_filter='Hydro.*'
release/sod/1d/tests/storageChecks --gtest_shuffle --gtest_repeat=10 --hpx:threads=2 --hpx:bind=none
```

HPX flags are needed only with HPX builds. Test discovery (`--gtest_list_tests`)
does not initialize HPX. Runtime cases initialize HPX once per process and run
GoogleTest on the console locality; worker localities serve the tested actions.
No change to production `main` is required. Temporary INI and Silo files use
unique directories so concurrent or repeated tests cannot overwrite one another.
The `gravityParallelChecks benchmark [level] [order]` command remains available.

To write a JUnit report:

```bash
ctest --test-dir release/sod/1d --output-on-failure --output-junit test-results.xml
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
| Gravity | Independent point-mass fields, self exclusion, vacuum, superposition, geometric scaling, multipole-order convergence against independent direct sums, worker/block reproducibility, interaction counts, failure/retry |
| Verification | Exact Sod/star states, shell-integrated sphere/Gaussian reference, streaming wraparound, injected error norms, accuracy gates, spatial convergence, reproducible sampling and uncertainty |
| Output and communication | Silo geometry/CGS fields/readback/overwrite, CLI/report checks, typed HPX archive roundtrip, zero-copy buffers, 2- and 3-locality execution |

## Build the complete matrix

The runner reads supported problem/dimension combinations from `bin` manifests;
it does not model unused dimensions as singleton axes. Each combination gets
its own build directory, configure/build/test failures stop the run, and each
build writes `test-results.xml`.

```bash
python3 tests/run_matrix.py --backend serial -j 4
python3 tests/run_matrix.py --backend hpx -j 4 --cmake-arg=-DHPX_DIR=/path/to/hpx/lib/cmake/HPX
python3 tests/run_matrix.py --backend serial --problem streaming --ndim 3 -j 4
```

`--build-type Debug`, `--unit-only`, `--build-root PATH`, repeated `--problem`
and repeated `--cmake-arg=-D...` options are supported. Existing compiler,
`CMAKE_PREFIX_PATH`, and dependency environment settings are inherited.
Hydro-only, radiation-only and gravity-only builds instantiate/link their own
physics tests. Collapse also checks hydro/gravity coupling. HPX-only tests are
registered only when HPX is enabled.

## Dependencies and distributed runs

GoogleTest 1.12+ is found with `find_package`; otherwise the default fallback
fetches GoogleTest 1.14.0 with a SHA-256 checked archive. No package manager or
administrator rights are required. Disable fallback with
`-DOCTOTIGERII_FETCH_GOOGLETEST=OFF`, select an installation with `GTest_ROOT`,
or use `-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=/path/to/googletest` offline.
`OCTOTIGERII_BUILD_TESTS=OFF` avoids GoogleTest discovery/download entirely.

For an HPX build:

```bash
ctest --test-dir release/collapse/3d --output-on-failure -L distributed
```

These tests launch independent HPX TCP localities on the same host, using two
threads per locality. They validate communication/correctness, not inter-node
cluster performance. `OCTOTIGERII_TEST_THREADS` changes thread count. Distributed
CTest entries run serially to avoid resource oversubscription. They require a
working HPX TCP parcelport and local socket access.
