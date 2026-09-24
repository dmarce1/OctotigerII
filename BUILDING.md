# Dimension builds {#problem_builds}

One CMake configuration builds `octoII-1d`, `octoII-2d`, and `octoII-3d`.
`octoII` is a relative symbolic link to `octoII-3d`, in both the build and
installation directories. Problems are selected at runtime with
`--problem.name=sod` on the command line. There is no default problem; an input
file alone does not select one.
Dimension is fixed by the executable; `--ndim` and `--mesh.ndim` are rejected.

From HOME, with an existing HPX installation:

```bash
cmake -S ~/workspace/OctotigerII -B ~/workspace/OctotigerII/release \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$HOME/local/hpx/Release;$HOME/local/silo;$HOME/local/hdf5"
cmake --build ~/workspace/OctotigerII/release -j 12
ctest --test-dir ~/workspace/OctotigerII/release --output-on-failure -j 2
~/workspace/OctotigerII/release/octoII-1d \
  --problem.name=sod --config="$HOME/workspace/OctotigerII/bin/hydro_tests/Sod/inputs" --hpx:threads=12
~/workspace/OctotigerII/release/octoII --problem.name=gravity-sphere --hpx:threads=12
```

To build HPX locally first, use `~/workspace/OctotigerII/build.sh release -j 12`.
The helper shares dependencies under `packages/TYPE` and builds all dimensions
in `TYPE/`. `--problem`, `--ndim`, and `--tests-only` are no longer build options.
Use a fresh directory when migrating from an old problem-specific build.
The old `OCTOTIGERII_PROBLEM` and `OCTOTIGERII_NDIM` cache choices are rejected.

## Physics modules

All three CMake options default to ON:

| CMake option | Module |
|---|---|
| `OCTOII_WITH_HYDRO` | Hydrodynamics |
| `OCTOII_WITH_RADIATION` | M1 radiation transport |
| `OCTOII_WITH_GRAVITY` | Self-gravity |

For example, append `-DOCTOII_WITH_GRAVITY=OFF` to CMake or `build.sh`.
The leading `-D` is CMake's option syntax, not part of the option name.
Disabled solver sources and problem implementations are omitted from the build.
The registry retains unavailable problem names so selecting one prints a clear
error and exits with status 1. `--help` lists problems and availability reasons.
An all-OFF build is supported; every production problem is then unavailable.
Gravity and Rayleigh–Taylor remain 3D-only. Kelvin–Helmholtz requires 2D or 3D.
Uniform external acceleration remains a hydro feature and does not require the
self-gravity module, but the existing 3D restriction remains in force.

Compiled availability is distinct from activity: Sod allocates and advances
hydro only, streaming radiation only, gravity fixtures gravity and density,
and collapse/polytrope hydro plus gravity. Problem identity is serialized with
Config for HPX localities. All localities must run the same dimensional binary.

## Problem manifests and implementation

CMake discovers `bin/*/*/CMakeLists.txt`, for example:

```cmake
octo_problem(sod DIMENSIONS 1 2 3 DEFAULT_DIMENSION 1 HYDRO)
```

`DIMENSIONS` and module requirements control registry availability. The legacy
`DEFAULT_DIMENSION` manifest field is descriptive; it does not select a build.
Each implementation uses its own namespace (`octotigerII::sod`,
`octotigerII::kelvin_helmholtz`, etc.) and implements the hooks declared in
`problems.hpp`: defaults, validation, initialization, analytic reference, and
analytic boundary data. CMake generates the dispatcher for each dimension.
Add a manifest, namespaced implementation, and inputs to add a problem.

The parser requires a command-line problem, reads all input files, applies the
selected problem's defaults once, then applies
input settings in the same precedence order. Problem-specific parameters remain
typed members of Config. Initializers continue receiving cell width through
Snapshot and may widen unresolved features during startup refinement probes.

## Compile-time dimension

Each dimension has its own generated header and libraries. Coordinates have
exactly `ndim` entries, blocks contain `cells^ndim` interiors, and there are
`2^ndim` children. Hydro stores `ndim+2` scalars and radiation `ndim+1`.
Lower-dimensional builds evolve only active momentum and flux components.
The M1 isotropic pressure remains E/3; CGS totals use unit transverse measure.
Do not mix generated dimension headers in an include search path.

## Dependencies and tests

CMake 3.22+, a C++20 compiler, Boost 1.71+, Silo, HDF5, and zlib are required.
HPX 1.11+ with distributed runtime is the default. Set
`-DOCTOTIGERII_WITH_HPX=OFF` for the serial backend, which uses Boost.Program_options.
The `release` and `serial` presets build all dimensions in separate directories.
GoogleTest 1.12+ is reused when installed, otherwise pinned 1.14.0 is fetched.
For offline builds set `OCTOTIGERII_FETCH_GOOGLETEST=OFF` and provide GTest_ROOT.
`OCTOTIGERII_BUILD_TESTS=OFF` omits tests. CTest names include dimension and suite.
`tests/run_matrix.py --module-matrix` also tests disabled-module combinations.

## Eclipse and profiling

The default compilation database is `release/compile_commands.json`. A source
has one compile command per dimension; select the desired variant in the IDE.
`build.sh` retains HPX APEX/PAPI setup. See [profiling](docs/profiling.md).
