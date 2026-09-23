# Problem-specific builds {#problem_builds}

The organization follows Castro's separation of common solver code from
problem-local setup, using CMake and OctotigerII's typed storage/runtime.
A build directory selects one problem and one dimension. All executable and
library translation units in that directory use the same generated
`include/octotigerII/buildConfig.hpp`; no dimension is sent as mutable runtime
metadata. Headers that depend on this file must be compiled through the CMake
library targets so they inherit the correct generated include directory.

## Problem manifests and implementation

The top-level CMake configuration discovers `bin/*/*/CMakeLists.txt` manifests.
For example, the Sod manifest contains:

```cmake
octo_problem(sod DIMENSIONS 1 2 3 DEFAULT_DIMENSION 1 HYDRO)
```

`HYDRO`, `RADIATION`, and `GRAVITY` select the required physics. CMake rejects
unsupported dimensions, including any gravity build outside 3D. The generated
header exposes `octotigerII::ndim` and `build::{problem,hydro,radiation,gravity}`.
Only the selected problem's `problem.cpp` is linked. The gravity FMM sources
are compiled only for gravity problems. Unneeded hydro/radiation adapters are
not application dependencies; the optional units/serialization checks may
still build both adapters to check the common physical interfaces.

Each problem implements the three hooks in `problems.hpp`:

- `problemDefaults(Config&)`: defaults before runtime inputs are applied.
- `validateProblem(Config const&)`: restrictions specific to that setup.
- `initializeProblem(Snapshot&, Config const&)`: initialize one interior block.

The mesh, storage and timestep scheduler do not branch on problem names.
The parser registers the selected problem's additional typed parameters;
Rayleigh–Taylor exposes its layer densities, interface pressure, and seed amplitude.
Add a directory, manifest, implementation, and inputs to introduce a
new problem. A new physics implementation still needs its own common-library
integration; manifests select existing modules rather than creating a solver.

Unlike Castro's generated `_prob_params`, this revision does not introduce a
parameter-code generator. The current fixtures use the existing typed Config
settings, including `Config::RayleighTaylorOptions`. Additional options should
remain typed rather than becoming an untyped map inside numerical kernels.

Uniform external gravity uses `Config::HydroOptions::acceleration` and can
be used by a `HYDRO` problem without selecting `GRAVITY` (which enables the
self-gravity FMM). Rayleigh–Taylor selects this hydro-only configuration in
2D and 3D, with downward acceleration along the last active axis.

## Compile-time dimension

`OCTOTIGERII_NDIM` sets the literal `ndim` constant. There is no runtime fallback
or dispatch between 1D/2D/3D solvers. Coordinates and geometry arrays contain
exactly `ndim` entries. Cartesian traversal carries over those entries directly.
A block has `cells^ndim` interiors and `(cells+2*ghostWidth)^ndim` padded temporary
cells. Child counts are `2^ndim`; each transport workspace has `ndim` directional
arrays. Stored hydro states have `ndim+2` scalars; radiation has `ndim+1`.
Silo writes only the corresponding coordinates and vector components.

Lower-dimensional hydro evolves only those momentum components; radiation
evolves only those flux components. This is not a 1D/2D calculation retaining
three-component transverse dynamics. The M1 closure still describes physical
radiation: its isotropic pressure is E/3 even in a 1D or 2D spatial calculation.
The 3D gravity algorithm and physical dimensional exponents retain their
mathematical constants. CGS diagnostic totals retain a unit transverse measure
in reduced dimensions without adding mesh axes or cells.

## Build and run choices

Build choices: problem, dimension, compiler, build type and HPX backend.
Runtime choices: block cells, level, domain bounds, supported boundaries,
physical parameters, timestep and output controls. Changing a runtime choice
never regenerates the build header. A gravity executable remains 3D even if
an input file attempts to say otherwise; runtime problem/dimension options
are rejected with a message identifying the CMake settings.

Use separate build directories such as `release/sod/1d`, `release/sod/2d`, and
`release/collapse/3d`. The executable names are `octotigerII-sod-1d`,
`octotigerII-sod-2d`, and `octotigerII-collapse-3d`. Do not mix executables from
different builds within an HPX run: their state shapes and action payloads differ.

The existing `build.sh` is noninteractive and passes these choices to CMake.
It shares its HPX dependency between problems of the same build type. The
interactive problem/dimension/level wizard is deliberately deferred; see
[the roadmap](ROADMAP.md).

## Eclipse and clangd

The checked-in default Eclipse build path and clangd compilation database point
to `release/sod/1d`, matching the default preset. For another variant, point
Eclipse's build directory and `.clangd`'s `CompilationDatabase` at that variant's
build directory. CMake writes `compile_commands.json` there, including the
matching generated include path. Do not combine dimension-specific generated
headers in one include search path.

## Tests

Tests use GoogleTest with individual CTest discovery. An installed GoogleTest
1.12+ is reused; otherwise CMake fetches the pinned 1.14.0 source without root.
`GTest_ROOT` or `CMAKE_PREFIX_PATH` selects an installation. For offline builds,
set `OCTOTIGERII_FETCH_GOOGLETEST=OFF` or provide an unpacked tree using
`FETCHCONTENT_SOURCE_DIR_GOOGLETEST`. Tests can be omitted with
`OCTOTIGERII_BUILD_TESTS=OFF`.

After `./build.sh release --problem sod --ndim 1`, run:

```bash
ctest --test-dir release/sod/1d --output-on-failure -j 2
ctest --test-dir release/sod/1d --output-on-failure -L unit
```

See [testing](docs/testing.md) for the full coverage list, direct GoogleTest
filters, distributed tests, and the problem/dimension matrix runner.

## Profiling dependencies

`build.sh` enables APEX and PAPI in HPX, enables PAPI in APEX, and lets HPX fetch
its matching APEX source. It reuses installed/module PAPI or builds PAPI 7.2.0
locally without root access. See [profiling](docs/profiling.md) for dependency
selection, annotated regions, and timing/hardware-counter commands.
