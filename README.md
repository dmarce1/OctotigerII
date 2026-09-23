# OctotigerII

A new, deliberately small CPU project built from the modular Octo-TIGER
numerical kernels and the separate diagonal Cartesian FMM. This is a new
application and CMake project, not a compatibility mode inside Octo-TIGER.
The source archive includes the project's Git history.

## Build and run

Each executable contains **one problem and one compile-time dimension**.
Select them with `OCTOTIGERII_PROBLEM` and `OCTOTIGERII_NDIM`. Keep a separate
build directory for each combination. Resolution and refinement level remain
runtime settings; changing either does not require recompilation.

From HOME, using an existing HPX installation:

```bash
cmake -S ~/workspace/OctotigerII -B ~/workspace/OctotigerII/release/sod/2d \
  -DCMAKE_BUILD_TYPE=Release \
  -DOCTOTIGERII_PROBLEM=sod -DOCTOTIGERII_NDIM=2 \
  -DCMAKE_PREFIX_PATH="$HOME/local/hpx/Release;$HOME/local/silo;$HOME/local/hdf5"
cmake --build ~/workspace/OctotigerII/release/sod/2d -j 12
ctest --test-dir ~/workspace/OctotigerII/release/sod/2d --output-on-failure

~/workspace/OctotigerII/release/sod/2d/octotigerII-sod-2d \
  --config="$HOME/workspace/OctotigerII/bin/hydro_tests/Sod/inputs" \
  --mesh.cells=16 --mesh.level=2 --hpx:threads=12
```

For a serial build, set `-DOCTOTIGERII_WITH_HPX=OFF` and choose another build
directory, such as `build-serial/sod/2d`. The numerical code and Silo output
are the same. This configuration requires the Boost.Program_options development
library; the HPX build uses HPX's own Program_options. `release` and `serial`
CMake presets select Sod 1D in separate
`release/sod/1d` and `build-serial/sod/1d` directories.

The existing noninteractive helper can build HPX locally first:

```bash
~/workspace/OctotigerII/build.sh release --problem kelvin-helmholtz --ndim 3 -j 12
```

It shares HPX under `packages/TYPE/hpx`, while application builds live in
`TYPE/PROBLEM/NDIMd`. Debug and RelWithDebInfo work the same way. It requires
installed CMake, a C++20 compiler, Git, Boost, hwloc, Silo, HDF5, and zlib.
Use `CC`, `CXX`, `CMAKE_PREFIX_PATH`, and `Silo_ROOT` for local installations.
HPX 1.11.0 is fetched by this helper and fetches its own Asio headers.
Ordinary CMake builds require CMake 3.22+, Boost 1.71+, and (unless disabled)
HPX 1.11+ with its distributed runtime. There is no accelerator build path yet.

All localities in an HPX run must use the same executable/build configuration.
Use your HPX installation's launcher and pass HPX options as `--hpx:name=value`.
Field partitions are distributed over those localities; workers prefer local
work and may steal tasks. The installed executable keeps its problem/dimension
name, and `cmake --install` installs the selected problem's inputs alongside it.

`OCTOTIGERII_BUILD_TESTS=OFF` disables the existing regression executables.
Google Test migration and the interactive builder are deferred in
[ROADMAP.md](ROADMAP.md). See [BUILDING.md](BUILDING.md) for the problem layout,
adding a problem, dimension semantics, and Eclipse configuration.

To archive source plus Git history, run `~/workspace/OctotigerII/arc.sh`.
Builds, dependencies, and output directories are excluded.

## Problems

| CMake problem | Directory under `bin/` | Supported `ndim` |
| --- | --- | --- |
| `sod` | `hydro_tests/Sod` | 1, 2, 3 |
| `kelvin-helmholtz` | `hydro_tests/KelvinHelmholtz` | 2, 3 |
| `streaming` | `radiation_tests/Streaming` | 1, 2, 3 |
| `radiation-pulse` | `radiation_tests/RadiationPulse` | 1, 2, 3 |
| `gravity-sphere` | `gravity_tests/Sphere` | 3 |
| `gravity-gaussian` | `gravity_tests/Gaussian` | 3 |
| `collapse` | `science/Collapse` | 3 |

Each directory owns its `CMakeLists.txt`, `problem.cpp`, and `inputs`.
Sod is planar along x in every dimension; Kelvin–Helmholtz uses x/y and is
extruded uniformly along z in 3D. Streaming propagates along the diagonal of
the compiled dimension. The pulse is initially isotropic. Static gravity
problems use `runtime.stopTime=0`; collapse evolves gas and gravity with the selected image boundaries.

Inputs contain dotted `key=value` entries; `--key=value` overrides them on
the command line, independent of argument order. For example, set
`gravity.multipoleOrder=5` and `gravity.openingAngle=0.5` in the gravity
problem's `inputs`, or override them with
`--gravity.multipoleOrder=4 --gravity.openingAngle=0.45`. The same input
can use INI sections, for example `[gravity]` followed by
`multipoleOrder=4`. The typed `Config` has corresponding groups such as
`config.gravity.multipoleOrder`. `--help` lists canonical names; the
previous snake_case spellings remain accepted as aliases. Specifying both
spellings of one setting in the same CLI or input file is an error. Repeated
`--config` files are applied in order, then command-line overrides apply.
The HPX build parses these settings with `hpx::program_options`; the
HPX-free serial build uses `boost::program_options` with the same settings.
`mesh.cells` is a power-of-two number of cells per block
per axis; `mesh.level` gives `2^level` blocks per axis. The current mesh is
fixed-level Cartesian. `problem.name` and `mesh.ndim` are no longer runtime
options: remove them from old input files and select their CMake equivalents.
Unknown settings and unsupported build combinations fail explicitly.

## Units and output

All inputs, stored physical quantities, and outputs use cgs: cm, s, g,
erg/cm³, cm²/s² for potential, and cm/s² for acceleration. Hydro momentum is
momentum density. Radiation flux is erg/(cm² s). The numerical radiation state
stores physical `F` in erg/(cm² s); calculations temporarily form `Q=F/c`
using physical `c=2.99792458e10 cm/s`.
`radiation.lightSpeedRatio` changes transport speed only.
Physical constants live in `octotigerII::constants`, declared in
`include/octotigerII/units/constants.hpp`: `G`, `c`, `atomicMassUnit` (`m_u`),
`boltzmann` (`k_B`), and `radiation` (`a_r`). Their values follow CODATA 2022;
the radiation constant is derived from exact `k_B`, Planck's constant, and `c`.

`octotigerII::units` defines Boost quantities using centimeters, grams,
seconds, and kelvin. Stored hydro and radiation fields, reconstructed states,
fluxes, mesh geometry, times, gravity fields, and integrated diagnostics retain
these types. Hydro stores a heterogeneous tuple because density, momentum
density, and energy density have different dimensions; radiation stores `ndim+1`
physical quantities `(E,F)`, with energy-density and energy-flux units. A flux has its own dimensions and must be
multiplied by `dt/dx` before it can update a state. Dimensionless coefficients,
ratios, and limiter weights remain plain numbers.

For example:

```cpp
using namespace octotigerII;
auto rho = units::Density::from_value(1e-7); // g/cm^3
auto temperature = units::Temperature::from_value(1e4); // K
units::Pressure p = rho / constants::m_u * constants::k_B * temperature;
// rho + p is a compile-time error.
```

Numeric values cross the application boundary at input parsing, HPX archives,
terminal output, and Silo. Serialization writes each quantity's CGS value and
restores its declared unit on load; it is not a checkpoint implementation.
The diagonal FMM's mathematical coefficients remain dimensionless, explicitly
normalized by one gram and one centimeter. The adapter restores physical units
before applying the shared `G`.

Each output frame is a Silo database with dimensional multimeshes,
zone-centered fields with CGS unit labels, coordinates labeled in cm,
time in seconds, and cycle. `frames.visit` lists the frame series.
Silo creation uses `DB_CLOBBER`.

Each example has its own relative `output.directory`; paths are relative to
the working directory. Frame zero and the final frame are always written;
`output.every` selects the intervening step cadence. Step summaries appear
on the terminal. Rerunning in the same directory overwrites matching filenames;
choose a new directory to retain an earlier run. Output is visualization,
not a restart checkpoint. `output.enabled=off` creates no output files.

The 1D/2D integrated diagnostics use a unit transverse measure: 1D totals
describe a 1 cm² area, and 2D totals describe a 1 cm thickness. There are no
inactive mesh axes, singleton cells, or stored transverse vector components. Mass is consequently in grams and energy in ergs
in every dimension. Numerical totals match the previous per-unit convention.
Open boundaries can change integrated quantities. No conserved gravitational
field energy is reported.

## Small architecture

- `storage::StoragePartition` owns typed columns shared by all enabled physics.
  `Field<T>` allocates a scalar column with configurable banks; transport uses two.
  Its range addresses are independent of any mesh or particle topology.
- `Subgrid` is geometry and an interior range; it owns no numerical arrays.
- `CartesianTopology` supplies block placement and reusable halo geometry.
- `Runtime` owns the field repository and one work scheduler per locality.
  Stages read immutable inputs and publish disjoint outputs after completion.
- Local kernels view field interiors directly. Only halos and solver work arrays
  occupy temporary worker storage. HPX handles remote transfers and coalescing.
- `simulation.cpp` coordinates timesteps and compiled physics; the selected
  `bin/.../problem.cpp` supplies initialization and problem-local configuration.
- `runtime.workerTasks=0` selects the HPX worker count. A positive value bounds
  active block tasks per locality. `runtime.workStealing=on/off` controls remote
  task execution; persistent ownership does not move when a task is stolen.

See [STORAGE.md](STORAGE.md) for ownership, stage safety, zero-copy behavior,
and how another topology or particle field can use the same storage API.

Hydro uses the existing modular unsplit MUSCL–Hancock integrator, primitive
PLM reconstruction, HLLC with HLL fallback, and positivity limiting. M1 uses
the same finite-volume scaffold, the existing Skinner–Ostriker closure and
HLL transport, and realizability limiting. Gas and radiation share machinery,
but there is currently no gas–radiation exchange solver or combined example.

Gravity uses a new cell-octree driver around the recovered compact Cartesian
plane-wave FMM operators. Moments and locals retain `(p+1)²` real coefficients
for `p=1..10`. An upward M2M pass, target-owned diagonal M2L interactions,
and a downward L2L pass produce potential and acceleration. The opening test
is `cell_width / separation < gravity.openingAngle`, with the angle strictly
below `1/sqrt(3)`. Near leaf pairs use Newtonian direct summation. Cells are
point masses at their centers and self contributions are omitted. No old
Cartesian Taylor gravity backend, torque correction, or alternate solver is
included. Direct summation in the tests is only an independent reference.

The collapse example uses a half gravity kick, a hydro transport step, a new
gravity solve, and a second half kick. The kick updates momentum and gas
kinetic energy while preserving internal energy. There are no gravitational
energy fluxes, `dphi/dt` evolution, or exact total-energy conservation claims.

## Scope of this first version

Included: CPU numerics, 1D/2D/3D fixed tiled meshes, per-face periodic/reflecting/outflow/analytic transport,
isolated 3D gravity, shared timesteps, distributed field storage and locality work queues, and Silo output.

Omitted: CUDA, HIP, Kokkos, Vc, CPPuddle, Unitiger, the old FMM, old problems
and test harnesses, SCF, binary-star setup, rotating frames, species/degenerate
EOS, radiation opacities/coupling/subcycling, dynamic AMR, shadow hierarchies,
dynamic repartitioning, checkpoints, and old command-line compatibility.

Transport fetches only the required halo ranges. The application FMM partitions
moments and local expansions at every tree level across HPX localities. Bounded
HPX workers process target cells within each locality; source moments and parent
locals are fetched in deduplicated batches. Density comes directly from field
storage, and gravity is written directly to the next field bank. No full density
or FMM hierarchy is gathered for the solve. The hierarchy persists across steps.

`--hpx:threads=N` sets worker threads per process; `runtime.workerTasks=0` uses
that count, and a positive value caps concurrent gravity workers per locality.
FMM cells can be divided across localities even when the mesh has only one block.
The console reports `cellsPerLocality` for the last solve and the total number of
bounded `workerTasks` dispatched. FMM ownership is static; remote transport work
stealing does not move FMM targets. In an HPX-free build the same field solver
runs serially. The standalone `gravity::solve` remains a serial numerical reference.

Interaction pairs and the opening criterion are preserved. Each target accumulates
its own direction of a pair in deterministic order; this avoids concurrent writes.
Pair counters still count unordered pairs once. Rounding can differ from the old
symmetric traversal. See [the FMM execution design](docs/parallel-fmm.md).

Snapshots, diagnostic reductions, direct-reference verification, and Silo output
still use the coordinating locality. There is a synchronization between tree
levels; coarse levels have limited concurrency. Dynamic AMR, dynamic FMM
repartitioning, and checkpoint/restart remain future work.

See [VALIDATION.md](VALIDATION.md) for the checks actually performed and
[PROVENANCE.md](PROVENANCE.md) for the source selection.

## Developer documentation

The source archive includes the generated manual at `docs/generated/html/index.html`.
Open that file directly in a browser, or regenerate it after changing the source.

Build and open the Doxygen reference from HOME:

```bash
sudo apt install doxygen
~/workspace/OctotigerII/docs.sh --open
```

Doxygen 1.9.8+ and Python 3 are sufficient; this does not compile the application.
The HTML includes source browsing, class documentation, storage ownership rules,
physical units, and author–year citations linked to [BIBLIOGRAPHY.md](BIBLIOGRAPHY.md).
An existing CMake build also has a `docs` target. See [docs/index.md](docs/index.md).

## Analytic comparisons

Sod, spherical gravity, and periodic streaming runs automatically report errors
against exact solutions and write exact/error Silo fields plus
`analytic-errors.json`. Gravity defaults to direct summation of the same
cell-center masses, using reproducible target sampling for large problems.
Use `--verification.directSamples=2048 --randomSeed=5489` to select a sample
size and global seed, or `--verification.gravityReference=continuum` for the
original smooth sphere/Gaussian comparison. Orders 1 through 10 are supported.
Use `--verification.analytic=on` to require a reference.
See [analytic verification](docs/analytic-verification.md) for available
solutions, accuracy gates, interpretation, and the standard problem interface.

Physical boundary options and analytic problem hooks are documented in [docs/boundaries.md](docs/boundaries.md).

Periodic and reflecting gravity, Ewald image sums, and examples are documented in
[docs/gravity-images.md](docs/gravity-images.md).
