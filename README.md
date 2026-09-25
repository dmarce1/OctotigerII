# OctotigerII

A new, deliberately small CPU project built from the modular Octo-TIGER
numerical kernels and the separate diagonal Cartesian FMM. This is a new
application and CMake project, not a compatibility mode inside Octo-TIGER.
The source archive includes the project's Git history.

## Build and run

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
in `TYPE/`. `--problem`, `--ndim`, and `` are no longer build options.
Use a fresh directory when migrating from an old problem-specific build.
The old `OCTOTIGERII_PROBLEM` and `OCTOTIGERII_NDIM` cache choices are rejected.

See [BUILDING.md](BUILDING.md) for module switches, manifests, dependencies, and tests.

## Problems

| Runtime problem | Directory under `bin/` | Supported `ndim` |
| --- | --- | --- |
| `sod` | `hydro_tests/Sod` | 1, 2, 3 |
| `kelvin-helmholtz` | `hydro_tests/KelvinHelmholtz` | 2, 3 |
| `rayleigh-taylor` | `hydro_tests/RayleighTaylor` | 3 |
| `streaming` | `radiation_tests/Streaming` | 1, 2, 3 |
| `radiation-pulse` | `radiation_tests/RadiationPulse` | 1, 2, 3 |
| `gravity-sphere` | `gravity_tests/Sphere` | 3 |
| `gravity-gaussian` | `gravity_tests/Gaussian` | 3 |
| `collapse` | `science/Collapse` | 3 |
| `polytrope` | `science/Polytrope` | 3 |

[Polytrope](docs/polytrope.md) initializes an isolated Lane–Emden star with hydro,
gravity, a configurable radius and center, and density-based AMR.

Each directory owns its `CMakeLists.txt`, `problem.cpp`, and `inputs`.
Sod is planar along x in every dimension; Kelvin–Helmholtz uses x/y and is
extruded uniformly along z in 3D. Streaming propagates along the diagonal of
the compiled dimension. The pulse is initially isotropic. Static gravity
problems use `runtime.stopTime=0`; collapse evolves gas and gravity with the selected image boundaries.

[Rayleigh–Taylor](docs/rayleigh-taylor.md) places heavy fluid above light fluid
under uniform downward gravity in 3D. The x and y faces are periodic and both
z faces reflect.

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
per axis; `mesh.level` gives the initial `2^level` blocks per axis. Enable
`amr.enabled=on` for mixed spatial levels, evolved shadows, and mass/error
refinement criteria. See [AMR configuration and design](docs/amr.md). `problem.name` selects the problem at runtime. Remove `mesh.ndim` from old input
files and choose the corresponding executable. Unknown settings, unavailable
modules, and unsupported problem dimensions fail explicitly.

## Units and output

All inputs, stored physical quantities, and outputs use cgs: cm, s, g,
erg/cm³, cm²/s² for potential, and cm/s² for acceleration. Hydro momentum is
momentum density. Radiation flux is erg/(cm² s). The numerical radiation state
stores physical `F` in erg/(cm² s); calculations temporarily form `Q=F/c`
using physical `c=2.99792458e10 cm/s`.
`radiation.lightSpeedRatio` changes transport speed only.
Physical constants live in `octotigerII::constants`, declared in
`octotigerII/units/constants.hpp`: `G`, `c`, `atomicMassUnit` (`m_u`),
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

Momentum, velocity, radiation flux, and gravitational acceleration are stored
as separate scalar components with `DBPutQuadvar1`, such as `momentumX`,
`momentumY`, and `momentumZ` in 3D. VisIt vector expressions named `momentum`,
`velocity`, `radiationFlux`, and `acceleration` combine the corresponding
components. These expressions are embedded with `DBPutDefvars` and loaded
automatically, without duplicating the field arrays. Select `momentumZ` for
a Pseudocolor plot or `momentum` for a Vector plot.

Available reference and signed-error components follow the same convention:
`momentumXExact` and `momentumXError` are scalars, while `momentumExact` and
`momentumError` are vector expressions. Only active components are stored;
in 1D/2D, the expressions supply zone-centered zeros in inactive directions.
AMR frames declare changing metadata and connectivity so VisIt refreshes its
domain list after regridding.

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
- `refinement::CellView` supplies typed fields, gradients, shadows, and signal
  speeds to common refinement criteria. `amr::Hierarchy` evolves the coarse
  shadows and supplies conservative transfers during regridding.
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

Hydro uses primitive PLM reconstruction, HLLC with HLL fallback, and positivity
limiting. Global and transport-only stepping retain the modular unsplit
MUSCL–Hancock integrator; coupled gravity uses an explicit midpoint stage built
from the numerical-flux divergence and gravity sources. M1 uses
the same finite-volume scaffold, the existing Skinner–Ostriker closure and
HLL transport, and realizability limiting. Gas and radiation share machinery,
but there is currently no gas–radiation exchange solver or combined example.

Gravity uses a new cell-octree driver around the recovered compact Cartesian
plane-wave FMM operators. Source moments and scalar-potential locals retain
`(p+1)²` real coefficients for `p=1..10`. Separate force locals extend through
degree p+1 so their gradients retain the same degree p as the source expansion;
see the [mutual-force proof](docs/parallel-fmm.md#scalar-reciprocity-and-mutual-force).
An upward M2M pass, target-owned diagonal M2L interactions,
and a downward L2L pass produce potential and acceleration. The opening test
is `cell_width / separation < gravity.openingAngle`, with the angle strictly
below `1/sqrt(3)`. Near leaf pairs use Newtonian direct summation. Cells are
point masses at their centers and self contributions are omitted. No old
Cartesian Taylor gravity backend, torque correction, or alternate solver is
included. Direct summation in the tests is only an independent reference.

Self-gravitating AMR hydro supports hierarchical or conventional gravity time
integration. The hierarchical default assigns each pair interaction to the
slower level's cadence, using physical gas midpoint predictors and accepted HOLD
endpoint impulses. Its energy update uses the actual numerical mass flux and each
interaction shell's endpoint-averaged potential. Globally synchronized runs
retain the two endpoint kicks and replace their temporary self-gravity work.
No `dphi/dt` evolution is needed. See [gravity energy](docs/gravity-energy.md)
for the sequences and boundary accounting. `gravity.energyTreatment=mullen` and
`gravity.conserveRegridEnergy=on` are independent defaults. Select `naive` to
retain kinetic-work kicks, or turn the regrid option off to transfer gas energy
alone. The enabled regrid path transfers E+rho*phi/2 and recovers gas energy after
the new gravity solve. The adaptive FMM uses a reciprocal interaction walk.

## Material fractions

The optional [composition module](docs/mass-fractions.md) stores material
partial densities and massless tracers. Active hydro density is derived from
the material sum, and all components use hydro's stored numerical mass flux.
Definitions accept explicit A,Z or mixtures of all 118 elements by mass.
The entropy auxiliary remains in hydro.

## Dual energy

[Dual-energy hydrodynamics](docs/dual-energy.md) is on by default. The auxiliary
is `A = rho (u/rho^gamma)^alpha` in CGS, with `alpha=1` by default. Pressure and
temperature use total-energy subtraction above a thermal fraction of `0.001`;
A is reset from that subtraction only above `0.1`, after each timestep's AMR
flux corrections. Set `hydro.dualEnergy.exponent` to any finite nonzero value
or `hydro.dualEnergy.enabled=off` to use total energy alone.

## Scope of this first version

Included: CPU numerics, 1D/2D/3D adaptive Cartesian meshes, per-face periodic/reflecting/outflow/inflow/analytic transport,
3D gravity with periodic/reflecting images, uniform external acceleration, shared or level-refined timesteps,
distributed field storage and locality work queues, evolved coarse shadows,
conservative prolongation/restriction and refluxing, Morton-ordered regridding,
an independent adaptive FMM octree, and Silo output.

Omitted: CUDA, HIP, Kokkos, Vc, CPPuddle, Unitiger, the old FMM, old problems
and test harnesses, SCF, binary-star setup, rotating frames, composition-dependent/degenerate
EOS, radiation opacities/coupling/subcycling, temporal AMR, checkpoints,
and old command-line compatibility.

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
The uniform solver counts unordered pairs once; the adaptive traversal counts
directed accepted interactions, as the image solver does. Rounding can differ from the old
symmetric traversal. See [the FMM execution design](docs/parallel-fmm.md).

Snapshots, diagnostic reductions, direct-reference verification, and Silo output
still use the coordinating locality. There is a synchronization between tree
levels; coarse levels have limited concurrency. AMR decisions and shadow storage
currently use the coordinator. Regridding reassigns contiguous Morton segments;
FMM ownership is rebuilt separately. Checkpoint/restart remains future work.

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

## Conservation diagnostics

Every invocation writes `conservation.csv` in `output.directory`, including with
`output.enabled=off`. It contains absolute grid sums, cumulative physical-boundary
transport, L1 norms and normalized drift, plus gas-plus-potential energy when
hydro and gravity are active. See [docs/conservation.md](docs/conservation.md).

Headers live directly under `octotigerII/`; the source root is the include path.

### Level-wise time refinement

AMR transport uses dyadic level subcycling by default. Every cell on a spatial
level shares a timestep; the runtime selects 1/2, 1/4, 1/8, etc. from the CFL
limit and rechecks it before each substep. Self-gravitating AMR runs use
`gravity.timeIntegration=hierarchical` by default; select `conventional` for
full-force updates at each active level's own cadence. Set
`--timestep.refinement=off` to use the globally synchronized reference path.
Runs with imposed uniform external acceleration still use global stepping,
including when self-gravity is also enabled.

The hierarchical option adapts HOLD's accepted interaction impulses to
finite-volume gas transport; it does not claim a symplectic gas integrator.
The implementation caches repeated nested endpoint fields but keeps an unpruned
partial-force upward pass and central coordination of the gravity ledgers, so
subcycling is not yet a performance guarantee. See [time refinement](docs/time-refinement.md) and the
[coupling derivation](docs/gravity-time-coupling-derivation.md).

Smooth fixed-mesh temporal tests observed orders approximately 2.04–2.09 in
both gravity modes. The [validation report](docs/validation/gravity-time-integration.txt)
records the tested scope; the result does not imply second-order behavior at
shocks or arbitrary changes of timestep groups.
