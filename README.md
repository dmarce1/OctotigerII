# OctotigerII

A new, deliberately small CPU project built from the modular Octo-TIGER
numerical kernels and the separate diagonal Cartesian FMM. This is a new
application and CMake project, not a compatibility mode inside Octo-TIGER.
The source archive includes a new local Git repository with one initial commit.

## Build and run

Requires a C++20 compiler and CMake 3.22 or later. The default build also
requires HPX 1.11 or later with distributed components, Silo with HDF5,
HDF5 development libraries, and zlib. There is no accelerator build path.
Set `CMAKE_PREFIX_PATH` to your installed dependencies as appropriate.
From HOME, with this directory at `~/OctotigerII`:

```bash
cmake -S ~/OctotigerII -B ~/OctotigerII/release \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$HOME/local/hpx/Release;$HOME/local/silo;$HOME/local/hdf5"
cmake --build ~/OctotigerII/release -j
ctest --test-dir ~/OctotigerII/release --output-on-failure

~/OctotigerII/release/octotigerII \
  --config="$HOME/OctotigerII/examples/sod.ini" \
  --output.format=silo --hpx:threads=12
```

`Silo_ROOT` can also identify a Silo installation explicitly. With HPX,
component placement is round-robin over the localities returned by HPX. Use
your HPX installation's normal launcher for multiple localities. Pass HPX
options in `--hpx:name=value` form.

A dependency-light build uses precisely the same numerical code:

```bash
cmake -S ~/OctotigerII -B ~/OctotigerII/build-serial \
  -DCMAKE_BUILD_TYPE=Release \
  -DOCTOTIGERII_WITH_HPX=OFF -DOCTOTIGERII_WITH_SILO=OFF
cmake --build ~/OctotigerII/build-serial -j
ctest --test-dir ~/OctotigerII/build-serial --output-on-failure
~/OctotigerII/build-serial/octotigerII --config="$HOME/OctotigerII/examples/collapse.ini"
```

The `release` and `serial` CMake presets express these two configurations.
An optional Python 3 interpreter enables one application/output check; it is
not a runtime dependency. `OCTOTIGERII_BUILD_TESTS=OFF` disables the test suite.
`cmake --install` installs the executable and example INIs.

## Examples

Every example runs as a normal problem through the same executable:

| INI / `problem.name` | Physics | Default global mesh |
| --- | --- | --- |
| `sod` | 1D ideal-gas shock tube | 64 |
| `kelvin-helmholtz` | 2D periodic shear instability | 32 × 32 |
| `gravity-sphere` | Isolated constant-density sphere, static solve | 8³ |
| `gravity-gaussian` | Isolated Gaussian density, static solve | 8³ |
| `streaming` | 1D periodic M1 streaming Gaussian | 64 |
| `radiation-pulse` | 2D initially isotropic M1 pulse | 16² |
| `collapse` | 3D gas and self-gravity | 8³ |

These are small examples, not resolution recommendations. Raise `mesh.cells`
or `mesh.level` to resolve the solution. `mesh.cells` is the number of cells
per block per active axis and must be a power of two. There are
`2^mesh.level` blocks along each active axis. All leaves have the same level;
there is no hidden extra refinement branch. Radiation examples accept
`mesh.ndim=1`, `2`, or `3`; gravity requires three dimensions.

Options are dotted `key=value` entries in an INI or `--key=value` arguments.
CLI values override INIs independently of argument order. `--help` lists
all supported options. Unknown options and unsupported combinations fail
explicitly; old option aliases are intentionally absent. Problem selection
sets the required physics, so no contradictory hydro/radiation/gravity enable
switches are needed. Static gravity examples use `runtime.stop_time=0`.

## Units and output

All inputs, stored physical quantities, and outputs use cgs: cm, s, g,
erg/cm³, cm²/s² for potential, and cm/s² for acceleration. Hydro momentum is
momentum density. Radiation flux is erg/(cm² s). The numerical radiation state
uses `Q=F/c`; output always converts with physical `c=2.99792458e10 cm/s`.
`radiation.light_speed_ratio` changes transport speed only.
`G=6.67430e-8 cm³/(g s²)` is fixed in the gravity module.

The examples default to CSV output with no visualization dependency. Each
snapshot includes cell centers, cell width, time, and all enabled fields.
`output.format=silo` writes dimensional multimeshes and zone-centered
variables plus `frames.visit`. Silo creation uses `DB_CLOBBER`.

Each example has its own relative `output.directory`; paths are relative to
the working directory. Frame zero and the final frame are always written;
`output.every` selects the intervening step cadence. `diagnostics.csv` records
every step. Rerunning in the same directory overwrites matching filenames;
choose a new directory to retain an earlier run. Output is visualization,
not a restart checkpoint. `output.enabled=off` creates no output files.

The 1D/2D integrated diagnostics are per unit inactive length/area. For
example, 1D mass is in g/cm²; three-dimensional mass is in g. Open boundaries
can change integrated quantities. No conserved gravitational field energy is
reported.

## Small architecture

- `Subgrid` owns the geometry, enabled fields, initial data, CFL calculation,
  halo-backed transport advance, and local gravity kick.
- `Runtime` owns the fixed directory, placement, and phase barriers. Its HPX
  `SubgridComponent` is a thin wrapper around the same `Subgrid` used by the
  serial build. There is no `node_server` or `node_client` and no old grid.
- `simulation.cpp` coordinates the shared timestep and the enabled physics.
- `problems.cpp` contains the seven initial-condition recipes.
- The numerical libraries are `OctotigerII::Mesh`, `Physics`, `Hydro`,
  `Radiation`, `Gravity`, and `Subgrid`. They do not include HPX headers.
- `OctotigerII::Application` supplies execution, orchestration, and output.

Hydro uses the existing modular unsplit MUSCL–Hancock integrator, primitive
PLM reconstruction, HLLC with HLL fallback, and positivity limiting. M1 uses
the same finite-volume scaffold, the existing Skinner–Ostriker closure and
HLL transport, and realizability limiting. Gas and radiation share machinery,
but there is currently no gas–radiation exchange solver or combined example.

Gravity uses a new cell-octree driver around the recovered compact Cartesian
plane-wave FMM operators. Moments and locals retain `(p+1)²` real coefficients
for `p=3,4,5`. An upward M2M pass, symmetric pair traversal with diagonal M2L,
and a downward L2L pass produce potential and acceleration. The opening test
is `cell_width / separation < gravity.opening_angle`, with the angle strictly
below `1/sqrt(3)`. Near leaf pairs use Newtonian direct summation. Cells are
point masses at their centers and self contributions are omitted. No old
Cartesian Taylor gravity backend, torque correction, or alternate solver is
included. Direct summation in the tests is only an independent reference.

The collapse example uses a half gravity kick, a hydro transport step, a new
gravity solve, and a second half kick. The kick updates momentum and gas
kinetic energy while preserving internal energy. There are no gravitational
energy fluxes, `dphi/dt` evolution, or exact total-energy conservation claims.

## Scope of this first version

Included: CPU numerics, 1D/2D/3D fixed tiled meshes, periodic/outflow transport,
isolated 3D gravity, shared timesteps, HPX component execution, CSV and Silo.

Omitted: CUDA, HIP, Kokkos, Vc, CPPuddle, Unitiger, the old FMM, old problems
and test harnesses, SCF, binary-star setup, rotating frames, species/degenerate
EOS, radiation opacities/coupling/subcycling, dynamic AMR, shadow hierarchies,
load balancing, checkpoints, and old command-line compatibility.

The runtime currently exchanges complete immutable snapshot directories. Its
communication grows quadratically with block count. The gravity octree is
assembled and solved on the coordinating locality. Thus this initial version
is a runnable foundation, not the final scalable distributed implementation.
Nearest-neighbor exchange, distributed FMM passes, and adaptive refinement
should be added to their respective modules rather than rebuilding a large
node class. Fixed-mesh state is the only supported runtime topology.

See [VALIDATION.md](VALIDATION.md) for the checks actually performed and
[PROVENANCE.md](PROVENANCE.md) for the source selection.
