# Source inventory — completed reorganization

This is the complete source inventory with the implemented changes highlighted
against the previous layout. **NEW** marks eight added files; **CHANGED** marks
files whose contents or responsibilities changed. Unmarked rows retain their
previous path and responsibility. The unmarked [current inventory](source-inventory.md)
is the navigation reference; the [reorganization plan](source-reorganization-plan.md)
records the implementation sequence and checks.

**REMOVED:** `octotigerII/subgrid/exchange.hpp`, the unused alternative snapshot
exchange header that also defined `FieldFluxPacket`. The retained packet header
has the same contents and is marked **CHANGED** for its now unique role. All
other previously listed source and input files remain. The current list contains
161 program source files, two build/documentation configuration files, and 24
runnable inputs: 187 file rows in all.

## Build definitions and standalone tools

| File | Contents |
|---|---|
| [`CMakeLists.txt`](../CMakeLists.txt) | Project options, dependencies, three dimension builds, install rules, and documentation target. |
| [`CMakePresets.json`](../CMakePresets.json) | Named HPX and serial configure, build, and test presets. |
| [`Doxyfile`](../Doxyfile) | Doxygen settings for the generated developer reference. |
| [`cmake/Problems.cmake`](../cmake/Problems.cmake) | Reads problem manifests, checks build compatibility, and generates the problem registry. |
| [`cmake/dimension/CMakeLists.txt`](../cmake/dimension/CMakeLists.txt) | **CHANGED:** Per-dimension libraries, module dependencies, executable, HPX linkage, and tests; includes the new AMR and runtime implementations in their existing targets. |
| [`cmake/buildConfig.hpp.in`](../cmake/buildConfig.hpp.in) | Template for compile-time dimension and enabled physics flags. |
| [`cmake/problems.cpp.in`](../cmake/problems.cpp.in) | Template for the compiled runtime problem dispatch table. |
| [`cmake/PatchApex.cmake`](../cmake/PatchApex.cmake) | Applies the pinned APEX header fix needed by the HPX dependency. |
| [`build.sh`](../build.sh) | Local build helper and dependency configuration. |
| [`arc.sh`](../arc.sh) | Creates a source archive including Git history. |
| [`docs.sh`](../docs.sh) | Builds the developer reference independently of the application. |
| [`profile.sh`](../profile.sh) | Runs a configured executable under profiling tools. |
| [`gravity_sphere_sweep.cpp`](../gravity_sphere_sweep.cpp) | Standalone gravity-sphere parameter sweep and plotting driver. |

## Problem modules

Every `problem.cpp` supplies some combination of defaults, validation, initial
conditions, boundary state, and analytic reference. Each adjacent
`CMakeLists.txt` declares the problem name, required modules, and dimensions.

| File | Contents |
|---|---|
| [`problem/gravity_tests/Gaussian/CMakeLists.txt`](../problem/gravity_tests/Gaussian/CMakeLists.txt) | Registers the 3D gravity Gaussian problem. |
| [`problem/gravity_tests/Gaussian/problem.cpp`](../problem/gravity_tests/Gaussian/problem.cpp) | Initializes a smooth Gaussian density for gravity accuracy tests. |
| [`problem/gravity_tests/Sphere/CMakeLists.txt`](../problem/gravity_tests/Sphere/CMakeLists.txt) | Registers the 3D gravity sphere problem. |
| [`problem/gravity_tests/Sphere/problem.cpp`](../problem/gravity_tests/Sphere/problem.cpp) | Initializes an isolated uniform sphere and its gravity reference. |
| [`problem/hydro_tests/KelvinHelmholtz/CMakeLists.txt`](../problem/hydro_tests/KelvinHelmholtz/CMakeLists.txt) | Registers the 2D/3D Kelvin–Helmholtz problem. |
| [`problem/hydro_tests/KelvinHelmholtz/problem.cpp`](../problem/hydro_tests/KelvinHelmholtz/problem.cpp) | Initializes shearing layers and their perturbation. |
| [`problem/hydro_tests/RayleighTaylor/CMakeLists.txt`](../problem/hydro_tests/RayleighTaylor/CMakeLists.txt) | Registers the 3D Rayleigh–Taylor problem. |
| [`problem/hydro_tests/RayleighTaylor/problem.cpp`](../problem/hydro_tests/RayleighTaylor/problem.cpp) | Initializes stratified layers under constant external acceleration. |
| [`problem/hydro_tests/Sod/CMakeLists.txt`](../problem/hydro_tests/Sod/CMakeLists.txt) | Registers the 1D/2D/3D Sod problem. |
| [`problem/hydro_tests/Sod/problem.cpp`](../problem/hydro_tests/Sod/problem.cpp) | Initializes the planar shock tube and reference data. |
| [`problem/radiation_tests/RadiationPulse/CMakeLists.txt`](../problem/radiation_tests/RadiationPulse/CMakeLists.txt) | Registers the 1D/2D/3D radiation pulse. |
| [`problem/radiation_tests/RadiationPulse/problem.cpp`](../problem/radiation_tests/RadiationPulse/problem.cpp) | Initializes a localized radiation pulse. |
| [`problem/radiation_tests/Streaming/CMakeLists.txt`](../problem/radiation_tests/Streaming/CMakeLists.txt) | Registers the 1D/2D/3D radiation streaming problem. |
| [`problem/radiation_tests/Streaming/problem.cpp`](../problem/radiation_tests/Streaming/problem.cpp) | Initializes directed radiation transport and reference state. |
| [`problem/science/Collapse/CMakeLists.txt`](../problem/science/Collapse/CMakeLists.txt) | Registers 3D self-gravitating collapse. |
| [`problem/science/Collapse/problem.cpp`](../problem/science/Collapse/problem.cpp) | Initializes a collapsing gas sphere. |
| [`problem/science/Polytrope/CMakeLists.txt`](../problem/science/Polytrope/CMakeLists.txt) | Registers the 3D self-gravitating polytrope. |
| [`problem/science/Polytrope/problem.cpp`](../problem/science/Polytrope/problem.cpp) | Scales a Lane–Emden sphere into hydrostatic initial data. |
| [`problem/science/RotatingStar/CMakeLists.txt`](../problem/science/RotatingStar/CMakeLists.txt) | Registers the 3D rotating star. |
| [`problem/science/RotatingStar/problem.cpp`](../problem/science/RotatingStar/problem.cpp) | Scales and interpolates the original oblate SCF equilibrium into conserved fields. |
| [`problem/science/RotatingStar/equilibrium.inc`](../problem/science/RotatingStar/equilibrium.inc) | Embedded 100×100 positive-quadrant density and internal-energy table. |
| [`problem/science/RotatingStar/import_equilibrium.py`](../problem/science/RotatingStar/import_equilibrium.py) | Recreates and checks the embedded table from the pinned source dataset. |

## Public headers: geometry, units, and configuration

| File | Contents |
|---|---|
| [`octotigerII/math/Real.hpp`](../octotigerII/math/Real.hpp) | Scalar floating-point type and mathematical constants. |
| [`octotigerII/math/Vector.hpp`](../octotigerII/math/Vector.hpp) | Fixed-size vector arithmetic. |
| [`octotigerII/math/FpeGuard.hpp`](../octotigerII/math/FpeGuard.hpp) | Scoped floating-point exception checks around numerical kernels. |
| [`octotigerII/units/cgs.hpp`](../octotigerII/units/cgs.hpp) | Compile-time CGS quantities and dimension-safe arithmetic. |
| [`octotigerII/units/constants.hpp`](../octotigerII/units/constants.hpp) | Physical constants in CGS units. |
| [`octotigerII/units/state.hpp`](../octotigerII/units/state.hpp) | Heterogeneous conserved-state tuple and componentwise operations. |
| [`octotigerII/mesh.hpp`](../octotigerII/mesh.hpp) | Cartesian indexing, block locations, patch geometry, and time metadata. |
| [`octotigerII/config.hpp`](../octotigerII/config.hpp) | Runtime options, parsing API, validation, and option serialization. |
| [`octotigerII/problems.hpp`](../octotigerII/problems.hpp) | Problem registry API and initialization/boundary/reference hooks. |
| [`octotigerII/problems/laneEmden.hpp`](../octotigerII/problems/laneEmden.hpp) | Lane–Emden solution and polytrope sampling API. |
| [`octotigerII/problems/rotatingStar.hpp`](../octotigerII/problems/rotatingStar.hpp) | Rotating-star SCF table access and physical scaling API. |
| [`octotigerII/profiling.hpp`](../octotigerII/profiling.hpp) | Timing regions and HPX/APEX profiling adapters. |

## Public headers: transport and physics

| File | Contents |
|---|---|
| [`octotigerII/physics/finiteVolume.hpp`](../octotigerII/physics/finiteVolume.hpp) | Generic unsplit MUSCL–Hancock reconstruction, flux integration, and slope limiting. |
| [`octotigerII/physics/boundary.hpp`](../octotigerII/physics/boundary.hpp) | Face boundary types, ghost-state transforms, and rotating free-boundary rules. |
| [`octotigerII/physics/frame.hpp`](../octotigerII/physics/frame.hpp) | Grid/inertial rotations, rigid mesh velocity, and rotation timestep limit. |
| [`octotigerII/hydro/hydroSystem.hpp`](../octotigerII/hydro/hydroSystem.hpp) | Euler conserved/primitive states, HLLC fluxes, and positivity safeguards. |
| [`octotigerII/hydro/dualEnergy.hpp`](../octotigerII/hydro/dualEnergy.hpp) | Entropy auxiliary conversion and synchronization rules. |
| [`octotigerII/radiation/m1.hpp`](../octotigerII/radiation/m1.hpp) | M1 closure and HLL radiation solver in scaled calculation variables. |
| [`octotigerII/radiation/radiationTransport.hpp`](../octotigerII/radiation/radiationTransport.hpp) | Physical radiation state, moving-face flux adapter, and timestep estimate. |
| [`octotigerII/composition/species.hpp`](../octotigerII/composition/species.hpp) | Material species, elements, mixtures, and fraction definitions. |
| [`octotigerII/composition/transport.hpp`](../octotigerII/composition/transport.hpp) | Species and tracer transport driven by the hydro mass flux. |

## Public headers: gravity

| File | Contents |
|---|---|
| [`octotigerII/gravity/diagonal/fmm.hpp`](../octotigerII/gravity/diagonal/fmm.hpp) | Harmonic-equivalent moments, local expansions, translations, and direct pair kernels. |
| [`octotigerII/gravity/solver.hpp`](../octotigerII/gravity/solver.hpp) | Standalone serial uniform-grid gravity solver and work statistics. |
| [`octotigerII/gravity/fieldSolver.hpp`](../octotigerII/gravity/fieldSolver.hpp) | Distributed field-solve requests and the uniform/adaptive solver facade. |
| [`octotigerII/gravity/adaptiveFieldSolver.hpp`](../octotigerII/gravity/adaptiveFieldSolver.hpp) | Adaptive sparse-octree solver interface. |
| [`octotigerII/gravity/gravityFields.hpp`](../octotigerII/gravity/gravityFields.hpp) | Potential and acceleration field state and component ordering. |
| [`octotigerII/gravity/boundary.hpp`](../octotigerII/gravity/boundary.hpp) | Valid gravity image-boundary combinations. |
| [`octotigerII/gravity/images.hpp`](../octotigerII/gravity/images.hpp) | Reflected/periodic image geometry and interaction lists. |
| [`octotigerII/gravity/ewald.hpp`](../octotigerII/gravity/ewald.hpp) | Ewald corrections for one-, two-, and three-axis periodicity. |
| [`octotigerII/gravity/fluxWork.hpp`](../octotigerII/gravity/fluxWork.hpp) | Gravitational energy work from accepted finite-volume mass fluxes. |
| [`octotigerII/gravity/rotationWork.hpp`](../octotigerII/gravity/rotationWork.hpp) | Balanced local gravity work from coordinate-weighted solves on a rotating grid. |

## Public headers: mesh hierarchy and storage

| File | Contents |
|---|---|
| [`octotigerII/subgrid/subgrid.hpp`](../octotigerII/subgrid/subgrid.hpp) | Temporary interior snapshots for initialization, diagnostics, and output. |
| [`octotigerII/subgrid/topology.hpp`](../octotigerII/subgrid/topology.hpp) | Block geometry, halo plans, reflux plans, and gravity-work face plans. |
| [`octotigerII/subgrid/view.hpp`](../octotigerII/subgrid/view.hpp) | Non-owning patch views joining interiors with compact halo data. |
| [`octotigerII/subgrid/fluxPacket.hpp`](../octotigerII/subgrid/fluxPacket.hpp) | **CHANGED:** Time-tagged face-flux packet values; now the sole FieldFluxPacket definition after removal of the unused exchange header. |
| [`octotigerII/refinement/criteria.hpp`](../octotigerII/refinement/criteria.hpp) | AMR criteria inputs, thresholds, and callback interface. |
| [`octotigerII/amr/hierarchy.hpp`](../octotigerII/amr/hierarchy.hpp) | **CHANGED:** Covered coarse-state evolution, reconstruction, conservative transfer, and startup mesh API. |
| [`octotigerII/amr/regridSelection.hpp`](../octotigerII/amr/regridSelection.hpp) | **NEW:** RegridResult and selectMesh declarations for refinement and coarsening decisions. |
| [`octotigerII/amr/interpolation.hpp`](../octotigerII/amr/interpolation.hpp) | Conservative interpolation and restriction helpers. |
| [`octotigerII/storage/layout.hpp`](../octotigerII/storage/layout.hpp) | Record ranges and partition layout independent of mesh geometry. |
| [`octotigerII/storage/buffer.hpp`](../octotigerII/storage/buffer.hpp) | Typed buffers and serial/HPX future abstractions. |
| [`octotigerII/storage/partition.hpp`](../octotigerII/storage/partition.hpp) | HPX storage component, range transfer, and serial fallback. |
| [`octotigerII/storage/field.hpp`](../octotigerII/storage/field.hpp) | Named fields, banked handles, and typed access to partitions. |
| [`octotigerII/storage/columns.hpp`](../octotigerII/storage/columns.hpp) | Structure-of-arrays wrappers for heterogeneous states. |
| [`octotigerII/storage/registry.hpp`](../octotigerII/storage/registry.hpp) | Application field directory/schema over generic storage. |
| [`octotigerII/storage/types.def`](../octotigerII/storage/types.def) | Distinct wire types registered for HPX serialization. |

## Public headers: runtime and verification

| File | Contents |
|---|---|
| [`octotigerII/runtime.hpp`](../octotigerII/runtime.hpp) | Published-state API, locality scheduling statistics, staged advance/regrid/gravity calls. |
| [`octotigerII/simulation.hpp`](../octotigerII/simulation.hpp) | Top-level run loop and conserved diagnostics interface. |
| [`octotigerII/conservation.hpp`](../octotigerII/conservation.hpp) | Volume-integral diagnostics and accumulated boundary transport budgets. |
| [`octotigerII/output.hpp`](../octotigerII/output.hpp) | Silo frame and conservation CSV output API. |
| [`octotigerII/verification/analytic.hpp`](../octotigerII/verification/analytic.hpp) | Problem reference solutions and volume-weighted error comparison. |
| [`octotigerII/verification/directGravity.hpp`](../octotigerII/verification/directGravity.hpp) | Independent point-pair gravity sampling interface. |

## Compiled implementations

| File | Contents |
|---|---|
| [`src/main.cpp`](../src/main.cpp) | CLI/HPX startup, run invocation, progress table, and completion report. |
| [`src/config.cpp`](../src/config.cpp) | Option definitions, INI/command-line merge, defaults, and cross-option validation. |
| [`src/mesh.cpp`](../src/mesh.cpp) | Mesh indexing, block ancestry, cell geometry, and time-state methods. |
| [`src/physics/finiteVolume.cpp`](../src/physics/finiteVolume.cpp) | Returns the displayed name of the finite-volume scheme. |
| [`src/hydro/hydroSystem.cpp`](../src/hydro/hydroSystem.cpp) | Hydro state conversions, HLLC fluxes, characteristic speeds, and limiting. |
| [`src/radiation/radiationTransport.cpp`](../src/radiation/radiationTransport.cpp) | Radiation M1 fluxes, HLL transport, and realizability limiting. |
| [`src/composition/species.cpp`](../src/composition/species.cpp) | Periodic-table data, species parser, validation, and initial fractions. |
| [`src/problems/laneEmden.cpp`](../src/problems/laneEmden.cpp) | Numerical Lane–Emden integration and interpolation. |
| [`src/subgrid.cpp`](../src/subgrid.cpp) | Snapshot allocation and initialization support. |
| [`src/topology.cpp`](../src/topology.cpp) | Leaf layout, halo source plans, fine/coarse reflux, and gravity-work face maps. |
| [`src/refinement/criteria.cpp`](../src/refinement/criteria.cpp) | Density/mass/error refinement decisions. |
| [`src/amr/hierarchy.cpp`](../src/amr/hierarchy.cpp) | **CHANGED:** Evolved coarse shadows, reconstruction, conservative state transfer, and combined gas/gravity energy transfer. |
| [`src/amr/regridSelection.cpp`](../src/amr/regridSelection.cpp) | **NEW:** Refinement/coarsening decisions, geometric buffers, periodic neighbors, and 2:1 balance. |
| [`src/amr/initialization.cpp`](../src/amr/initialization.cpp) | **CHANGED:** Startup AMR refinement and lookahead loop, initial-condition resampling, and timestep/signal-speed estimates; uses the separate mesh-selection API. |
| [`src/gravity/diagonal/fmm.cpp`](../src/gravity/diagonal/fmm.cpp) | Harmonic FMM coefficient algebra, translations, and interaction operators. |
| [`src/gravity/solver.cpp`](../src/gravity/solver.cpp) | Standalone serial uniform-grid FMM reference solve. |
| [`src/gravity/fieldSolver.cpp`](../src/gravity/fieldSolver.cpp) | Distributed uniform-grid FMM plus facade selecting uniform or adaptive solve. |
| [`src/gravity/adaptiveFieldSolver.cpp`](../src/gravity/adaptiveFieldSolver.cpp) | Distributed sparse adaptive FMM traversal and field publication. |
| [`src/gravity/gravityFields.cpp`](../src/gravity/gravityFields.cpp) | Translation unit for the gravity field type. |
| [`src/gravity/images.cpp`](../src/gravity/images.cpp) | Periodic/reflected gravity image offsets and pair lists. |
| [`src/gravity/ewald.cpp`](../src/gravity/ewald.cpp) | Periodic Newtonian Ewald sums and corrections. |
| [`src/storage.cpp`](../src/storage.cpp) | Storage component registration and typed HPX actions. |
| [`src/runtime.cpp`](../src/runtime.cpp) | **CHANGED:** Runtime construction/destruction, snapshot gathering, backend identity, and public bookkeeping queries. |
| [`src/runtime/internal.hpp`](../src/runtime/internal.hpp) | **NEW:** Shared private stage types, Runtime implementation state, LocalExecutor declarations, numerical/storage templates, and HPX action declarations. |
| [`src/runtime/localExecutor.cpp`](../src/runtime/localExecutor.cpp) | **NEW:** Locality worker queues, halo assembly, block transport and source operations, and HPX component/action registration. |
| [`src/runtime/stages.cpp`](../src/runtime/stages.cpp) | **NEW:** Stage dispatch, draining of task results and errors, completion checks, and scheduling statistics. |
| [`src/runtime/transport.cpp`](../src/runtime/transport.cpp) | **NEW:** Transport timestep selection, AMR subcycling, accepted-bank publication, and failed-step restoration. |
| [`src/runtime/gravity.cpp`](../src/runtime/gravity.cpp) | **NEW:** Gravity rungs, force stages, endpoint energy work, gravity-field publication, and coupled-step rollback. |
| [`src/runtime/regrid.cpp`](../src/runtime/regrid.cpp) | **NEW:** Mesh-selection loop, replacement of topology/field directories and executors, and combined-energy recovery after regridding. |
| [`src/simulation.cpp`](../src/simulation.cpp) | Run loop, timestep selection, global diagnostics, and advance composition. |
| [`src/output.cpp`](../src/output.cpp) | Silo mesh/field writing and conservation CSV calculations. |
| [`src/verification/solutions.cpp`](../src/verification/solutions.cpp) | Analytic reference formulas supplied to problem hooks. |
| [`src/verification/analytic.cpp`](../src/verification/analytic.cpp) | Reference sampling, numerical comparison, and error reports. |
| [`src/verification/directGravity.cpp`](../src/verification/directGravity.cpp) | Independent direct gravity samples and error statistics. |

## Test source and harnesses

| File | Contents |
|---|---|
| [`tests/CMakeLists.txt`](../tests/CMakeLists.txt) | GoogleTest targets, discovery, labels, and distributed test registration. |
| [`tests/googleTestMain.cpp`](../tests/googleTestMain.cpp) | Small shared GoogleTest executable entry point. |
| [`tests/googleTestMain.hpp`](../tests/googleTestMain.hpp) | GoogleTest startup adapter. |
| [`tests/runtimeMain.hpp`](../tests/runtimeMain.hpp) | Starts tests under HPX or the serial runtime. |
| [`tests/testSupport.hpp`](../tests/testSupport.hpp) | Common fixtures, option parsing, temp paths, and assertions. |
| [`tests/test_problem/problem.cpp`](../tests/test_problem/problem.cpp) | Simple combined-physics problem used by application tests. |
| [`tests/application.py`](../tests/application.py) | Launches registered problems and checks CLI/module contracts. |
| [`tests/distributed.py`](../tests/distributed.py) | Launches test executables across local TCP HPX localities. |
| [`tests/run_matrix.py`](../tests/run_matrix.py) | Configures/builds/tests dimensions and optional physics module combinations. |
| [`tests/amrChecks.cpp`](../tests/amrChecks.cpp) | **CHANGED:** AMR hierarchy, conservative transfer, startup refinement, mesh-selection behavior, and regrid conservation checks; uses the separate selection API. |
| [`tests/analyticChecks.cpp`](../tests/analyticChecks.cpp) | Analytic solutions and comparison metrics. |
| [`tests/boundaryChecks.cpp`](../tests/boundaryChecks.cpp) | Physical ghost, image, and rotating outflow boundary behavior. |
| [`tests/conservationChecks.cpp`](../tests/conservationChecks.cpp) | Integrated mass/energy/momentum budgets and output records. |
| [`tests/directGravityChecks.cpp`](../tests/directGravityChecks.cpp) | Sampled direct-pair gravity verification. |
| [`tests/dualEnergyChecks.cpp`](../tests/dualEnergyChecks.cpp) | Entropy auxiliary, synchronization, and AMR transport. |
| [`tests/ewaldChecks.cpp`](../tests/ewaldChecks.cpp) | One-, two-, and three-periodic-axis Ewald kernels. |
| [`tests/finiteVolumeChecks.cpp`](../tests/finiteVolumeChecks.cpp) | Hydro/radiation updates and finite-volume invariants. |
| [`tests/gravityChecks.cpp`](../tests/gravityChecks.cpp) | Standalone serial FMM behavior and validation. |
| [`tests/gravityEnergyChecks.cpp`](../tests/gravityEnergyChecks.cpp) | Gravitational energy bookkeeping and reciprocity. |
| [`tests/gravityFluxWorkChecks.cpp`](../tests/gravityFluxWorkChecks.cpp) | Gravity work from numerical face mass fluxes. |
| [`tests/gravityParallelChecks.cpp`](../tests/gravityParallelChecks.cpp) | Distributed field-solver agreement, work partitioning, and image cases. |
| [`tests/gravityTimeIntegrationChecks.cpp`](../tests/gravityTimeIntegrationChecks.cpp) | Hierarchical/conventional gravity rungs and conservation under subcycling. |
| [`tests/hydroChecks.cpp`](../tests/hydroChecks.cpp) | Euler state, flux, and admissibility tests. |
| [`tests/laneEmdenChecks.cpp`](../tests/laneEmdenChecks.cpp) | Lane–Emden values and polytrope interpolation. |
| [`tests/limiterChecks.cpp`](../tests/limiterChecks.cpp) | Spatial slope limiter behavior. |
| [`tests/mathChecks.cpp`](../tests/mathChecks.cpp) | Fixed-size vectors and state arithmetic. |
| [`tests/meshChecks.cpp`](../tests/meshChecks.cpp) | Indexing, block geometry, and time metadata. |
| [`tests/numericalChecks.cpp`](../tests/numericalChecks.cpp) | Numerical convergence and physics reference comparisons. |
| [`tests/optionsChecks.cpp`](../tests/optionsChecks.cpp) | INI/CLI parsing, defaults, and invalid option rejection. |
| [`tests/partialGravityChecks.cpp`](../tests/partialGravityChecks.cpp) | Selected-source/target gravity solves for timestep rungs. |
| [`tests/polytropeChecks.cpp`](../tests/polytropeChecks.cpp) | Polytrope initialization and equilibrium properties. |
| [`tests/radiationChecks.cpp`](../tests/radiationChecks.cpp) | M1 radiation states, fluxes, and bounds. |
| [`tests/rayleighTaylorChecks.cpp`](../tests/rayleighTaylorChecks.cpp) | Rayleigh–Taylor setup and evolution properties. |
| [`tests/rotatingGravityChecks.cpp`](../tests/rotatingGravityChecks.cpp) | Balanced rotation work, mode conservation, AMR transitions, and temporal order. |
| [`tests/rotatingStarChecks.cpp`](../tests/rotatingStarChecks.cpp) | SCF data, scaling, shape, initialization, and startup refinement. |
| [`tests/rotatingTransportChecks.cpp`](../tests/rotatingTransportChecks.cpp) | Moving hydro/radiation faces, relative CFL, and temporal order. |
| [`tests/serializationChecks.cpp`](../tests/serializationChecks.cpp) | HPX wire round trips for options, geometry, and states. |
| [`tests/siloChecks.cpp`](../tests/siloChecks.cpp) | Silo output variables, geometry, and rotating mesh coordinates. |
| [`tests/speciesChecks.cpp`](../tests/speciesChecks.cpp) | Species parser, partial-density transport, and field storage. |
| [`tests/storageChecks.cpp`](../tests/storageChecks.cpp) | Field partitions, banks, handles, and transfer semantics. |
| [`tests/timeRefinementChecks.cpp`](../tests/timeRefinementChecks.cpp) | AMR transport subcycling, coarse shadows, and flux registers. |
| [`tests/unitsChecks.cpp`](../tests/unitsChecks.cpp) | CGS quantity dimensions and arithmetic. |

## Runnable input files

These are configuration, not implementation, but belong in a code-organization
map because they define the examples and defaults users run. Every problem
module has a sibling `inputs` file. These ten files mirror the ten problems
above and set their runnable defaults:

| File | Contents |
|---|---|
| [`problem/gravity_tests/Gaussian/inputs`](../problem/gravity_tests/Gaussian/inputs) | Gravity Gaussian defaults. |
| [`problem/gravity_tests/Sphere/inputs`](../problem/gravity_tests/Sphere/inputs) | Uniform gravity sphere defaults. |
| [`problem/hydro_tests/KelvinHelmholtz/inputs`](../problem/hydro_tests/KelvinHelmholtz/inputs) | Kelvin–Helmholtz setup. |
| [`problem/hydro_tests/RayleighTaylor/inputs`](../problem/hydro_tests/RayleighTaylor/inputs) | Rayleigh–Taylor setup. |
| [`problem/hydro_tests/Sod/inputs`](../problem/hydro_tests/Sod/inputs) | Sod shock tube setup. |
| [`problem/radiation_tests/RadiationPulse/inputs`](../problem/radiation_tests/RadiationPulse/inputs) | Radiation pulse setup. |
| [`problem/radiation_tests/Streaming/inputs`](../problem/radiation_tests/Streaming/inputs) | Radiation streaming setup. |
| [`problem/science/Collapse/inputs`](../problem/science/Collapse/inputs) | Collapse setup. |
| [`problem/science/Polytrope/inputs`](../problem/science/Polytrope/inputs) | Polytrope setup. |
| [`problem/science/RotatingStar/inputs`](../problem/science/RotatingStar/inputs) | Oblate rotating-star setup. |
| [`tests/test_problem/inputs`](../tests/test_problem/inputs) | Small combined-physics application test setup. |
| [`examples/collapse.ini`](../examples/collapse.ini) | User-facing collapse run. |
| [`examples/gravity-gaussian.ini`](../examples/gravity-gaussian.ini) | User-facing Gaussian gravity run. |
| [`examples/gravity-sphere.ini`](../examples/gravity-sphere.ini) | User-facing uniform sphere gravity run. |
| [`examples/kelvin-helmholtz.ini`](../examples/kelvin-helmholtz.ini) | Kelvin–Helmholtz run. |
| [`examples/kelvin-helmholtz-time-refinement.ini`](../examples/kelvin-helmholtz-time-refinement.ini) | Kelvin–Helmholtz run with AMR time refinement. |
| [`examples/mass-fractions.ini`](../examples/mass-fractions.ini) | Material fraction transport example. |
| [`examples/polytrope.ini`](../examples/polytrope.ini) | Self-gravitating polytrope example. |
| [`examples/radiation-pulse.ini`](../examples/radiation-pulse.ini) | Radiation pulse example. |
| [`examples/rayleigh-taylor-amr.ini`](../examples/rayleigh-taylor-amr.ini) | Refined Rayleigh–Taylor example. |
| [`examples/rotating-star.ini`](../examples/rotating-star.ini) | Corotating star with hierarchical gravity rungs. |
| [`examples/rotating-star-inertial.ini`](../examples/rotating-star-inertial.ini) | Same star on a stationary grid. |
| [`examples/sod.ini`](../examples/sod.ini) | Sod shock tube example. |
| [`examples/streaming.ini`](../examples/streaming.ini) | Streaming radiation example. |

