# Source file inventory

This is a map of the maintained program source as of 2026-09-26, after the
source reorganization and radiation-coupling work. It covers maintained
C++/Python/shell/CMake source files in the root, `bin/`, `cmake/`,
`octotigerII/`, `src/`, and `tests/`, plus build/documentation configuration
and runnable input files. It excludes prose documentation, generated build trees, vendored
packages, output, caches, and binaries. Paths are relative to the repository
root. `octotigerII/` holds public declarations and inline implementations;
`src/` holds compiled implementations and the private runtime header.

The [highlighted change map](source-inventory-proposed.md) records the completed
changes against the previous layout. The [reorganization plan](source-reorganization-plan.md)
explains the chosen boundaries and validation sequence.

## Build definitions and standalone tools

| File | Contents |
|---|---|
| [`CMakeLists.txt`](../CMakeLists.txt) | Project options, dependencies, three dimension builds, install rules, and documentation target. |
| [`CMakePresets.json`](../CMakePresets.json) | Named HPX and serial configure, build, and test presets. |
| [`Doxyfile`](../Doxyfile) | Doxygen settings for the generated developer reference. |
| [`cmake/Problems.cmake`](../cmake/Problems.cmake) | Reads problem manifests, checks build compatibility, and generates the problem registry. |
| [`cmake/dimension/CMakeLists.txt`](../cmake/dimension/CMakeLists.txt) | Per-dimension libraries, module dependencies, executable, HPX linkage, and tests. |
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
| [`bin/gravity_tests/Gaussian/CMakeLists.txt`](../bin/gravity_tests/Gaussian/CMakeLists.txt) | Registers the 3D gravity Gaussian problem. |
| [`bin/gravity_tests/Gaussian/problem.cpp`](../bin/gravity_tests/Gaussian/problem.cpp) | Initializes a smooth Gaussian density for gravity accuracy tests. |
| [`bin/gravity_tests/Sphere/CMakeLists.txt`](../bin/gravity_tests/Sphere/CMakeLists.txt) | Registers the 3D gravity sphere problem. |
| [`bin/gravity_tests/Sphere/problem.cpp`](../bin/gravity_tests/Sphere/problem.cpp) | Initializes an isolated uniform sphere and its gravity reference. |
| [`bin/hydro_tests/KelvinHelmholtz/CMakeLists.txt`](../bin/hydro_tests/KelvinHelmholtz/CMakeLists.txt) | Registers the 2D/3D Kelvin–Helmholtz problem. |
| [`bin/hydro_tests/KelvinHelmholtz/problem.cpp`](../bin/hydro_tests/KelvinHelmholtz/problem.cpp) | Initializes shearing layers and their perturbation. |
| [`bin/hydro_tests/RayleighTaylor/CMakeLists.txt`](../bin/hydro_tests/RayleighTaylor/CMakeLists.txt) | Registers the 3D Rayleigh–Taylor problem. |
| [`bin/hydro_tests/RayleighTaylor/problem.cpp`](../bin/hydro_tests/RayleighTaylor/problem.cpp) | Initializes stratified layers under constant external acceleration. |
| [`bin/hydro_tests/Sod/CMakeLists.txt`](../bin/hydro_tests/Sod/CMakeLists.txt) | Registers the 1D/2D/3D Sod problem. |
| [`bin/hydro_tests/Sod/problem.cpp`](../bin/hydro_tests/Sod/problem.cpp) | Initializes the planar shock tube and reference data. |
| [`bin/radiation_tests/RadiationPulse/CMakeLists.txt`](../bin/radiation_tests/RadiationPulse/CMakeLists.txt) | Registers the 1D/2D/3D radiation pulse. |
| [`bin/radiation_tests/MatterCoupling/CMakeLists.txt`](../bin/radiation_tests/MatterCoupling/CMakeLists.txt) | Registers the combined gas/radiation relaxation fixture. |
| [`bin/radiation_tests/MatterCoupling/problem.cpp`](../bin/radiation_tests/MatterCoupling/problem.cpp) | Initializes smooth moving gas with radiation supplied by the selected source equilibrium. |
| [`bin/radiation_tests/RadiationPulse/problem.cpp`](../bin/radiation_tests/RadiationPulse/problem.cpp) | Initializes a localized radiation pulse. |
| [`bin/radiation_tests/Streaming/CMakeLists.txt`](../bin/radiation_tests/Streaming/CMakeLists.txt) | Registers the 1D/2D/3D radiation streaming problem. |
| [`bin/radiation_tests/Streaming/problem.cpp`](../bin/radiation_tests/Streaming/problem.cpp) | Initializes directed radiation transport and reference state. |
| [`bin/science/Collapse/CMakeLists.txt`](../bin/science/Collapse/CMakeLists.txt) | Registers 3D self-gravitating collapse. |
| [`bin/science/Collapse/problem.cpp`](../bin/science/Collapse/problem.cpp) | Initializes a collapsing gas sphere. |
| [`bin/science/Polytrope/CMakeLists.txt`](../bin/science/Polytrope/CMakeLists.txt) | Registers the 3D self-gravitating polytrope. |
| [`bin/science/Polytrope/problem.cpp`](../bin/science/Polytrope/problem.cpp) | Scales a Lane–Emden sphere into hydrostatic initial data. |
| [`bin/science/RotatingStar/CMakeLists.txt`](../bin/science/RotatingStar/CMakeLists.txt) | Registers the 3D rotating star. |
| [`bin/science/RotatingStar/problem.cpp`](../bin/science/RotatingStar/problem.cpp) | Scales and interpolates the original oblate SCF equilibrium into conserved fields. |
| [`bin/science/RotatingStar/equilibrium.inc`](../bin/science/RotatingStar/equilibrium.inc) | Embedded 100×100 positive-quadrant density and internal-energy table. |
| [`bin/science/RotatingStar/import_equilibrium.py`](../bin/science/RotatingStar/import_equilibrium.py) | Recreates and checks the embedded table from the pinned source dataset. |

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
| [`octotigerII/radiation/matterCoupling.hpp`](../octotigerII/radiation/matterCoupling.hpp) | Typed local and transport-forced matter exchange interfaces. |
| [`octotigerII/radiation/coupledPatch.hpp`](../octotigerII/radiation/coupledPatch.hpp) | Source-aware midpoint patch update used by shadows and diffusion regressions. |
| [`octotigerII/radiation/diffusionFlux.hpp`](../octotigerII/radiation/diffusionFlux.hpp) | Thick-cell face flux interpolation and moving material source equilibrium. |
| [`octotigerII/radiation/couplingDiagnostics.hpp`](../octotigerII/radiation/couplingDiagnostics.hpp) | Nonfatal optical-depth, trapping, and reduced-speed estimates. |
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
| [`octotigerII/subgrid/fluxPacket.hpp`](../octotigerII/subgrid/fluxPacket.hpp) | Time-tagged face-flux packet values and the sole FieldFluxPacket definition. |
| [`octotigerII/refinement/criteria.hpp`](../octotigerII/refinement/criteria.hpp) | AMR criteria inputs, thresholds, and callback interface. |
| [`octotigerII/amr/hierarchy.hpp`](../octotigerII/amr/hierarchy.hpp) | Covered coarse-state evolution, reconstruction, conservative transfer, and startup mesh API. |
| [`octotigerII/amr/regridSelection.hpp`](../octotigerII/amr/regridSelection.hpp) | RegridResult and selectMesh declarations for refinement and coarsening decisions. |
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
| [`src/main.cpp`](../src/main.cpp) | CLI/HPX startup, run invocation, progress table with total/leaf subgrid counts by level, and completion report. |
| [`src/config.cpp`](../src/config.cpp) | Option definitions, INI/command-line merge, defaults, and cross-option validation. |
| [`src/mesh.cpp`](../src/mesh.cpp) | Mesh indexing, block ancestry, cell geometry, and time-state methods. |
| [`src/physics/finiteVolume.cpp`](../src/physics/finiteVolume.cpp) | Returns the displayed name of the finite-volume scheme. |
| [`src/hydro/hydroSystem.cpp`](../src/hydro/hydroSystem.cpp) | Hydro state conversions, HLLC fluxes, characteristic speeds, and limiting. |
| [`src/radiation/radiationTransport.cpp`](../src/radiation/radiationTransport.cpp) | Radiation M1 fluxes, HLL transport, and realizability limiting. |
| [`src/radiation/matterCoupling.cpp`](../src/radiation/matterCoupling.cpp) | Conservative implicit gray source integration with safeguarded nonlinear solves. |
| [`src/runtime/radiation.cpp`](../src/runtime/radiation.cpp) | Distributed coupled midpoint/source stages and interval rollback. |
| [`src/composition/species.cpp`](../src/composition/species.cpp) | Periodic-table data, species parser, validation, and initial fractions. |
| [`src/problems/laneEmden.cpp`](../src/problems/laneEmden.cpp) | Numerical Lane–Emden integration and interpolation. |
| [`src/subgrid.cpp`](../src/subgrid.cpp) | Snapshot allocation and initialization support. |
| [`src/topology.cpp`](../src/topology.cpp) | Leaf layout, halo source plans, fine/coarse reflux, and gravity-work face maps. |
| [`src/refinement/criteria.cpp`](../src/refinement/criteria.cpp) | Density/mass/error refinement decisions. |
| [`src/amr/hierarchy.cpp`](../src/amr/hierarchy.cpp) | Evolved coarse shadows, reconstruction, conservative state transfer, and combined gas/gravity energy transfer. |
| [`src/amr/regridSelection.cpp`](../src/amr/regridSelection.cpp) | Refinement/coarsening decisions, geometric buffers, periodic neighbors, and 2:1 balance. |
| [`src/amr/initialization.cpp`](../src/amr/initialization.cpp) | Startup AMR refinement and lookahead loop, initial-condition resampling, and timestep/signal-speed estimates. |
| [`src/gravity/diagonal/fmm.cpp`](../src/gravity/diagonal/fmm.cpp) | Harmonic FMM coefficient algebra, translations, and interaction operators. |
| [`src/gravity/solver.cpp`](../src/gravity/solver.cpp) | Standalone serial uniform-grid FMM reference solve. |
| [`src/gravity/fieldSolver.cpp`](../src/gravity/fieldSolver.cpp) | Distributed uniform-grid FMM plus facade selecting uniform or adaptive solve. |
| [`src/gravity/adaptiveFieldSolver.cpp`](../src/gravity/adaptiveFieldSolver.cpp) | Distributed sparse adaptive FMM traversal and field publication. |
| [`src/gravity/gravityFields.cpp`](../src/gravity/gravityFields.cpp) | Translation unit for the gravity field type. |
| [`src/gravity/images.cpp`](../src/gravity/images.cpp) | Periodic/reflected gravity image offsets and pair lists. |
| [`src/gravity/ewald.cpp`](../src/gravity/ewald.cpp) | Periodic Newtonian Ewald sums and corrections. |
| [`src/storage.cpp`](../src/storage.cpp) | Storage component registration and typed HPX actions. |
| [`src/runtime.cpp`](../src/runtime.cpp) | Runtime construction/destruction, snapshot gathering, backend identity, and public bookkeeping queries. |
| [`src/runtime/internal.hpp`](../src/runtime/internal.hpp) | Shared private stage types, Runtime implementation state, LocalExecutor declarations, numerical/storage templates, and HPX action declarations. |
| [`src/runtime/localExecutor.cpp`](../src/runtime/localExecutor.cpp) | Locality worker queues, halo assembly, block transport and source operations, and HPX component/action registration. |
| [`src/runtime/stages.cpp`](../src/runtime/stages.cpp) | Stage dispatch, draining of task results and errors, completion checks, and scheduling statistics. |
| [`src/runtime/transport.cpp`](../src/runtime/transport.cpp) | Transport timestep selection, AMR subcycling, accepted-bank publication, and failed-step restoration. |
| [`src/runtime/gravity.cpp`](../src/runtime/gravity.cpp) | Gravity rungs, force stages, endpoint energy work, gravity-field publication, and coupled-step rollback. |
| [`src/runtime/regrid.cpp`](../src/runtime/regrid.cpp) | Mesh-selection loop, replacement of topology/field directories and executors, and combined-energy recovery after regridding. |
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
| [`tests/amrChecks.cpp`](../tests/amrChecks.cpp) | AMR hierarchy, conservative transfer, startup refinement, mesh-selection behavior, and regrid conservation checks. |
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
| [`tests/radiationCouplingChecks.cpp`](../tests/radiationCouplingChecks.cpp) | Independent thermal/momentum relaxation, source order, dual energy, and stiff pressure balance. |
| [`tests/radiationDiffusionChecks.cpp`](../tests/radiationDiffusionChecks.cpp) | Transparent/thick flux limits and evolved Fourier diffusion. |
| [`tests/radiationConservationChecks.cpp`](../tests/radiationConservationChecks.cpp) | Combined physical/weighted budgets, CSV ledgers, and diagnostic Silo roundtrip. |
| [`tests/radiationIntegrationChecks.cpp`](../tests/radiationIntegrationChecks.cpp) | Production coupling, temporal order, rotating/AMR/gravity conservation, and atomic failure checks. |
| [`tests/radiationDepthChecks.cpp`](../tests/radiationDepthChecks.cpp) | Coupled predictor stability and conservation on three refinement levels. |
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
| [`bin/gravity_tests/Gaussian/inputs`](../bin/gravity_tests/Gaussian/inputs) | Gravity Gaussian defaults. |
| [`bin/gravity_tests/Sphere/inputs`](../bin/gravity_tests/Sphere/inputs) | Uniform gravity sphere defaults. |
| [`bin/hydro_tests/KelvinHelmholtz/inputs`](../bin/hydro_tests/KelvinHelmholtz/inputs) | Kelvin–Helmholtz setup. |
| [`bin/hydro_tests/RayleighTaylor/inputs`](../bin/hydro_tests/RayleighTaylor/inputs) | Rayleigh–Taylor setup. |
| [`bin/hydro_tests/Sod/inputs`](../bin/hydro_tests/Sod/inputs) | Sod shock tube setup. |
| [`bin/radiation_tests/RadiationPulse/inputs`](../bin/radiation_tests/RadiationPulse/inputs) | Radiation pulse setup. |
| [`bin/radiation_tests/MatterCoupling/inputs`](../bin/radiation_tests/MatterCoupling/inputs) | Combined gas/radiation relaxation setup. |
| [`bin/radiation_tests/Streaming/inputs`](../bin/radiation_tests/Streaming/inputs) | Radiation streaming setup. |
| [`bin/science/Collapse/inputs`](../bin/science/Collapse/inputs) | Collapse setup. |
| [`bin/science/Polytrope/inputs`](../bin/science/Polytrope/inputs) | Polytrope setup. |
| [`bin/science/RotatingStar/inputs`](../bin/science/RotatingStar/inputs) | Oblate rotating-star setup. |
| [`tests/test_problem/inputs`](../tests/test_problem/inputs) | Small combined-physics application test setup. |
| [`examples/collapse.ini`](../examples/collapse.ini) | User-facing collapse run. |
| [`examples/gravity-gaussian.ini`](../examples/gravity-gaussian.ini) | User-facing Gaussian gravity run. |
| [`examples/gravity-sphere.ini`](../examples/gravity-sphere.ini) | User-facing uniform sphere gravity run. |
| [`examples/kelvin-helmholtz.ini`](../examples/kelvin-helmholtz.ini) | Kelvin–Helmholtz run. |
| [`examples/kelvin-helmholtz-time-refinement.ini`](../examples/kelvin-helmholtz-time-refinement.ini) | Kelvin–Helmholtz run with AMR time refinement. |
| [`examples/mass-fractions.ini`](../examples/mass-fractions.ini) | Material fraction transport example. |
| [`examples/polytrope.ini`](../examples/polytrope.ini) | Self-gravitating polytrope example. |
| [`examples/radiation-pulse.ini`](../examples/radiation-pulse.ini) | Radiation pulse example. |
| [`examples/radiation-matter.ini`](../examples/radiation-matter.ini) | Combined gray radiation/matter example. |
| [`examples/rayleigh-taylor-amr.ini`](../examples/rayleigh-taylor-amr.ini) | Refined Rayleigh–Taylor example. |
| [`examples/rotating-star.ini`](../examples/rotating-star.ini) | Corotating star with hierarchical gravity rungs. |
| [`examples/rotating-star-inertial.ini`](../examples/rotating-star-inertial.ini) | Same star on a stationary grid. |
| [`examples/sod.ini`](../examples/sod.ini) | Sod shock tube example. |
| [`examples/streaming.ini`](../examples/streaming.ini) | Streaming radiation example. |

## Architecture notes

- The public `Runtime` API remains in `octotigerII/runtime.hpp`. Its lifecycle
  and bookkeeping implementation is in `src/runtime.cpp`; private runtime
  responsibilities are grouped under `src/runtime/`. `internal.hpp` supplies
  shared types, state, templates, and action declarations. `localExecutor.cpp`
  owns the locality queue and the single set of HPX action registrations.
- `stages.cpp` dispatches work, drains every task before returning an error,
  verifies completion, and accumulates scheduling statistics. The transport,
  gravity, and regrid callers publish completed field banks or replacement
  mesh/field directories after their required work succeeds. Generic field
  storage remains in `octotigerII/storage/` and `src/storage.cpp`.
- AMR state evolution and conservative transfer belong to `amr/hierarchy.hpp`
  and `src/amr/hierarchy.cpp`. Mesh-selection decisions belong to
  `amr/regridSelection.hpp` and `src/amr/regridSelection.cpp`; startup refinement
  and runtime regridding both use that interface.
- The three gravity implementations have distinct roles: `solver.cpp` is the
  standalone serial reference, `fieldSolver.cpp` implements the distributed
  uniform-grid solver and selects the solver, and `adaptiveFieldSolver.cpp`
  implements the distributed adaptive solver. `FieldSolver` selects the
  adaptive solver whenever AMR is enabled.
- Each problem's one-line manifest registers its name and requirements.
  `cmake/Problems.cmake` discovers these manifests and generates the registry
  from `cmake/problems.cpp.in`; problem-specific physics stays in its module.
- Most generic transport behavior is implemented in
  `octotigerII/physics/finiteVolume.hpp`; its `.cpp` returns the scheme name.
  `src/gravity/gravityFields.cpp` likewise supplies a small translation unit
  for an inline field type. Header and implementation responsibilities are
  described individually above.
