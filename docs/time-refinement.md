# Level-wise time refinement

`timestep.refinement = on` is the default. Set `--timestep.refinement=off`
to retain globally synchronized timesteps and the existing gravity
kick/transport/solve/kick sequence.

This implementation applies to mixed-level AMR hydro and/or radiation transport,
including material partial densities, tracers, and self-gravitating gas. All
active cells at a given time level advance together. By default time levels
equal spatial levels. Uniform meshes have one
timestep and do not gain temporal subcycling. Runs with imposed uniform external
acceleration retain global stepping, both with and without self-gravity.

Self-gravitating AMR runs select their momentum schedule independently of the
energy and regrid controls:

```ini
[timestep]
refinement = on
coarseLevel = 0

[gravity]
timeIntegration = hierarchical
energyTreatment = mullen
conserveRegridEnergy = on
```

`gravity.timeIntegration=hierarchical` is the default. `conventional` evaluates
the full source force for each active level at that level's own cadence, using
time-interpolated inactive sources. The hierarchical schedule assigns slow-slow
and slow-fast interactions to the slower level's interval and recurses on
fast-fast interactions. The bins are shared by time level; this is not an
independent timestep choice for every cell. These schedules have one cadence
on a uniform mesh. Disabling AMR or timestep refinement selects the global
reference path regardless of `gravity.timeIntegration`.

## Timesteps

`timestep.coarseLevel = n` groups spatial levels `0..n` into a common time
level `n`. Finer spatial levels retain their own time levels. Equivalently,
`timeLevel = max(spatialLevel, n)`. The default `0` preserves the original
schedule. The option accepts `0..16` and is inactive when timestep refinement
is disabled. The selected level need not contain leaves; selecting at or above
the finest occupied spatial level gives one common timestep on the AMR mesh.

`Runtime::stableTimestep()` returns the minimum CFL/acceleration limit over
all blocks in the coarsest occupied time level
when time refinement applies. Signal-speed maxima still include every level for
the AMR travel/buffer checks. The caller may shorten the synchronization interval,
for example to reach the requested output or stopping time.

For example, with `coarseLevel=5`, spatial levels `0..5` advance together,
level 6 starts at half their step, and level 7 starts at a quarter of it.
The common step respects the tightest limit anywhere in levels `0..5`.
Grouping can reduce the number of HOLD shells, gravity solves, and whole-mesh
source-work passes, at the cost of more frequent coarse-cell transport.
Its benefit depends on the occupied levels and their CFL limits; it must be
measured for the workload.

For each parent interval, the next finer time level starts with half its
parent's step (or 1/2^k for a gap of k time levels). Before **every** substep,
the runtime reduces the CFL limit over all blocks on that level. It halves the
step repeatedly until the limit is met. Thus 1/2, 1/4, 1/8, and smaller ratios are
supported. A level may shrink its step during the parent interval; it does not
grow it again until the next parent interval. A binary fraction clock preserves
alignment and ends exactly at the parent's synchronization time.

The coarsest level is also checked before each accepted step. If a caller
supplies an interval longer than its current CFL limit, that interval is itself
subdivided dyadically. A step that cannot advance floating-point physical time
is rejected.

Self-gravitating timesteps include the acceleration constraint as well as the
hydrodynamic CFL limit. Changing a step ratio at an aligned interval boundary
does not change the level masks used by open gravity interaction shells.

## States and fluxes

Transport-only stepping retains a level's old and provisional transport
endpoints while its children advance. Fine halos interpolate those conserved
endpoints at the fine step's starting time, then apply spatial prolongation
and the MUSCL-Hancock predictor.

Coupled gravity instead supplies a separate coarse forecast from the old state
plus the coarse duration times the initial full numerical right-hand side.
That initial derivative uses canonical fine-face fluxes at coarse/fine
boundaries, including a probe-time reflux correction. This forecast is stored
in a fourth state bank (bank 3); the actual coarse transport endpoint remains
separate until its accepted fine fluxes have refluxed. Both fine hydro halos
and conventional gravity's inactive source densities interpolate the dedicated
forecast. Interpolating the unrefluxed transport endpoint would leave an O(H)
state error at a coarse/fine interface even with midpoint fluxes. Same-level
halo reads use the same time and bank on every block.

Each fine boundary-flux register accumulates `sum(dt_f * F_f) / dt_parent`.
Multiplication by `dt_parent` during reflux therefore compares the coarse flux
integral against the sum of fine flux integrals over exactly the same interval.
Area averaging is retained at coarse/fine faces. Hydro, radiation, and species
use the same interval accounting. Coarse corrections are applied only after
all finer levels have reached that coarse endpoint.
At a spatial coarse/fine interface whose two sides share one time level,
reflux uses the fine side's current substep flux. The accumulated interval
register is used only across different time levels. This distinction also
applies when the common coarse group subdivides an oversized caller interval.

A third state bank (bank 2) retains the published state during an entire
synchronization interval. If any phase fails, all tasks/transfers are drained and the checkpoint
is restored; public time, boundary totals, and accepted-step counts do not
advance. Numerical failures propagate to the caller after restoration. The
runtime does not automatically retry an inadmissible coarse update. A newly
tighter fine CFL limit is handled by further subdivision before its next step.

Regridding, snapshots, and output remain synchronized operations. Shadow states
used for AMR error estimation are separately CFL-subcycled over the interval.
`Runtime::statistics().levelSteps[level]` counts accepted transport steps for
each spatial level on the refinement path. Occupied spatial levels in one
time group receive equal counts. Snapshot step counters describe global
synchronization intervals.

## Progress inside a synchronization interval

Set `runtime.progressLevel = n` to print progress for time levels `n` and
coarser. The default `-1` disables these lines. For the QueenBee4 level-7
case, `--runtime.progressLevel=7` exposes every substep; `5` hides finer
substeps. Choose a value at least as large as `timestep.coarseLevel` to see
the grouped coarse steps.

Progress is flushed to the coordinator's standard error stream, so it appears
in task 0's Slurm `.err` log with the supplied batch script. Lines identify the
coarse step, time level, active stage, substep start time and duration in seconds,
and completed fraction of the containing parent interval. `wall_s` is elapsed
wall time since runtime construction; `substep_wall_s` includes waits and finer
descendants. `stage=done` confirms that the substep's reflux and gravity closure
have finished. Earlier stage lines identify work about to begin, rather than
claiming its completion. The `cfl` line carries a trial duration; `begin` carries
the accepted duration after applying the CFL limit.

The gravity path also reports preparation before the first substep, including
scratch allocation, the initial flux probe, and gravity source-rate assembly.
During gravity source-rate assembly, it reports `blocks_done` and
`blocks_total` after a completed block at least every 30 seconds, and once
at the end. If a block waits indefinitely on a future, the completed count
stops at the preceding block; the line does not identify the specific future.
These are progress messages, not additional Silo frames or global conservation
records: intermediate levels can be at different physical times. A failure
can roll back substeps already reported as done within the current coarse step.

## Finite-volume gravity coupling

`Runtime::advanceGravity(dt)` owns a complete self-gravitating gas interval,
including transport, force solves, source updates, reflux, and energy closure.
The simulation driver selects it automatically. It requires synchronized initial
gravity and returns with the gas and gravity synchronized again; callers must
not add a second set of endpoint kicks around it. The manual global-step APIs
remain available for reference calculations.

The gas states used by reconstruction and halo interpolation represent physical
momentum and energy. The coupled path uses a finite-volume explicit midpoint
stage: a zero-duration Riemann-flux probe supplies the full discrete hydro
right-hand side, to which gravity adds momentum forcing and discrete mass-flux
work. A half-step of that combined derivative predicts cell-centered midpoint
states. Spatial reconstruction and the shared Riemann solver then evaluate the
midpoint flux, with the Hancock time predictor disabled to avoid advancing the
stage twice. The conservative update retains the original state as its base.

Midpoint and coarse forecasts also receive a predictor-only convex kinetic
energy remainder of O(dt²), helping keep cold accelerating states physical.
This term is not added to the accepted conservative gas-energy update.

For halos, form midpoint donor states before spatial prolongation, so coarse
donor slopes include the same evolution as interior states. Initial derivatives
are populated on all levels; active derivatives are refreshed before their next
prediction. Inactive derivatives may be old by a coarse interval, which is
consistent with the conditional second-order argument in the
[derivation](gravity-time-coupling-derivation.md). Coarse forecasts use the
canonical probe-time derivative described above. The reconstruction, Riemann
solver, and flux limiters remain shared with the global transport path.

The communication optimization caches the initial two-layer halo from the
probe for reuse by the corrector. Midpoint stage data still cross block
boundaries, also using two ghost layers. Reuse requires the same physical
snapshot, donor time interpolation, and field version; a kick or other change
to the conservative update base requires fresh data. Raw AMR donor values
are cached before prolongation and physical-boundary transforms, so the
midpoint donors can be formed and limited at their own time. This changes
data reuse without replacing the numerical midpoint derivative.

Each hierarchical shell keeps its opening field while descendants advance,
recording provisional impulses. After reflux, its endpoint solve replaces those
impulses by the paired opening/closing quadrature. This accepts the same HOLD
interaction impulse quadrature, but adapts the intermediate states to Eulerian
mass transport; it is not the particle Hamiltonian map and has no symplectic
claim. Conventional mode instead closes the full-force impulse on the active
level alone, so opposite interaction directions can use different time
quadratures.

The FMM preserves matching source and target force degrees with a separate
auxiliary local expansion of degree p+1. This is distinct from potential
reciprocity; see the [mutual-force derivation](parallel-fmm.md#scalar-reciprocity-and-mutual-force).

With `gravity.energyTreatment=mullen`, every open shell accumulates mass fluxes
on all faces throughout its interval. Closure uses the canonical fine-subface
flux on both sides of coarse/fine interfaces. A level's initial energy forecast
includes its nested fast set; descendants' work on inactive coarse neighbors is
deferred until that coarse level closes. Shell closure replaces recorded
provisional work, including these transfers, with endpoint work. All gravity
ledgers close before regridding or publication. See
[gravity energy](gravity-energy.md) for the accounting and independent options.

This implementation establishes the coupling structure before optimizing its
execution. Nested endpoint fields are reused at the same physical time to avoid
duplicate partial solves. Partial source solves still perform an unpruned
upward tree pass, and a central coordinator manages the source and work ledgers. Added field
storage, remote reads, and coordination can outweigh saved force work. Neither
mode currently carries a speedup or scaling claim; performance must be
evaluated separately from conservation and temporal accuracy.

## Gravity coupling validation

### Coarse time grouping, 2026-09-30

The grouped schedule passed the following local checks:

| Build and scope | Passed |
|---|---:|
| Serial, 1D option parsing and validation | 44 |
| Serial, 1D time refinement | 11 |
| Serial, 1D grouped and existing mixed-level radiation integration | 2 |
| Serial, full 3D gravity time integration suite | 15 |
| HPX, two localities, 1D time refinement | 11 |
| HPX, two localities, 1D grouped and existing mixed-level radiation integration | 2 |
| HPX, two localities, 3D grouped and nested gravity conservation | 3 |

Grouping tests cover common CFL selection, partial grouping with finer
subcycling, a selected coarse time level above the finest spatial level,
subdivision of an oversized common interval, spatial-level step counts,
conservation, and synchronized regridding. Both hierarchical and conventional
gravity schedules are exercised. In the smooth fixed-mesh grouped HOLD test,
spatial levels 1 and 2 share time level 2 while spatial level 3 subcycles.
The observed orders for two successive timestep halvings were:

| Field | First halving | Second halving |
|---|---:|---:|
| Density | 2.04328 | 2.08299 |
| Momentum | 2.04390 | 2.08380 |
| Gas energy | 2.04349 | 2.08304 |

These support second-order time accuracy for this smooth fixed-mesh test.
They do not establish a speedup or cluster scaling. A serial Polytrope
application smoke run also completed a coarse step and emitted flushed stage
lines with `runtime.progressLevel=2` and `timestep.coarseLevel=2`.

The distributed checks can be reproduced after building the named test targets:

```bash
python3 tests/distributed.py release/1d/tests/timeRefinementChecks-1d
python3 tests/distributed.py release/1d/tests/radiationIntegrationChecks-1d \
  --gtest_filter='RadiationIntegration.GroupedCoarseLevelsSubcycleAndConserve:RadiationIntegration.MixedLevelSubcyclingAndRegriddingConserve'
python3 tests/distributed.py release/3d/tests/gravityTimeIntegrationChecks-3d \
  --gtest_filter='GravityTimeIntegration.GroupedCoarseLevelsSubcycleAndConserve:GravityTimeIntegration.GroupingAboveFinestLevelUsesOneStableConservativeStep:GravityTimeIntegration.NestedThreeLevelShellAndFluxRegistersConserve'
```

The serial temporal gate is
`GravityTimeIntegration.GroupedHoldSmoothTemporalConvergence` in
`gravityTimeIntegrationChecks-3d`.

### Previous coupling validation

Smooth fixed-mesh temporal regressions give observed orders of approximately
2.04–2.09 for density, momentum, and gas energy in both integration modes.
These support second-order accuracy for the tested smooth problem; they do not
extend the claim to shocks, changing masks, or unbounded step ratios. The
coupled checks also exercise conservation, nested level registers, CFL ratios,
global-step compatibility, and synchronized regridding. Results and execution
scope are recorded in the [gravity integration validation report](validation/gravity-time-integration.txt).

The local physical-flux predictor experiment did not meet these fixed-mesh
second-order requirements. Its approximately first-order temporal results are
preserved as [historical evidence](validation/local-gravity-predictor.txt),
with the [original failed assertions](validation/local-gravity-predictor-migration.txt).
The original second-order gates have been restored. The
[cached-midpoint 3D serial validation](validation/cached-midpoint/gravity-serial-3d.md)
passes the nonrotating and rotating temporal gates with component orders
approximately 2.04–2.09. The report distinguishes these new runs from historical
midpoint results and scopes the separate halo request measurements.

## Initial hydro/radiation validation

`timeRefinementChecks` covers the default/off option, a hot fine region needing
at least four substeps, hydro/radiation conservation, temporal halo interpolation,
regridding, and nested registers across three occupied levels. A 1D streaming
Gaussian crossing the coarse/fine interface checks second-order convergence.
Existing AMR, species, storage, option, and gravity-energy suites are also run.

The original validation used the HPX-free serial build in 1D, 2D, and 3D.
An HPX installation was unavailable for a distributed execution test at that
stage. These results predate the gravity coupling described above and do not
establish its convergence, stability, or distributed performance.

Validation result: **109 tests passed across 10 test executables**. The streaming
L1 errors for 8, 16, and 32 cells per block were 1.541870e+09, 4.186460e+08,
and 1.023420e+08 (energy density integrated over x in cgs); the final observed
order was 2.032330. Full output is in
[the validation log](validation/time-refinement.txt).
