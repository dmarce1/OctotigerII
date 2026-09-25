# Level-wise time refinement

`timestep.refinement = on` is the default. Set `--timestep.refinement=off`
to retain globally synchronized timesteps and the existing gravity
kick/transport/solve/kick sequence.

This implementation applies to mixed-level AMR hydro and/or radiation transport,
including material partial densities, tracers, and self-gravitating gas. All
active cells at a given spatial level advance together. Uniform meshes have one
timestep and do not gain temporal subcycling. Runs with imposed uniform external
acceleration retain global stepping, both with and without self-gravity.

Self-gravitating AMR runs select their momentum schedule independently of the
energy and regrid controls:

```ini
[timestep]
refinement = on

[gravity]
timeIntegration = hierarchical
energyTreatment = mullen
conserveRegridEnergy = on
```

`gravity.timeIntegration=hierarchical` is the default. `conventional` evaluates
the full source force for each active level at that level's own cadence, using
time-interpolated inactive sources. The hierarchical schedule assigns slow-slow
and slow-fast interactions to the slower level's interval and recurses on
fast-fast interactions. The bins are shared by spatial level; this is not an
independent timestep choice for every cell. These schedules have one cadence
on a uniform mesh. Disabling AMR or timestep refinement selects the global
reference path regardless of `gravity.timeIntegration`.

## Timesteps

`Runtime::stableTimestep()` returns the CFL limit of the coarsest occupied level
when time refinement applies. Signal-speed maxima still include every level for
the AMR travel/buffer checks. The caller may shorten the synchronization interval,
for example to reach the requested output or stopping time.

For each parent interval, the next finer spatial level starts with half its
parent's step (or 1/2^k for a gap of k spatial levels). Before **every** substep,
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

A third state bank (bank 2) retains the published state during an entire
synchronization interval. If any phase fails, all tasks/transfers are drained and the checkpoint
is restored; public time, boundary totals, and accepted-step counts do not
advance. Numerical failures propagate to the caller after restoration. The
runtime does not automatically retry an inadmissible coarse update. A newly
tighter fine CFL limit is handled by further subdivision before its next step.

Regridding, snapshots, and output remain synchronized operations. Shadow states
used for AMR error estimation are separately CFL-subcycled over the interval.
`Runtime::statistics().levelSteps[level]` counts accepted transport steps for
each level on the refinement path. Snapshot step counters describe global
synchronization intervals.

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

Smooth fixed-mesh temporal regressions give observed orders of approximately
2.04–2.09 for density, momentum, and gas energy in both integration modes.
These support second-order accuracy for the tested smooth problem; they do not
extend the claim to shocks, changing masks, or unbounded step ratios. The
coupled checks also exercise conservation, nested level registers, CFL ratios,
global-step compatibility, and synchronized regridding. Results and execution
scope are recorded in the [gravity integration validation report](validation/gravity-time-integration.txt).

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
