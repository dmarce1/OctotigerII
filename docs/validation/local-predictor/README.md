# Local source predictor validation

Results collected on 2026-09-28 from serial Release and HPX builds. The method
uses two initial ghost layers and restores physical exterior traces at half
time. These records establish numerical behavior and requested halo data;
they do not establish whole-application speedup or measured network bandwidth.

## Completed numerical checks

| Suite | Dimension | Result | Evidence |
|---|---:|---:|---|
| Local radiation predictor | 1 | 8/8 | [Log](local-predictor-serial-1d.log) |
| Local radiation predictor | 3 | 6/6 | [Log](local-predictor-serial-3d.log) |
| Radiation diffusion | 1 | 8/8 | [Log](diffusion-serial-1d.log) |
| Coupled radiation integration | 3 | 10/10 | [Log](integration-serial-3d.log) |
| Three-level coupled depth | 3 | 2/2 | [Log](depth-serial-3d.log) |
| Radiation boundaries | 3 | 4/4 | [Log](boundary-serial-3d.log) |
| Finite volume | 3 | 13/13 | [Log](finite-volume-serial-3d.log) |
| Time refinement | 3 | 8/8 | [Log](time-refinement-serial-3d.log) |
| Photon heating | 3 | 8/8 | [Log](photon-heating-serial-3d.log) |
| HPX radiation integration, one locality | 1 | 8/8 | [Log](radiation-integration-hpx-single-1d.log) |

The [gravity report](../local-gravity-predictor.txt) records 12/12 gravity
integration checks, 5/5 rotating-gravity checks, analytic local-predictor tests,
a coarse-forecast positivity regression, and a final runtime smoke check.
The [historical failures](../local-gravity-predictor-migration.txt) preserve the
old fixed-mesh second-order gates rather than hiding the accuracy change.

The long production runs preceded removal of unused gas-midpoint/gas-rate
scratch and the added full-forecast positivity guard. Subsequent runtime
rebuilds, HPX 1D integration/accounting, gravity helper checks, and a hierarchical
subcycling smoke check cover those final changes. The guard reduces slopes and
their physical driving together when the full coarse forecast is inadmissible;
it does not clip the forecast independently.

## What the checks establish

Two-layer patches give exactly the same updates and accepted fluxes as patches
with four layers whose outer two contain NaNs. A periodic domain split preserves
shared fluxes and updates exactly. Homogeneous exchange matches one conservative
source solve through optical intervals of 1e6. Sharp opacity transitions and
local positivity fallback remain admissible and conservative.

Physical boundary checks sample analytic center and face values at half time,
verify reflecting and insulating normal-face swaps, and retain tangential
variation at free boundaries. Production tests cover rotating boundaries,
gas/radiation energy and momentum, all three gravity schedules, subcycling,
regridding, photon heating, and failed-step rollback. The AMR shadow hierarchy
is exercised even when the shadow refinement criterion is disabled.

## Independent smooth-wave solution

For a small perturbation around isotropic radiation in massive stationary
pure-scattering gas, let q=F/c. The continuum equations are

```text
E_t = -chat q_x
q_t = -(chat/3) E_x - chat chi q.
```

On the unit periodic domain, with k=2*pi and s=chat*t, initial
E=1+epsilon*cos(k*x), q=0 has solution

```text
alpha = chi/2, omega = sqrt(k*k/3 - alpha*alpha)
A(s) = epsilon exp(-alpha*s) [cos(omega*s) + alpha/omega sin(omega*s)]
B(s) = epsilon k/(3*omega) exp(-alpha*s) sin(omega*s)
E = 1 + A(s) cos(k*x), q = B(s) sin(k*x).
```

The test uses cell averages, epsilon=1e-7, chi=4, chat/c=0.2, final s=0.4,
and 2*N timesteps. Nonlinear radiation corrections relative to the perturbation
are O(epsilon); gas feedback is suppressed by E/(rho*c*c). Both are much smaller
than the measured errors. Space and time refine together.

| Cells | Normalized L1 error in E and F/c | Error ratio |
|---:|---:|---:|
| 16 | 0.0280237 | — |
| 32 | 0.00673455 | 4.16 |
| 64 | 0.00154647 | 4.35 |

The independent accelerated Euler wave gives density orders 1.86 and 2.01 with
matching momentum and energy behavior. This establishes the tested joint
space–time order, not joint order for every adaptive self-gravitating flow.

The unchanged coupled Fourier test measures decay divided by the analytic
discrete equilibrium-diffusion rate of 1.00607, 1.00067, and 1.00009 at cell
optical depths 10, 100, and 1000. Full and reduced light speed agree. These
results apply to a prepared smooth mode and the tested limiter regime. They
do not establish uniform accuracy through an unresolved initial relaxation
layer or wherever positivity limiting changes the face flux.

## Temporal accuracy and timestep changes

The existing coupled radiation fixed-mesh temporal test passes unchanged in
3D. Its radiation-energy errors are

```text
uniform:  0.00788578, 0.00136858, 0.000312252
adaptive: 0.00079268, 0.000184844, 0.0000451919.
```

Homogeneous external-acceleration work also retains ratios near four. These
particular cases do not establish second-order time integration of an arbitrary
fixed spatial operator. The local physical derivative differs from the
numerical-flux derivative. Gravity's fixed-mesh tests now show approximately
first-order temporal convergence; their renamed diagnostic gates explicitly
check that contract, alongside the independent joint-refinement gate.

The former three-level test required the root interval to be no more than twice
the finest CFL interval because every leaf supplied a coarse half-step
prediction. That operation is gone. The [old-cap failure log](depth-obsolete-midpoint-cap-serial-3d.log)
records only obsolete timestep assertions; physical evolutions and conservation
gates completed successfully. The revised test requires the coarse CFL, checks
that the initial coarse interval is fully recovered, and retains level-step
counts, conservation, and admissibility. Each level's actual substeps still
obey its current local CFL.

## Communication

The [accounting report](halo-README.md) preserves the full fixture and raw logs.
The measured 1D requested halo values decrease from 288 to 64 per step, a
77.78% reduction. The 3D after count is 43,008; its 76.923% reduction is relative
to a calculated old baseline, not a measured old 3D run. Single-locality HPX
accounting passed with profiling disabled in 1D and enabled in 3D.

The production multi-process test remains unverified: the sandbox rejected
local TCP sockets and permission escalation did not complete. Requested
numerical halo data excludes message overhead, interior transfers, task
migration, gravity solves, and reflux. No application speedup is claimed.

## Reproduction

Configure following [BUILDING.md](../../../BUILDING.md), then build the named
targets, for example:

```sh
cmake --build BUILD --target localPredictorChecks-1d localPredictorChecks-3d \
  localGravityPredictorChecks-1d radiationDiffusionChecks-1d \
  radiationIntegrationChecks-3d radiationDepthChecks-3d radiationBoundaryChecks-3d \
  finiteVolumeChecks-3d timeRefinementChecks-3d photonHeatingChecks-3d \
  gravityTimeIntegrationChecks-3d rotatingGravityChecks-3d haloTrafficChecks-1d
```

Serial executables are in BUILD/1d/tests or BUILD/3d/tests and run without
arguments. HPX tests use `--hpx:threads=2 --hpx:bind=none`. The serial coupled
integration suite took approximately 271 seconds and the rotating-gravity
suite approximately 451 seconds; these are validation timings, not application
performance benchmarks.
