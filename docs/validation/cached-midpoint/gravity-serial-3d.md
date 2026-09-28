# Cached numerical midpoint: 3D serial validation

Date: 2026-09-28. Working-tree build based on commit
`c45735658a684ab5c3dc4c7fe716e5d8ede46d22`.

Build directory: `/tmp/octoii-gravity-local-serial`. Configuration: Release,
`/usr/problem/c++`, HPX disabled, profiling enabled. The existing second-order
gravity and rotating-gravity assertions are unchanged.

## Reproduction

Run from the repository root:

```sh
cmake --build /tmp/octoii-gravity-local-serial --target gravityTimeIntegrationChecks-3d rotatingGravityChecks-3d haloTrafficChecks-3d radiationIntegrationChecks-3d -j2
/tmp/octoii-gravity-local-serial/3d/tests/gravityTimeIntegrationChecks-3d --gtest_filter='GravityTimeIntegration.*SmoothTemporalConvergence'
/tmp/octoii-gravity-local-serial/3d/tests/rotatingGravityChecks-3d --gtest_filter='RotatingGravityIntegration.TemporalConvergence'
/tmp/octoii-gravity-local-serial/3d/tests/haloTrafficChecks-3d
/tmp/octoii-gravity-local-serial/3d/tests/radiationIntegrationChecks-3d --gtest_filter='RadiationIntegration.PhysicalAndReducedSpeedConserveOnPeriodicGrid:RadiationIntegration.CoupledTransportConvergesAtSecondOrderInTime:RadiationIntegration.ExternalAccelerationIncludesRadiationMomentumInItsWork:RadiationIntegration.MixedLevelSubcyclingAndRegriddingConserve:RadiationIntegration.RotatingGridConservesInertialCombinedBudgets:RadiationIntegration.NumericalFailureRollsBackAndNextStepMatchesCleanRun:RadiationIntegration.AddedEquilibriumRadiationHasNoLocalExchange:RadiationIntegration.GasRadiationGravityConserveAcrossAllSchedules'
```

## Results

The four requested targets built successfully. Every test invocation below
completed with exit status zero; no acceptance threshold was changed.

| Invocation | Result | Scope |
| --- | --- | --- |
| Gravity temporal filter | 2/2 passed | Original hierarchical and conventional fixed-mesh temporal gates |
| Rotating temporal filter | 1/1 passed | Original rotating hierarchical and conventional temporal gate |
| Halo traffic suite | 6/6 passed | Request accounting, midpoint reuse, opening-kick refresh, compact/full-halo equivalence |
| Focused radiation filter | 8/8 passed | Temporal order, source work, budgets, AMR/regrid, rotation, rollback, equilibrium, self-gravity |

The remaining non-temporal tests in the gravity and rotating executables were
not rerun in this validation. The radiation filter excludes the zero-opacity
equivalence and invalid-argument tests; the numerical-failure rollback test
is included.

| Gravity schedule | Density orders | Momentum orders | Gas-energy orders |
| --- | --- | --- | --- |
| Hierarchical | 2.05254, 2.08750 | 2.04190, 2.07969 | 2.05264, 2.08761 |
| Conventional | 2.05254, 2.08749 | 2.04181, 2.07968 | 2.05264, 2.08760 |
| Rotating hierarchical | 2.04609, 2.08513 | 2.04974, 2.08132 | 2.04627, 2.08524 |
| Rotating conventional | 2.04608, 2.08513 | 2.04973, 2.08132 | 2.04627, 2.08524 |

These are the original smooth mixed-mesh experiment: 2, 4, and 8 steps
compared with a 32-step reference at fixed spatial resolution. The original
error-ratio threshold below 0.35 is unchanged. This verifies temporal order for
the tested smooth case, not shock accuracy or arbitrary refinement ratios.

### Coupled radiation checks

The original radiation temporal gate uses 1, 2, and 4 steps against a 64-step
reference. Its required successive error reduction greater than 3.5 is
unchanged.

| Mesh | Radiation-energy relative errors |
| --- | --- |
| Uniform | 0.00788585, 0.00136863, 0.000312280 |
| Adaptive | 0.000786152, 0.000182610, 0.0000441965 |

The independent accelerated radiation/matter work comparison also passed its
original temporal-order and absolute-error gates at physical and reduced
light speeds. Combined budgets passed after mixed-level subcycling,
coarsening/refinement, rotation, and numerical failure followed by retry.
The self-gravity case covers global, hierarchical, conventional, and rotating
hierarchical schedules, with the physical/reduced light-speed cases selected
by the existing test. These are the existing benchmark scopes, not a proof of
uniform accuracy in every stiff limit.

### Halo request accounting

Each row is one complete coupled step on the test fixture. Both consecutive
steps reproduced these counts for configured worker counts one and two, with
work stealing disabled/enabled where selected by the test. This serial build
does not execute the HPX worker or cross-locality stealing paths.

| Case | Halo calls | Scalar range requests | Requested values | Calculated legacy values |
| --- | ---: | ---: | ---: | ---: |
| No opening kick | 32 | 336 | 86,016 | 186,368 |
| Opening acceleration kick | 40 | 448 | 114,688 | 186,368 |

The measured counts are 53.85% and 38.46% below the calculated legacy
five-gas/three-radiation-halo baseline, respectively. The opening kick must
refresh the changed gas update base. The legacy baseline was calculated from
the same donor plans, not measured in a pre-change 3D run.

Logs in this directory:

- `build-gravity-serial-3d.log`
- `gravity-temporal-serial-3d.log`
- `rotating-gravity-temporal-serial-3d.log`
- `halo-traffic-serial-3d.log`
- `radiation-focused-serial-3d.log`

Serial halo counters measure requested values and contain no remote payload.
They do not establish distributed bandwidth or application speedup.

The halo suite was rebuilt and rerun after adding the direct compact/full-halo
comparison and correcting execution-locality statistics. That comparison covers
nonuniform AMR donors, reflecting and time-dependent analytic boundaries, and
radiation flux limiting; all face fluxes and updated cells agree exactly. The
gravity/radiation integration logs above precede that diagnostics-only change;
their numerical code is unchanged.
