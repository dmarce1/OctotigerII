# Opaque rotating-star startup validation — 2026-09-27

The governing construction and its approximations are described in
[the design](../../radiating-star-design.md). The earlier leaking equilibrium
work is preserved in [the research record](../../research/radiating-star-2026-09/README.md).

The implemented `radiatingStar` problem uses physical light speed, constant
gray absorption opacity 0.34 cm²/g, an isolated rotating gas/radiation SCF star,
zero prescribed photon heating, and an insulating radiation-energy box. Gas
retains outflow boundaries. No field is reset toward the initial profile.

## Actual startup evolution

Both meshes cover the same box and run to the same physical time, 6.49913 s,
using the runtime's full reported stable timestep (with a final shortened step).
This is only **2.70e-4 dynamical times**. Serial and single-locality HPX results
agree at the printed precision. The raw values are in
[star-checks.log](star-checks.log); HPX's four-test result is in
[hpx-star.log](hpx-star.log).

| Diagnostic | 8³ cells | 16³ cells |
| --- | ---: | ---: |
| Accepted steps | 3 | 6 |
| Mass-normalized density L1 drift | 1.43074e-3 | 1.32959e-3 |
| Meridional RMS Mach (rotation removed) | 1.96068e-3 | 1.43041e-3 |
| Relative radiation-energy change | +1.98793e-2 | -1.18049e-4 |
| Scaled gas+radiation+gravity energy-budget residual | 1.71377e-16 | 1.46426e-16 |
| Scaled mass-budget residual | 0 | -1.17206e-16 |
| Maximum `|F|/(cE)` | 0.837711 | 0.951343 |

Photon heating and both inward/outward radiation-energy boundary ledgers are
exactly zero. Gas and gravitational boundary-energy transport are included in
the total-energy budget. Gas temperature, pressure, and radiation realizability
remain admissible.

The coarsest radiation-energy change is appreciable despite the long continuum
thermal time. These very coarse Cartesian samples underresolve the initial
stratification; the continuum diffusion estimate does not bound numerical
redistribution in an unresolved model. The finer result improves, but two
startup points do not establish an asymptotic convergence order or long-term
stability. The density-drift reduction is modest and should not be oversold.

The conservative gas-energy variable becomes negative in a small amount of
dilute surface material after transport/gravitational energy accounting. The
existing dual-energy method still provides positive selected thermal energy.
The summed magnitude of negative gas energy, divided by the total energy norm,
is `5.36512e-11` and `2.54128e-12`; affected mass fractions are `1.71416e-9` and
`1.94141e-10`. These quantities are explicitly reported rather than hidden by
raising the atmosphere floor. Their further mesh dependence and any long-term
effect remain part of a resolved stability study.

The supplied example also completed two steps to 1 s with its maximum AMR level
temporarily restricted to 2; see [example-startup.log](example-startup.log).
That startup used the original source starting guess and happened not to
encounter the particular coarse-grid failure described below. The final source
kernel is verified by the full-CFL two-mesh tests above.

## Focused checks and one startup defect fixed

- All 16 source tests pass, including optically stiff moving balance through
  `c chi dt=1e12`, temporal-order references, and a captured bright-inflow
  atmosphere cell. See [source-checks.log](source-checks.log).
- All four stellar tests pass in serial and HPX. They verify full-c/no-source
  defaults, independent pressure-gradient and M1-force consistency, supported
  configuration, and the evolution above.
- All ten serial coupled-runtime integration tests pass after the source fix,
  including rotation, AMR/subcycling, gravity schedules, imposed acceleration
  work, and failure rollback. These and the source/star checks are recorded in
  [serial-coupling-star-integration.log](serial-coupling-star-integration.log).
- The insulating-boundary tests pass in 1D/3D serial and 3D HPX. They cover
  actual physical faces inside padded predictor patches, AP advection, limiter
  fallback, mixed-level evolution/regridding, and accepted energy/momentum
  ledgers. See [serial-boundary-diffusion.log](serial-boundary-diffusion.log)
  and [hpx-boundary.log](hpx-boundary.log).
- The six 3D and eight 1D diffusion checks pass, including the 1D coupled Fourier
  decay comparison. The final 1D source-kernel relink is recorded separately in
  [serial-diffusion-1d.log](serial-diffusion-1d.log).

The first coarse startup exposed a nonlinear initial-guess problem in a dilute
cell. Its accepted transport drive increased radiation energy from about 2.8
to 1.65e5 erg/cm³. Holding the radiation at its old value when starting an
implicit stage assigned that transported energy to trial gas heat, causing
enormous trial emission and exhausting cone-constrained retries. The captured
failure is preserved in [atmosphere-source-before.log](atmosphere-source-before.log).

The solver now considers the admissible transported state and selects it only
when its implicit residual is smaller than the existing guess. This preserves
already-balanced stiff states and resolves the bright-inflow failure. The
source equations, integration coefficients, conservation pairing, atmosphere,
timestep, and tolerances are unchanged. The exact failure is a permanent test.

A separate HPX test failure was only normalization: two opposing pressure
impulses near 4.33e17 differed by one ULP (64), while cell momentum stayed zero.
Normalizing by their tiny net difference gave a false order-one error. Using
the magnitudes of the accepted inward and outward impulses gives a residual
about 6.04e-17. The test tolerance was unchanged; three independent repeats and
the integrated HPX test pass with this physically meaningful normalization.

## Reproduction and remaining evidence

Configure the ordinary full-physics Release build with tests enabled; serial
uses `OCTOTIGERII_WITH_HPX=OFF`, HPX uses `ON`. Build these targets:

```sh
cmake --build BUILD --target radiatingStarChecks-3d radiationBoundaryChecks-3d \
  radiationCouplingChecks-3d radiationIntegrationChecks-3d -j2
ctest --test-dir BUILD --output-on-failure \
  -R '^3d\.(radiatingStarChecks|radiationBoundaryChecks|radiationCouplingChecks|radiationIntegrationChecks)\.'
```

The example is [examples/radiating-star.ini](../../../examples/radiating-star.ini).
For the recorded small AMR startup:

```sh
BUILD/octoII-3d --config examples/radiating-star.ini --amr.maxLevel=2 \
  --runtime.stopTime=1 --output.enabled=off --verification.analytic=off \
  --verification.directSamples=0
```

Long evolutions on meshes that resolve the core and surface, perturbation
response, and temporal/spatial refinement remain necessary before presenting
this as a stable star maintained for many dynamical times. These results
establish a runnable source-free construction, local balances, insulating
boundaries, and conservation in startup; they do not establish that later claim.
