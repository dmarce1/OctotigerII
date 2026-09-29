# Binary SCF validation, 2026-09-29

These measurements apply to the binary SCF implementation added to the working
tree based on `fc683f9a562d772f2a5fa971ab65514d5828b0d3`. They establish the
initializer's algebra, transfer, and first-step integration, not a validated
mass-transfer experiment or a performance/scaling result.

## Independent checks

* The isolated convolution is compared with direct all-pairs summation at every
  target of an 8^3 asymmetric density distribution, including edge targets.
  The bound is 3e-14 relative. A single corner mass checks self-term omission
  and absence of wrapped periodic forces.
* Single-polytrope pressure and enthalpy are compared with their analytic power
  laws for n=.5, 1, 1.5, 3, and 5. Finite differences independently check
  `dH=dP/rho`. Unequal core/envelope indices and a density discontinuity test
  interface matching and inversion on both sides.
* Detached and semi-detached binary tests check the actual convergence
  residual. Roche filling is the equality of the donor's Bernoulli constant
  with the directly evaluated L1 effective potential.
* Integrated material fields recover the specified component masses after
  conservative transfer, within 1e-11 relative. Material sums equal hydro
  density, including the atmosphere exactly once. Total and auxiliary energy
  recover pressure with evolution gamma=5/3.
* Physical rescaling checks density, pressure, velocity, and angular velocity
  independently of structural indices and grid rotation.
* Failures are tested for invalid input, exhausted iterations, a converged but
  unresolved reference, a cropped hydro domain, and an overly coarse hydro mesh.

The integration test advances through rotating-frame gravity and hydro for
`0.001/omega` seconds. Its mass budget closes within 2e-12 relative and its total
gas-plus-gravity energy budget within 2e-9 of the initial energy norm, including
boundary transport. This is a deliberately coarse **wiring/conservation test**
with a relaxed virial gate. It does not measure physical long-term equilibrium.

The serial test build disables HPX, radiation, mass fractions, and profiling.
The HPX build includes the normal modules but selects gravity/hydro only for this
problem. The same initializer and integration tests are exercised in both.
Configuration serialization is checked in the HPX build. The integration suite
also passes with two and three HPX localities on one host. This validates process
handoff; it is not a distributed performance measurement.

## Reference resolution study

Equal component masses of 1.98847e33 g, separation 1e11 cm, n=1.5 throughout,
primary filling factor .8 and donor filling factor 1. Density tolerance 1e-5;
the virial gate is explicitly relaxed to 1 for the convergence experiment.

| Cells per axis | Iterations | Density L1 residual | Bernoulli residual | Virial residual | Omega (rad/s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 32 | 34 | 8.6303e-6 | 1.6632e-6 | 0.138275 | 5.1917197e-4 |
| 64 | 78 | 7.2912e-6 | 1.5830e-6 | 0.043014 | 5.1858249e-4 |
| 128 | 27 | 7.9293e-6 | 1.8923e-6 | 0.011894 | 5.1841007e-4 |

Full measured values are in [resolution.json](binary-scf/resolution.json).
The virial error decreases by factors of approximately 3.21 and 3.62 under
doubling, approaching second-order behavior in this case. Omega changes by
about 0.033% between the two finest resolutions. The 128^3 virial error is still
**1.19%**. An iteration residual below 1e-5 must not be reported as 1e-5 physical
equilibrium accuracy.

This is consistent with finite-resolution errors in the cell-mass gravity and
stellar surface/core sampling. It does not isolate every error source or prove
the same convergence rate for density jumps and more condensed cores. Core
resolution, hydro discretization, and subsequent evolution require their own
checks.

## Complete hydro handoff

The supplied `problem/science/BinaryScf/inputs` example uses mass ratio .7,
primary filling factor .9, donor filling factor 1, and 64^3 reference and
initial hydro meshes. It converged in 139 iterations:

| Quantity | Measured value |
| --- | ---: |
| Undamped density residual | 9.5757e-6 |
| Bernoulli residual | 2.2606e-6 |
| Reference virial residual | 0.0316679 |
| Hydro-mesh virial residual after the ordinary FMM solve | 0.0314803 |
| Angular velocity | 4.7867565e-4 rad/s |
| Total hydro mass, including atmosphere | 3.3803990091e33 g |

The requested stellar mass is 3.380399e33 g. The difference is the explicitly
added atmosphere. See [hydro-handoff-64.json](binary-scf/hydro-handoff-64.json).
This run passed both default 5% virial gates and the mass-retention check.

A stronger bipolytropic test with n_c=3, n_e=1.5 and density jumps 2 in both
stars, interface density fractions .5, mass ratio .7, and fills .95/1 also
converged on a 64^3 reference. Its 5.52% virial residual exceeds the default
gate. This is a concrete example of an algebraically converged model that
should be refined before use.

The supplied 128^3 bipolytrope example subsequently completed the full ordinary
FMM/hydro handoff on 4,096 leaf blocks (2,097,152 cells):

| Quantity | Measured value |
| --- | ---: |
| Iterations | 95 |
| Density residual | 9.4869e-6 |
| Maximum per-star Bernoulli residual | 1.7583e-6 |
| Reference virial residual | 0.0133714 |
| Hydro-mesh virial residual | 0.0129424 |
| Primary/donor core mass fractions | 0.120353 / 0.124483 |

See [bipolytrope-handoff-128.json](binary-scf/bipolytrope-handoff-128.json).
At this resolution, direct inversion across the density jump initially stalled
with alternating interface cells. The accepted implementation uses a proximal
inversion of the continuous enthalpy graph for interface-crossing updates,
including the common-pressure mixed state. Both the density update and each
star's actual Bernoulli balance must converge; the latter prevents a small
preconditioned update from concealing an enthalpy error. No gradual alteration
of the requested structural parameters was needed for the recorded run.

## Shared-K q=0.7 DWD model

The current DWD inputs use `scf.commonPolytropicK=on`: both stars have
`P=K*rho^(5/3)`, the donor surface passes through L1, and the primary
Bernoulli constant is solved from its mass at the common K. This replaces
the independent pressure normalizations used in the historical model below.
The physical scaling is 0.6 + 0.42 solar masses at separation 3e9 cm. The
common K is an output of that scaling, not a prescribed degeneracy constant.

The [reference resolution record](binary-scf/dwd-common-k-resolution.json)
was produced by:

```sh
./release/3d/tests/binaryScfChecks-3d \
  --gtest_also_run_disabled_tests --gtest_filter=BinaryScf.DISABLED_SharedKResolutionStudy \
  --hpx:threads=2 --hpx:bind=none
```

| Reference cells per axis | Virial residual | Accretor diameter (cm) | Donor diameter (cm) | Common K (CGS) |
| ---: | ---: | ---: | ---: | ---: |
| 32 | 0.0670069 | 1.81526e9 | 2.07389e9 | 2.40277e12 |
| 64 | 0.0174796 | 1.79367e9 | 2.06633e9 | 2.58170e12 |
| 128 | 0.00442889 | 1.78967e9 | 2.06910e9 | 2.62709e12 |

The diameters are volume-equivalent. On the 128^3 mesh, the primary/donor
x extents are 1.89804e9 / 2.51490e9 cm, y extents are 1.80313e9 / 1.99294e9
cm, and z extents are 1.70823e9 / 1.89804e9 cm. These cell-based surface
measurements are quantized at the reference spacing. Geometry diagnostics
exclude residual density tails outside the converged Bernoulli surface;
the historical diameters below included those tails.

At 128^3 the donor is 15.6% larger in volume-equivalent diameter, the
primary enthalpy fill is 0.84065, and the donor fill is 1. The maximum
component density residual is 6.895e-6 and the mass-weighted Bernoulli
residual is 1.466e-6. Between 64^3 and 128^3, the primary diameter changes
by 0.22%, the donor diameter by 0.13%, and K by 1.76%. The virial error
falls by about four under each doubling. These are discrete equilibrium
checks, not evidence of converged long-term mass transfer.

The shared-K numerical test checks `P/rho^(5/3)` throughout both reference
stars, their mass ratio, the larger donor, and the L1 surface condition.
The gravity/hydro integration test now evolves a coarse shared-K q=0.7
binary through its first step and checks mass and energy budgets.
The numerical suite passes 12 tests in both HPX and serial builds; the HPX
integration suite passes four and serialization passes three. The distributed
integration suite also passes with two and three HPX localities on one host.

The QueenBee4 preset uses a twelve-separation evolution box, block size 8,
base level 2, maximum level 7, and density refinement only. Its finest
spacing gives 50.9 / 58.9 cells across the primary/donor volume-equivalent
diameters. The compact reference remains two initial separations wide,
with 37.7 / 43.6 cells across those diameters. Hydro refinement cannot
replace reference resolution; both are reported explicitly.

The full preset's [time-zero AMR handoff](binary-scf/dwd-common-k-amr.json)
passed on one local HPX locality with 1,523,712 active hydro cells. Native
AMR gravity gives a virial residual of 0.00412104, compared with 0.00442889
on the reference. The initial period is 88.1957 seconds. This run stopped
at time zero; full-resolution evolution and QueenBee4 execution remain
unverified. The stop-time-dependent mesh padding can change the initial
leaf count when evolution is enabled.

```sh
./release/octoII-3d --config=problem/science/BinaryScf/qb4-dwd.ini \
  --scf.evolveOrbits=0 --output.enabled=off \
  --output.directory=/tmp/scf-common-k-amr --verification.analytic=off \
  --hpx:threads=8 --hpx:bind=none
```

## Historical independent-K q=0.7 model and density AMR

An earlier version of the QueenBee4 input used independent pressure
normalizations for its two n=1.5 stars. It therefore did not represent the
requested common-K DWD model: its lower-mass donor was slightly smaller.
The measurements below are retained as a historical handoff check only.
That version used a Roche-filling donor, a uniform
128^3 SCF box two separations wide, and an evolution box eight separations
wide. Its initial AMR grid has eight cells per block, base level 2, finest
level 6, and density refinement only. The initial hydro handoff was run here
on one HPX locality with evolution stopped at time zero; it has **not** been
timed or evolved on QueenBee4.

| Quantity | Measured value |
| --- | ---: |
| Initial orbital period | 88.09114 s |
| SCF iterations | 19 |
| Accretor diameter | 2.25434e9 cm |
| Accretor width on SCF / finest hydro mesh | 47.4 / 48.1 cells |
| Donor width on SCF / finest hydro mesh | 45.5 / 46.2 cells |
| Reference / hydro virial residual | 0.003621 / 0.002526 |
| Active hydro cells | 1,294,336 |

See [dwd-compact-amr.json](binary-scf/dwd-compact-amr.json) for the reference
and native gravity diagnostics. The two-separation reference contained both
lobes without hitting its boundary in this case. These are volume-equivalent
diameters measured from occupied reference cells; the Roche distortion means
the width along a particular axis differs.

## Reproduce

```sh
cmake --build release --target binaryScfChecks-3d binaryScfIntegrationChecks-3d serializationChecks-3d octoII-3d -j 4
./release/3d/tests/binaryScfChecks-3d --hpx:threads=2 --hpx:bind=none
./release/3d/tests/binaryScfIntegrationChecks-3d --hpx:threads=2 --hpx:bind=none
./release/3d/tests/serializationChecks-3d --hpx:threads=2 --hpx:bind=none
ctest --test-dir release -R '^distributed.binaryScfIntegrationChecks-3d' --output-on-failure

./release/3d/tests/binaryScfChecks-3d --gtest_also_run_disabled_tests \
  --gtest_filter=BinaryScf.DISABLED_ResolutionStudy --hpx:threads=2 --hpx:bind=none

./release/octoII-3d --config=problem/science/BinaryScf/inputs \
  --output.enabled=off --output.directory=/tmp/scf-production64 \
  --verification.analytic=off --hpx:threads=4 --hpx:bind=none

./release/octoII-3d --config=problem/science/BinaryScf/bipolytrope.ini \
  --output.enabled=off --output.directory=/tmp/scf-bipolytrope128 \
  --verification.analytic=off --hpx:threads=8 --hpx:bind=none
```

Serial minimum-physics configuration:

```sh
cmake -S . -B build-scf-serial -DCMAKE_BUILD_TYPE=Release \
  -DOCTOTIGERII_WITH_HPX=OFF -DOCTOTIGERII_WITH_PROFILING=OFF \
  -DOCTOII_WITH_RADIATION=OFF -DOCTOII_WITH_MASS_FRACTIONS=OFF
cmake --build build-scf-serial --target binaryScfChecks-3d binaryScfIntegrationChecks-3d -j 4
./build-scf-serial/3d/tests/binaryScfChecks-3d
./build-scf-serial/3d/tests/binaryScfIntegrationChecks-3d
```

GCC 13.3 was used; the existing HPX release build has `-march=native`, while
the serial build uses ordinary Release flags. Before explicit y/z reflection
projection was added, these builds could follow different transverse drift
modes during accelerated iteration. Both now pass; bitwise identity across
compilers is not required. Distributed tests need local network sockets, which
were disallowed by the initial sandbox and then tested with permission outside
it.

## Not yet established

No multi-orbit quiet-binary evolution, mass-transfer convergence, cold DWD
mass-radius relation, reproduction of a named historical binary, or scientific
stability classification is claimed here. The uniform reference construction
has not been benchmarked at machine scale. In particular, these results do not
support a Gordon Bell performance claim.
