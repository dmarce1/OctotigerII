# Rotating-grid validation — 2026-09-26

These results cover the rigidly rotating grid, inertial hydro/radiation states,
balanced gravitational rotation work, and the original oblate `rotatingStar`
model. See [the derivation](../../rotating-frame.md) and
[the benchmark documentation](../../rotating-star.md).

## Build and execution scope

Release serial executables were built in
`/tmp/octotigerII-conservation-serial`, with `OCTOTIGERII_WITH_HPX=OFF`.
The existing `release` build uses HPX, profiling, hydro, radiation, gravity, and
material fractions. Serial 1D/2D/3D executables and the HPX 3D executable built.
HPX application tests used one locality and two threads unless specified below.
The distributed check used two independent TCP localities on the same machine.
This is not a cluster performance measurement.

| Selection | Result |
|---|---:|
| 3D serial rotating-gravity suite | 5 passed |
| 3D serial rotating-transport suite | 9 passed |
| 3D serial gravity energy/time integration, finite-volume, conservation, and time-refinement regressions | 47 passed |
| 3D serial configuration / boundaries / Silo | 42 / 9 / 4 passed |
| 1D and 2D serial configuration/transport plus 3D stellar-model CTest selection | 110 passed |
| 3D HPX rotating transport, one locality | 9 passed |
| 3D HPX stellar model / boundaries / Silo / serialization | 24 passed |
| 3D HPX rotating gravity, two localities, excluding the convergence reference | 4 passed |

Nonzero rotation is rejected in 1D; the 1D selection verifies that configuration
contract and the applicable moving-face/zero-rotation tests. This was a targeted
regression run, not a rerun of every test in the project. `git diff --check`
also passed.

Key raw logs are retained here:
[serial gravity](gravity-serial.log),
[two-locality gravity](gravity-two-localities.log), and
[serial transport](transport-serial.log).

To repeat the gravity tests after building their targets:

```bash
/tmp/octotigerII-conservation-serial/3d/tests/rotatingGravityChecks-3d

OCTOTIGERII_TEST_LOCALITIES=2 \
GTEST_FILTER=-RotatingGravityIntegration.TemporalConvergence \
python3 tests/distributed.py release/3d/tests/rotatingGravityChecks-3d
```

The serial gravity suite took about 508 seconds on this host, including about
420 seconds for the convergence reference. Its CTest timeout is 900 seconds.
The two-locality selection took about 127 seconds. Timings include concurrent
work and are not performance benchmarks.

## Conservation and convergence checks

The rotation-work tests exercise full, nested, and rung operators on uniform
and adaptive trees with both direct and multipole interactions. Their largest
normalized balanced-power residual was `1.46e-16`; direct forces agree with
ordinary local `rho w dot g` work to `4.73e-16`. The multipole cases deliberately
have nonzero raw torque power, so they test the cancellation that motivated
the balanced work rather than only central-force cases.

Global, hierarchical, and conventional integration were exercised with both
signs of rotation. The tests check boundary-corrected mass and total inertial
gas-plus-gravity energy, actual fine substeps, applicable momentum balance,
coarsening/refinement followed by continued evolution, and the independent
naive energy option. Mass and total-energy tolerances are `8e-13` and `1e-12`.

The smooth off-axis gravity test holds its adaptive spatial mesh fixed, uses
`omega=0.2`, and compares 2, 4, and 8 coarse steps against 32 steps:

| Integrator | Density orders | Momentum orders | Gas-energy orders |
|---|---|---|---|
| Hierarchical | 2.04609, 2.08513 | 2.04974, 2.08132 | 2.04627, 2.08524 |
| Conventional | 2.04608, 2.08513 | 2.04973, 2.08132 | 2.04627, 2.08524 |

The separate rotating entropy-pulse transport test measured temporal orders
2.02476 and 2.02539. Its normalized temporal errors were `1.08197e-4`,
`2.65889e-5`, and `6.53127e-6` against a 128-step reference. These measurements
test temporal convergence, not spatial resolution of the stellar surface.
Fixed-spatial-grid temporal order was not established here for the legacy
global-step MUSCL–Hancock source/transport/source path; the gravity convergence
table concerns the two refined numerical-RHS midpoint paths.

Transport tests also verify zero-rotation equality with the previous update,
inertial pressure work, moving-contact sampling, relative CFL, uniform inertial
states under changing normals, and radiation faces moving faster than the
reduced light speed. Boundary tests cover the relative gas diode and a vacuum
radiation exterior. Silo checks inspect rotating-grid coordinates and verify
that field values are preserved.

## Short rotating-star runs

The example uses the original Octo-Tiger SCF table, not a spun-up spherical
polytrope. Its independent density-quadrature Bernoulli check has spread
`2.23e-4` in the original code units (test tolerance `2e-3`). The table's source
revision, checksum, and reproducible importer are recorded in
[DATA.md](../../../bin/science/RotatingStar/DATA.md).

The matched cases use the default example: box `[-2e9,2e9]` cm, 4 cells per
subgrid direction, maximum level 3, and a final time of 0.25 s. They begin with
512 level-3 subgrids and regrid to 56 level-2 plus 64 level-3 subgrids. A full
stellar rotation is about 47.1741 s, so these are short execution checks.

The serial corotating hierarchical run is reproduced by:

```bash
/tmp/octotigerII-conservation-serial/octoII-3d \
  --config=examples/rotating-star.ini --runtime.stopTime=0.25 \
  --verification.analytic=off --output.enabled=off \
  --output.directory=/tmp/octo-rotating-star-serial-smoke
```

The matched HPX stationary-grid and conventional runs use:

```bash
./release/octoII-3d --config=examples/rotating-star.ini \
  --runtime.stopTime=0.25 --frame.omega=0 \
  --verification.analytic=off --output.enabled=off \
  --output.directory=/tmp/octo-rotating-star-inertial-validation \
  --hpx:localities=1 --hpx:threads=2 --hpx:bind=none

./release/octoII-3d --config=examples/rotating-star.ini \
  --runtime.stopTime=0.25 --gravity.timeIntegration=conventional \
  --verification.analytic=off --output.enabled=off \
  --output.directory=/tmp/octo-rotating-star-conventional-validation \
  --hpx:localities=1 --hpx:threads=2 --hpx:bind=none
```

The result table below is generated from the retained conservation CSV files.
Energy drift and budget residual use the application's gas-plus-gravity norm;
mass drift includes boundary fluxes. The angular-momentum change is a grid
integral without an escaped-angular-momentum correction.

| Frame / integrator / backend | Steps | Mass drift | Total-energy drift | Energy-budget residual | Peak-density change | Grid Jz change |
|---|---:|---:|---:|---:|---:|---:|
| [Corotating / hierarchical / serial](star-corotating-hierarchical-serial.csv) | 4 | 0.000e+00 | -9.272e-17 | -4.845e-17 | -3.08958% | +0.28757% |
| [Stationary / hierarchical / HPX](star-inertial-hierarchical-hpx.csv) | 3 | 1.470e-16 | -1.391e-16 | -6.543e-17 | -3.08043% | +0.28178% |
| [Corotating / conventional / HPX](star-corotating-conventional-hpx.csv) | 4 | 0.000e+00 | -9.272e-17 | -5.276e-17 | -3.08958% | +0.28757% |

The roughly 3% peak-density decline is also present on the stationary grid.
The default grid samples the initial peak at only `0.919712 rho0`, and the
Cartesian finite-volume state is not a discrete hydrostatic equilibrium of its
cell-centered gravity operator. This comparison identifies a shared spatial
equilibrium/transport error; it does not establish an accurate long-duration
star. Resolution and full-period studies remain necessary for that claim.
Angular momentum is diagnosed, not conserved to roundoff. Korobkin's
force-operator changes are not part of this implementation.

## Radiation and material-fraction application checks

A 2D radiation pulse completed seven steps with `omega=1 rad/s` and reduced
light speed `c_hat=1 cm/s`. Boundary-corrected radiation energy drift was zero
at printed precision; the two integrated flux drifts were below `1.3e-17`.
The [CSV](radiation-rotating-serial.csv) is retained.

```bash
/tmp/octotigerII-conservation-serial/octoII-2d \
  --problem.name=radiation-pulse --mesh.lower=-1 --mesh.upper=1 \
  --mesh.cells=8 --mesh.level=1 --mesh.periodic=off --frame.omega=1 \
  --radiation.lightSpeedRatio=3.3356409519815204e-11 \
  --runtime.stopTime=.1 --verification.analytic=off --output.enabled=off \
  --output.directory=/tmp/octo-rotating-radiation-smoke
```

A one-locality HPX corotating star also completed a step with material fractions
enabled. Boundary-corrected mass and total-energy drift were zero at printed
precision, and the energy-budget residual was `-8.75e-17`.
This is a configuration/execution smoke test, not a separate species-convergence
study. Its [CSV](star-species-hpx.csv) is retained.

```bash
./release/octoII-3d --config=examples/rotating-star.ini \
  --mesh.lower=-1.5e9 --mesh.upper=1.5e9 --amr.maxLevel=2 \
  --runtime.stopTime=.05 --massFractions.enabled=on \
  '--massFractions.species=fuel:0.7:He;ash:0.3:O' \
  --verification.analytic=off --output.enabled=off \
  --output.directory=/tmp/octo-rotating-species-smoke \
  --hpx:localities=1 --hpx:threads=2 --hpx:bind=none
```

Nonzero rotation currently requires two or three dimensions and outflow/free
boundaries on all faces. Gravity uses isolated boundaries and remains 3D. Two
extra coordinate-weighted field solves are required for each relevant gravity
field evaluation; zero rotation skips that work.
