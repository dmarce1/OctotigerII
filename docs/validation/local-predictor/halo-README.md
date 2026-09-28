# Local predictor validation records

These records were collected on 2026-09-28 while converting coupled transport to
the local predictor with two ghost layers. They establish the scope of the halo
request measurement and selected numerical regressions. They are not network
bandwidth or application performance benchmarks.

## Communication fixture

`HaloTraffic.CoupledStepReadsOneInitialGasAndRadiationHalo` uses the homogeneous,
periodic `photon-source` problem, four cells per block, spatial level one, and
two complete coupled timesteps. AMR, gravity, species transport, and work
stealing are disabled. The test derives the expected requested ranges and remote
ownership from the two-layer donor plans. The numerical timestep and physical
interval are the same before and after conversion. Each row below gives counts
for one step; the second step repeats them.

| Case | Halo calls | Scalar range requests | Scalar values | Remote payload bytes |
| --- | ---: | ---: | ---: | ---: |
| Serial 1D, measured before | 16 | 72 | 288 | 0 |
| Serial 1D, measured after | 4 | 16 | 64 | 0 |
| HPX 1D, one locality, measured after | 4 | 16 | 64 | 0 |
| HPX 3D, one locality, measured after | 16 | 168 | 43,008 | 0 |

The one-dimensional measured reduction in requested scalar values is **77.78%**.
The three-dimensional pre-conversion count, **186,368**, is calculated from the
old sequence of five gas and three radiation halo reads and the unchanged donor
plans. The corresponding **76.923%** reduction uses this calculated baseline;
there is no successful measured three-dimensional baseline in this record.

Stored gas has `ndim + 5` scalar components and radiation has `ndim + 1`. Gas
includes the auxiliary, nuclei, and electron columns in addition to density,
momentum, and energy. Counts follow every actual scalar source read, including
derived density columns, repeated ranges, and AMR time endpoints. The fixture
excludes derived species columns by disabling species transport.

- [Before, serial 1D](halo-before-serial-1d.log): the new one-halo assertions are
  deliberately run against the old runtime and fail with the measured baseline.
- [After, serial 1D](halo-after-serial-1d.log): all four accounting tests pass.
- [After, HPX 1D](halo-after-hpx-single-1d.log): all four accounting tests pass.
- [After, HPX 3D](halo-after-hpx-single-3d.log): all four accounting tests pass.

The production two-locality run is **unverified**. The sandbox rejects the local
TCP sockets used by `tests/distributed.py`, and escalation did not complete.
All production runs in this record therefore have zero remote bytes. A
preliminary two-locality accounting unit run did complete earlier and counted
eight remote scalar values, or 64 numerical buffer bytes, in its endpoint test;
that is not a distributed production comparison.

On a host permitting local TCP sockets, the registered production test can be
repeated with:

```sh
python3 tests/distributed.py /path/to/build/1d/tests/haloTrafficChecks-1d
```

See [the accounting definition](../../halo-communication-accounting.md) for
counter semantics and exclusions. In particular, numerical halo payload omits
message overhead, interior transfers, reflux, gravity solves, and task migration.

## Numerical regressions

[HPX one-dimensional radiation integration](radiation-integration-hpx-single-1d.log)
records eight passing tests: periodic conservation at physical and reduced light
speed, temporal refinement, acceleration and radiation momentum work, AMR
subcycling and regridding, zero-opacity equivalence, invalid-step rejection,
rollback, and equilibrium exchange. The one-dimensional filter omits the test
whose rotating-grid setup requires higher dimensions.

The temporal refinement test's fine-pair error ratios are approximately 4.28 on
the fixed mesh and 4.36 with AMR. These establish the tested regime's temporal
behavior; they do not establish uniform second order for arbitrarily stiff
sources or prove spatial diffusion accuracy. Those questions require separate
local predictor, diffusion, and physical-boundary tests.
