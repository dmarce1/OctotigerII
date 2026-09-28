# Cached numerical midpoint validation

Date: 2026-09-28. This implementation restores the original numerical midpoint
stages and caches repeated halo inputs. The original time-order assertions are
unchanged. See [3D results and commands](gravity-serial-3d.md) and
[halo accounting](../../halo-communication-accounting.md).

## Storage and validity

- Coupled radiation: the initial two-layer gas/radiation halos are fetched once.
  After boundary processing, their first physical ghost layer is retained for
  the corrector's original-state limiter and update base.
- Midpoint gas/radiation halos are each fetched once into reusable worker
  buffers. Hydro and radiation's material evaluation share the gas halo.
- Gravity without coupled radiation retains the complete raw initial halo,
  including AMR donor slots. Raw initial states and numerical rates are combined
  before nonlinear prolongation, preserving the original midpoint construction.
- Each executor preallocates a buffer pool for its owned blocks plus its worker
  count. Stolen work can grow the pool; returned buffers retain their capacity.
  Regridding constructs a new pool for the new topology. Worker scratch vectors
  retain capacity after first use.
- Probe work may be stolen. Each corresponding corrector is assigned to the
  locality holding its initial cache, without migrating cached state.
- The cache key includes physical time, source-rate identity and interval, input
  bank, and donor time interpolation. New probes, rollback, and mesh replacement
  invalidate it. A global opening force kick refreshes initial gas data while
  reusing unchanged radiation across the bank copy.

Only halo storage and scheduling change: accepted numerical fluxes, source
integration, gravity work, AMR reflux, and endpoint forecasts retain the
original method. The owning AMR shadow still uses the original four-layer local
construction; it does not exchange distributed halos.

## Serial checks

| Check | Result |
| --- | --- |
| Original gravity and rotating-gravity temporal gates, 3D | All pass; orders 2.04–2.09 |
| Radiation integration, 1D | 8/8 pass, including temporal order, AMR, regrid, rollback |
| Focused radiation integration, 3D | 8/8 pass, including coupled self-gravity and rotating budgets |
| Radiation diffusion, 1D | 8/8 pass, including the forced opaque equilibrium |
| Halo requests and full/compact-halo equivalence, 1D and 3D | 6/6 pass in each dimension |

The full/compact comparison uses identical midpoint data and compares every
accepted face flux and cell update exactly, including nonuniform AMR,
reflecting/analytic boundaries, and radiation realizability fallback.

Additional reproduction commands from the repository root:

```sh
cmake --build /tmp/octoii-radiation-serial --target haloTrafficChecks-1d radiationIntegrationChecks-1d radiationDiffusionChecks-1d -j2
/tmp/octoii-radiation-serial/1d/tests/haloTrafficChecks-1d
/tmp/octoii-radiation-serial/1d/tests/radiationIntegrationChecks-1d
/tmp/octoii-radiation-serial/1d/tests/radiationDiffusionChecks-1d
```

Logs: `halo-1d.txt`, `radiation-integration-1d.txt`,
`radiation-diffusion-1d.txt`, and the 3D logs linked above. The serial build is
Release with HPX disabled and profiling enabled.

## HPX checks

`haloTrafficChecks-1d` passes all six tests with two HPX threads on one locality
and on two independent TCP localities on this host. Application profiling is
disabled in this Release build, exercising the HPX compile guards in that
configuration. Successful logs are `halo-hpx-1d.txt` and
`halo-hpx-2localities-final-1d.txt`.

```sh
cmake --build /tmp/octoii-radiation-hpx --target haloTrafficChecks-1d -j2
/tmp/octoii-radiation-hpx/1d/tests/haloTrafficChecks-1d --hpx:threads=2 --hpx:bind=none
python3 tests/distributed.py /tmp/octoii-radiation-hpx/1d/tests/haloTrafficChecks-1d
```

The distributed launcher requires permission to open local TCP sockets. With
stealing disabled, the coupled fixture requests 1,024 remote numerical payload
bytes per step, increasing to 1,408 bytes after an opening force kick. The
all-remote uncached baseline would be 2,304 bytes from its measured scalar count;
that remote baseline was not run. These are same-host correctness measurements.
Stealing is enabled in one schedule but actual probe migration is not forced.

An earlier run used gtest `ScopedTrace` across HPX suspension points. Its
thread-local trace stack is incompatible with coroutine migration and crashed
on destruction. That test-only trace was removed. The superseded distributed
log `halo-hpx-2localities-1d.txt` records that failed attempt; the `final` log
above records the corrected successful run.

## Request reduction

| Fixture, scalar values per unforced coupled step | Original | Cached | Reduction |
| --- | ---: | ---: | ---: |
| 1D | 288 measured | 128 measured | 55.56% |
| 3D | 186,368 calculated | 86,016 measured | 53.85% |

An opening force kick increases cached counts to 176 and 114,688 respectively.
The rejected first-order local predictor used half the unforced cached count;
the midpoint exchange is retained to meet the original accuracy requirement.
These tests measure requested numerical halo values, not whole-application
network bytes, achieved bandwidth, or time to solution.
