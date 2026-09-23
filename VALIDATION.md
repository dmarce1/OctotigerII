# Validation of OctotigerII

## Threaded and distributed FMM (2026-09-23)

The application now uses a persistent partitioned FMM hierarchy through
`Runtime::solveGravity()`. Validation uses GCC 13.3.0, Release builds, and the
existing HPX 1.11.0 installation. The standalone serial solver remains available
as an independent traversal reference.

The dedicated `gravity.partitioned` test covers:

- Identical unordered M2L and P2P pair counts relative to the old symmetric solver.
- Potential and all acceleration components within `3e-12` of each component's
  maximum absolute reference value, including zeros and cancellation points.
- Orders 1, 3, and 10; opening angles 0.1, 0.4, 0.5, and 0.57; dense asymmetric
  masses, one isolated mass, and a completely empty mass distribution.
- 4³, 8³, and 16³ meshes; one and multiple field blocks; FMM partition boundaries
  crossing field blocks; empty coarse-level partitions.
- One versus four worker tasks, repeated solves on alternating banks, exact
  preservation of the other state fields, invalid density rejection, draining
  failed work, and a corrected retry with unchanged published gravity.

Two- and three-process TCP CTests run the same checks, with two HPX threads per
process. A separate three-process run with one HPX thread per process also passed.
Every process reported its owned evaluated leaves. For 4,096 leaves, the
three-process split was 1,366 / 1,365 / 1,365. Reported field digests were identical
across one, two, and three processes for every dedicated case, and across one and
four worker tasks. These checks use independent processes on one host.

A separate three-process collapse application run completed a transport step,
both gravity kicks, and a new FMM solve. Its final gravity fields passed the
existing direct-reference checks. The complete suites passed **49/49 CTests**:

| Build | Checks passed |
| --- | ---: |
| HPX gravity sphere 3D | 14/14 |
| HPX collapse 3D | 14/14 |
| Serial gravity sphere 3D | 11/11 |
| HPX Sod 1D (gravity excluded) | 10/10 |

The serial partitioned test was rerun after the final input-validation guard.
Build logs contained no compiler warnings. Doxygen regenerated successfully with
warnings treated as errors. Test and timing logs are in
`docs/validation/parallel-fmm/`.

### Local timing sample

An asymmetric dense 32³ mesh, order 3, theta 0.5, with eight cells per field-block
axis was timed inside `FieldSolver::solve()`. The first solve warms operator
caches; the table is the median of the following three solves. Field transfers,
FMM exchanges, and publication writes are included; initialization of the fixture,
snapshots, output, and direct verification are excluded. HPX affinity was disabled.
The environment provides an eight-CPU quota. All configurations ran on one host.

| HPX processes | Threads per process | Median solve time (s) | Relative to one thread |
| ---: | ---: | ---: | ---: |
| 1 | 1 | 1.058092e+00 | 1.00× |
| 1 | 4 | 3.322818e-01 | 3.18× |
| 2 | 2 | 3.422519e-01 | 3.09× |

Reproduce with `tests/gravityParallelChecks benchmark 2 3` and the desired HPX
thread/launcher settings. The positional benchmark arguments are block level and
multipole order; block width is eight cells. This is a small timing sample, not
cluster scaling evidence. Separate physical-node testing has not been performed.
Snapshots, diagnostics, Silo output, and direct-reference verification remain
centralized. See [parallel FMM](docs/parallel-fmm.md) for execution details.


## Direct gravity verification and orders 1–10 (2026-09-23)

This revision passed 32/32 CTest checks across HPX Release builds: gravity-sphere
3D (11), Sod 1D (10), and collapse 3D (11). Existing installed dependencies were
reused. The gravity build's application and direct-reference checks were rerun
after the final reference mass calculation was aligned exactly with the solver's
`density * cellVolume` evaluation order. No multi-locality test was rerun for
this revision; canonical source ordering and reversed snapshot ordering are
covered by the direct-reference test.

Checks cover order validation (1–10 accepted; 0/11 rejected), independent
point-mass potential/acceleration and self exclusion, vacuum targets and
sources, automatic work limits and explicit count overrides, repeated/changed
seeds, target uniqueness, snapshot-order independence, sample uncertainty and
warnings, JSON output, accuracy gates, and HPX option serialization. Continuum
reference identities and mesh convergence remain tested explicitly in continuum
mode. The 8³ dense-source solver check runs all ten orders against its independent
full direct sum. Relative acceleration RMS decreases from 5.809166e-02 at p=1
to 1.356410e-08 at p=10; potential RMS at p=10 is 6.506875e-11.

The user's 64³ uniform-sphere setup (`mesh.cells=8`, `mesh.level=3`, theta=0.5)
was also run with the default direct reference, output disabled, and four HPX
threads. Both orders used the same 1,024 targets, seed 5489, and all 17,256
nonzero source masses. The sphere mass was 5.266113e+30 g.

| Order | Potential relative L2 | Acceleration-X relative L2 |
| ---: | ---: | ---: |
| 3 | 3.330820e-04 | 5.764977e-03 |
| 5 | 1.966291e-05 | 3.846277e-04 |

At p=5 the approximate 95% sampling half-widths for these L2 values were
4.343574e-06 and 5.961419e-05 respectively. These are estimates over sampled
targets, not bounds on unsampled errors. The direct comparison removes the
continuous-sphere discretization floor from the FMM accuracy measurement.

## Problem and compile-time-dimension revision

Checked with GCC 13.3.0, C++20, Release builds on Linux x86-64. HPX builds use
real HPX 1.11.0, with no substituted runtime or physics. Numerical libraries
compile with `-Wall -Wextra -Wpedantic`; the final matrix builds emitted no
compiler warnings.

Every supported problem/dimension combination was configured, built, and tested:

| Problem | Dimensions checked | CTest checks per build |
| --- | --- | ---: |
| Sod | 1, 2, 3 | 6 |
| Kelvin–Helmholtz | 2, 3 | 6 |
| Streaming | 1, 2, 3 | 6 |
| Radiation pulse | 1, 2, 3 | 6 |
| Gravity sphere | 3 | 6 |
| Gravity Gaussian | 3 | 6 |
| Collapse | 3 | 7 |

The serial matrix passed **85/85 checks across 14 builds**. Three additional
HPX builds—streaming 1D, streaming 2D, and collapse 3D—passed **22/22 CTest
checks**, using two worker threads per locality. Each HPX build also passed
storage and transport checks across two independent TCP processes, with parcel
coalescing enabled and a one-byte zero-copy threshold. Streaming 1D transport
also passed across three TCP processes, including fewer blocks than localities.
These are separate processes on one host, not multi-node or scaling measurements.

CMake rejects KH 1D, gravity sphere 2D, collapse 1D, Sod 4D, and unknown problems.
A fresh KH configuration without an explicit dimension selects its manifest's
2D default. The noninteractive helper validates the same manifests before doing
dependency work. CMake presets, builder argument/help syntax, and standalone
manifest validation were checked. The existing HPX installation was reused;
this revision did not repeat a clean HPX dependency bootstrap.

## What the current checks establish

CTest launches the existing assertion-based regression programs and the Python
application check. Google Test conversion is **deferred**, as is the interactive
builder; see [ROADMAP.md](ROADMAP.md).

- Coordinates contain exactly `ndim` elements; hydro and radiation states occupy
  `ndim+2` and `ndim+1` scalar quantities respectively. There are no inactive axes.
  Cartesian traversal, flattening, child slots/parents, and invalid-axis rejection
  are checked for every compiled dimension.
- The selected problem's initial and evolved values agree between a tiled mesh
  and a single-block reference mesh. Periodic runs conserve each evolved field;
  hydro positivity and M1 realizability are checked. Sod remains planar, and KH
  remains uniform along its extruded z direction in 3D.
- Conservation tolerances are relative to each component's L1 norm. This allows
  floating-point cancellation of opposite physical CGS fluxes whose signed total
  is zero, without imposing an inappropriate dimensionful absolute tolerance.
- The independent 1D streaming profile translates at ĉ while every cell retains
  physical `F=cE`, tested at full and quarter light speed.
- The collapse gravity kick applies its momentum impulse and preserves gas
  internal energy. Gravity compares the discrete point-mass solution with an
  independent direct sum, checks near/far work, and verifies order convergence
  and physical length scaling.
- Generic storage checks include irregular capacities, multiple partitions per
  locality, independent layouts and same-type fields, configurable banks, local
  views, remote transfers, quantity-type rejection, retained-buffer retirement,
  fresh identities, and pointer serialization chunks for density and physical F.
  Sod retains the failed numerical-stage rollback/retry check. Other transport
  problems reject an invalid step and check publication and task accounting.
- HPX archives roundtrip dimension-sized geometry, typed snapshots, and flux
  packets. Runtime problem/dimension overrides are rejected. Application output
  completes its requested time and writes the expected time series.
- Silo readback checks the actual compiled geometry, coordinates, physical units,
  field values, cycle/time, absence of inactive vector components, and repeated
  writes using `DB_CLOBBER`.

Streaming relative L1 energy errors, 128 cells at t=0.2 s:

| ĉ/c | Relative L1 |
| --- | ---: |
| 1 | 1.418577e-03 |
| 0.25 | 4.121192e-04 |

Gravity relative RMS errors, opening angle 0.5:

| Order | Potential | Acceleration |
| --- | ---: | ---: |
| 3 | 2.025507e-05 | 7.413578e-04 |
| 4 | 3.701440e-06 | 1.836427e-04 |
| 5 | 1.558004e-06 | 4.481547e-05 |

## Reproduce

Select one problem/dimension as described in [BUILDING.md](BUILDING.md), build,
and run CTest in that build directory. Repeat for the supported combinations
to cover the whole matrix. For example, from HOME:

```bash
ctest --test-dir ~/workspace/OctotigerII/release/streaming/1d --output-on-failure
python3 ~/workspace/OctotigerII/tests/distributed.py \
  ~/workspace/OctotigerII/release/streaming/1d/tests/storageChecks
python3 ~/workspace/OctotigerII/tests/distributed.py \
  ~/workspace/OctotigerII/release/streaming/1d/tests/numericalChecks transport
OCTOTIGERII_TEST_LOCALITIES=3 python3 ~/workspace/OctotigerII/tests/distributed.py \
  ~/workspace/OctotigerII/release/streaming/1d/tests/numericalChecks transport
```

The Doxygen manual generates with 1.9.8 and warnings treated as errors. All 12
bibliography entries have resolving citation keys. Generated internal HTML
file/anchor links were checked. Both standalone `docs.sh` and the CMake `docs`
target were exercised.

Sanitizer, leak, accelerator, multi-node, throughput, and memory-performance
measurements were not performed for this revision. The tests do not establish
production astrophysical accuracy or mesh convergence of every problem.
Gravity tests compare discrete point masses, not an exact continuum sphere;
global gas-plus-gravity energy conservation is not claimed.

## Earlier revisions

The preceding custom-storage/physical-flux revision passed 13 serial and 14 HPX
CTest checks in its runtime-selected multi-problem executable, plus two- and
three-locality checks. The original implementation also passed historical serial
AddressSanitizer/UndefinedBehaviorSanitizer checks with leak detection disabled.
Those historical sanitizer results do not validate the present revision.

## Standard analytic verification

The `analytic.references` and `analytic.convergence` CTests validate independent
exact solutions and numerical convergence. Supported problems are Sod, uniform
sphere gravity, spherical-cutoff Gaussian gravity, and periodic streaming.
See [analytic verification](docs/analytic-verification.md) for definitions and
limitations. Existing direct-summation gravity checks are retained.

### Verified on 2026-09-23

All 63 CTests passed across seven Release builds: Sod 1D/3D, streaming 1D/3D,
uniform sphere 3D, Gaussian gravity 3D with HPX, and radiation-pulse 2D (the
unavailable-reference path). The other six builds used the serial backend.
HPX tests included configuration serialization. These are single-locality HPX
checks; no new multi-locality claim is made.

Representative relative L1 results from `analytic.convergence`:

| Build | Field | Coarse N | Coarse error | Fine N | Fine error |
| --- | --- | ---: | ---: | ---: | ---: |
| Sod 1D | density | 32 | 2.552316e-02 | 64 | 1.280111e-02 |
| Sod 3D | density | 16 | 4.090404e-02 | 32 | 2.591659e-02 |
| Gaussian gravity 3D (HPX) | potential | 16 | 4.876158e-03 | 32 | 1.221247e-03 |
| Sphere gravity 3D | potential | 16 | 3.850622e-02 | 32 | 1.304695e-02 |
| Streaming 1D | radiationEnergy | 32 | 2.377071e-02 | 64 | 6.118829e-03 |
| Streaming 3D | radiationEnergy | 16 | 1.318741e-01 | 32 | 3.710528e-02 |
