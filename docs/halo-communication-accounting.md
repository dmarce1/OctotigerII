# Measuring halo requests

`haloReadStatistics()` in `octotigerII/subgrid/view.hpp` returns cumulative
counts on the calling locality. Take a snapshot immediately before an operation
and another after its tasks have completed, then subtract the snapshots. Counters
are available in serial builds and HPX builds, including builds with profiling
disabled. They are updated once per completed halo call with relaxed atomics;
there are no blocking locks or counter resets.

| Counter | Meaning |
| --- | --- |
| `haloCalls` | Calls to assemble a field or state halo, including empty plans |
| `fieldReads` | Requested scalar-field ranges, local and remote |
| `values` | Scalar elements in those ranges, local and remote |
| `remoteFieldReads` | Scalar-field ranges whose owner differs from the requesting locality |
| `remoteValues` | Scalar elements requested from another locality |
| `remotePayloadBytes` | Numerical buffer bytes requested from another locality |

The counters follow issued donor reads. They include duplicate requests, the full
contiguous requested ranges even when only part of a range is copied, both time
endpoints when an AMR donor is interpolated, and every source read used to form a
derived column. Local reads normally retain existing storage; `values` therefore
measures requested work, not network traffic. Multiple partitions on the same
locality do not contribute to remote counts. An asynchronous read that later fails
is still counted as a request. Totals include unsuccessful intervals and predictor
work; measure successful steps separately when comparing production algorithms.

When APEX is enabled, the same quantities are recorded as samples named
`transport.halo.calls`, `transport.halo.field_reads`, `transport.halo.values`,
`transport.halo.remote_field_reads`, `transport.halo.remote_values`, and
`transport.halo.remote_payload_bytes`. Sum the samples to obtain totals. These
are sampled counts, not timers.

Distributed measurements must sum the snapshot differences over all participating
localities, because requests are counted where tasks execute. A snapshot taken
while tasks are active can contain updates from different points in time; take
measurements at synchronization points. Task stealing can change the fraction of
requested values that is remote even when the halo algorithm is unchanged.

These counts exclude request messages, serialization and transport metadata,
message aggregation, reflux, gravity solves, interior transfers, task migration,
and every field read made outside halo assembly. Thus `remotePayloadBytes` is a
count of requested numerical payload, not a measurement of total network bytes
or achieved bandwidth. Network tracing is needed for whole-application byte
counts. Comparing algorithms also requires the same mesh, physical interval,
error tolerance, physics, and placement policy.

The `haloTrafficChecks` test target verifies scalar and multi-column reads,
repeated AMR time endpoints, derived source columns, local versus remote
placement, and empty analytic plans. Its distributed tests use independent HPX
processes on one host to exercise real remote reads; they do not establish
multi-node performance.

## Cached numerical midpoint measurements

The production test
`HaloTraffic.CoupledStepReusesInitialAndMidpointHalosAcrossWorkerSchedules`
measures two complete coupled steps on the homogeneous, periodic `photon-source`
problem with four cells per block, spatial level one, no AMR, no species, and no
gravity. It checks one initial gas halo, one initial radiation halo, one midpoint
gas halo, and one midpoint radiation halo per block per step. Reusing those
values during the remaining predictor and corrector work must not issue more
donor reads. The plans still describe two ghost layers.

The test repeats this measurement with one and two worker tasks, then with work
stealing enabled. Total requested ranges and values must be independent of these
scheduling choices, and the accepted cell states must be identical. When stealing
is disabled, it also checks the exact remote fraction from block ownership and
donor plans. When stealing is enabled, the execution locality can change that
fraction. The test does not guarantee that a task is actually stolen; a run with
no stolen tasks cannot establish cache reuse across an actual locality change.
The second accepted step verifies that a previous step's cached values are not
reused without fetching fresh initial and midpoint data.

Let `G = ndim + 5` and `R = ndim + 1` be the stored gas and radiation component
counts. Gas includes auxiliary, nuclei, and electron columns in addition to
density, momentum, and energy, even with species transport disabled. The earlier
uncached numerical midpoint sequence requests `5 G + 3 R` scalar values per donor
cell. The cached implementation requests `2 G + 2 R`. This retains the initial and midpoint
exchanges while removing repeated requests within those stages.

The cached implementation produced the following counts in 1D and 3D serial
runs on 2026-09-28. The original 1D baseline is measured; the original 3D
baseline is calculated. Logs are in `validation/cached-midpoint/`.

| Fixture, per step | Uncached baseline | Cached measured | Reduction |
| --- | ---: | ---: | ---: |
| 1D halo calls | 16 measured | 8 | 50% |
| 1D scalar range requests | 72 measured | 32 | 55.56% |
| 1D requested scalar values | 288 measured | 128 | 55.56% |
| 3D halo calls | 64 calculated | 32 | 50% |
| 3D scalar range requests | 728 calculated | 336 | 53.846% |
| 3D requested scalar values | 186,368 calculated | 86,016 | 53.846% |

The opening-kick test verifies the necessary gas refresh (`3 G + 2 R`),
uniformity, and the imposed combined gas/radiation momentum change. It measures
176 requested values in 1D and 114,688 in 3D, reductions of 38.89% and 38.46%
against those same uncached baselines. Cache memory is distinct from read
volume: only the first physical layer of the original coupled state is
retained, while each fresh exchange still supplies two layers for reconstruction.

The rejected local method's counts below are half the unforced cached counts.
That extra exchange is retained to recover the required numerical time order.

The one-dimensional uncached baseline was measured on 2026-09-28 and is retained
in `validation/local-predictor/halo-before-serial-1d.log`. Its assertions expected
the earlier local-predictor experiment's smaller count and intentionally failed;
the logged uncached counts remain the comparison baseline. The three-dimensional
uncached counts are calculated from its sequence and donor plans. Neither
baseline establishes whole-application network bytes.

## Superseded local-predictor experiment

The earlier local-predictor implementation requested only the initial halos. It
measured 64 requested scalar values in the one-dimensional fixture and 43,008 in
the three-dimensional fixture on one locality. Those figures, retained in
`validation/local-predictor/halo-after-serial-1d.log` and
`validation/local-predictor/halo-after-hpx-single-3d.log`, describe that experiment.
They are **not** results or targets for the cached numerical midpoint method.

## HPX remote requests

The cached implementation passes all six halo tests with two HPX worker threads
on one locality and on two independent local TCP localities (application
profiling disabled). The successful distributed log is
`validation/cached-midpoint/halo-hpx-2localities-final-1d.txt`.
With stealing disabled, the 1D coupled fixture requests 128 remote values
(1,024 numerical payload bytes) per step, or 176 values (1,408 bytes) with an
opening force kick. The same totals were observed with stealing enabled; that
does not prove a probe was actually stolen.

The historical uncached baseline was measured only locally. Its 288 requested
values imply 2,304 remote payload bytes on this all-remote donor placement,
but that remote baseline has not been measured. No multi-node bandwidth,
whole-application network-byte, or time-to-solution improvement is claimed.
