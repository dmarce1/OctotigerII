# Profiling OctotigerII

`build.sh` enables HPX task descriptions, APEX, HPX's PAPI counters, APEX's
PAPI measurements, and HPX worker idle-rate counters. HPX 1.11.0 fetches its
matching APEX source using `HPX_WITH_FETCH_APEX=ON`; the script deliberately
keeps HPX's default APEX repository and tag.
After fetching, `cmake/PatchApex.cmake` supplies a missing `<cstdint>` include
in that pinned APEX version's `gzstream.hpp`. The fix is idempotent and does
not modify an already-fixed header.

PAPI is a separate development library, not part of a standard Linux install.
The default dependency policy is:

1. Check installed headers and libraries, `CMAKE_PREFIX_PATH`, pkg-config,
   and `Papi_ROOT` / `PAPI_ROOT` (also EasyBuild's `EBROOTPAPI`).
2. Try the `papi` environment module, or `OCTOTIGERII_PAPI_MODULE` if supplied.
3. Reuse a previous local installation, or clone and build PAPI 7.2.0 under
   `packages/TYPE/papi`. This needs `make` and a C compiler but no root access.

The probe compiles and links against PAPI; it does not require access to
hardware counters on a cluster login node. The same checked include/library
paths are passed to HPX and APEX. Both `HPX_WITH_PAPI=ON` and
`APEX_WITH_PAPI=ON` are needed for their respective integrations.

Set `OCTOTIGERII_PAPI=system` to require an installed/module copy, or `build`
to use the local source build. `auto` is the default. Load your compiler and
site dependency modules before running the script on a cluster. The script
never invokes a system package manager or changes kernel counter permissions.

## Build and collect a first profile

Commands below start from HOME, assuming the checkout is in `~/workspace/OctotigerII`.

```bash
cd ~/workspace/OctotigerII
./build.sh release  -j 12
./profile.sh ./release/octoII-3d --problem.name=gravity-sphere \
  --problem.name=gravity-sphere --config=bin/gravity_tests/Sphere/inputs \
  --mesh.cells=4 --mesh.level=1 --output.enabled=off \
  --hpx:threads=12
```

`profile.sh` creates a unique directory for each process under
`OCTOTIGERII_PROFILE_DIR` (default: `./profiles`), then runs the supplied command
without changing its working directory. It enables APEX screen output, CSV,
and TAU-format text profiles unless you have supplied the respective
`APEX_SCREEN_OUTPUT`, `APEX_CSV_OUTPUT`, or `APEX_PROFILE_OUTPUT` settings already.
Writing the TAU-format files does not require installing TAU.

For a distributed launch, wrap **each locality's executable** with `profile.sh`,
not the launcher itself. HPX's TCP transport does not use APEX's MPI CSV gather:
each process writes `apex_profiles.csv`, so a shared output directory loses
reports. This wrapper prevents that collision. In our two-locality test the
pinned APEX version also produced empty CSVs on worker localities, while its
TAU-format `profile.LOCALITY.0.0` files contained the timers and samples on both
localities. Use those files for distributed results; the filename identifies
the HPX locality. Locality zero's CSV is convenient for spreadsheet inspection.
Alternatively set a different `APEX_OUTPUT_FILE_PATH` yourself for each locality
and enable `APEX_PROFILE_OUTPUT=1`.

These variables select report formats; HPX's APEX instrumentation
is compiled into this HPX build regardless of whether screen output is enabled.

Hardware events are opt-in at runtime; begin with timing only, then try:

```bash
export APEX_PAPI_METRICS='PAPI_TOT_INS PAPI_TOT_CYC'
```

Run the same executable again. These are requested events, not a guarantee
that a particular CPU, virtual machine, or cluster policy exposes them. Use
`papi_avail` and `papi_native_avail` from your PAPI installation to inspect the
machine. If counter access is denied, unset `APEX_PAPI_METRICS` and collect
ordinary timing data; on a managed cluster consult the administrators rather
than changing permissions in the build script. Do not request the same PMU
events simultaneously through HPX's `/papi/...` counters and APEX unless you
have checked the platform's event-set limits.

List HPX's counters with the executable's `--hpx:list-counters` option. Thread
idle rates are available from this build. Queue wait-time measurements need
additional HPX instrumentation (`HPX_WITH_THREAD_QUEUE_WAITTIME=ON`); that
higher-overhead option is not enabled by this script.

## Instrumentation and interpretation

| Names | Meaning |
|---|---|
| `gravity.clear`, `gravity.p2m`, `gravity.m2m` | Scheduled worker batches for initialization and the upward pass |
| `gravity.m2l_l2l`, `gravity.p2p_l2l` | Current fused downward worker batches, including interaction traversal and L2L |
| `gravity.publish` | Scheduled field-publication batches, including leaf evaluation and storage access |
| `gravity.exchange.plan`, `gravity.exchange.pack`, `gravity.exchange.unpack` | Remote source planning and coefficient packing/unpacking |
| `gravity.operator_setup` | M2L operator construction on cache misses |
| `runtime.*.worker` | Transport, timestep, and gravity-kick worker tasks |
| `hydro.advance`, `radiation.advance` | Synchronous per-block numerical transport |
| `transport.reconstruct_predict`, `transport.fluxes`, `transport.update` | Shared finite-volume kernels; use their parent region to distinguish physics |
| `transport.signal_speed`, `gravity.kick` | Synchronous timestep and gravity-source work |
| `verification.*`, `diagnostics`, `output.*` | Reference calculations, error norms, diagnostics, Silo and JSON output |
| `gravity.serial.*` | Separate serial reference FMM used by verification tests |
| `*.wall_ns` | Sampled elapsed wall time in nanoseconds, including waits and child tasks |
| `gravity.multipole_pairs`, `gravity.direct_pairs`, `gravity.worker_tasks` | Root-locality samples of global work counts per completed solve |

The current distributed FMM fuses interactions and L2L within each target's
work. The labels reflect that implementation; they do not claim isolated M2L
arithmetic timing. The existing pair counters count an unordered source/target
pair once, whereas the partitioned solver evaluates the two directed
contributions separately. Account for that when normalizing instructions.

Task names use `hpx::annotated_function`, allowing HPX/APEX to observe task
suspension and resumption. Nested synchronous regions use exception-safe
`apex::start` / `apex::stop` guards. HPX 1.11's `scoped_annotation` changes
the APEX task name and is not a substitute for independent nested timers.

No explicit APEX region spans an HPX future wait. Larger operations instead
sample a steady-clock duration after completion, with names ending in
`.wall_ns`. Those are latency distributions, not CPU-time totals. Nanoseconds preserve
short durations in this APEX version, whose CSV writer rounds samples to integers.
Nested timers overlap and worker times accumulate across threads: adding
all rows does not yield the program's elapsed time. Profiler overhead and
cold operator caches also affect the first run.

## Optional builds

`OCTOTIGERII_WITH_PROFILING` defaults to `ON`. With an ordinary HPX installation
without APEX, task names still work and APEX regions/samples become no-ops.
An HPX-free build also compiles without APEX or PAPI. Calls made outside an
active HPX task do not initialize APEX on their own.

Set `-DOCTOTIGERII_WITH_PROFILING=OFF` when configuring directly with CMake to
remove OctotigerII's hooks. HPX's own built-in APEX instrumentation remains
unless you also use an HPX built without APEX.

CUDA profiling is not enabled here. When CUDA kernels are added, CUPTI/NVML
support can provide device timing and transfers alongside these CPU labels.
CPU launch-region time must remain separate from actual GPU execution time.
