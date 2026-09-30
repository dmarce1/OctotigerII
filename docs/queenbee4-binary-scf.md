# Four-node QueenBee4 DWD proxy run

Use [the q=0.7 input](../problem/science/BinaryScf/qb4-dwd.ini). It constructs
the donor at L1 and evolves the ideal-gas n=1.5 polytropic binary for 0.2 of
the **initial** orbital period. `scf.framesPerOrbit=100` writes at time zero
and at each 0.01 initial orbit, producing 21 Silo frames if the run completes.
The 128^3 uniform SCF box is two separations wide. The density-refined AMR
evolution box is twelve separations wide, with a 32^3 base mesh and 1024^3
finest equivalent spacing around the stars. Mass and shadow refinement are
disabled; proper nesting and one buffer cell remain active.

The shared-K constraint solves the accretor surface while retaining the donor
at L1. Both stars use `P=K*rho^(5/3)`. The common dimensional K follows from
the specified 0.6-solar-mass accretor and 3e9 cm separation; the input does not
prescribe a composition-dependent electron-degeneracy constant.

The box and maximum AMR level give 50.9 finest cells across the
accretor's volume-equivalent diameter and 58.9 across the donor. The uniform
reference has 37.7 and 43.6 cells across them; refining the hydro mesh does not recover missing reference detail.
See [the validation record](validation/binary-scf.md) for both resolutions
and the current handoff diagnostics.

## Update and build in the current compute session

```bash
cd /work/dmarce1/OctotigerII
git pull --ff-only
./build.sh release -j 12
```

`build.sh release` sets `HPX_WITH_APEX=ON` and
`OCTOTIGERII_WITH_PROFILING=ON`, and disables mmap-backed HPX coroutine stacks
after a QueenBee4 stack-mapping failure. It reuses the existing Release
dependency builds where possible. Confirm the settings if desired:

```bash
grep '^HPX_WITH_APEX:BOOL=ON' packages/release/hpx/build/CMakeCache.txt
grep '^HPX_WITH_THREAD_STACK_MMAP:BOOL=OFF' packages/release/hpx/build/CMakeCache.txt
grep '^OCTOTIGERII_WITH_PROFILING:BOOL=ON' release/CMakeCache.txt
```

Log out of that session after the build, then request four CPU nodes for two
hours. Supply your own LONI allocation name:

```bash
ssh qbd.loni.org
salloc -A YOUR_ALLOCATION -p workq -N 4 -n 4 -c 64 -t 02:00:00
```

Inside the allocation, launch one HPX locality per node. HPX detects the
Slurm hosts; there is no need to provide a nodefile or AGAS address.

```bash
cd /work/dmarce1/OctotigerII
export OCTOTIGERII_PROFILE_DIR="$PWD/profiles/qb4-dwd-$SLURM_JOB_ID"
srun -u -N 4 -n 4 --ntasks-per-node=1 -c 64 --cpu-bind=cores \
  ./profile.sh ./release/octoII-3d \
  --config=problem/science/BinaryScf/qb4-dwd.ini \
  --output.directory="output/qb4-dwd-$SLURM_JOB_ID" \
  --verification.analytic=off --hpx:threads=64
```

Each locality constructs the same compact SCF reference independently; the
evolution uses distributed AMR gravity and hydro. SCF finishes before the first timestep; the positive `scf.evolveOrbits=0.2`
in the input starts evolution immediately afterward. Increase it on a later
run if the first run leaves enough wall time. Use a distinct output directory
for each run. The current code has no restart from Silo, so a job terminated
at the allocation limit leaves its completed frames but cannot continue from
the last frame. `profile.sh` now saves cumulative APEX snapshots every 60
wall-clock seconds on every locality, including during SCF. A normal exit
still produces the final reports; a walltime kill leaves the last completed
snapshots. Change the interval with `OCTOTIGERII_PROFILE_INTERVAL_SECONDS`
(0 disables periodic snapshots). See [profiling snapshots](profiling.md#profiling-snapshots-while-running).

The output directory contains `scf.json`, `conservation.csv`, `frames.visit`,
and numbered `frame_*.silo` files. `profile.sh` creates a separate APEX directory
for each locality. Inspect the TAU-format `profile.*` files on all four nodes;
the pinned APEX version can leave worker CSV files empty. SCF reports
`scf.solve`, `scf.gravity_fft`, `scf.l1_search`, `scf.common_k_constraint`, `scf.density_update`,
`scf.mixing`, and `scf.handoff_block`. Existing regions cover FMM workers,
hydro transport, diagnostics, and Silo output. The names are useful timing
regions, not a timer around every individual function.

## Progress diagnostics and coarse time grouping

Add `--runtime.progressLevel=7` to the executable arguments in the batch script
to see substep progress through every occupied time level. These flushed lines
appear in `logs/qb4-dwd-JOBID/srun-JOBID.0-task-0.err` and include the active
gravity/transport stage, physical substep time and duration, parent interval
fraction, and elapsed wall time. A smaller value reports only that time level
and coarser. The default `-1` disables progress lines.

`--timestep.coarseLevel=5` groups spatial levels 0 through 5 at the smallest
CFL/acceleration limit among them, with levels 6 and 7 subcycling above that
group. The default `0` retains the original schedule. Compare equal physical
intervals when measuring this option: fewer HOLD shells may help, but the
grouped coarse cells advance more frequently. See [time refinement](time-refinement.md).

## Postmortem: job 1060438, 2026-09-29

The four-node run used commit `1eb13ce0756fdb08520288b5cf76b98f39f26e94`
with `HPX_WITH_THREAD_STACK_MMAP:BOOL=OFF`, following successful rebuild job
1060378. Slurm records 19:14:31–23:15:01 America/Chicago, elapsed 04:00:30,
state `TIMEOUT`. Step 1060438.0 was cancelled with signal 15; its reported
maximum task RSS was 23,443,608 KiB (about 22.36 GiB).

SCF converged after 26 iterations. Initialization populated 3,312 leaf blocks
on four localities and wrote `frame_000000.silo` at about 19:23. At termination,
`conservation.csv` contained only step 0 at time 0. This file is flushed at
each completed global step, so no complete evolution step was recorded.
The application logs contain no recurrence of the earlier stack-mapping
failure. This run therefore establishes that disabling mmap avoided that
observed failure for four hours, but it did not resolve the lack of completed
evolution steps.

The last two retained APEX snapshots on each locality show the following
approximate increases over a nominal 60-second interval. APEX prints large
counts with limited precision; these are differences of rounded counters.

| Locality | Counter | Increase |
|---|---|---:|
| 0 | `schedule_parcel` | 275,700 |
| 1 | `potentialRead` | 4,100 |
| 2 | `potentialRead` | 4,100 |
| 3 | `potentialRead` | 112,150 |
| 3 | `momentumRead` | 64,250 |
| 3 | `primary_namespace_increment_credit_action` | 101,060 |
| 3 | `densityWrite` / `momentumWrite` / `energyWrite` | 628 / 471 / 157 |

No recorded increase occurred in `gravity.adaptive.*` or the hydro
`transport.reconstruct_predict`, `transport.fluxes`, and `transport.update`
regions on any locality during those intervals. Communication and field
updates were nevertheless continuing. `schedule_parcel` is HPX receive-side
parcel scheduling, not a send counter. Action names describe wire types:
`MassFlux` aliases `MomentumDensity`, so `momentumRead` includes mass-flux
transfers. HPX task timer counts can include yield/resume accounting and must
not be interpreted as a count of completed timesteps or gravity solves.

The leading code-level explanation is costly serialized field assembly in
the HOLD coordinator. `gravitySourceRate`, `provisionalGravity`, and
`closeGravityFrame` iterate over blocks and call `gravity::fluxWork`.
For each nonphysical interface face, `fluxWork` issues two potential reads;
canonical work also reads a mass flux per face. It waits for one block's
transfers before its caller moves to the next block. A single-endpoint call
currently fetches the same potential twice. Other coordinator routines also
read and write whole-mesh fields one block at a time. This mechanism fits the
observed traffic and absent kernel activity, but the retained profiles do not
identify the live caller or prove that it was progressing through distinct
blocks. A localized stall, repeated work, or extremely slow finite work cannot
yet be distinguished. The cluster performance impact of coarse time grouping
has not been measured.

The next diagnostic run should retain Release optimization and the mmap-off
build, enable progress through level 7, and use a shorter wall limit. First
keep `timestep.coarseLevel=0` to locate the original slow stage. Then compare
grouping on the same mesh over the same physical interval. If a stage remains
active too long, collect its current block index and outstanding transfers or
sample process stacks while the job is alive. If the field-assembly hypothesis
is confirmed, prioritize owned-block execution, batched neighbor reads, and
eliminating duplicate endpoint reads. Changing the time schedule alone does
not establish that the communication cost has been fixed.

The evidence remains under `/work/dmarce1/OctotigerII/{logs,profiles,output}/qb4-dwd-1060438`
on QueenBee4. The compared snapshot pairs are 239/240 on localities 0, 2, and 3,
and 238/239 on locality 1. Accounting can be reproduced with:

```bash
sacct -j 1060438 --parsable2 --format=JobID,Start,End,Elapsed,State,ExitCode,MaxRSS
```

QueenBee4 has 64 CPU cores per ordinary compute node according to
[LONI's system description](https://hpc.loni.org/resources/hpc/). The
`workq` allocation and `salloc` syntax follow
[LONI's current Slurm guide](https://hpc.loni.org/documentation/job-submission/).
The one-locality-per-node `srun` launch follows
[HPX's Slurm documentation](https://docs.hpx.dev/latest/html/manual/running_on_batch_systems.html).

These are nonrelativistic white-dwarf **polytropic proxies** with shared K.
They do not include relativistic degeneracy or a thermal/composition-dependent
white-dwarf EOS. Evolution uses ideal-gas gamma=5/3.
