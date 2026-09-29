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
cd ~/workspace/OctotigerII
git pull --ff-only
./build.sh release -j 12
```

`build.sh release` sets `HPX_WITH_APEX=ON` and
`OCTOTIGERII_WITH_PROFILING=ON`; it reuses the existing Release dependency
builds where possible. Confirm the settings if desired:

```bash
grep '^HPX_WITH_APEX:BOOL=ON' packages/release/hpx/build/CMakeCache.txt
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
cd ~/workspace/OctotigerII
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
the last frame. Let the process exit normally to collect complete APEX reports.

The output directory contains `scf.json`, `conservation.csv`, `frames.visit`,
and numbered `frame_*.silo` files. `profile.sh` creates a separate APEX directory
for each locality. Inspect the TAU-format `profile.*` files on all four nodes;
the pinned APEX version can leave worker CSV files empty. SCF reports
`scf.solve`, `scf.gravity_fft`, `scf.l1_search`, `scf.common_k_constraint`, `scf.density_update`,
`scf.mixing`, and `scf.handoff_block`. Existing regions cover FMM workers,
hydro transport, diagnostics, and Silo output. The names are useful timing
regions, not a timer around every individual function.

QueenBee4 has 64 CPU cores per ordinary compute node according to
[LONI's system description](https://hpc.loni.org/resources/hpc/). The
`workq` allocation and `salloc` syntax follow
[LONI's current Slurm guide](https://hpc.loni.org/documentation/job-submission/).
The one-locality-per-node `srun` launch follows
[HPX's Slurm documentation](https://docs.hpx.dev/latest/html/manual/running_on_batch_systems.html).

These are nonrelativistic white-dwarf **polytropic proxies** with shared K.
They do not include relativistic degeneracy or a thermal/composition-dependent
white-dwarf EOS. Evolution uses ideal-gas gamma=5/3.
