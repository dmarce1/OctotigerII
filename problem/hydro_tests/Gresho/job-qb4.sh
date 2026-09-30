#!/usr/bin/env bash
#SBATCH --job-name=qb4-gresho-hpxdiag
#SBATCH --account=loni_graphshpx
#SBATCH --partition=workq
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=00:30:00
#SBATCH --output=slurm-gresho-%j.out
#SBATCH --error=slurm-gresho-%j.err
set -euo pipefail
cd /work/dmarce1/OctotigerII
job_id="${SLURM_JOB_ID:?Submit with sbatch}"
log_dir="$PWD/logs/qb4-gresho-$job_id"
executable="$PWD/relwithdebinfo-hpxdebug/octoII-3d"
[[ -x "$executable" ]]
mkdir -p "$log_dir"
cp problem/hydro_tests/Gresho/qb4.ini "$log_dir/input.ini"
cp problem/hydro_tests/Gresho/job-qb4.sh "$log_dir/job.sh"
git rev-parse HEAD > "$log_dir/commit.txt"
sha256sum "$executable" > "$log_dir/executable.sha256"
export OCTOTIGERII_PROFILE_DIR="$PWD/profiles/qb4-gresho-$job_id"
export OCTOTIGERII_PROFILE_INTERVAL_SECONDS=60
printf 'Job %s started at %s\n' "$job_id" "$(date --iso-8601=seconds)"
trap 'status=$?; printf "Job finished at %s, status=%s\n" "$(date --iso-8601=seconds)" "$status"' EXIT
srun --unbuffered --label --kill-on-bad-exit=1 \
    --output="$log_dir/srun-%J-task-%t.out" --error="$log_dir/srun-%J-task-%t.err" \
    -N 4 -n 4 --ntasks-per-node=1 -c 64 --cpu-bind=cores \
    ./profile.sh "$executable" --config="$log_dir/input.ini" \
    --output.directory="$PWD/output/qb4-gresho-$job_id" --hpx:threads=64 \
    --hpx:ini=hpx.stacks.use_guard_pages=0 --hpx:ini=hpx.trace_depth=0 \
    --hpx:ini=hpx.lock_detection=1 --hpx:ini=hpx.throw_on_held_lock=1 \
    --hpx:ini=hpx.minimal_deadlock_detection=1 --hpx:ini=hpx.spinlock_deadlock_detection=1 \
    --hpx:ini=hpx.spinlock_deadlock_detection_limit=10000000 \
    --hpx:ini=hpx.logging.level=3 --hpx:ini=hpx.logging.destination=cerr
