# Threaded and distributed FMM

With AMR enabled, the independent sparse cell octree uses a dual tree walk,
including different source and target sizes. See [adaptive gravity](amr.md).
The uniform path described below remains available when AMR is disabled.

## Execution and ownership

The application calls `Runtime::solveGravity()`. A persistent `gravity::FieldSolver`
creates one `FmmPartition` component on every HPX locality. At each cell-octree
level, the global x-contiguous cell indices are divided into nearly equal,
contiguous ranges. Each locality stores only its owned moments and local
expansions. Empty ranges at coarse levels are supported. Leaf moments contain
only mass; higher levels contain the existing `(p+1)^2` coefficients.

This partition is independent of the field-block placement. Contiguous leaf
segments are mapped to field ranges once. Density is fetched directly from those
ranges, including remote ranges when an FMM partition cuts a field block. There
is no full-density gather or replicated global hierarchy. Small coarse levels
can naturally be needed by several localities.

Within a locality, bounded HPX tasks claim batches of 16 targets. By default the
number of tasks is bounded by `hpx::get_os_thread_count()`; `runtime.workerTasks`
can select another bound. No per-cell HPX task is created. One target expansion
has exactly one writer. Static FMM ownership does not use transport work stealing.

## Stages

1. Clear the local hierarchy and read masses from the published field bank.
2. Traverse upward, level by level. Owners fetch required child moments and
   independently translate and sum eight children in fixed slot order.
3. Traverse downward, level by level. Owners fetch the same-level source moments
   needed for M2L (or leaf masses for P2P), and their parent local expansions.
   Each target sums its interactions and then its translated parent local.
4. Evaluate potential and acceleration and write the other field bank. Copy the
   other enabled state fields into that bank without modifying them.
5. Publish the new bank and generation after every locality has completed.

Each level has a completion barrier. Source requests are deduplicated by owner
and split into bounded parcels; independent reads launch before waiting. Local
sources are referenced directly, and remote halos are discarded after the stage.
The persistent hierarchy is reused at the next solve.

All outstanding tasks and transfers are drained before an ordinary numerical or
remote-action exception escapes. The runtime does not change its published bank,
generation, gravity time, or physical time on a failed solve. A later solve clears
the hierarchy before retrying. As with the existing runtime, process loss and
out-of-memory recovery are not fault-tolerant execution modes.

## Interaction equivalence

The opening condition is still `distance * theta > 1` in cell-width units.
For a target at a given level, candidates are the children of source parents
whose separation from the target parent fails that opening condition. A
nonleaf pair is accepted if its current separation passes. Leaves evaluate all
remaining candidates directly, excluding self terms.

For `theta < 1/sqrt(3)`, descendants of an accepted pair remain separated:
`d_child >= 2*d_parent - sqrt(3) > 1/theta`. Consequently an unaccepted parent
also has unaccepted ancestors, and this interaction list reproduces the original
recursive pair traversal. The test suite checks identical unordered pair counts.

The original traversal updates both targets of a pair together. The parallel
solver evaluates the two directed contributions on their respective target
owners. This performs the same two translations/direct contributions and
requires no atomics on coefficients. Pair statistics count only `target < source`,
preserving their previous meaning. Accumulation order can differ from the old
solver. Within the new solver it is independent of thread and locality count.

`gravity::solve` remains the standalone serial reference used by numerical tests.
The application uses `FieldSolver` in both HPX and serial builds; an HPX-free
build executes the same target traversal sequentially.

## Running and testing

For an HPX build, append `--hpx:threads=12` to the existing problem command.
The same executable must be launched on every locality using the installed HPX
launcher and parcelport. The existing `tests/distributed.py` is a local-host TCP
correctness launcher; it is not a cluster submission script.

From HOME, replacing the project path if needed:

```bash
ctest --test-dir ~/workspace/OctotigerII/release \
  --output-on-failure -R 'gravity\.(partitioned|distributed)'

OCTOTIGERII_TEST_LOCALITIES=3 OCTOTIGERII_TEST_THREADS=2 \
python3 ~/workspace/OctotigerII/tests/distributed.py \
  ~/workspace/OctotigerII/release/tests/gravityParallelChecks
```

The console's `cellsPerLocality` reports the number of leaves evaluated on each
locality in the last solve. `workerTasks` counts bounded tasks dispatched across
all FMM stages, including exchange planning; it is not a timing or utilization
measurement. See `VALIDATION.md` for measured test results.

## Remaining limits

Tree-level barriers and sparse coarse levels limit strong scaling. Ownership is
static and balanced by cell count rather than measured interaction cost. Geometry
metadata is still replicated by the existing runtime. Snapshot collection,
diagnostics, Silo output, and the sampled direct-reference calculation remain
centralized and can limit full-application scaling. No GPU path is introduced.
