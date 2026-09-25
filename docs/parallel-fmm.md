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
only mass; higher levels contain `(p+1)^2` source coefficients. Each target has
an order-p scalar local and a separate order-(p+1) auxiliary force local. Image
boundaries add one background/trace coefficient to each expansion.

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
   needed for M2L (or leaf masses for P2P), and both parent local expansions.
   Each target sums its interactions and then its translated parent locals.
4. Evaluate potential and acceleration and write the other field bank. Copy the
   other enabled state fields into that bank without modifying them.
5. Publish the new bank and generation after every locality has completed.

Each level has a completion barrier. Source requests are deduplicated by owner
and split into bounded parcels; independent reads launch before waiting. Local
sources are referenced directly, and remote halos are discarded after the stage.
The persistent hierarchy is reused at the next solve.

`FieldSolver::solve(FieldSolveRequest)` uses the same hierarchy for an
independent partial field. It accepts a density handle, per-block weighted bank
selections, an output handle/bank, and a target mask. Zero source weights exclude
a block. Two banks can interpolate a density or form a signed density difference;
the latter requires `allowSignedDensity=true`. Target subtrees with no selected
leaves are skipped. This overload writes only its selected gravity output and
does not publish the runtime bank or copy other state fields.
Partial-field pair statistics count directed interactions because a source and
target mask need not be symmetric.

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

## Scalar reciprocity and mutual force

These are separate properties. A reciprocal scalar potential operator gives
the energy coupling its discrete binding-energy identity. A force operator must
also exchange equal and opposite total forces between two source groups.

Let A and B be an accepted node pair, R the separation of their geometric
centers, and M_A^alpha and M_B^beta their mass moments, including the usual
multi-index factorial normalization. With the even pair kernel g(R), the
order-p scalar interaction uses

```
U_AB = sum[|alpha|<=p, |beta|<=p]
       (-1)^|beta| M_A^alpha M_B^beta D^(alpha+beta) g(R).
```

Exchanging A and B changes R to -R and preserves this expression. The corrected
adaptive traversal presents the same accepted node pairs in both directions,
so the scalar operator retains this reciprocity on a fixed mesh.

Differentiating only an order-p target local leaves target force degree p-1
while the source still has degree p. Those unequal ranges do not preserve
action/reaction. The force calculation therefore uses the same source moments
through degree p with an auxiliary local through degree p+1. Its gradient has
target degree p, giving

```
F_(A<-B),d = -sum[|alpha|<=p, |beta|<=p]
             (-1)^|beta| M_A^alpha M_B^beta D^(alpha+beta+e_d) g(R).
```

The derivative has parity (-1)^(|alpha|+|beta|+1). Reversing the pair and
exchanging the two multi-indices therefore negates the force term by term.
This applies to each accepted pair and to any symmetric collection of pairs,
including a hierarchical interaction shell. For a kick, the target masses must
match the masses used when those targets serve as reciprocal sources. Reflecting
image boundaries can exert a net physical wall force; zero total force is the
appropriate test for isolated or fully periodic boundaries.

The scalar locals and their operators are unchanged. Only accelerations use the
auxiliary locals. Their Newtonian plane-wave quadrature integrates source-plus-
local degree 2p+1, using p+1 radial points and 2p+2 angular points. Ewald auxiliary
M2L uses derivatives through 2p+1 and retains the same explicit background
trace. `gravity.multipoleOrder=1..10` still selects source/scalar order p;
internal order 11 is used only for the highest-order auxiliary local. The
separate force expansion adds local storage, parent-local communication, and a
second M2L/L2L evaluation without increasing source moments or changing pair
acceptance.

`partialGravityChecks` independently tests scalar superposition/reciprocity,
partial-output isolation, signed sources, and full/shell/cross-only force
balance. On its fixed p=5, theta=0.5 fixture, the normalized shell-force defect
fell from 9.03e-7 to 1.42e-16 on the uniform mesh and from 2.07e-7 to 4.13e-17 on
the adaptive mesh. The normalization is the maximum component of net force
divided by the sum of absolute cell-force components. Direct-force cases also
balance to roundoff. A separate small-node regression covers source orders 1,
5, and 10 with both Newtonian and Ewald auxiliary operators.

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
