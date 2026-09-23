# Distributed fields and execution {#storage_design}

`storage::StoragePartition` owns named typed scalar columns in an HPX component.
`storage::Field<T>` allocates a field across a shared `PartitionSet`; a
`FieldHandle<T>` carries its identity and address directory. No partitioned-vector
container is used. The serial backend uses the same component logic through local
shared pointers. The store has no mesh type, level, neighbor, parent, or child.

Hydro uses `ndim+2` columns and radiation `ndim+1`. Gravity uses four in its
3D-only builds. Prescribed gravity
density is a separate column. Radiation stores **(E, F)** with exactly `ndim` flux components: E is
`EnergyDensity`, and every F component is `EnergyFlux`, in erg/(cm² s).
All stored and transmitted scalars remain Boost CGS quantities.

## Addressing and ownership

`storage::Layout` assigns arbitrary positive record counts to indivisible
`Range{partition, offset, count}` allocations. Partition capacities are their
actual occupied totals, so an uneven decomposition needs no maximum-capacity
padding. Each field has its own layout and configurable bank count. The transport
runtime currently requests two separately allocated banks per scalar field.
A range can represent mesh cells or particle attributes without changing the store.

The application shares one partition component per locality across all enabled
fields. The generic `PartitionSet` also permits several partitions on one locality.
Fields with related indexing share a layout; unrelated layouts may use the same
component directory. Two fields with the same quantity type get different IDs
and different arrays. Registration of a C++ type never aliases its physical fields.

`FieldHandle<T>` copies only metadata and component handles. Lookups validate the
field's actual C++ quantity type, bank, and element range. The component holds a
short metadata lock during lookup; kernels operate outside it. Callers enforce
immutable inputs and disjoint output ownership. The low-level API intentionally
does not choose a mesh traversal or a timestep policy.

Each acquired local buffer retains the selected bank's shared allocation.
Retirement removes a field from future lookup without invalidating existing
buffers. IDs are never reused within a partition set. `Field::retire()` must be
called after outstanding writers/transfers drain. New allocation/layout versions
use a fresh field, then retire the old one at a synchronized boundary. In-place
resizing, automatic compaction, and data migration are not implemented.
Dropping a `Field` alone does not retire it: handles may still reference it.
Storage is ultimately released with the last component handle or explicit retirement.

One `storage/types.def` list declares supported wire quantity types; `src/storage.cpp`
registers the component and its actions. Add a new distinct scalar type there to
make it remotely transferable. New fields of an existing type need no new action
registration or source file. The serial backend can instantiate other scalar types
directly. Do not opt arbitrary nontrivial types into bitwise serialization.

`Subgrid` now contains only an ID, geometry and an interior range. The fixed
Cartesian topology adapter owns the block directory and maps coordinates to
storage ranges. It is the only part of execution that interprets `mesh.level`.
The old owning-subgrid components and whole-directory halo-exchange routines
have been removed. No per-block HPX component owns physics data. A different mesh adapter can
reuse the store and supply its own geometry and halo plan. Refinement level
is deliberately absent from the storage key; an adapter can group allocations
by level without requiring a different field implementation.

The current application still supplies uniform, fixed-level Cartesian meshes.
Dynamic refinement, mixed-level transport with refluxing, particle dynamics,
and ownership migration are not implemented by this storage change. Regridding
would build new allocation and topology generations at a completed-stage
boundary, transfer valid interiors, and retire the old handles after readers
finish. It need not change the generic field representation.

## Checkout and publication

Each runtime operation holds a stage lock. A transport stage reads the published
bank and reserves disjoint ranges in the other bank. Each block is claimed once
from a locality's queue. A worker:

1. Acquires views of its input columns.
2. Starts every independent halo read before awaiting the halo results.
3. Constructs a `PatchView` from the interior columns and compact ghost buffers.
4. Runs the unsplit MUSCL-Hancock kernel into the output columns.
5. Waits for any remote writeback before reporting completion.

Local interior inputs and outputs are direct views into field partitions.
`PatchView` resolves interior and ghost indices without copying the interior
into a padded array. Predictor and flux work arrays, and compact halo buffers,
are reused by each worker across blocks and stages. Their count is bounded by
`runtime.worker_tasks` per locality (zero selects the HPX worker-thread count,
limited by total block count). There are no persistent per-block ghost values.

A stage publishes its output bank only after every locality finishes, every
writeback completes, and completed task counts cover the directory. Failures
drain outstanding tasks and reads, leave the published bank and time unchanged,
and permit retry. A separate monotonically increasing dispatch token prevents
stale queue claims. Gravity kicks use the same transaction. Gravity assignments
validate every input, fill the unpublished bank, and then publish it.

Raw field handles deliberately expose low-level reads and writes. Their caller
must enforce immutable readers and exclusive writers. `Runtime` provides that
policy; independent users of the generic store must provide their own policy.
Checkpoint readers would select a published bank under this same stage lock.
This revision does not add a checkpoint file format or restart implementation.

## Locality scheduling

The coordinator creates one executor component per locality, initializes the
field repository, and starts each phase. Each executor owns a queue of block
IDs, immutable geometry, and a bounded set of worker workspaces. Workers consume
local work first, then ask other executors for unclaimed work when
`runtime.work_stealing=on` (the default). Queue claims are brief critical sections;
no lock is held across a network operation or numerical kernel.

Stealing executes a task remotely and sends its result back to the persistent
owner. It does not migrate storage. Set `runtime.work_stealing=off` to measure
whether remote execution pays for a particular workload. Stage scheduling
requires a bounded number of coordinator actions per locality, rather than one
coordinator action for each block. Local workers schedule individual blocks.

Snapshots contain compact interiors only. Transport never broadcasts the global
snapshot directory. Output/diagnostics and the present gravity driver still
gather snapshots on the coordinator. The mathematical tree inside the FMM is
independent of field ownership and remains part of that gravity algorithm.
Distributed gravity and distributed diagnostic reductions are separate work.

## HPX communication and zero-copy

A topology halo plan lists exact contiguous source runs and their destination
ghost indices. It includes periodic wrapping, outflow clamping, edge and corner
cells. It is built once for owned blocks; a stolen block's plan is temporary.
There is no search through the global field snapshots during an update.

Read actions return `hpx::serialization::serialize_buffer<T>` referencing the
source bank. A custom deleter retains its shared allocation until serialization
and transmission release the buffer. Eligible CGS quantity arrays are explicitly
marked bitwise serializable. HPX can send their storage as pointer chunks and
receive into the final transfer buffer without intermediate archive copies.
The implementation never assumes that zero-copy eliminates network traffic or
that a remote pointer can be dereferenced on the caller's locality.

Field read/write actions opt into HPX's coalescing handler when available. To
load it, supply `--hpx:ini=hpx.parcel.message_handlers=1`. HPX controls the
coalescing thresholds and parcelport settings. An HPX build without that plugin
uses its normal parcel path. This code adds no independent message batching
layer. Existing input buffers remain immutable for the full transfer lifetime.

Ghost assembly still copies received values into compact ghost workspaces.
A remotely executed task also copies its result into the owner's preallocated
output range. Those are explicit layout/writeback operations; the local interior
path has neither copy. Many short halo runs may be below HPX's default zero-copy
threshold, and packing could outperform them on some hardware; no unmeasured
performance advantage is claimed.

Small inline geometry arrays are serialized element by element. HPX's deferred
zero-copy receive can otherwise retain an address into a temporary object that
is subsequently moved while deserializing a vector of descriptors or snapshots.
Field payloads remain eligible for zero-copy. The distributed test deliberately
lowers the threshold to one byte to exercise this distinction.

## Validation commands

Run the ordinary CTest suite for the selected backend. With an HPX TCP build:

```bash
cd ~/workspace/OctotigerII
python3 tests/distributed.py release/sod/1d/tests/storageChecks
python3 tests/distributed.py release/sod/1d/tests/numericalChecks transport
```

These launch two independent local processes with coalescing enabled and a
one-byte zero-copy threshold. They require available loopback sockets.
`OCTOTIGERII_TEST_LOCALITIES=3` also exercises uneven placement and empty
localities. `OCTOTIGERII_TEST_COALESCING=0` and
`OCTOTIGERII_TEST_CHUNK_THRESHOLD=4096` select alternative transport settings.
The one-byte threshold is a stress-test setting, not a performance recommendation.
