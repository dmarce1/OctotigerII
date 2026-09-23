# Adaptive mesh refinement

AMR is optional and uses a separate `OctotigerII::Refinement` module for
criteria, `OctotigerII::Amr` for conservative transfers and shadow evolution,
and the mesh/runtime layers for topology, communication, and publication.
The build-time dimension remains 1, 2, or 3; gravity requires 3D.

## Inputs

All levels below are **block** levels. Each leaf block retains `mesh.cells`
cells per axis. The cell width is `(mesh.upper-mesh.lower)/(mesh.cells*2^level)`.
Reduced-dimensional cell masses retain the existing unit transverse measure.

| Setting | Default | Meaning |
| --- | --- | --- |
| `amr.enabled` | `off` | Enable initial and subsequent adaptation |
| `mesh.level` | problem default | Initial uniform block level |
| `amr.minLevel` | `-1` | Minimum level; `-1` uses `mesh.level` |
| `amr.maxLevel` | `6` | Maximum block level, at most 16 |
| `amr.regridEvery` | `4` | Maximum timesteps between refinement checks |
| `amr.maxCellMass` | `0` | Maximum cell mass in grams; zero disables this criterion |
| `amr.shadowTolerance` | `0.05` | Relative fine/shadow difference; zero disables this criterion |
| `amr.shadowFloor` | `1e-8` | Field normalization floor as a fraction of its maximum magnitude |
| `amr.coarsenFactor` | `0.25` | Coarsening threshold relative to the refinement threshold |
| `amr.signalBuffer` | `1` | Safety factor for signal travel, at least one |
| `amr.bufferCells` | `1` | Additional cells around tags |
| `amr.hydro` | `on` | Include hydro components in the shadow criterion |
| `amr.radiation` | `on` | Include radiation components in the shadow criterion |

Example settings for the Rayleigh–Taylor problem:

```ini
amr.enabled=on
mesh.level=1
amr.minLevel=1
amr.maxLevel=3
amr.regridEvery=4
amr.maxCellMass=0
amr.shadowTolerance=0.05
amr.signalBuffer=1.25
amr.bufferCells=1
```

The complete input is `examples/rayleigh-taylor-amr.ini`. After building:

```bash
release/rayleigh-taylor/3d/octotigerII-rayleigh-taylor-3d \
  --config=examples/rayleigh-taylor-amr.ini
```

Open `output/rayleigh-taylor-amr/frames.visit` in VisIt. Rayleigh–Taylor requires
a 3D build. Every AMR frame includes `MetadataIsTimeVarying=1` and
`ConnectivityIsTimeVarying=1`, allowing VisIt to refresh its block list and
connectivity as blocks split or merge. After rebuilding an older output writer,
write a new series and close/reopen the database in VisIt to discard its cache.

Criteria combine by OR. A cell exceeding the mass limit is refined until the
limit is met or `amr.maxLevel` is reached. An entire block is split when any
cell requests refinement. Complete sibling families can merge when both the
fine cells and proposed coarse cells satisfy the coarsening threshold.
Face, edge, and corner neighbors, including periodic neighbors, remain 2:1
balanced. Only active leaves count toward conserved totals and gravity.

## Shared criterion interface

A criterion is `Real(refinement::CellView const&)`. The view supplies position,
cell width/volume, time, level, cell mass, directional signal speeds and
accelerations, and typed
hydro/radiation samples. Each field sample provides its conserved values,
physical gradients, coarse shadow values at the cell center, and normalization
scales. Hydro and radiation gradients preserve their field/length dimensions.
Return a nonnegative score greater than one to request refinement. Criteria
cannot mutate fields or mesh topology. Additional problem-specific criteria can
be passed to the `Runtime` constructor using the same interface.

The shadow score is the maximum over selected field components of

```
abs(fine-shadow) /
  (shadowTolerance * (max(abs(fine),abs(shadow)) + shadowFloor*fieldScale))
```

This combines spatial detail and coarse/fine evolution error; it is not a
calibrated Richardson estimate with a guaranteed error constant.

## Independently evolved shadows

The hierarchy starts with conservative volume averages of active cells.
Coarse states used by the refinement estimator are advanced independently with the same hydro/radiation
MUSCL–Hancock kernels and synchronized physical timestep. Shadow patches exist
on leaf blocks and their ancestors, with half as many cells per axis; this also
provides a coarse comparison inside an unrefined root block. Neighboring shadow
states supply their boundaries. Where that level has no shadow coverage,
active/coarser states supply boundary data.

Gravity kicks also act on shadow hydro states. Their acceleration comes from
the restricted physical gravity field plus the configured external acceleration.
Active leaf entries are refreshed after a transport step, while covered coarse
states retain their independent evolution. Refinement checks compare against
those evolved states **before** replacing them with fresh conservative averages
for the next comparison interval. Shadows are never additional gravity masses,
physical transport cells, or contributions to output integrals.

## Refinement frequency and signal travel

The physical padding in direction d is

```
signalBuffer * maximumPredictedSpeed[d] * nextDt * regridEvery
    + bufferCells * cellWidth
```

The maximum is taken over the active domain, including the current signal speed
plus the magnitude of acceleration times the lookahead interval. Radiation uses
the configured transport speed. This covers faster propagation outside the
initial tagged cell. Tag boxes wrap around periodic boundaries. Signals from
any selected criterion can request this buffering, including cells already at
the finest level. The runtime tracks consumed travel and checks the next step
against the remaining budget; an increased speed or timestep can trigger an
earlier regrid. The step-count interval is an upper bound, not a reason to
ignore a depleted travel budget.

## Conservative transport and transactional regridding

The transport timestep is shared by all spatial levels. Fine-to-coarse ghost
values are volume averages. Coarse-to-fine ghosts use limited linear
reconstruction in conserved variables with a common admissibility limiter.
Fine boundary fluxes replace the corresponding coarse flux by an area-weighted
reflux correction before the output bank is published. The correction preserves
mass, momentum, gas energy, radiation energy, and radiation flux integrals to
roundoff for closed/periodic transport. Physical boundary fluxes and gravity
source terms retain their usual effects.

Regridding freezes one immutable pre-regrid source structure. It computes the
complete balanced destination mesh, sorts leaves by the **Morton space-filling
curve**, and assigns contiguous curve segments to storage localities. All
unchanged, split, and merged cells transfer from the original source, using
conservative limited prolongation or volume restriction. New storage, executors,
and halo plans are constructed before one final publication. Failure leaves
the previous published mesh and values available. Time is unchanged by regridding.

## Independent adaptive gravity octree

Gravity builds its own sparse cell octree from active hydro/density cell
centers and masses `rho*cellVolume`. It extends below hydro blocks, holds one
multipole/local expansion per node, and partitions numerical expansions
independently of hydro field ownership. It is rebuilt after topology changes;
subsequent solves update moments while reusing geometry and interaction lists.

A dual tree traversal handles differing cell sizes. Fixed geometric centers
and dyadic scales give integer relative offsets on a common lattice, so
translation transforms remain memoized across density updates. Rescaling
moments and locals accounts for source and target widths. Periodic Ewald
acceptance uses the next image distance, independently of the nearest-image
Newtonian calculation. Reflection uses the actual image geometry. Expansions
crossing a nearest-image branch are opened. Ewald M2Ls retain the dense
translation and background trace term.

The sampled direct reference now accepts mixed-level cells and uses each
cell's actual volume. Target selection is reproducible in physical x-fast
center order. Direct-reference error norms are cell-sampled, not volume-weighted
continuum norms; the continuum comparison remains volume-weighted.

## Current limits

There is no time subcycling or restart format. Refinement decisions, transfer
planning, and shadow storage use the coordinating locality and temporary
snapshots; they are not yet a distributed persistent mesh directory. Active
transport, flux reconciliation, and gravity run through the distributed field
store. Hydro self gravity retains the existing kinetic-energy kick and is not
an exactly energy-conserving gravity discretization.

VisIt files contain active leaves only and include `refinementLevel` when AMR
is enabled. This supports mixed cell widths without overlapping coarse output.
