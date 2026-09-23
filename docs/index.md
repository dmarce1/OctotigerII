# OctotigerII developer reference {#mainpage}

OctotigerII separates numerical fields, mesh geometry, and distributed execution.
The current application supports adaptive Cartesian meshes, Euler hydro,
uncoupled M1 radiation transport, and Newtonian gravity with open, periodic,
and reflecting boundaries.

Start with [problem builds](../BUILDING.md), [storage and execution](../STORAGE.md),
[the numerical conventions](numerics.md), [parallel gravity](parallel-fmm.md),
[profiling](profiling.md), or the
[bibliography](../BIBLIOGRAPHY.md). The Modules, Classes, and Files pages provide
API documentation and links to the source. Search accepts names such as
`StoragePartition`, `RadiationSystem`, or `Runtime`.

## Build and view

From HOME, assuming the project is at `~/workspace/OctotigerII`:

```bash
sudo apt install doxygen
~/workspace/OctotigerII/docs.sh --open
```

Doxygen 1.9.8 or newer and Python 3 are required. Graphviz and LaTeX are optional
and are not required by this configuration. No application build, HPX installation,
or network service is needed to generate or read the HTML. Open
`~/workspace/OctotigerII/docs/generated/html/index.html` in a browser later.
The generated site includes its search assets and works from local files.

For an existing CMake build:

```bash
cmake --build ~/workspace/OctotigerII/release/sod/1d --target docs
```

That target writes `release/sod/1d/docs/html/index.html`. In Eclipse, create an External
Tools entry for `${workspace_loc:/OctotigerII/docs.sh}` with argument `--open`,
or open the generated `index.html` directly in a browser.

## Architecture

| Layer | Owns | Does not interpret |
| --- | --- | --- |
| `StoragePartition` | Named typed columns, banks, allocation lifetimes | Mesh cells, tree relationships, particles |
| `FieldHandle<T>` | Field identity and partition/range directory | Halo geometry or update order |
| `CartesianTopology` | Blocks, coordinates, field ranges, halo plans | Numerical field values |
| `refinement` | Read-only criterion interface, mass and shadow criteria | Topology changes and storage placement |
| `amr::Hierarchy` | Conservative transfers, evolved coarse shadows | Physical leaf ownership |
| `LocalExecutor` | Ready-work queue and bounded reusable workspaces | Field allocation policy |
| `Runtime` | Stage synchronization and publication | Physics-specific flux formulas |

A local checkout references the stored interior. A remote checkout receives a
typed buffer. The temporary `PatchView` joins those values to compact halos.
All persistent fields, including radiation flux, contain physical CGS quantities.

A transport stage reads bank A, writes bank B, waits for all workers and writebacks,
then publishes B. A failed stage drains work and preserves A. See @ref storage
and @ref runtime for the contracts that make this safe.

## Scope and validation

The storage component accepts arbitrary range lengths and several layouts in
the same partitions. The adaptive mesh supplies conservative transfers, evolved
shadows, coarse/fine refluxing, and Morton-ordered placement. See [AMR](amr.md).
Particle integration and restart checkpoints remain future work. The FMM partitions its numerical expansions across localities,
reads density from distributed field ranges, and publishes gravity through the
same field store. Diagnostics and output still gather on the coordinator.

See [validation](../VALIDATION.md) for reproducible tests and their limits.
See [GoogleTest testing](testing.md) for unit/integration suites and the matrix runner.
See [documenting code](documenting.md) for the citation and formatting conventions.

- [Per-face physical boundaries](boundaries.md)
