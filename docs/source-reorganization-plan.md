# Source reorganization

The approved layout is implemented. The [current inventory](source-inventory.md)
lists every maintained source and input file; the
[highlighted change map](source-inventory-proposed.md) identifies the moved and
retired responsibilities. Validation is recorded in
[the reorganization checks](validation/source-reorganization.md).

## Goal and boundary

Give each file one understandable job while preserving the numerical method,
public `Runtime` API, field layout, HPX action names, timestep order, and output
formats. A source snapshot was taken before the moves, including the existing
uncommitted rotating-grid and gravity work. The separate user edit to
`examples/polytrope.ini` is preserved byte for byte.

## Sequence

1. **Resolve the duplicate packet definition.** Removed the unused
   `octotigerII/subgrid/exchange.hpp` after checking maintained includes and
   installation rules. `octotigerII/subgrid/fluxPacket.hpp` is the single
   `FieldFluxPacket` definition, exercised by serialization tests. Runtime
   halo plans are unchanged.

2. **Separate AMR state from AMR decisions.** Kept the covered coarse states,
   reconstruction, and conservative transfers in `amr/hierarchy.hpp` and
   `src/amr/hierarchy.cpp`. Moved `RegridResult`, `selectMesh`, and their
   geometric selection helpers to `amr/regridSelection.hpp` and
   `src/amr/regridSelection.cpp`. Updated the three direct callers in startup
   initialization, runtime regridding, and AMR tests. The existing 2:1 balance
   and periodic-neighbor rules are unchanged.

3. **Split the runtime without changing its public face.** Added a private
   `src/runtime/internal.hpp` for shared stage types, implementation state,
   templates, and HPX action declarations. Kept
   `octotigerII/runtime.hpp` and the public class names intact. Moved the
   `LocalExecutor` methods and all of its HPX component/action registration to
   `src/runtime/localExecutor.cpp`, with existing action identifiers unchanged.
   Stage dispatch, task completion, and error draining moved to
   `src/runtime/stages.cpp`, AMR transport subcycling to
   `src/runtime/transport.cpp`, gravity rungs and endpoint energy work to
   `src/runtime/gravity.cpp`, and mesh replacement/regrid recovery to
   `src/runtime/regrid.cpp`. `src/runtime.cpp` remains the short lifecycle and
   public API facade. Successful field-bank publication remains in each
   transport/gravity/regrid caller, after its work completes.
   `cmake/dimension/CMakeLists.txt` adds the new implementation files to the
   same Application and Amr targets, with the same HPX definitions and linkage.

The three gravity solvers retain their separate roles: `solver.cpp` is the
standalone serial reference, `fieldSolver.cpp` handles the uniform distributed
mesh and selects the implementation, and `adaptiveFieldSolver.cpp` handles
adaptive distributed meshes. Combining them would mix distinct data layouts
and complicate the reorganization. The one-line problem manifests also remain
local to their problems; centralizing them would make adding a problem touch a
shared registry by hand again.

## Checks at each step

- Build the 1D, 2D, and 3D serial targets and the HPX 3D target. The runtime
  split must compile with profiling both enabled and disabled, because both
  configurations need the HPX compile guards and dependency.
- After the packet change, run serialization and application checks. After the
  AMR split, run AMR, time-refinement, gravity-energy, and regrid tests.
- After each runtime move, run the matching transport or gravity suite; at the
  end run the full configured CTest suites plus the existing two-locality
  gravity tests. Compare a fixed rotating-star run at `frame.omega=0` and at
  its corotating rate with the saved baseline CSVs. Conservation, timestep
  counts, and hierarchy changes must remain within the existing tolerances.
- Keep HPX component registration in exactly one translation unit. Ensure all
  dispatched work is drained before a failed stage returns and that field-bank
  publication remains the final successful step.

The reorganization adds eight files, retires one, and changes the contents
of several existing files. That is 161 program source files in the completed
layout, plus the same two configuration files and 24 runnable inputs. File
counts are a navigation aid; the goal is clearer ownership of behavior.
