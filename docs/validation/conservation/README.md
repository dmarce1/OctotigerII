# Header layout, CLI selection, and conservation validation

Validated with GCC 13.3.0, Release, HPX disabled, all physics modules enabled.
Built `octoII-1d`, `octoII-2d`, `octoII-3d` and the `octoII` symlink, plus
options, AMR and conservation GoogleTest executables for each dimension.

142 selected checks passed: 139 GoogleTest cases and three application suites.
The application suite exercises every available problem per dimension and
checks unavailable dimensions/modules, missing problem selection, config-file
precedence, Silo output and conservation-file columns. The first application
run exposed an overly specific test expectation for an empty option value;
Boost rejects that syntax before semantic validation. The test was corrected
and all three application suites passed on rerun.

Numerical coverage includes:

- Deep 2:1 refinement across faces, edges, corners and periodic neighbors.
- Conservative refinement/coarsening and mixed-level periodic refluxing.
- Hydro and radiation boundary-flux budgets on uniform and mixed-level meshes.
- Analytic radiation inflow/outflow and reflecting-wall momentum traction.
- Cumulative transport surviving regridding and rejected negative timesteps.
- A nonzero L1 momentum norm when signed momentum cancels to zero.
- Zero transport through periodic domain boundaries.
- Kinetic + thermal = gas energy, and self-gravitational energy with 1/2 rho phi.
- Per-step diagnostic output with Silo disabled.

A separate streaming run checked CSV numeric values for finiteness, zero initial
normalized drift, and normalized boundary-corrected drift below 1e-12.

HPX compilation and distributed execution were not available in this environment.
Boundary reductions follow the existing serialized PhaseResult and locality
reduction paths; no new standard threading or synchronization was introduced.

See ../../conservation.md for units, norms and the scope of gravity diagnostics.
