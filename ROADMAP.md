# Problem-build roadmap

Requested sequence:

1. **Implemented:** Castro-style problem directories. Each executable contains
   one problem implementation and is identified by problem and dimension.
   Sod: 1D/2D/3D; Kelvin–Helmholtz: 2D/3D; streaming and radiation pulse:
   1D/2D/3D; gravity sphere, gravity Gaussian and collapse: 3D only.
   Keep resolution, refinement level, stopping time and other run settings in inputs.
2. **Implemented with step 1:** CMake selects compile-time `ndim`.
   Coordinates, vectors, state tuples, traversal, halos and workspaces use exactly
   that dimension. Remove runtime dimension dispatch and padded inactive axes.
   Physical coefficients and the intrinsically 3D gravity algorithm retain their
   mathematical constants. Lower-dimensional CGS totals retain the documented
   unit transverse measure; it does not create inactive mesh axes or cells.
3. **Implemented:** GoogleTest unit tests, fixtures, parameterized cases and
   typed finite-volume tests, individually discovered by CTest. Existing numerical,
   storage, serialization, Silo and distributed regressions are retained as named
   tests. See [testing](docs/testing.md) for coverage and matrix commands.
4. **Deferred by request:** plain-text interactive builder. Prompt for problem,
   supported dimension, build mode/backend, block cells, refinement level and run
   controls. Use the same problem manifests and CMake options as command-line
   builds; save runtime choices to an inputs file. Changing level must not rebuild
   the executable. No prompt-driven builder is implemented in this step.
