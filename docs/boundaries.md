# Physical boundary conditions

Each active coordinate has independent lower and upper faces. In 3D the options are:

```ini
[mesh.boundary]
xLower = periodic
xUpper = periodic
yLower = reflecting
yUpper = outflow
zLower = analytic
zUpper = outflow
```

Use the same dotted names on the command line, for example
`--mesh.boundary.yLower=reflecting`. A 1D build accepts only x faces; a 2D
build accepts x and y faces. Inactive-axis options and unknown boundary names
are errors. Analytic faces require a problem evaluator, as described below.

The five values are:

| Value | Ghost-cell treatment |
| --- | --- |
| `periodic` | Read the corresponding cell on the opposite side of the domain. Both faces of that axis must be periodic. |
| `reflecting` | Mirror each ghost layer across the face and negate only the normal hydro momentum or normal radiation flux. Density, energy, and tangential components retain their values. |
| `outflow` | Copy the outermost interior cell, then zero inward normal hydro momentum and inward normal radiation flux independently. Lower faces clamp positive components; upper faces clamp negative components. Other values, including total energy, are copied unchanged. |
| `inflow` | Copy the outermost interior cell without clamping momentum or radiation flux. Both inward and outward flow are permitted. |
| `analytic` | Evaluate the problem's prescribed physical state at the actual ghost-cell center and current stage time. |

Different axes and nonperiodic lower/upper faces can use different rules.
Reflections compose at edges and corners, reversing each crossed normal
component. If any crossed face is analytic, the prescribed **complete state**
at the original ghost position takes precedence over wrapping, clamping, and
reflection. A single problem evaluator defines this extension consistently;
users implementing custom problems must choose compatible data at face intersections.

`inflow` names constant extrapolation here; it does not prescribe a separate
inlet state. Use `analytic` for prescribed data. Older configurations that need
the previous unclamped `outflow` behavior should select `inflow`.

## Defaults and precedence

Existing problem defaults are retained: Kelvin–Helmholtz, radiation-pulse, and
streaming default to all periodic; Sod, collapse, and gravity tests default to
all outflow. Sod and Kelvin–Helmholtz can now override their defaults.
Rayleigh–Taylor defaults to periodic transverse faces and reflecting walls on
both ends of the last active axis (z in 3D, y in 2D).

`mesh.periodic=on` remains a shorthand for setting every active face to
periodic; `mesh.periodic=off` sets every face to outflow. For each input source,
the shorthand is applied first, then explicit face settings, independent of
line or argument order. Later INI files override earlier files, and command-line
settings override all INI files. A later shorthand resets all earlier face
settings. Pair validation runs after the complete configuration is assembled.
A single periodic face is rejected; its partner is never silently changed.

For example, a 3D Sod shock tube with reflecting x walls, periodic y faces,
and outflow z faces can use:

```sh
./release/octoII-3d --problem.name=sod --mesh.periodic=off \
  --mesh.boundary.xLower=reflecting --mesh.boundary.xUpper=reflecting \
  --mesh.boundary.yLower=periodic --mesh.boundary.yUpper=periodic
```

## Analytic problem hooks

`problemBoundary(Config const&)` returns a `ProblemBoundary` callback taking
CGS position and physical time. It supplies hydro primitives and/or the
radiation conserved state for the physics enabled by the problem. Hydro is
converted to conserved quantities using the configured equation of state.
Invalid states and missing evaluators are errors, not fallback outflow states.
Boundary evaluation is independent of `verification.analytic`; disabling error
reports does not disable prescribed boundary data.

Implemented evaluators:

* **Sod:** the exact infinite-domain Riemann solution, including its evolution
  after waves reach the simulation boundary.
* **Streaming:** the freely translating Gaussian beam, moving at the configured
  transport speed and storing the physical flux F = cE in its propagation
  direction. Only axes with paired periodic faces wrap their coordinates.

Kelvin–Helmholtz, Rayleigh–Taylor, radiation-pulse, and the gravity problems do not currently
provide analytic boundary data. They reject analytic faces during configuration.
Add a callback in the selected `bin/.../problem.cpp` to support another problem;
an analytic boundary prescription does not require a full-domain exact solution.

Callbacks are constructed independently on each locality from the serialized
configuration. Cached halo plans store geometry, reflection masks, and lower/upper outflow masks;
analytic states are evaluated anew at the current stage time, including for
stolen work. MUSCL–Hancock fills at the step's starting time and predicts face
states to the half step. The owning-patch interface also refreshes ghosts at
the end time. The application rebuilds temporary halo values at each stage.

The Sod verification reference is unavailable for periodic x faces and remains
time limited for reflecting/outflow x faces. With analytic data on both x
faces it has no boundary-crossing time limit. Streaming verification accepts
axes with paired periodic or paired analytic faces; other valid mixed-face
runs still execute but do not claim that full-domain reference is exact.

## Gravity

Gravity supports outflow, inflow, paired periodic faces, and reflecting faces.
Outflow and inflow both have an isolated gravitational exterior. The
periodic source lattice may span one, two, or three axes. Reflecting faces add
same-sign image masses; paired reflecting faces give an infinite even extension
with twice the domain period. When all axes repeat, a uniform compensating
background removes the volume-averaged density. The 1P/2P line/sheet zero modes
remain. Analytic gravity is not implemented.

See [periodic and reflecting gravity](gravity-images.md) for the image geometry,
acceptance tests, Ewald kernels, potential conventions, and verification.

## Tests

`boundaryChecks` covers face parsing, precedence, invalid pairs, inactive axes,
ghost layers, mixed faces/edges/corners, analytic precedence and physical time,
inflow copying, directional outflow clamping, reflecting conservation,
decomposition agreement, and the supported gravity boundaries.
It is registered with CTest and with the existing two-/three-locality HPX test
launcher. HPX serialization tests cover the per-face configuration and halo
metadata. Run the focused checks from a configured build directory with:

```sh
ctest --output-on-failure -R 'Boundary|GravityBoundaries|FiniteVolume'
```
