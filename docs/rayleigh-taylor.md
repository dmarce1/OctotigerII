# Rayleigh–Taylor instability

Select `--problem.name=rayleigh-taylor` in `octoII-3d` or `octoII`.
Rayleigh–Taylor requires hydro and is available only in 3D, with z vertical.

```sh
cmake -S . -B release -DCMAKE_BUILD_TYPE=Release -DOCTOTIGERII_WITH_HPX=OFF
cmake --build release -j 8
./release/octoII-3d --problem.name=rayleigh-taylor --config=bin/hydro_tests/RayleighTaylor/inputs
```

The default 3D box is `[-0.5,0.5]` cm on each axis. Both x faces and both y
faces are periodic. Both z faces are reflecting. These defaults can be
overridden through the usual per-face options. The shared input file relies
on these defaults.

The [linked Wikipedia article](https://en.wikipedia.org/wiki/Rayleigh%E2%80%93Taylor_instability)
uses vertically unbounded fluids with perturbations vanishing at infinity;
it does not specify a finite upper numerical boundary. Reflecting top and
bottom walls are our finite-box choice, also used by
[Stone & Gardiner (2007), section III](https://arxiv.org/abs/0707.1022).
This setup uses a cubic domain and a deterministic single-mode seed, rather
than reproducing that paper's aspect ratio and random perturbations.

## Initial state and controls

All values are CGS. Let `L=mesh.upper-mesh.lower`, `z0` be the domain midpoint,
and `z` denote the vertical coordinate. Defaults are:

| Input | Default | Meaning |
| --- | --- | --- |
| `rayleighTaylor.densityLower` | 1 g/cm³ | Density below the interface |
| `rayleighTaylor.densityUpper` | 2 g/cm³ | Density above the interface |
| `rayleighTaylor.interfacePressure` | 2.5 dyn/cm² | Continuous pressure at `z0` |
| `rayleighTaylor.perturbation` | 0.01 cm/s | Vertical velocity seed amplitude |
| `hydro.acceleration.z` | −0.1 cm/s² | Uniform external gravity |
| `hydro.gamma` | 1.4 | Ideal-gas adiabatic index |
| `runtime.stopTime` | 10 s | End time |

Each layer begins with hydrostatic pressure
`P(z)=P0+rho*a_z*(z-z0)`. The vertical seed is
`v_z=A*cos(pi*(z-z0)/L)^2*product_d cos(2*pi*(x_d-z0)/L)`, with the product
over transverse axes. Transverse velocities vanish. The seed is periodic
transversely, has zero transverse mean, and vanishes at both vertical walls.
The interface lies on cell faces for every supported even resolution.

The upper density must exceed the positive lower density. Gravity must be
vertical and nonpositive, and the hydrostatic pressure must stay positive
through the upper wall. Set the perturbation or gravity to zero for control
runs. Runtime inputs reject nonfinite values and incompatible profiles.

The build enables hydrodynamics and prescribed gravity; it does not enable
self gravity. Uniform acceleration is supported by 3D hydro builds through
`hydro.acceleration.x/y/z`. The runtime applies symmetric
half kicks around each transport step, updating momentum and the kinetic
part of total energy while preserving internal energy during each kick.
The timestep includes an acceleration limit. If self gravity is compiled,
the external acceleration is added to its field. The source update follows
the same distributed ownership and work-stealing rules as transport.

This is the existing second-order finite-volume scheme with split gravity,
not an exactly well-balanced hydrostatic discretization. Reflecting ghosts
mirror pressure as well as density; finite-resolution wall and interface
transients therefore remain. Use resolution studies for quantitative growth
rates. The nonlinear compressible problem has no implemented exact solution;
`verification.analytic=on` correctly rejects it.

## Verification

`rayleighTaylorChecks` checks the default walls, parameter validation,
hydrostatic initialization, momentum/energy kicks, acceleration timestep,
mass conservation, agreement across block decompositions, and growth of
the seeded transverse mode against a zero-gravity control. HPX builds also
run source and decomposition checks on two and three localities.

```sh
ctest --test-dir build-rt --output-on-failure -R 'RayleighTaylor|Boundary|FiniteVolume'
```
