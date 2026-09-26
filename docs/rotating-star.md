# Original Octo-Tiger rotating-star benchmark

`rotatingStar` initializes the oblate, uniformly rotating n=3/2 equilibrium from
Octo-Tiger's original rotating-star test. Its axisymmetric density and internal
energy come from the original self-consistent-field (SCF) table. The model has
approximately a 2/3 polar-to-equatorial radius ratio and uses an ideal gas with
gamma=5/3. Both hydrodynamics and self-gravity are enabled; radiation is absent.

The table is compiled into the executable. There is no runtime input-data file
to find or download. [Data provenance and reproducible extraction](../bin/science/RotatingStar/DATA.md)
record the original revision, binary hash, generator, and license. The source
benchmark is described in [section 3.7 of the Octo-Tiger paper](https://doi.org/10.1093/mnras/stab937).

## Run on one locality

From the repository directory, the short corotating example is:

```bash
./release/octoII-3d --config=examples/rotating-star.ini \
  --hpx:localities=1 --hpx:threads=8 --hpx:bind=none
```

The comparison with a stationary grid uses the same physical initial state:

```bash
./release/octoII-3d --config=examples/rotating-star-inertial.ini \
  --hpx:localities=1 --hpx:threads=8 --hpx:bind=none
```

These examples end at 0.25 seconds for an initial check. A complete stellar
rotation takes 47.17408127431507 seconds at their density scale; set
`runtime.stopTime=47.17408127431507` for that duration. Increase `amr.maxLevel`
for a spatial convergence comparison. The initial resolution and duration are
not a substitute for the long-duration, higher-resolution published benchmark.

Open each output directory's `frames.visit` in VisIt to compare the density
profile, flattening, and velocity. Gas energy and momentum are inertial
quantities even in the corotating calculation. Conservation records include
boundary fluxes, which matter with the free (outflow) boundaries used here.

## Scaling and controls

The original table has G=1, dimensionless spin
`Omega*=0.5155532816213834`, and a density normalization near one. Choose the
length unit `L0=star.radius` and density unit `rho0=star.centralDensity` in CGS.
The physical units are

```
t0 = 1/sqrt(G rho0)
v0 = L0/t0
e0 = G rho0^2 L0^2
Omega_star = Omega*/t0
```

For this fixed model, `star.radius` is the **length unit**, whereas the spherical
`polytrope` problem uses that setting as its surface radius. The tabulated
equatorial and polar surfaces are approximately `0.905 L0` and `0.605 L0`.
The table's maximum stored density is `1.0000118570968204 rho0`, owing to its
finite SCF grid. `star.centralDensity` therefore specifies its original
normalization rather than imposing an additional renormalization.

With the supplied `L0=1e9 cm` and `rho0=1e6 g/cm^3`, stellar spin is
`0.13319147161856015 rad/s`. `frame.omega` independently selects grid rotation;
zero gives a stationary grid. Changing it does not change the stellar spin or
any initial inertial conserved quantity. If the density scale changes, update
the example's explicit `frame.omega` to the derived stellar spin when corotation
is desired. There is no arbitrary spin multiplier: it would invalidate this
fixed equilibrium.

`star.center.x/y/z` shift the star's center and spin axis. Omitted centers
default to the box midpoint; the example centers both star and grid rotation
on the origin. Omitted box bounds are `-2 L0` and `+2 L0`. An off-axis center
does not remain at a fixed grid coordinate on a rotating grid. The initializer
specifies its physical state at time zero.

The original initializer floors dimensionless density and internal energy
independently to `star.atmosphereFraction` (default `1e-10`). Thus the physical
floors are `f rho0` and `f e0`. Momentum is calculated from the unfloored stellar
density: the vacuum atmosphere is initially stationary in the inertial frame,
as in the original problem. The atmosphere and interpolation introduce small
departures from the continuum equilibrium. Density and energy are interpolated
independently using the original four-point bicubic construction, with exact
Lagrange coefficients replacing rounded decimal coefficients.

The original table fixes n=1.5 and gamma=5/3; incompatible overrides are rejected.
Only outflow boundaries are supported for this isolated-star problem. Startup
refinement temporarily widens an unresolved core so it can be detected; final
stored values use the physical table scale.

## Verification

The original table contains no gravitational potential or acceleration. It is
therefore not advertised as an exact gravity reference. The default
`verification.gravityReference=direct` compares the gravity solver with the
discrete mass distribution. `verification.gravityReference=continuum` reports
the unavailable reference explicitly.

The model checks verify original table entries and flattening, dimensional
scaling and inertial kinetic energy, floor behavior, independence of stellar
spin from grid rotation, INI-only selection, and startup refinement of an
off-center core. An independent cylindrical integration of the density checks
the SCF Bernoulli relation `h + Phi - Omega_star^2 R^2/2 = constant`; its
tolerance accounts for the finite table and integration resolution.

The Cartesian finite-volume method is not a discretely balanced SCF solver.
Small pulsations and surface diffusion are expected. Assess their convergence
and the angular-velocity profile alongside total mass, inertial energy, and
angular momentum. Energy conservation alone does not show that the equilibrium
has been accurately maintained, and gravity energy bookkeeping does not imply
exact conservation of total angular momentum.
