# Spherical polytropic star

`polytrope` is a 3D hydro + self-gravity problem. Radiation is absent. The
continuum stellar interior is in hydrostatic equilibrium, with zero velocity.
A tabulated Lane–Emden solution supplies the radial profile; generic indices
such as 1.5 require numerical ODE integration rather than an elementary formula.

## Build and run

From HOME:

```bash
cd ~/workspace/OctotigerII
./build.sh release  -j 12
./release/octoII-3d --problem.name=polytrope \
  --problem.name=polytrope --config=bin/science/Polytrope/inputs \
  --hpx:threads=12 --hpx:bind=none
```

Open `output/polytrope/frames.visit` in VisIt. Use `density` for the structure,
`refinementLevel` for AMR, and the `velocity` vector expression for departures
from equilibrium. `runtime.stopTime=0` writes only the initialized star.

## Parameters

All dimensional values are CGS. CLI settings override INI values.

| Option | Default | Meaning |
| --- | ---: | --- |
| `star.radius` | `1e9` | First-zero surface radius, cm |
| `star.centralDensity` | `1e6` | Central density, g/cm³ |
| `star.polytropicIndex` | `1.5` | Lane–Emden index, 0 < n < 5 |
| `star.center.x`, `.y`, `.z` | box midpoint | Independent center coordinates, cm |
| `star.atmosphereFraction` | `1e-8` | Minimum density divided by central density |
| `hydro.gamma` | `1+1/n` | Evolution adiabatic index; independently overridable |
| `mesh.lower`, `mesh.upper` | `-2R`, `+2R` | Cube bounds; explicit settings take precedence |
| `amr.refineDensity` | `0.01*rho_c` | Density above which cells request refinement, g/cm³ |

The supplied input enables AMR, sets base/minimum level 1 and maximum level 3,
and disables mass/shadow tagging. It runs for 0.1 s and samples 32 cells for the
direct gravity comparison. For example, modify the INI with:

```ini
star.radius=8e8
star.center.x=2e8
star.center.y=-1e8
star.center.z=0
```

The example input explicitly sets all three center coordinates to zero; delete
those entries to use the midpoint defaults. Moving the center does not translate
the box. Ensure the entire star remains inside the chosen bounds if an isolated,
complete star is intended. Increasing `amr.maxLevel` improves spatial resolution.

## Model

The regular Lane–Emden initial conditions are theta(0)=1 and theta'(0)=0:

```
dtheta/dxi = -mu/xi²
 dmu/dxi   = xi² theta^n
```

The integrator starts with the regular central series, uses fourth-order
Runge–Kutta steps, brackets the first zero xi1, and interpolates with cubic
Hermite polynomials. Tables are cached by n. With a=R/xi1:

```
rho = rho_c theta^n
Pc  = 4 pi G a² rho_c² / (n+1)
P   = Pc theta^(n+1)
M(r)= 4 pi a³ rho_c mu
Phi = -G M(R)/R - (n+1) Pc/rho_c theta       inside
Phi = -G M(R)/r                             outside
```

Thus radius and central density determine the pressure normalization and total
mass. The reference acceleration is -G M(r) (x-center)/r³, with zero at the center.
These scalings follow [Princeton's polytrope notes](https://www.astro.princeton.edu/~gk/A403/polytrop.pdf).

During startup tagging, the core scale a=R/xi1 is widened to at least one cell
width while retaining the central density. Using the core scale also handles
centrally concentrated indices whose cores are much smaller than their surface
radii. This exposes an off-center star to the default density criterion. Final
leaf values use the requested physical radius. An AMR maximum level with cell
width greater than a is rejected. Higher density thresholds may require finer
initial sampling. The final profile is sampled at cell centers, consistently
with the other initializers and analytic comparisons.

## Equilibrium and verification

The hydro solver requires positive density and pressure. The atmosphere uses
rho_floor=f*rho_c and P_floor=Pc*f^(1+1/n). It has a small mass and is not an exact
Lane–Emden equilibrium. The Cartesian finite-volume scheme and discrete gravity
are also not a discretely balanced equilibrium method. Expect residual velocities;
measure their convergence with resolution rather than assuming exact stationarity.

Direct gravity comparisons remain the default. Set
`verification.gravityReference=continuum` for density, pressure, velocity, and
gravity comparisons with the equilibrium target and corresponding Silo fields.
The continuum gravity target omits the atmosphere's gravity. For unstable choices
of n/gamma the evolved model can legitimately depart from the equilibrium.

Checks cover the n=1 closed form, n=1.5 surface constants, radius scaling,
dP/dr=-rho G M/r², INI radius/center controls, detection of a narrow displaced star,
and a short hydro/gravity evolution. The standalone all-physics test build also
runs the Lane–Emden unit checks; problem-specific equilibrium tests run in the
polytrope build.

See the [validation record](validation/startup-polytrope/README.md) for the
serial regression results and the limits of those checks.
