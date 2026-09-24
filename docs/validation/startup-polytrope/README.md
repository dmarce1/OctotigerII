# Startup AMR and Lane–Emden star validation

Validated with GCC 13.3, C++20, Release optimization, Boost 1.83, GoogleTest
1.14, and `OCTOTIGERII_WITH_HPX=OFF` in three dimensions.

| Build/check | Result |
| --- | --- |
| Collapse full CTest suite | 135/135 passed |
| Polytrope full CTest suite | 138/138 passed |
| All-physics standalone AMR suite (hydro, gravity, radiation) | 13/13 passed |
| Final polytrope-specific tests, including missed coarse peaks | 3/3 passed |
| Collapse initial Silo density compared at every cell center | Maximum absolute error 1.819e-12 g/cm³ |

The full logs accompany this file. The polytrope-specific tests were rerun
after adding assertions that the physical, unbroadened root grid entirely
misses the displaced star for both tested indices, 1.5 and 3. Startup still
reaches the finest level, captures the core, and leaves every stored density
equal to the requested physical profile.

## Collapse output regression

The supplied Collapse input was run with `runtime.stopTime=0` and analytic
verification disabled. An independent Silo reader checked each output cell
against

```
rho(x,y,z) = 1e4 * (0.01 + exp(-(x*x+y*y+z*z)/(2*(3e8)^2)))
```

The frame contains 120 blocks and 7,680 cells. Its maximum density is
9,938.55697400955 g/cm³. This is the Gaussian sampled on the fine mesh; the
old coarse peak is no longer prolonged into a flat central plateau.
The run and check output are in `collapse-initial.log` and
`collapse-initial-check.log`; the standalone reader is included as
`check-collapse-silo.cpp`.

## What the tests establish

- Startup visits the root and each newly refined level before the next
  criterion check, with freshly evaluated initial conditions.
- Final snapshots match initialization at each leaf's own coordinates and
  width, including the combined hydro/radiation build.
- Later refine/coarsen operations conserve the evolved integrals.
- Density tagging is independent of cell volume, unlike mass tagging.
- Lane–Emden integration matches the n=1 closed form and n=1.5 surface
  constants, radius scalings, and hydrostatic pressure-gradient identity.
- INI radius and center settings, CLI overrides, radial symmetry, and zero
  central acceleration are checked.
- A short evolution conserves stellar mass within 1e-7 relative error and
  keeps core velocities below one percent of the central sound speed.

## Limits

An HPX installation was unavailable in the validation environment, so HPX
compilation, serialization, and distributed execution were not run. HPX paths
use HPX mutexes, including the gravity and Lane–Emden caches and storage
metadata. Gravity/Subgrid tests now start HPX when that backend is enabled;
test discovery remains independent of runtime startup.

The star is a continuum equilibrium with a tenuous atmosphere. The short
evolution check does not establish long-term numerical equilibrium or
stability. The discretization is not an exactly balanced hydrostatic scheme;
see [the scenario documentation](../../polytrope.md).
