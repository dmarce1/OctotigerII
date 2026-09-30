# Stationary Gresho–Chan vortex

This is the compressible Euler test described by Springel (2010), section 8.5,
following Gresho & Chan (1990) and Liska & Wendroff (2003):
https://academic.oup.com/mnras/article/401/2/791/1147356

The square domain has side length L and the vortex is centered in the box.
With dimensionless cylindrical radius r = sqrt(x*x + y*y)/L measured from that
center, density is 1, gamma is 5/3, and the azimuthal speed is

- r < 0.2: v_phi = 5 r
- 0.2 <= r < 0.4: v_phi = 2 - 5 r
- r >= 0.4: v_phi = 0

Pressure is respectively `5 + 12.5 r^2`,
`9 + 12.5 r^2 - 20 r + 4 log(r/0.2)`, and `3 + 4 log(2)`.
The code uses CGS quantities with unit density, speed, and length scales in the
supplied inputs. The radial momentum equation is dp/dR = rho v_phi^2/R.
Density, pressure, and velocity are time independent in the exact solution.
There are no gravitational, radiation, or rotating-frame source terms.

`octoII-2d` runs the planar test. `octoII-3d` extends the same state uniformly
along z, sets v_z=0, and uses periodic boundaries in all three axes. This is
an exact 3D vortex-column solution, but does not exercise intrinsically 3D
flow dynamics. The reference reports cell-center errors; these are not errors
against exact cell-volume averages. The vorticity jumps at r=0.2 and r=0.4,
so smooth-solution convergence rates should not be assumed.

## Four-node QueenBee4 diagnostic

Build the current revision on a compute node with the existing diagnostic
HPX prefix, then submit from the repository:

```sh
sbatch problem/hydro_tests/Gresho/job-qb4.sh
```

The script uses four localities, one per node, with 64 HPX threads each and the
`loni_graphshpx` allocation. It retains lock verification and deadlock
instrumentation. `hpx.trace_depth=0` avoids the independently reproduced HPX
1.11/GCC stack-unwinding crash; see `docs/queenbee4-binary-scf.md`.

`qb4.ini` uses 8^3 cells per block, spatial levels 2–3 (equivalent to 32^3 and
64^3 cell widths), and independent time levels with subcycling enabled.
`amr.refineSpeed=0.9` tags the fastest part of the vortex. This general optional
criterion tags the magnitude of grid-frame gas velocity; zero disables it.
Density tagging cannot distinguish the vortex because its density is uniform.
The shadow criterion is disabled. Regridding is checked every coarse step;
additional cell padding is zero, while the mandatory signal-travel buffer
remains enabled. The initial 3D mesh has 48 level-2 leaves and 128 level-3
leaves, totaling 90,112 cells. The first tested coarse interval takes four
fine intervals; CFL limits can require more than two.

This first run stops at t=0.02, to diagnose progress and conservation cheaply.
It is not the literature's t=3 accuracy benchmark and not a scaling result.
Per-level stages appear in task-0 stderr, and the normal table reports leaf
counts and synchronized times. Each job has separate logs, profiles, Silo
frames, `conservation.csv`, and `analytic-errors.json`. The submission script
also records its inputs, repository commit, and executable hash.

## Checks

`greshoChecks-2d` and `greshoChecks-3d` check centrifugal balance by numerical
pressure differentiation, continuity, positive states, distinct AMR levels,
fine-step counts, and conservation after transport and regridding. The 2D
suite also checks decreasing velocity error with mesh refinement at t=0.05.
These tests use the standard distributed test launcher when HPX is enabled.
A longer run and resolution study are still needed for a full accuracy claim.
