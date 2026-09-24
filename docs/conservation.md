# Conservation diagnostics

`output.directory/conservation.csv` is truncated at the start of an invocation
and flushed at initialization and every completed timestep. `output.every` and
`output.enabled` control Silo output only. Conserved sums are not printed to the
terminal. The CSV uses CGS and 17 digits after the decimal
point to retain small conservation errors. In 1D/2D, the existing unit transverse
area/thickness convention is retained.

For each evolved conserved component (mass, active momentum components, gas
total energy, radiation energy, and active radiation-flux components), columns
provide:

- `grid`: volume integral over active leaf interiors; ghosts and covered coarse
  cells are excluded.
- `in`, `out`: cumulative negative/positive parts of outward-oriented numerical
  transport, split independently per component. Both columns are nonnegative.
- `corrected`: `grid + out - in`.
- `l1`: current sum of absolute cell contributions, `sum(abs(U_i) * V_i)`.
- `norm`: `max(initial_l1, current_l1, in + out)`.
- `drift_scaled`: `(corrected - initial_grid) / norm`, or zero when the norm is
  zero (all contributing values are then zero).

The L1 scale avoids dividing by a cancelling signed sum such as net momentum.
Including transported magnitude handles components initially zero. The supplied
absolute columns also permit a different normalization in postprocessing.
Radiation-flux integrals have units erg cm/s; divide by physical `c^2` for the
usual radiation-momentum integral (the reduced-light-speed equation can change
its relationship to exchange/source conservation).

Boundary transport uses the actual time-centered MUSCL-Hancock numerical flux,
multiplied by face area and the accepted timestep. Only leaf faces on physical
nonperiodic boundaries contribute; periodic and internal faces contribute zero.
AMR boundary faces use their own cell size, and regridding preserves accumulated
transport. Failed transport/reflux phases do not publish boundary sums. Worker
and locality reductions include stolen tasks exactly once.

For signed vector components, `in` and `out` are the negative/positive parts of
the outward flux of that component; they are not classifications by the sign of
mass flow. Momentum transport includes pressure traction at reflecting walls.
Analytic boundaries can supply inflow. Source terms (self gravity or imposed
acceleration) can change boundary-corrected gas energy and momentum.

## Hydro plus gravity

Additional columns report kinetic and thermal energy, and when self gravity is
active:

    potential_energy = sum(0.5 * rho * phi * V)
    gas_gravity_energy = kinetic_energy + thermal_energy + potential_energy

The gas-gravity normalization is the larger of the initial and current values
of `sum((abs(Egas) + abs(0.5*rho*phi)) * V)`. This avoids cancellation between
positive gas energy and negative binding energy. Its normalized drift measures
change of the on-grid gas-plus-gravity energy.

The gas-gravity drift is deliberately an **on-grid** diagnostic. It is not an
open-boundary corrected invariant: gas energy flux alone does not account for
nonlocal gravitational energy transport or work by imposed acceleration. No
boundary flux of gravitational energy is implemented here, and no exact energy
conservation is imposed on the existing split gravity kicks or regridding. For
an isolated closed system this diagnostic measures their energy error. The
factor 1/2 applies to the self-gravitational potential stored by the solver.

## Existing AMR safeguards

Mesh selection already enforces a maximum one-level jump across faces, edges,
corners and periodic neighbors after refinement and coarsening. Hydro and
radiation already apply a separate reflux phase: for each coarse face adjacent
to fine leaves, the coarse-cell update is corrected with the area average of
fine-face fluxes. With the current synchronized timestep, this is algebraically
equivalent to replacing that coarse flux before taking its divergence. No
subcycling is used; a future subcycled integrator must also time-integrate the
fine flux register.

## Current block-size constraints

`mesh.cells` remains a power of two from 4 through 128. The shadow solver evolves
a half-sized patch, and `MeshLayout` requires an even count of at least two.
The combined cell/block hierarchy uses binary cell locations with a level offset
of `log2(mesh.cells)`. The gravity tree also repeatedly halves the grid. These
implementation assumptions require changes before admitting 2-cell blocks or
non-power-of-two block sizes; the finite-volume equations do not require them.
