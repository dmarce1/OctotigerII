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
physical radiation-momentum integral. Radiation-matter coupling with a reduced
light speed instead conserves a weighted combined momentum, described below.

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
acceleration, or radiation-matter exchange) can change boundary-corrected gas
energy and momentum. Radiation energy and flux also exchange with the gas;
their separate drift columns remain useful diagnostics, but are not invariants
when coupling is active.

## Gas, radiation, and gravity together

When both hydro and radiation are enabled, new columns are appended after the
existing columns. Write `w = c/chat = 1/radiation.lightSpeedRatio`. The on-grid
combined quantities are

    physical_total_energy = gas_energy + potential_energy + radiation_energy
    rsla_total_energy = gas_energy + potential_energy + w * radiation_energy
    physical_total_momentum = gas_momentum + radiation_flux_integral / c^2
    rsla_total_momentum = gas_momentum + w * radiation_flux_integral / c^2

Potential energy is zero when self gravity is disabled. The `physical_` values
are the physical gas-plus-radiation quantities; the `rsla_` values are the
invariants of the implemented reduced-speed exchange equations. They coincide
when `chat=c`. With `chat<c`, conservation of the weighted quantities does not
mean that physical energy and momentum are conserved by that approximation.

Each combined energy/momentum has the same `grid`, `in`, `out`, `corrected`,
`l1`, `norm`, and `drift_scaled` suffixes as the separate ledgers. Boundary
contributions use the same sums and radiation weights as the corresponding
grid quantity. Their `in` and `out` are sums of the existing signed-component
ledgers, rather than a fresh sign split of the combined face flux. This leaves
`corrected = grid + out - in` exact while providing a conservative magnitude
scale even when components oppose one another. The `l1` scale similarly sums
the absolute gas, radiation, and binding contributions before cancellation.
The binding contribution retains its existing factor of one half on the grid
and its full potential-transport weight at the boundary.

The existing `angular_momentum_z_g_cm2_s_grid` remains gas-only. Additional
grid columns report

    radiation_angular_momentum_z = integral [x_inertial cross F/c^2]_z dV
    physical_total_angular_momentum_z = gas_angular_momentum_z + radiation_angular_momentum_z
    rsla_total_angular_momentum_z = gas_angular_momentum_z + w * radiation_angular_momentum_z

These are inertial moments about the origin, including on a rotating mesh.
They exclude angular momentum transported through the boundary. No angular
momentum boundary ledger or claim of exact global angular-momentum conservation
is implied by these grid integrals. Local paired radiation exchange conserves
the weighted angular moment at the cell position; spatial transport and gravity
must be assessed separately.

The maximum cell optical depth, radiation trapping parameter, and RSLA
criterion are also appended. Their definitions, Silo fields, and limitations
are given in [radiation-coupling.md](radiation-coupling.md).

## Hydro plus gravity

Additional columns report kinetic and thermal energy, and when self gravity is
active:

    potential_energy = sum(0.5 * rho * phi * V)
    gas_gravity_energy = gas_energy + potential_energy

With dual energy enabled, `thermal_energy` integrates the internal energy chosen
by the pressure threshold. It can differ from `gas_energy - kinetic_energy` in
kinetically dominated cells. The independently conserved `gas_energy` is used
for conservation and gas-gravity accounting. The entropy auxiliary is excluded
from this ledger because synchronization is a source for it; see
[dual-energy.md](dual-energy.md).

The gas-gravity normalization is the larger of the initial and current values
of `sum((abs(Egas) + abs(0.5*rho*phi)) * V)`. This avoids cancellation between
positive gas energy and negative binding energy.

The potential-energy `in` and `out` columns record the signed split of
`dt*A*F_out*phibar_face` on physical boundaries. Combined energy has `grid`,
`corrected`, `norm`, and `drift_scaled` columns; the corrected value adds both
gas and potential outward transport and subtracts inward transport. It is the
energy invariant for a fixed mesh and reciprocal self-potential operator when
there is no radiation-matter exchange. With coupling, use the combined RSLA
energy ledger above.
Imposed acceleration, a nonreciprocal gravity operator, and changes of mesh
with `gravity.conserveRegridEnergy=off` can cause drift. See [gravity-energy.md](gravity-energy.md) for the exact
work formula and the distinction between this coupling and momentum kicks.
The factor 1/2 applies only to the on-grid self-potential integral, not the
potential advected through a physical boundary.

## Existing AMR safeguards

Mesh selection already enforces a maximum one-level jump across faces, edges,
corners and periodic neighbors after refinement and coarsening. Hydro and
radiation already apply a separate reflux phase: for each coarse face adjacent
to fine leaves, the coarse-cell update is corrected with the area average of
fine-face fluxes. With a synchronized timestep, this is algebraically
equivalent to replacing that coarse flux before taking its divergence. With
time refinement, accepted fine fluxes are accumulated over the coarse interval
before reflux. Gas and radiation have separate registers, and their combined
budget uses both accepted transport contributions.

## Current block-size constraints

`mesh.cells` remains a power of two from 4 through 128. The shadow solver evolves
a half-sized patch, and `MeshLayout` requires an even count of at least two.
The combined cell/block hierarchy uses binary cell locations with a level offset
of `log2(mesh.cells)`. The gravity tree also repeatedly halves the grid. These
implementation assumptions require changes before admitting 2-cell blocks or
non-power-of-two block sizes; the finite-volume equations do not require them.

## Separating energy error sources

The cumulative columns `gravity_reciprocity_defect_erg` and
`gravity_regrid_energy_change_erg` measure the endpoint potential-operator defect
and actual total-energy changes during remapping, respectively.
`gravity_energy_budget_residual_scaled` is the scaled combined-energy drift
minus those two contributions. All use actual evolved fields, not a fit to the
observed drift. The original physical-drift column remains unchanged in meaning;
no measured error is removed from the cell solution or that column.

These decompositions are filled by the normal `run()` driver. A custom caller
that constructs `Diagnostics` directly must supply cumulative step and regrid
measurements itself; `diagnose()` alone cannot infer temporal history.
Compensated summation is used for global grid integrals and norms, reducing
order-dependent summation noise. Boundary transport retains its existing ledger.
