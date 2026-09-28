# Helmholtz hydrodynamic closure

Select `hydro.eos=helmholtz`. The default remains `ideal`, and the existing
`white-dwarf` closure is retained. The implementation calls our C++ translation
of F. X. Timmes's Helmholtz evaluator with his unchanged, checksum-verified
precomputed table. See [source credit](../lib/helmholtz/NOTICE.md).
There is no Fortran runtime and no table generation during a simulation.

A runnable example is [helmholtz-advection.ini](../examples/helmholtz-advection.ini):

```sh
/path/to/octoII-1d --config=examples/helmholtz-advection.ini
/path/to/octoII-1d --config=examples/helmholtz-advection.ini \
  --radiation.enabled=on --radiation.opacity=1
```

The fixture advects a smooth helium/iron composition profile at initially
uniform density, pressure, and velocity. It is a conservation and integration
fixture, not an exact nonlinear EOS contact or stellar equilibrium solution.
Its two material species can be changed with the existing mass-fraction options.
Existing problem initializers continue to prescribe their own CGS primitives;
turning on Helmholtz does not make an ideal-gas stellar profile an equilibrium
of the new EOS. Existing ideal-gas analytic comparisons are disabled for this
closure rather than presented as Helmholtz reference solutions.

## Composition and photons

Helmholtz requires `massFractions.enabled=on` and material species. At every
material field read, density and two composition moments are derived from the
actual partial densities, excluding massless tracers:

```
rho = sum rho_i
nuclei = sum rho_i/A_i          electrons = sum rho_i Z_i/A_i
Abar = rho/nuclei              Zbar = electrons/nuclei
```

Mixture species use the existing mass-weighted reciprocal-A and Z/A definitions.
These moments are included in temporary hydro states for reconstruction,
prediction, halo exchange, radiation source solves, and AMR transfer. Their
face transport uses the same donor-cell concentration rule as species transport.
Published material fields derive the moments from species, so they cannot drift
as independent conserved composition estimates. Species contacts remain first
order under the existing composition transport scheme.

With separately evolved radiation, the gas uses `pgas`, `egas`, `sgas`, their
fixed-composition derivatives, and `gam1_gas`. Photons enter exclusively through
the radiation equations. Without separate radiation, the equilibrium photon
contribution is included using the total EOS outputs. The local radiation
exchange API rejects a Helmholtz closure that includes photons.

The Newtonian acoustic speed is `sqrt(Gamma1 * P/rho)` for the selected
components. Timmes's relativistic `cs` is not used by the Newtonian hydro solver.
The radiation source evaluates `T(rho,e,Abar,Zbar)` for each nonlinear iterate;
its emission derivative uses `dT/du = 1/(rho*cv)`. It does not reuse an ideal-gas
`T proportional to u` relation. Paired gas/photon changes still conserve
`Egas + Er/(chat/c)` and the corresponding momentum invariant before floor
corrections.

## Inversion and dual energy

Temperature inversions for energy, pressure, and entropy use bracketed
Newton/bisection on the supplied table's temperature interval. Inputs above the
maximum temperature or outside the density/composition domain fail explicitly.
The density coordinate is `rho*Zbar/Abar`, with limits 1e-12 to 1e15 g/cm^3.
The upper temperature bound is 1e13 K. No density extrapolation is introduced.

For Helmholtz, the auxiliary is `A = rho*s/s0`, where `s0=1e8 erg/(g K)`.
Entropy may be signed; its normalization only supplies density units and a
reasonable numerical magnitude. Smooth adiabatic flow transports this auxiliary
conservatively. Reliable total energy resets it at shocks. Radiation heating
updates it using the selected EOS and the actual energy exchange. The existing
constant-gamma entropy power is retained only for the older closures;
`hydro.dualEnergy.exponent` must be 1 with Helmholtz.

The reliability numerator is internal energy **above the energy at the
configured temperature floor**, so a large degenerate ground-state energy does
not falsely establish that the small thermal contribution is resolved. The
existing lower pressure-selection and upper synchronization thresholds apply
relative to conserved total energy.

## Temperature floor and accounting

`hydro.temperatureFloor` defaults to 1000 K and must lie in [1000,1e13) K.
At fixed density and composition, the selected EOS determines `e_floor`.
Thermodynamic trial evaluations below that energy use the floor temperature,
but the standalone EOS evaluator continues to reject out-of-table temperatures.

Accepted updates first try the independently advected entropy when total-minus-
kinetic energy is unreliable. If both estimates are below the supported floor,
the accepted-state pass raises stored energy to `K+rho*e_floor` and resets the
entropy auxiliary. It records the actual representable energy increase.
Source iterations may temporarily cool below the floor while emitting at the
floor temperature; the accepted-state correction accounts for that artificial
heat supply. This is a numerical thermostat, not additional physical heating.

`Runtime::eosFloorEnergy()` and `Runtime::eosFloorCells()` report cumulative
energy added and the number of cell correction events. CSV columns are
`eos_floor_energy_erg` and `eos_floor_cell_corrections`. Gas and combined
radiation energy budgets subtract this explicit source in corrected columns;
raw on-grid energies remain available. Predictor, covered shadow cells, and
source subdivision attempts do not contribute. The complete coupled-step
rollback restores the ledger. Regridding preserves cumulative totals. The
owning patch solver returns its floor ledger in `StepResult`.

The floor pass runs at the accepted global synchronization boundary. Intermediate
AMR substeps use bounded thermodynamics until that pass; the counts are accepted
cell corrections, not the number of all temporary visits below the floor.

## Validation scope

`helmholtzHydroChecks` exercises EOS round trips across gas/pair/degenerate
states, photon exclusion, the constant-entropy pressure derivative against
sound speed, entropy fallback, explicit floor energy changes, composition and
tracer handling, radiation LTE and relaxation, and coupled species transport
through AMR/regridding. See the validation record for measured runs.

The supplied table is not uniformly reliable at arbitrarily cold, dense states:
the electron thermal derivative issue recorded in
[the library validation](../lib/helmholtz/VALIDATION.md) remains. The hydro adapter
rejects states with nonpositive/nonfinite selected total heat capacity, pressure,
or sound-speed squared. A positive total heat capacity does not certify the
accuracy of its smaller electron contribution. These tests do not establish
merger accuracy, stellar equilibrium, shock convergence, or performance.
