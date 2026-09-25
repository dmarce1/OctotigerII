# Dual-energy hydrodynamics

Dual energy is enabled by default in each hydro build. It adds one scalar to
`hydro::ConservedState`, its flux, field storage, halos, and HPX serialization.
It uses the same conservative transport, restriction, prolongation, and AMR
reflux machinery as the Euler fields.

## Variable and evolution

Let rho be mass density, u internal energy **density**, gamma the constant
ideal-gas adiabatic index, and alpha a finite, nonzero configurable exponent.
The requested variable is

    A = rho (u / rho^gamma)^alpha
    u = rho^gamma (A / rho)^(1/alpha).

For smooth adiabatic flow, D(u/rho^gamma)/Dt = 0. Combining this with mass
continuity gives

    dA/dt + div(A v) = 0.

Thus any nonzero alpha works mathematically, including negative values, when
rho and u are positive. Alpha=0 would leave only A=rho and lose all thermal
information. Alpha=1 is the entropy-density form (up to the constant gamma-1
when entropy is defined using pressure); alpha=1/gamma gives u^(1/gamma).
These transformations are equivalent for smooth continuum solutions but need
not give identical finite-volume errors, mixing, or shock profiles.

At shocks the passive equation alone does not generate the physical entropy
increase. The reliable total-energy solution resets A using the upper
threshold below. Additional nonadiabatic heating/cooling would require its own
consistent auxiliary update; the current code has no gas-radiation exchange
solver.

Because both gamma and alpha are runtime values, the implementation takes powers
of a dimensionless ratio using fixed references rho0=1 g/cm^3 and u0=1 erg/cm^3:

    A = rho [(u/u0)/(rho/rho0)^gamma]^alpha.

This has exactly the requested numerical value when the inputs are expressed in
CGS. The stored A has density units, its flux has mass-flux units, and its
primitive A/rho is dimensionless. There is no selectable code-unit scale.
Logarithmic conversion avoids intermediate power overflow. Exponents that make
a required value fall outside the positive normal floating-point range fail
explicitly instead of silently clipping the entropy. Extremely small exponents
also lose thermal information through roundoff even though the continuum inverse
exists; order-unity exponents are the useful numerical choices.

## Two thresholds

Define the local resolved thermal energy uE=E-K, with
K=sum(momentum_i^2)/(2 rho). E is gas total energy density; gravitational
potential energy is not included in this denominator.

| Condition | Pressure and temperature | Auxiliary update |
| --- | --- | --- |
| uE <= eta1 E | Use u recovered from A | Preserve transported A |
| eta1 E < uE <= eta2 E | Use uE | Preserve transported A |
| uE > eta2 E | Use uE | Reset A from rho and uE |

The defaults are eta1=0.001 and eta2=0.1. The comparisons are strict `>`;
validation requires 0 <= eta1 < eta2 < 1. This uses each cell's own E, without a
neighbor-maximum criterion. Synchronization changes only A: E, density, and
momentum keep their conservative updates. Valid auxiliary pressure can handle
zero or negative E-K due to kinetic-energy cancellation; the full state must
still be finite with positive A and valid density and pressure floors.
Conservative E itself may be nonpositive when the auxiliary supplies valid
thermodynamics (for example after conservative gravity work in an atmosphere).
It is retained without clipping; nonpositive E is never used for pressure or
auxiliary synchronization. With dual energy disabled, E must remain positive.

Every transport timestep finishes with synchronization after all AMR reflux
corrections and before publication of the next field bank. The owning-patch
solver and coarse shadow evolution also synchronize completed steps. Gravity
kicks in a conservative self-gravity step retain the evolved A until the final
mass-flux energy correction, then apply the upper criterion. Standalone kicks
perform the check immediately after their momentum and kinetic-energy updates.
Reconstruction carries A/rho independently of pressure, including between the
two thresholds; HLLC star states preserve this specific auxiliary value.

## Options

An INI file may contain:

```ini
[hydro.dualEnergy]
enabled = on
exponent = 1
pressureThreshold = 0.001
syncThreshold = 0.1

[hydro]
meanMolecularWeight = 1
```

Equivalent command-line options are
`--hydro.dualEnergy.enabled=off`, `--hydro.dualEnergy.exponent=-0.5`,
`--hydro.dualEnergy.pressureThreshold=0.001`, and
`--hydro.dualEnergy.syncThreshold=0.1`. Disabling dual energy restores pressure
from E-K and suppresses auxiliary initialization and synchronization. The extra
storage slot remains zero in normally initialized disabled runs.

`hydro.meanMolecularWeight` is the positive mean particle mass in atomic mass
units, default 1. Temperature is T=P mu m_u/(rho k_B); changing mu affects
reported temperature, not pressure or dynamics. Both P and T use the same
selected internal energy.

A newly constructed primitive has auxiliary=0, which tells `conservedState` to
initialize A from its density and pressure. A positive primitive auxiliary is
an explicitly supplied A/rho and is preserved. When reusing a reconstructed
primitive to prescribe a new thermodynamic initial state, clear its auxiliary
if A should be reinitialized from the new pressure.

## Output and verification

Silo adds `dualEnergy` (normalized A, g/cm^3, enabled runs only),
`internalEnergy` (selected u, erg/cm^3), and `temperature` (K).
`pressure` and `minimumPressure` use the same dual-energy selection.
The thermal-energy integral now integrates selected u. In cells using the
auxiliary, kinetic plus reported thermal energy can differ from the conserved
gas-energy integral. Conservation and gas-plus-gravity diagnostics continue to
use the independently conserved E. A is intentionally absent from the global
conservation ledger because synchronization supplies entropy at shocks.

`dualEnergyChecks` runs in 1D, 2D, and 3D and covers nonzero positive/negative
exponents, compression, both thresholds and equality, complete thermal-energy
cancellation, disabled behavior, independent reconstruction, contact fluxes,
cold and intermediate-threshold periodic transport, conservative AMR transfer,
and synchronization at every runtime step after mixed-level refluxing.
Options, output, and HPX serialization tests cover the new fields as well.

```bash
cd "$HOME/workspace/OctotigerII"
ctest --test-dir release -R dualEnergyChecks --output-on-failure
```

Use the actual configured build directory if it differs from `release`.
