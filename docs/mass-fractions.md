# Material fractions and massless tracers

Enable at build time with `OCTOII_WITH_MASS_FRACTIONS=ON` (the default), and
at runtime with the following INI section. Runtime default is off.

```ini
[massFractions]
enabled = on
species = fuel:0.7:He=70%,O=30%;ash:0.3:A=56,Z=26;dye:0.2:A=0,Z=0
```

Each semicolon-separated definition is `name:initialFraction:composition`.
For a pure element use, for example, `helium:1:He` or `helium:1:HELIUM`.
For an explicit isotope use `carbon12:1:A=12,Z=6`. Atomic masses are numerical
values in atomic mass units. Names must be identifiers and unique ignoring case.
CLI values need shell quoting, for example:

```bash
./octoII-3d --problem.name=collapse --massFractions.enabled=on \
  '--massFractions.species=fuel:70%:helium=70%,oxygen=30%;ash:30%:A=56,Z=26;dye:0.2:A=0,Z=0'
```

Both initial material fractions and each elemental mixture must sum to one.
Mixture entries accept decimal mass fractions or percentages, never number
fractions. All 118 element symbols and full English names are case insensitive;
aluminum/aluminium, cesium/caesium and sulfur/sulphur are accepted. The table
uses [CIAAW abridged standard atomic weights](https://ciaaw.org/abridged-atomic-weights.htm)
(2024); entries without a standard weight have representative isotope mass
numbers from the [IUPAC table](https://iupac.org/wp-content/uploads/2022/05/IUPAC_Periodic_Table-04May22.pdf).
Explicit A,Z values override elemental defaults. Element mixtures retain their
constituents and expose effective values satisfying

    1/Aeff = sum(Xe/Ae)
    Zeff/Aeff = sum(Xe*Ze/Ae).

The existing ideal-gas EOS and configured `hydro.meanMolecularWeight` remain in
use. Composition metadata does not assume an ionization model or implement
nuclear reactions.

## Stored state and coupling

The composition module owns partial densities `qs = rho*Xs`. With the module
enabled, the persistent active hydro field banks allocate no total-density
column. Their density handle sums the material columns on demand, including
halo reads and gravity-source reads. Hydro working states, snapshots, and AMR
shadow estimates contain reconstructed rho; these are not an independently
evolved total-density field. With the module disabled, ordinary hydro density
storage and transport are used.

Hydro stores its final, positivity-limited numerical mass flux on every face.
Composition partitions this flux using donor-cell concentrations selected by
its sign. The largest material component receives the residual flux so that
material fluxes sum to the hydro mass flux to roundoff. The same stored mass
flux supplies gravitational energy work, even when composition is disabled.

Concentration transport is currently first order at contacts. Hydro retains
its existing MUSCL-Hancock reconstruction. No clipping or post-update
renormalization of species is used. Negative partial densities fail explicitly.
AMR restricts partial densities conservatively, prolongs density with injected
concentrations, and refluxes each partial density with fine-face fluxes.

The entropy auxiliary stays exclusively in hydro. It is never included in
species normalization or in the sum defining mass density.

## Massless tracers and initial conditions

`A=0,Z=0` declares a passive massless tracer. It stores `rho*C` and advects
with flux `C*F_rho`, but contributes to neither inertial density nor gravity.
Tracer concentrations are finite and nonnegative, need not sum to one, and
may exceed one. Material A,Z must satisfy `A >= Z > 0`; setting only one to
zero is invalid. At least one material species is required.

Built-in problems initialize every component using its configured initial
fraction or tracer concentration. A custom `initializeProblem` may fill
`Snapshot::species` in definition order to prescribe spatially varying partial
densities; initialSnapshot validates size and positivity and derives gas rho
from those fields. The custom initializer must supply compatible momentum,
energy and entropy. Analytic inflow uses configured initial concentrations
multiplied by the analytic hydro density. Other boundary types copy or average
scalar partial densities through the existing halo map.

Silo exports `partialDensity_NAME` and `massFraction_NAME` for materials,
`tracerDensity_NAME` and `tracer_NAME` for tracers. Tracer density means carrier
rho times concentration, not gravitating mass density. Density-unit quantities
are in g/cm^3; fractions and concentrations are dimensionless.
