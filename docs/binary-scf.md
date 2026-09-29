# Binary self-consistent field initialization

`binary-scf` constructs a circular, synchronously rotating, isolated binary in
3D, with gravity and ideal-gas hydro. Each star has independent core and envelope
structural indices. Evolution uses **gamma = 5/3**, irrespective of those indices.
Radiation, Helmholtz, asynchronous rotation, and a shared envelope are excluded.

The initial implementation uses a **uniform reference mesh**, followed by
conservative transfer to the evolution mesh. `scf.referenceWidth` sets the
reference box width in units of the separation, independently of the larger
evolution box. The reference must still contain both stellar lobes; an occupied
boundary cell is rejected. SCF does not iterate on distributed AMR leaves.
Construction runs once per process and is cached across block initializations;
multiple HPX localities currently repeat the reference solve. Evolution uses
the ordinary distributed gravity and hydro machinery, including AMR when
requested.

## Run

```sh
cmake --build release --target octoII-3d -j 4
./release/octoII-3d --config=problem/science/BinaryScf/inputs \
  --hpx:threads=4 --hpx:bind=none
```

The supplied input constructs an unequal-mass semi-detached binary and stops
at time zero. `scf.json` reports the reference solution and the populated hydro
mesh; `conservation.csv` records the usual budgets. Silo output is unchanged.
Set `scf.evolveOrbits=1` to evolve for one initial orbital period directly
after construction. This sets the stop time from the converged
`2*pi/omega`; alternatively set a positive `runtime.stopTime` in seconds. Do not
specify both positive controls. `scf.framesPerOrbit=100` writes 100 evenly spaced Silo frames per initial
orbit, plus the frame at time zero. The evolution timestep is clipped to each
output time, so this is a physical-time cadence rather than a step count.

Other input files in the same directory:

* `detached.ini`: both stars underfill their Roche lobes.
* `bipolytrope.ini`: core indices 3, envelope indices 1.5, density jumps 2,
  and a Roche-filling donor. These structural choices follow the family studied
  by [Kadam et al. (2018), section 2](https://arxiv.org/abs/1809.04884).
  Its masses and interface locations are illustrative, not a reproduction of
  that paper's particular binaries.
* `dwd-polytropes.ini`: n=3/2 polytropes with a shared K and Roche-filling
  donor. The accretor filling factor is solved from its mass at that same K.
  The dimensional K follows from the chosen mass and separation; it is not
  fixed to the electron-degeneracy value for a specified composition.
* `qb4-dwd.ini`: q=0.7 DWD proxy with a Roche-filling donor, compact SCF
  reference, density AMR, and timed orbital output. See the
  [QueenBee4 instructions](queenbee4-binary-scf.md).

The HPX-free executable accepts the same inputs without the `--hpx:` options.
Radiation and mass fractions can both be disabled at build time.

## Physical controls

All dimensional inputs use CGS.

| Input | Meaning |
| --- | --- |
| `scf.primaryMass` | Primary mass, excluding atmosphere (default 1.98847e33 g) |
| `scf.massRatio` | Donor mass / primary mass (default 1) |
| `scf.separation` | Final separation of the two mass centers (default 1e11 cm) |
| `scf.commonPolytropicK` | Share K between identical single polytropes; solve primary fill (default off) |
| `scf.referenceWidth` | Uniform SCF box width / separation, in [1.5,10]; default 3 |
| `scf.primary.coreIndex`, `scf.donor.coreIndex` | Positive, finite core structural indices |
| `scf.primary.envelopeIndex`, `scf.donor.envelopeIndex` | Positive, finite envelope structural indices |
| `scf.primary.interfaceFraction`, `scf.donor.interfaceFraction` | Core-side interface density / central density, in (0,1] |
| `scf.primary.densityJump`, `scf.donor.densityJump` | Core-side / envelope-side interface density, at least 1 |
| `scf.primary.fill`, `scf.donor.fill` | Enthalpy filling factors, in (0,1]; default .8 and 1 |
| `scf.evolveOrbits` | Stop after this many initial orbital periods; default 0 |
| `scf.framesPerOrbit` | Timed Silo frames per initial orbit; 0 uses `output.every` steps |

**`scf.donor.fill=1` puts the donor surface at the inner Lagrange saddle.**
For star s the surface Bernoulli constant is

```
Phi_eff = Phi - (omega^2/2)*((x-x_rotation)^2+y^2)
C_s = min_s(Phi_eff) + fill_s*(Phi_eff(L1)-min_s(Phi_eff))
H_s = C_s - Phi_eff
```

Thus `fill` is neither a radius fraction nor a volume filling fraction. Both
fill values below one give a detached binary. Both equal to one give marginal
contact at L1; overcontact/shared-envelope models are not supported.

With `scf.commonPolytropicK=on`, both stars must have the same structural
index and no density jump. The donor fill remains prescribed. At each gravity
iteration, the solver uses `rho=A*H^n`, fixes A from the donor mass, and solves
the primary surface constant to recover the primary mass at that same A.
Thus both stars have `K=A^(-1/n)/(n+1)`. The primary fill is an output;
`scf.primary.fill` only seeds the starting density in this mode. A primary
that cannot fit within its lobe is rejected. For n=3/2, isolated stars with
this barotrope recover `R proportional to K M^(-1/3)`; binary tides modify
the radii. Ideal-gas evolution still uses gamma=5/3.

The core mass fractions are **outputs**, rather than independently imposed
constraints. Adjust the interface density fractions to select core sizes.
Equal core/envelope indices with density jump 1 recover a single polytrope.
The defaults use n=1.5 throughout. Positive indices above 5 are not excluded
algebraically for composite stars, but the requested model must converge to
bounded, resolved lobes. The incompressible n=0 limit is not implemented.
Hydro stability is a separate question: allowing an index does not assert
convective, dynamical, or mass-transfer stability with gamma=5/3.

## Core/envelope matching and energy

With interface densities `rho_ci=f*rho_0`, `rho_ei=rho_ci/jump`, the common
interface pressure is determined by the central enthalpy. The structural
enthalpy is the continuous integral of `dP/rho` from vacuum:

```
Envelope: H = (n_e+1)*P_i/rho_ei * (rho/rho_ei)^(1/n_e)
Core:     H = H_i + (n_c+1)*P_i/rho_ci * ((rho/rho_ci)^(1/n_c)-1)
```

Pressure and this enthalpy match at the density jump. Intermediate densities
during iteration/remapping represent an interface mixture at the common pressure.
The update inverts the continuous enthalpy graph with a proximal step when a
cell crosses the interface. This admits a mixed state exactly at H=H_i and
avoids alternating between the two pure densities. Each star must separately
pass an enthalpy-balance tolerance as well as the density-update tolerance.
Temperature continuity is not imposed; a density jump can represent an entropy
or composition contrast. The structural enthalpy need not equal the ideal-gas
thermal enthalpy when the structural indices differ from 1.5.

Initialization constructs primitive density, pressure, and inertial rigid-body
velocity, then uses `HydroSystem::conservedState` to set **both** total energy and
the dual-energy auxiliary. In particular,

```
internal energy density = pressure/(5/3 - 1)
total energy density    = internal energy density + rho*|v|^2/2
```

An explicitly specified `frame.omega` is respected. Otherwise the grid rotates
at the converged binary angular velocity. Hydro always stores inertial momentum.
No subsequent relaxation, external driving, or angular-momentum removal is
applied during evolution.

When mass fractions are built, initialization uses five material fields in this
order: `primary_core`, `primary_envelope`, `donor_core`, `donor_envelope`,
`atmosphere`. Their sum is exactly the hydro density. The default material labels
are helium; custom material compositions may be supplied under those same names
and in the same order. Those labels do not supply a degeneracy EOS.

## Numerical method and acceptance

The reference solve uses G=1 and total stellar mass 1. Target component masses
are normalized each iteration; the initial orbital angular momentum fixes the
otherwise free length scale. The final model is rescaled to the requested mass
and separation. A seed based on approximate Roche radii is only a starting guess.
The potential is recomputed from the current density before updating enthalpy.

Repeated gravity uses a zero-padded isolated FFT convolution of the cell-center
`-1/r` kernel, omitting the self term exactly as `gravity::solve` does. It has no
periodic images and introduces no gravitational softening. The kernel is
verified against independent all-pairs summation, including boundary cells.
L1 is located using direct evaluation of the same point masses on the symmetry
axis. A connected-component restriction excludes the unphysical exterior
centrifugal branch. Reflection symmetry about y=0 and z=0 is imposed.

Bounded Anderson mixing accelerates the density iteration; `scf.history=0`
selects ordinary damped iteration. Negative extrapolations and excessive steps
fall back to ordinary mixing. Acceptance uses the **undamped**, pre-mass-
normalization proposed density change, separately normalized by each star's
mass, so damping and mass rescaling cannot manufacture convergence.
The separate per-star Bernoulli check prevents the interface preconditioner
from concealing a stalled enthalpy balance, including in a low-mass donor.

| Control | Default |
| --- | --- |
| `scf.cells` | 64 per axis; power of two, 16 through 128 |
| `scf.tolerance` | 1e-5, maximum component L1 density residual |
| `scf.maxIterations` | 1000; exceeding it is an error |
| `scf.relaxation` | 0.4 |
| `scf.history` | 4, allowable range 0 through 6 |
| `scf.virialTolerance` | 0.05 |
| `scf.atmosphereFraction` | 1e-10 in reference density and pressure units |

The scalar virial diagnostic is `abs(2*T+W+3*integral(P dV))/abs(W)`, where W
includes the factor 1/2. It is evaluated with the **final density and its own
potential**, before any handoff. A converged iteration that fails the virial
gate is rejected. The default 5% gate is a coarse initialization safeguard,
not a publication-accuracy criterion; tighten it for the intended experiment.

By default the hydro box and initial level follow the reference scale and
resolution. Transfer integrates overlaps of piecewise constant reference cells,
preserving stellar mass and pressure integrals across uniform meshes and AMR.
Translation to put the barycenter at the coordinate origin can mix neighboring
reference cells. Rotation is evaluated at hydro cell centers. These operations
can change kinetic/binding energies and local force balance, so the runtime
recomputes gravity and checks **hydro-mesh mass and virial error** before evolution.
An explicitly cropped box is rejected if it loses more than `scf.tolerance` of
the expected total mass. The uniform atmosphere is added once, is separately
tracked, and is included in the hydro gravity solve and subsequent budgets.

`scf.json` separates the reference diagnostics from the hydro-mesh diagnostics.
It reports each star's volume-equivalent diameter, x/y/z extents of occupied
cells inside its Bernoulli surface, solved filling factor, and (for single
polytropes) K in CGS. Surface widths are resolved only to the reference cell
spacing. These distinguish tidal elongation from volume-equivalent size.
Conservation of mass under transfer does not imply conservation of its binding
energy or preservation of local hydrostatic balance. A finer AMR mesh cannot
recover structure missing from the uniform reference. Check both resolutions.

## Verification and remaining work

See [the validation record](validation/binary-scf.md) for commands and measured
results. Unit tests cover analytic single-polytrope limits, interface matching,
`dH=dP/rho`, isolated gravity, independent indices, Roche filling, material mass
closure, energy initialization, failure gates, and configuration serialization.
Integration tests exercise the first gravity/hydro step and its budgets.
The orbital-cadence integration test checks three Silo frames at exact initial
period fractions even when `output.every` would suppress intermediate frames.

When profiling is enabled, APEX names the synchronous `scf.solve`,
`scf.gravity_fft`, `scf.l1_search`, `scf.common_k_constraint`,
`scf.density_update`, `scf.mixing`, and
`scf.handoff_block`
regions. The existing runtime regions cover FMM workers, hydro transport,
diagnostics, and Silo writing. See [profiling](profiling.md); wrap each HPX
locality with `profile.sh` to keep its report separate.

The resolution study is explicitly opt-in:

```sh
./release/3d/tests/binaryScfChecks-3d \
  --gtest_also_run_disabled_tests --gtest_filter=BinaryScf.DISABLED_ResolutionStudy \
  --hpx:threads=2 --hpx:bind=none
```

This implementation does not establish quiet evolution for many orbits,
converged mass-transfer rates, agreement with a particular published binary, or
large-machine SCF scaling. Those require separate experiments. The current
uniform 128^3 limit also restricts how well very condensed cores can be resolved.
A distributed AMR SCF solve and cold-white-dwarf structural EOS are subsequent
extensions, as is the V1309 common-envelope branch.
