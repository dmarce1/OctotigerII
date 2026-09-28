# Rotating radiation star: research record, September 2026

This records the equations, decisions, experiments, evidence, and open questions
from the radiating-star discussion. It is intended to survive changes to the
active implementation and to make the more ambitious equilibrium work resumable.
It is a technical research record, not a verbatim conversation transcript.

## Decision at the pause point

On 2026-09-27 the user replaced the fully radiating equilibrium target with a
simpler benchmark: an optically thick rotating gas/radiation star, physical
light speed, no photon heater, no radiation energy escaping the computational
domain, and a bulk thermal time much longer than its dynamical time. Surface
adjustment and slow thermal evolution are acceptable. Escaping photon angular
momentum had already been excluded as a reason to delay the benchmark.

The active design is [../../radiating-star-design.md](../../radiating-star-design.md).
The earlier [leaking equilibrium equations](leaking-equilibrium-design.md) are
preserved here. They describe a different problem, including artificial spatial
opacity and a distributed photon source, and must not be mistaken for the
current source-free setup.

## What is saved

- [research-snapshot.tar.gz](research-snapshot.tar.gz) preserves 145 files:
  experimental Python solvers and audits, intermediate and accepted reference
  arrays, logs, the exported C++ table, and a snapshot of relevant repository
  source. It includes rejected experiments; its contents are not all validated.
- [snapshot-manifest.json](snapshot-manifest.json) records each file's original
  byte count and SHA-256, plus the archive's checksum.
- [prototypes/table-problem.cpp](prototypes/table-problem.cpp) preserves the
  unfinished table-based OctoII problem wrapper. Its table reader was not
  integrated or tested in the application. Do not treat it as a runnable module.
- [opaque_star_scales.cpp](opaque_star_scales.cpp) and its
  [output](opaque_star_scales.txt) preserve the calculation supporting the new
  simple benchmark's optical depths and thermal/dynamical scale separation.
- The maintained offline generator and instructions are under
  [../../../tools/science](../../../tools/science/radiating_star_reference.md).
  A dated copy is also inside the archive.

These are local repository files. This operation does not create a Git commit
or a remote backup. The conversation remains useful context, but the equations
and numerical evidence can now be recovered without reconstructing it.

## Coupling and integration recovered before the star work

An earlier OctoI `dominic-wip` branch contained gas/radiation coupling that had
not carried into the OctoII transport-only path. The restored equal-gray-opacity
terms use physical laboratory moments and inertial velocity:

\[
 B=aT^4,\quad\chi=\rho\kappa,\qquad
 Q=\chi\left[c(E-B)-\frac{\mathbf v\cdot\mathbf F}{c}\right],\qquad
 \mathbf G=\frac\chi c(\mathbf F-\mathsf P\mathbf v-B\mathbf v).
\]

Here gas gains total energy at rate `Q` and momentum at rate `G`; radiation
receives the opposite exchange at physical light speed. Gas internal heating
is `Q-v·G`, not `Q` alone. The accepted representable radiation change is paired
with the gas change to preserve the combined discrete budgets. Stable kinetic
energy differences are used for the dual-energy auxiliary variable.

The source solver uses local SDIRK2 stages, with `gamma=1-1/sqrt(2)`, Newton
iteration, positivity/realizability checks, and transactional retries. Transport
uses source-aware midpoint states. The usual Hancock predictor is disabled in
that coupled path because the midpoint construction has already predicted the
state. Gravity, rotating-grid kinematics, subcycling, coarse/fine flux registers,
and source/boundary ledgers have to use the same accepted stage data. See
[coupling](../../radiation-coupling.md) and
[validation](../../radiation-validation.md) for implementation and evidence.

### Why an implicit spatial pressure solve was unnecessary

For stationary gas,

\[
 \partial_t\mathbf F=-c\widehat c\nabla\cdot\mathsf P
                        -\widehat c\chi\mathbf F.
\]

The diffusion limit is `F=-c div(P)/chi`, or `-c grad(E)/(3 chi)` for isotropic
pressure. Applying an explicit pressure kick and then independently damping it
exponentially erases the continually replenished diffusion flux when the source
is stiff. An explicit gradient with backward-Euler damping instead gives

\[
 \mathbf F^{n+1}=
 \frac{\mathbf F^n-c\widehat c\Delta t(\nabla\cdot\mathsf P)^*}
 {1+\widehat c\chi\Delta t}.
\]

It has the correct stiff limit without a global implicit spatial solve.
Higher-order stages must retain the analogous forcing/damping balance. A
separate correction of the numerical radiation **energy face flux** is needed
so HLL dissipation does not overwhelm physical diffusion in thick cells. Passing
the local source test alone does not establish this asymptotic property.

The reduced-speed approximation also has a dynamic-diffusion restriction:
roughly `chat >> v max(1,tau)` on the relevant physical scale. Nonfatal optical
depth/trapping diagnostics were added, but the active star uses `chat=c` and
does not require an RSLA design decision.

## Structural choices and why exact radiative balance became complicated

For a monatomic ideal gas, `pg=rho Rgas T`, `ug=3 pg/2`; in thick LTE,
`E0=aT^4`, `pr=E0/3`. A constant-entropy mixture is not generally a single
polytrope, and a mechanical barotrope does not automatically satisfy a
stationary luminosity equation. The user explicitly relaxed global constant
entropy in favor of stability.

A useful stable structural family is

\[
 P=K\rho^{1+1/n}=\rho\mathcal RT+\frac{aT^4}{3},\quad n=3.5,
 \qquad H=\int\frac{dP}{\rho}=(n+1)K\rho^{1/n}.
\]

`H` is a barotropic integral, not thermodynamic enthalpy when entropy varies.
Uniformly rotating hydrostatic structure obeys
`H+Phi-Omega² R²/2=C`, with `laplacian(Phi)=4 pi G rho`. The implemented
isolated multipole self-consistent-field solve constructs this model. Its
structural exponent `9/7` is below the mixture's adiabatic `Gamma1>=4/3`, giving
stable adiabatic stratification. This does not prove stability under arbitrary
radiative, thermal, or nonaxisymmetric perturbations.

For constant-beta Eddington structure, uniform rotation changes the heating
needed to balance diffusion to

\[
 j=\frac{c(1-\beta)}\kappa(4\pi G\rho-2\Omega^2),
\]

which is negative in sufficiently tenuous layers. In the variable-beta family,
write `Psi=Phi-Omega²R²/2`, `b=dpr/dP`. Then

\[
 \mathbf F_0=\frac c\kappa b\nabla\Psi,\qquad
 j=\frac c\kappa\left[b(4\pi G\rho-2\Omega^2)
                       -\frac{db}{dH}|\nabla H|^2\right].
\]

The second term helps positivity but the formal zero-temperature surface makes
the required heater singular for this family. A finite outgoing luminosity also
cannot coexist with `E->0` while satisfying `|F|<=cE`. This motivated a separate
transport envelope, not a density-floor patch to an alleged exact equilibrium.

The actual axisymmetric leaking construction needed Poisson gravity,
meridional gas force balance, full M1 momentum transport, moving-gas thermal
balance `Q-v·G=0`, and `j=div(F)+Q`. The photon heater was a prescribed spatial
function; it was never a feedback force that reset evolved fields. No artificial
mechanical torque was introduced. Photon spin loss remained physical, so even
this construction did not promise exact stationarity of every gas component.

## Spherical full-M1 control: implemented, with limited evolution evidence

The `radiating-sphere` module uses an artificial fixed spatial opacity
`kappa(r)=kappa0 max(0,1-rho_cut/rho_reference(r))²` and a transparent gas
envelope. The heater is inferred from the same continuous flux interpolant.
It is distributed, not a compact Ensman-style central lamp. Its differential
balance and reference-tolerance checks pass.

At default optical coefficient 240 and cutoff density fraction 0.016, the
material radius is `6.3941644 alpha`, opacity transition `3.5708999 alpha`,
and central optical depth about 286.41. Over 10.8652 seconds, the coarse 8³ and
16³ production meshes gave density drifts 0.0029424 and 0.0011921, and corrected
energy errors at most about 1.2e-15. This is a plumbing/conservation check over
only a small fraction of a dynamical time, **not** a long-term stability result.

## Rotating leaking reference: candidate converged, application unfinished

The saved spin-0.1 candidate uses `n=3.5`, central gas fraction 0.8,
dimensionless opacity coefficient 100 and cutoff density 0.015. Units are
`4 pi G=rho_c=Rgas Tc=1`, with `F/c` stored as the flux variable. The spherical
reference has radius 15.207814, mass 23.592019, and opacity transition 8.568680.
The physical dimensionless light speed is 520.3917149 for `rho_c=1`, `mu=0.6`.

The equatorial rotation rate is prescribed. Away from the equator, the solved
force balance gives differential rotation. The accepted graded grid is 193×65
in `(r,mu)`, with outer radius four times the opacity-transition radius.
The final file is `scratch/radiating_star_reference_193_65_0.1_moving.npz` inside
the archive; the raw exported include is `scratch/radiating_star_reference.inc`.

The moving-equation algebraic residual reached 8.19e-10. A dense near-transition
check found positive heating with minimum 6.834e-6; angular velocity squared was
positive throughout the checked material. Independent continuum checks for the
static meridional reference showed approximately second-order convergence:

| RMS residual | 97×33 | 193×65 |
| --- | ---: | ---: |
| Radial radiation momentum | 1.313e-5 | 3.430e-6 |
| Polar radiation momentum | 5.464e-6 | 1.355e-6 |
| Vertical gas force | 4.592e-7 | 1.162e-7 |

The preserved logs distinguish the static convergence audit from the final
moving-equation solve. There was no completed OctoII evolution of the rotating
table reference. Adiabatic Rayleigh/Solberg–Høiland checks passed in the thick
LTE and transparent-gas regions, but vertical shear remained:
`R partial_z(log Omega)` was about -0.038 in a thick subset. Thermal diffusion
can allow shear instabilities excluded by the adiabatic criteria. The earlier
claim of adiabatic stability must not be promoted to a demonstrated stable star.

### Numerical lessons worth retaining

1. Centered Cartesian collocation could yield tiny nonlinear residuals while
   producing negative inferred rotation squared and nonconvergent fields.
2. An integrated gas-enthalpy variable improved conditioning but did not fix
   the centered radiation-momentum mode.
3. Spherical coordinates align the opacity interface. Representing
   `F_theta/c=-sqrt(1-mu²) E q` avoids a polar-axis singularity.
4. Second-order outward radial differentiation removed the checkerboard mode;
   grading around the narrow M1 critical layer was necessary.
5. Keep both the spherical background and rotating correction interpolants
   piecewise at the opacity cutoff. A global cubic interpolant smeared a
   derivative change and produced roughly 4% force error despite small nodal
   residuals. Enforce opaque thermal balance after interpolation as well.
6. Source-free exterior `j=0` is a prescribed equation. Small signed residuals
   of an interpolated flux divergence should not be clipped into a new heater.
7. Positivity and independent continuum derivatives matter more than the
   nonlinear solver's success flag. The spin-0.2 branch lacked the same
   near-cutoff positivity margin; spin 0.1 was the accepted candidate.

The toroidal radiation equation and retained moving-gas thermal correction were
implemented separately and checked against Cartesian manufactured derivatives.
The thermal equation used was
`E-B-2 beta Fphi/c+beta²(Pphiphi+B)=0`. Terms tiny compared with `E` could still
be appreciable compared with the very small required heater; this motivated
retaining the complete balance instead of adding isolated second-order terms.

## Restarting the deferred investigation

Extract the archive into an empty directory, retaining `scratch/` so the old
audit scripts can import their neighboring modules. Some scripts contain old
`/tmp` paths; use the maintained generator for a fresh calculation. NumPy 1.26.4
and SciPy 1.11.4 were used. The reproduction command is:

```sh
python3 tools/science/generate_radiating_star_reference.py \
  --nr 193 --nm 65 --extent 4 --spin 0.1 \
  --output /tmp/radiating-star-reference.npz \
  --include /tmp/radiating-star-reference.inc
```

Before reviving the table module: validate its reader and interpolation against
the saved arrays, repeat moving-equation continuum residuals and heating
positivity, measure initial force/thermal errors on Cartesian meshes, and run
resolved perturbation/evolution tests with conservation ledgers. The main open
physics issue is diffusive/shear stability; the main integration issue is an
accurate transparent-layer interpolation and boundary treatment. Neither was
resolved by the decision to simplify the first benchmark.

The original source papers and the distinction between their results and our
derived benchmark equations are listed in
[the archived design](leaking-equilibrium-design.md#references-and-scope).
