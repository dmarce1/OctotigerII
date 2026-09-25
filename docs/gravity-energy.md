# Gravity energy from the numerical mass flux

## Independent controls

```ini
[gravity]
timeIntegration = hierarchical
energyTreatment = mullen
conserveRegridEnergy = on
```

`timeIntegration=hierarchical|conventional` selects the momentum schedule for
AMR time refinement. `hierarchical` is the default; `conventional` evaluates
the full force at each active spatial level's own cadence.
`energyTreatment=mullen|naive` independently selects the timestep energy source.
`conserveRegridEnergy=on|off` independently selects the regridding invariant.
The energy and regrid options do not select the momentum schedule.
The naive mode retains kinetic-work increments. It still records the same
endpoint potential boundary flux, so
its combined-energy diagnostic is directly comparable to the Mullen mode.

For example, append `--gravity.energyTreatment=naive` to an existing run to
change only the timestep treatment. Append `--gravity.conserveRegridEnergy=off`
to change only remapping. Use separate `output.directory` values for comparisons.

Set `--timestep.refinement=off` for the globally synchronized reference sequence
below. Disabling AMR also selects that sequence. Imposed uniform external
acceleration currently requires global stepping, including when self-gravity
is enabled. See [time refinement](time-refinement.md) for the level schedule.

## Globally synchronized sequence

The global momentum integrator retains its endpoint kicks:

1. Save potential and acceleration at time n, and zero the kick-work register.
2. Apply the first half momentum kick with the field at n.
3. Advance hydro and species, storing the actual numerical density flux; reflux.
4. Solve gravity for the updated density at n+1.
5. Apply the second half momentum kick with the field at n+1.
6. In Mullen mode, replace accumulated self-gravity kick work with endpoint-averaged
   mass-flux work, then synchronize hydro's entropy auxiliary. Naive mode retains
   the kick work.

Potential and acceleration remain synchronized at each endpoint. There is no
midpoint-density gravity solve and no separately evolved dphi/dt field.
`simulation::run` performs the sequence automatically. Manual Runtime clients
bracket the usual kick/advance/solve/kick sequence with `beginGravityEnergy()`
and `finishGravityEnergy(dt)`. Regridding is disallowed inside that bracket.
Calling a kick outside it retains the standalone kinetic-work behavior.

## Time-refined sequence

`Runtime::advanceGravity(dt)` performs the complete refined interval; the
simulation driver uses it automatically. The gas states remain physical while
an explicit midpoint predictor uses the discrete numerical-flux divergence,
momentum forcing, and gravitational mass-flux work rate. The midpoint donor
states are formed before halo prolongation and face reconstruction; the shared
Riemann solver evaluates those states with the Hancock time predictor disabled.
The conservative update starts from the original gas state. Each hierarchical
shell records provisional
impulses from its frozen opening field, then replaces them after reflux by the
paired endpoint impulse. This follows HOLD's accepted impulse quadrature; the
intermediate gas map accounts for Eulerian mass transport and is not claimed to
be symplectic. Conventional mode uses each active level's full-force endpoint
impulse instead.

Fine halos and conventional inactive-source densities use a separate coarse
forecast: the old state plus the coarse duration times the initial full
numerical derivative, with canonical fine fluxes in its probe-time reflux.
It occupies state bank 3, separate from the unrefluxed transport endpoint and
the checkpoint in bank 2. A predictor-only O(dt²) convex kinetic remainder is
included in this forecast and midpoint states; it does not enter the accepted
conservative energy update.

Energy closure uses interaction shells in both momentum modes. A shell records
the time-integrated mass flux across every face for its entire interval, then
evaluates the work formula below with its two endpoint potentials. Reflux and
work use the same canonical fine subface transfers. A coarse level forecasts
the work of its nested active set; descendants defer their exact work on coarse
neighbors until that level closes. The closure replaces this nested forecast
by the shell work plus descendant transfers. Faster cells replace only the
provisional work they recorded for the closing ancestor shell.

Thus the accepted Mullen source is the sum of all shell endpoint work, including
both sides of every transfer between timestep groups. No density-rate solve or
global redistribution of an energy error is needed. Every shell and deferred
transfer closes before publication or regridding. The full equations and
conditional smooth-time accuracy argument are in the
[coupling derivation](gravity-time-coupling-derivation.md). Smooth fixed-mesh
regressions observed second-order temporal convergence in both modes; see the
[validation report](validation/gravity-time-integration.txt) for scope and
results. No performance or general shock-convergence claim follows from the
energy identity.

## Discrete energy source

Let F be the accepted interval-average numerical mass flux, oriented from left
to right, and let
phibar = (phi^n + phi^(n+1))/2. On an equal-width interior face,

    gbar_face = -(phibar_R - phibar_L)/dx
    delta_E_L = delta_E_R = dt * F * gbar_face / 2.

These are energy-density increments. This is the mass-flux energy source of
[Mullen, Hanawa & Gammie (2021), equations 57–59](https://arxiv.org/abs/2012.01340)
on a uniform Cartesian mesh. Face gravity here is a potential difference;
the existing FMM acceleration continues to drive momentum.

The equivalent finite-volume expression for cell i is

    delta_E_i = -dt/V_i * sum_faces[A_f * F_out * (phibar_f - phibar_i)].

At a 2:1 interface, use each fine subface's flux and area on both sides and
one shared face potential, interpolated in proportion to the two cell widths.
Thus the coarse energy work uses the same effective flux as mass reflux.
Periodic faces pair normally. At physical boundaries, extrapolate the face
potential using the endpoint-averaged normal acceleration and include the
outward transport `dt*A_f*F_out*phibar_f` in the energy ledger.

## Why the global-step kick energy is replaced

Each momentum kick temporarily adds its kinetic-energy change, preserving
thermal energy during that substep. This is needed for a consistent pressure
in the subsequent hydro predictor. The retained self-gravity contribution is

    delta_E_kick,self = dt * g_self dot (p_before + impulse_total/2).

After both kicks, the final gas energy is

    E_final = E_after_kicks - sum(delta_E_kick,self) + delta_E_mass_flux.

This avoids counting gravitational work twice. Neither momentum kick is
undone. Work by an imposed external acceleration remains in the energy update.
The correction changes only gas energy and the final entropy synchronization.

The refined Mullen path adds provisional flux work directly and replaces its
recorded value at shell closure. Its momentum corrections do not add a separate
kinetic-work energy source. Naive mode retains kinetic-work energy and does not
have the conservative source identity of the Mullen update.

## Conservative regridding

With `conserveRegridEnergy=on`, construct the scalar

    Q = Egas + rho*phi/2

on the old synchronized mesh. Restrict Q by volume averaging and prolong Q
with conservative minmod slopes. Once density has been transferred and gravity
has been solved on the destination mesh, recover

    Egas_new = Q_transferred - rho_new*phi_new/2.

Q is held in a transfer register, separate from the gas-energy field, so a
negative binding-plus-gas energy is never passed to an EOS as thermal/gas energy.
Mesh selection and its shadow estimates continue to use actual gas states.
The sum of Q and the mass are preserved through both refinement and coarsening;
the correction is local to the remapped cells, not a global energy redistribution.
Entropy remains in hydro and is synchronized after recovery. Dual energy can
supply pressure if the recovered conservative gas energy is locally inadmissible;
without that fallback an inadmissible state fails rather than clipping energy.

Manual Runtime callers must have synchronized old gravity before a changed
energy-conserving regrid. After `regrid()` returns true, call `solveGravity()` (or
`setGravity()` for prescribed test fields) before another timestep or regrid.
The next successful field publication recovers gas energy once and rebuilds the
shadow hierarchy from the recovered leaves. An unchanged mesh does not apply
any correction. The ordinary simulation driver performs this sequence.

With the option off, gas energy itself is conservatively transferred, reproducing
the former remap. The resulting change in binding energy remains visible as a
jump in total energy.

## Adaptive gravity reciprocity

The previous adaptive tree walk split the target first when two equal-size nodes
failed the opening test. Reversing source and target could then accept different
child/parent groups, so the two directions used different approximate kernels.
The walk now descends both equal-sized nonleaf nodes together, and otherwise
descends the larger nonleaf node (or the available nonleaf when its partner is a
leaf). This is traversal of the existing gravity tree, not AMR mesh refinement.
It makes the interaction partition symmetric under exchange of source and target.
No multipole order increase or global energy correction is needed for this fix.

Force balance additionally requires equal source and target force degrees.
The scalar potential retains order p while an auxiliary force local retains
degree p+1, so its differentiated target terms match source degree p. The
[kernel proof](parallel-fmm.md#scalar-reciprocity-and-mutual-force) explains why
this requirement is separate from scalar reciprocity.

The independent bilinear regression uses unrelated densities on a mixed-level
mesh at orders 2 and 5; the potential must satisfy
`sum(m_a*phi_b) = sum(m_b*phi_a)` to roundoff. Isolated and periodic mass-flux
energy regressions also require roundoff conservation on fixed adaptive meshes.

## Conservation diagnostics and limits

On a fixed mesh, the coupling balances gas energy against `sum(rho*phi*V/2)`
to roundoff when the potential operator is reciprocal under volume weighting.
For a general approximate potential solve, its independently measurable defect
over one globally synchronized timestep is

    R = sum[V*(rho^n*phi^(n+1) - rho^(n+1)*phi^n)/2].

Tests separate this FMM reciprocity defect from the work discretization error;
no global energy redistribution is used to hide it. The corrected adaptive walk
reduces R to roundoff in the tested cases. A higher multipole order improves
field accuracy but is not a substitute for reciprocity. External imposed forces
can change the self-gravity energy invariant.

For time refinement the corresponding identity uses the sum of this contraction
for every shell interval, with that shell's masked endpoint potentials. The
root endpoint contraction alone does not generally equal the summed defect of
nonreciprocal partial solves. The distinction and boundary terms are derived in
[the coupling document](gravity-time-coupling-derivation.md).

`conservation.csv` now includes cumulative `gravity_reciprocity_defect_erg`,
`gravity_regrid_energy_change_erg`, and `gravity_energy_budget_residual_scaled`.
The last is the boundary-corrected total-energy drift minus these two measured
contributions, divided by the usual energy norm. These are explanatory columns;
the physical `gas_gravity_energy_drift_scaled` column still includes all errors.
The current reciprocity column measures the root endpoint contraction. In a
time-refined run with a nonreciprocal solver, it does not isolate all partial
solve defects; the physical total-energy drift remains the direct budget check.
For naive mode the residual also contains its nonconservative timestep work.
Global diagnostic integrals use compensated summation to reduce accumulation
roundoff; neither the cell solution nor the physical drift is forced to zero.

In dilute atmospheres, conservative work can make E or E-K nonpositive even
when the entropy auxiliary gives valid pressure. Dual energy then retains
finite conservative E and uses the auxiliary for thermodynamics; E is not
clipped. This may create a local discrepancy between conservative and
thermodynamic energies, which remains visible in diagnostics. With dual
energy disabled, an inadmissible result still fails explicitly.
