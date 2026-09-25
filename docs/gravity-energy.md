# Gravity energy from the numerical mass flux

The momentum integrator retains its existing endpoint kicks:

1. Save potential and acceleration at time n, and zero the kick-work register.
2. Apply the first half momentum kick with the field at n.
3. Advance hydro and species, storing the actual numerical density flux; reflux.
4. Solve gravity for the updated density at n+1.
5. Apply the second half momentum kick with the field at n+1.
6. Replace accumulated self-gravity kick work with endpoint-averaged mass-flux
   work, then synchronize hydro's entropy auxiliary.

Potential and acceleration remain synchronized at each endpoint. There is no
midpoint-density gravity solve and no separately evolved dphi/dt field.
`simulation::run` performs the sequence automatically. Manual Runtime clients
bracket the usual kick/advance/solve/kick sequence with `beginGravityEnergy()`
and `finishGravityEnergy(dt)`. Regridding is disallowed inside that bracket.
Calling a kick outside it retains the standalone kinetic-work behavior.

## Discrete energy source

Let F be the final numerical mass flux, oriented from left to right, and let
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

## Why the kick energy is replaced

Each momentum kick temporarily adds its kinetic-energy change, preserving
thermal energy during that substep. This is needed for a consistent pressure
in the subsequent hydro predictor. The retained self-gravity contribution is

    delta_E_kick,self = dt * g_self dot (p_before + impulse_total/2).

After both kicks, the final gas energy is

    E_final = E_after_kicks - sum(delta_E_kick,self) + delta_E_mass_flux.

This avoids counting gravitational work twice. Neither momentum kick is
undone. Work by an imposed external acceleration remains in the energy update.
The correction changes only gas energy and the final entropy synchronization.

## Conservation limits

On a fixed mesh, the coupling balances gas energy against `sum(rho*phi*V/2)`
to roundoff when the potential operator is reciprocal under volume weighting.
For a general approximate potential solve, its independently measurable defect
per timestep is

    R = sum[V*(rho^n*phi^(n+1) - rho^(n+1)*phi^n)/2].

Tests separate this FMM reciprocity defect from the work discretization error;
no global energy redistribution is used to hide it. The current adaptive FMM
can produce a nonzero R. Regridding conserves gas and species integrals but
changes the discrete gravitational binding energy; this patch does not apply
a compensating thermal adjustment at a mesh change. External imposed forces
can also change the self-gravity energy invariant.

In dilute atmospheres, conservative work can make E or E-K nonpositive even
when the entropy auxiliary gives valid pressure. Dual energy then retains
finite conservative E and uses the auxiliary for thermodynamics; E is not
clipped. This may create a local discrepancy between conservative and
thermodynamic energies, which remains visible in diagnostics. With dual
energy disabled, an inadmissible result still fails explicitly.
