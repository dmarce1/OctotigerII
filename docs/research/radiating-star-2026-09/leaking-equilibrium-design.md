> Archived design, 2026-09-27. The active benchmark now uses an opaque, source-free star with physical light speed. This is the earlier leaking/heated design; see README.md for the final research status, which supersedes the provisional status below.

# Rotating radiative star: design equations

Status: the structural EOS, rotating diffusion-interior SCF solver, and
nonrotating full-M1 control are implemented and tested. The rotating full-M1
reference is still under construction; no stable rotating result is claimed.
The target is a rotating configuration with balanced gravity, pressure, and
radiative heating, and stable stratification. Global constant entropy is no
longer required. The user explicitly excluded escaping photon angular momentum
as a concern on 2026-09-27: no compensating external torque is to be added, and
photon spin loss is not a blocker for this benchmark. The radiation momentum
equations and gas coupling remain intact. Their physical spin drift is allowed
and should be reported separately from numerical drift; its size has not yet
been measured for a completed model. Thermal and meridional force balance
remain requirements. This reference is not an exactly stationary solution of
every component of the full evolution equations when photon torque is nonzero.

## Thermodynamics and a spherical reference

For a monatomic ideal gas of fixed mean molecular weight, write
`Rgas = kB/(mu*mu_atomic)` and

    pg = rho Rgas T,       ug = (3/2) rho Rgas T.

In the optically thick LTE limit only,

    E0 = a T^4,            pr = E0/3,       P = pg + pr.

The evolution retains separate gas and radiation energies and the full M1
pressure tensor. A structural polytropic exponent must not replace the gas
adiabatic index, 5/3, or impose LTE in the transport atmosphere.

A useful nonrotating diffusion reference is the Eddington constant-beta model:

    beta = pg/P,
    P = K rho^(4/3),
    K = [3(1-beta)/a]^(1/3) [Rgas/beta]^(4/3),
    j = epsilon0 rho,
    epsilon0 = 4 pi c G (1-beta)/kappa,
    L(r) = epsilon0 m(r).

For 0 < beta < 1 its temperature gradient is `nabla = 1/4`, whereas

    nabla_ad = (8-6 beta)/(32-24 beta-3 beta^2),
    Gamma1 = (32-24 beta-3 beta^2)/(24-21 beta).

Thus it is convectively stable in this approximation and has `Gamma1 > 4/3`.
Those statements do not establish thermal or rotating stability. In particular,
constant opacity and constant heating per mass admit a thermally neutral
homology in the idealized spherical model.

This reference cannot simply be spun up. For uniform rotation let

    Psi = Phi - Omega^2 R^2/2,
    A = c(1-beta)/kappa.

The diffusion flux and the required volumetric heating are then

    F0 = A grad Psi,
    j_required = A (4 pi G rho - 2 Omega^2).

The centrifugal contribution makes `j_required` negative in sufficiently
tenuous layers. Heating proportional to density alone does not maintain this
rotating model. Relaxing constant entropy does not remove this constraint.

An alternative diffusion seed permits a nonconstant gas-pressure fraction:

    P = K rho^gamma_s,       gamma_s = 1 + 1/n,       3 < n < 5,
    rho Rgas T + a T^4/3 = K rho^gamma_s.

The temperature equation has a unique positive root for positive density.
For example, `n=3.5` gives `gamma_s=9/7 < 4/3 <= Gamma1`. It is therefore
stably stratified against local adiabatic convection in the LTE interior.
Uniform rotation also avoids differential shear in this seed. These are
limited stability statements, not a thermal or global stability proof.

Define the barotropic integral `H=integral dP/rho`; because entropy varies,
this is not the thermodynamic specific enthalpy. Mechanical balance gives

    H = (n+1) K rho^(1/n),       H + Psi = constant.

At constant opacity, set

    beta = pg/P,
    b = dpr/dP = 4(1-beta)(gamma_s-beta)/[gamma_s(4-3 beta)],
    F0 = (c/kappa) b grad Psi,
    j_required = (c/kappa) [b(4 pi G rho-2 Omega^2)
                           - (db/dH) |grad H|^2].

Here `db/dH < 0`, so the last term is positive and can offset the centrifugal
term. Positivity must still be checked throughout the deformed star. The
source is inferred from the structure, has a distributed outer tail, and
need not be a compact central heater. For `3<n<5`, extrapolation to the formal
zero-pressure surface makes it singular: `j_required` scales as
`H^((n-7)/4)`. This is another reason to match at finite optical depth to a
transport atmosphere rather than extend the diffusion polytrope to vacuum.
This construction is a candidate interior seed, not a full-M1 equilibrium.

## Steady equations of the actual transport model

Use inertial radiation energy density `E`, physical flux `F`, radiation pressure
tensor `Pr`, and inertial gas velocity `v`. Let `chi = rho kappa`; `Q` and `G`
are the gas energy and momentum gained per unit volume and time:

    Q = chi [c(E-a T^4) - v dot F/c],
    G = chi/c [F - Pr v - a T^4 v].

These are the retained equal-opacity mixed-frame terms documented in
[radiation-coupling.md](../../radiation-coupling.md). They are not a fully relativistic
model. Use the same M1 closure as the evolution code:

    f = |F|/(c E),
    D = (3+4 f^2)/(5+2 sqrt(4-3 f^2)),
    Pr = E [(1-D) I/2 + (3D-1) n n/2],       n = F/|F|.

At `F=0`, use `Pr=E I/3` without evaluating `n`.

Let `j` be externally supplied photon energy per volume and time, and `Sgamma`
the corresponding photon momentum source. The proposed benchmark uses a
prescribed photon source isotropic in the laboratory frame, so `Sgamma=0`.
A source isotropic in a moving source frame would instead have
`Sgamma=j vsource/c^2` to first order and would supply angular momentum. No such
angular-momentum supply or external mechanical force is required here.

For an axisymmetric reference with purely azimuthal velocity
`v=R Omega(R,z) e_phi`, let subscript `p` denote meridional components. The
required initial force and thermal balances are

    laplacian Phi = 4 pi G rho,
    grad_p pg = -rho grad_p Phi + rho R Omega^2 e_R + G_p,
    div Pr = Sgamma - G,
    div F = j - Q,
    Q - v dot G = 0.

The azimuthal gas equation is evolved rather than forced to vanish:

    rho d(v_phi)/dt = G_phi,
    dK/dt = v_phi G_phi.

Consequently stationary internal energy requires `Q-v dot G=0`, while total
gas energy can change by the rotational work `Q=v dot G`. For photon heating
only, the reference source is

    j_star = div F_star + Q_star
           = div F_star + v_star dot G_star.

This reduces to the interior diffusion heating formula at leading order.
Prescribe the reference source during evolution; do not recompute it to cancel
instantaneous radiation losses or force the evolved fields back to equilibrium.
The continuity equation holds automatically for axisymmetric pure rotation.
If meridional circulation develops, the evolution retains its continuity,
inertia, and energy advection.

The opacity, source prescription, rotation or flow selection, and boundary
conditions must close this system. Merely counting equations is insufficient:
allowing arbitrary heating retains structural freedom in the spherical and
diffusion limits. Conversely, prescribing arbitrary entropy, rotation, and
heating profiles together can overconstrain it. No full-M1 solution or
existence result is claimed here.

## Angular momentum: retained physics, no compensating torque

For a closed material boundary, isolated self-gravity, and full light speed,

    dJtotal/dt = Jdot_source - Jdot_out,
    Jdot_source = integral [r cross Sgamma]_z dV,
    Jdot_out = surface_integral R (Pr n)_phi dA.

Ordinary outward emission from rotating material carries angular momentum.
To first order `Pr_nphi = v_phi F_n/c^2`; for a uniformly bright spherical
surface this gives `Jdot_out = (2/3) Omega Rstar^2 L/c^2`. A point photon source
at the origin supplies no orbital angular momentum. An extended co-rotating
source supplies `integral j Omega R^2/c^2 dV`, generally insufficient for a
compact heater. These facts do not require torque compensation for the selected
benchmark: angular-momentum loss is an explicitly accepted exception to strict
stationarity. Do not remove `G_phi`, suppress azimuthal radiation stress at the
boundary, or reset the spin.

For the reference with balanced radiation energy and gas internal energy,

    Lout = integral j dV - integral Q dV
         = Lsource - dKrotation/dt.

Thus radiation drag can supply a small part of the escaping luminosity from
rotation. The full evolving energy ledger must include gas, radiation, gravity,
the prescribed photon input, and boundary fluxes. Angular-momentum accounting
should distinguish physical outgoing flux from numerical error; exact spin
stationarity is not an acceptance criterion.

With `eta=chat/c`, inject photon sources consistently with the reduced-speed
equations:

    dE/dt + eta div F = eta (j-Q),
    d(F/c^2)/dt + eta div Pr = eta (Sgamma-G).

The radiation reference balances above are unchanged. The exactly paired global
invariants instead contain `E/eta` and `F/(eta c^2)`. Physical and weighted
budgets must both be reported. Initial equilibrium verification should use
full light speed; agreement at reduced light speed requires a separate check.

## Boundary and stability requirements

The interior must connect to a radiation atmosphere with an outgoing outer
condition, isolated gravity, regular axis/center, and the selected material
boundary condition. Extending an LTE polytrope to zero temperature while
carrying finite luminosity violates `|F| <= c E`. An arbitrary density floor is
also not a stationary atmosphere. Neither is an acceptable surface closure.

Evidence required before calling the benchmark successful:

- Vanishing continuum residuals of meridional gas force, radiation momentum,
  thermal balance, and Poisson gravity under refinement of the reference solve.
  Report the physical azimuthal gas acceleration instead of cancelling it.
- Explicit energy and angular-momentum input/output ledgers, including photon
  momentum. No compensating mechanical torque or its work is added.
- Mesh/time convergence of equilibrium drift using OctoII's actual transport,
  coupling, gravity, atmosphere, and rotating-frame treatment. A continuum
  equilibrium is not automatically an exactly preserved discrete solution.
- Stable response to resolved perturbations, with convection, radial modes,
  thermal modes, and rotational/shear modes distinguished. Outward-increasing
  entropy alone does not guarantee stability with radiative diffusion and
  differential rotation.
- Separate optical-depth/reduced-speed checks, and recovery of the diffusion
  flux where applicable. A long thermal time alone does not establish heating
  balance. Photon spin-down is an accepted physical effect, to be distinguished
  from numerical loss.

## References and scope

- [Skinner & Ostriker (2013)](https://arxiv.org/abs/1306.0010): moment equations,
  mixed-frame coupling, and reduced-speed approximation. The balances above
  are derived from the retained equations; the paper is not a reference
  solution for this rotating benchmark.
- [Espinosa Lara & Rieutord (2013)](https://arxiv.org/abs/1212.0778): simultaneous
  two-dimensional stellar structure, thermal balance, and flow construction.
- [Caleo, Balbus & Potter (2015)](https://arxiv.org/abs/1502.06479): differential
  rotation compatible with radiative equilibrium.
- [Caleo, Balbus & Tognelli (2016)](https://arxiv.org/abs/1605.03750): rotational
  stability with thermal diffusion and the importance of viscosity.
- [Princeton Eddington-model notes](https://www.astro.princeton.edu/~burrows/classes/403/eddington.pdf):
  spherical constant-beta structural reference.

## Implemented spherical control and opacity assumption

The `radiating-sphere` problem implements a full-M1 nonrotating control, with
the source and opacity supplied through the same spatial material callback as
the evolving radiation solver. See its
[problem notes](../../../bin/science/RadiatingSphere/README.md) and
[example configuration](../../../examples/radiating-sphere.ini).

The chosen finite-surface construction has a transparent gas envelope. During
reference generation it uses

    kappa_star(r) = kappa0 max(0,1-rho_cut/rho_ref(r))^2.

The resulting opacity is frozen as a function of position in the evolution,
and is not recomputed from the evolving density. This is an artificial opacity
law for a controlled benchmark, not a realistic stellar opacity prescription.
At the transparent transition, the gas continues onto a polytrope matched in
pressure and pressure derivative; the radiation follows the outgoing vacuum
M1 branch independently of gas temperature. The heater is the divergence of
the continuous reference flux and is nonnegative, with a distributed outer
core tail. It is not a compact central emitter. Neither opacity nor heating
feeds back on the evolving solution to maintain equilibrium.

For `n=3.5`, central gas fraction `0.8`, `kappa0 rho_c alpha=240`, and
`rho_cut/rho_c=0.016`, independent differential checks recover hydrostatic
balance, M1 momentum balance, gravity, and luminosity closure. Reference
integration tolerance and the central series starting radius are varied
independently. The default material surface is at `6.3941644 alpha`, outside
the opacity transition at `3.5708999 alpha`; central optical depth is about
`286.41`.

The initial production evolution comparison uses full light speed and the
same physical interval of 10.8652 s on uniform meshes:

| Cells per axis | Mass-weighted density L1 drift | RMS Mach | Radiation energy change | Source/boundary-corrected energy error |
| --- | ---: | ---: | ---: | ---: |
| 8 | 0.0029424 | 0.0027012 | +0.096564 | -1.39e-16 |
| 16 | 0.0011921 | 0.0018710 | -0.015870 | -1.18e-15 |

These coarse meshes and this interval, only about `0.0028/sqrt(G rho_c)`,
test plumbing, conservation, and improving discretization error. They do not
establish long-term stability or convergence of all stellar observables. The
numerical gas floor and finite outflow boundary remain approximations to vary.

The full rotating transport boundary-value problem, a stable solution, and
its OctoII initializer remain to be completed. Small algebraic residuals in
experimental rotating reference solves have not been accepted as sufficient:
the angular velocity must be real everywhere in the material, heating
nonnegative, and continuum residuals and fields must converge with resolution.
