# Radiation-matter coupling and diagnostics

Radiation stores physical inertial energy density `Er` and flux `F`, with
`|F| <= c Er` using physical light speed `c`. The configured reduced speed
`chat = radiation.lightSpeedRatio * c` changes the transport and exchange time
scales, not the units of the stored flux. Gas energy includes inertial kinetic
energy. A rotating mesh changes integration surfaces and vector components;
it does not replace the velocity in radiation-matter exchange by velocity
relative to the mesh.

The default `radiation.opacityModel=constant` uses `radiation.opacity` as true
gray absorption in cm²/g. Optional `radiation.scatteringOpacity` adds elastic
scattering to momentum coupling and flux extinction, without adding thermal
emission or absorption. The gas temperature uses the selected dual-energy
thermal state and the configured mean molecular weight. The separate
[opaque rotating-star construction](radiating-star-design.md) supplies a
structural EOS and initial mechanical balance; coupling alone does not impose
radiative equilibrium.

## Governing exchange equations

Write `eta = chat/c`, gas velocity `v = momentum/rho`, and let `P` be the full
M1 radiation-pressure tensor. With Planck absorption `kappa_P` and flux
extinction `kappa_F`, the implemented first-order mixed-frame source rates are

    Q = rho c kappa_P (Er - a T^4) + rho (kappa_F - 2 kappa_P) v dot F/c
    G = rho/c [kappa_F F + kappa_P (Er - a T^4) v - kappa_F (Er v + P v)].

The local equations are

    dEgas/dt = Q,            d(momentum)/dt = G,
    dEr/dt   = -eta Q,       dF/dt          = -eta c^2 G.

These follow the first-order terms of equations (5) in
[Skinner and Ostriker (2013)](https://arxiv.org/html/1306.0010), approximating
their energy-weighted absorption mean by the Planck mean and their flux-weighted
extinction mean by the Rosseland total mean in the diffusion limit. Setting
`kappa_P=kappa_F=kappa` recovers the former equal-opacity source exactly.
Their additional terms quadratic in `v/c`
are not included, and this is not a relativistic hydrodynamics model.
Their reduced equation further replaces `a T^4` by `Er` in the
velocity-dependent momentum source. The joint implicit solve already has
the gas temperature, so the implementation retains the temperature term.
Keeping the full M1 tensor avoids replacing anisotropic radiation pressure by
`Er I/3` outside diffusion. Large material speeds require revisiting the
physical approximation, not merely reducing the numerical timestep.

## Analytic ionized-gas opacity

Set `radiation.opacityModel=ionized-gas` to evaluate separate gray means at
the current density and gas temperature. For hydrogen mass fraction `X` and
metal mass fraction `Z`, the provisional fully ionized free-free approximation
is

    kappa_R,ff = 3.68e22 (1+X)(1-Z) rho T^(-7/2) cm^2/g
    kappa_P,ff = 37 kappa_R,ff
    kappa_es   = 0.4 Ye cm^2/g
    kappa_P    = kappa_P,ff
    kappa_R,total = kappa_R,ff + kappa_es.

`rho` is in g/cm³ and `T` in K. `Ye` is the evolved electron composition
moment divided by gas density when available; otherwise it is `(1+X)/2`.
`radiation.hydrogenFraction` and `radiation.metalFraction` are spatially
uniform input parameters, defaulting to 0.7 and 0.02. The free-free prefactor
is a Kramers approximation; the approximate factor of 37 between Planck and
Rosseland free-free means is used in [gray atmosphere modeling](https://www.aanda.org/articles/aa/full_html/2020/01/aa35033-19/aa35033-19.html).
These values are **not calibrated stellar opacities**. This model omits
bound-free and line absorption, ionization transitions, dust, Compton energy
exchange, and electron conduction. Its free-free expression assumes hot,
nondegenerate, nearly fully ionized gas and fixed bulk `X,Z`; the evolved
electron moment only adjusts scattering. Adding the two Rosseland-scale
terms is another approximation: the Rosseland mean of their frequency-wise
sum is not generally the sum of their separate means. Do not apply it to a
degenerate white-dwarf interior or a cool photosphere.

The nonlinear source solve evaluates these coefficients at every trial gas
temperature, including their temperature derivatives. The opacity-aware
face flux and optical-depth diagnostics use `kappa_R,total`; radiation
emission and thermalization use `kappa_P`. The source flux mean tends to
`kappa_R,total` for nearly isotropic radiation and rises continuously to at
least `kappa_P` as `|F|/(c Er)` approaches one. This last interpolation
protects the streaming cone when the Planck mean exceeds the Rosseland mean;
it is a gray approximation, not a measured spectral flux mean. For radiation
whose spectrum differs strongly from local thermal emission, multigroup
transport is needed to obtain the actual energy and flux means. The analytic
model cannot be combined with `radiation.opacity`,
`radiation.scatteringOpacity`, or a problem's fixed spatial material profile.
Existing constant-opacity configurations retain their previous behavior
when `scatteringOpacity=0`.

The retained temperature term also matters at the streaming boundary in the
equal-opacity limit. For
`F=c Er n` and `beta_n=v dot n/c`, the source satisfies

    d(Er - n dot F/c)/dt = chat chi a T^4 (1-beta_n) >= 0

for subluminal material velocity. Replacing `a T^4` by `Er` in the momentum
source instead gives `chat chi (a T^4-beta_n Er)`, which can point outside
the radiation cone for a cold absorbing gas moving along the beam. This
is a property of the retained equations, separate from the numerical
admissibility checks on implicit stages and transport fluxes.

The source solver uses a two-stage L-stable singly diagonally implicit
Runge–Kutta method with diagonal coefficient `1-1/sqrt(2)`. A safeguarded
local nonlinear solve enforces positive selected gas thermal energy and
`|F| <= c Er` to scaled floating-point roundoff. It uses analytic derivatives
of the gas feedback and M1 tensor; its residual tolerance scales with the
radiation equation terms so weak radiation is not discarded relative to the
gas energy. Failed local solves subdivide the interval; a failure after
the bounded retry budget raises an error without modifying the caller's
states. Clipping the radiation flux would silently lose paired momentum
exchange and is not used to repair source iterates.

For a large incoming transport increment, the Newton solver also tries the
transported radiation state as its starting guess, when admissible and when
its stage residual is smaller. Otherwise it retains the existing guess. This
avoids spuriously assigning incoming radiation to trial gas heat in a dilute
surface cell, while preserving accurate already-balanced stiff states. Only
the nonlinear initial iterate changes; the stage equations and paired accepted
exchange are unchanged.

## Prescribed material and photon sources

A problem can provide a `RADIATION_MATERIAL` manifest hook that constructs an
immutable callback for opacity `kappa(x,t)` and prescribed photon power
`H(x,t)`, in `cm^2/g` and `erg/(cm^3 s)` respectively. Coordinates are the
physical grid coordinates used by the initializer; a rotating axisymmetric
model uses the same reference profile. The callback depends on the problem
configuration, position, and time, not the evolving temperature or radiation
state. Both returned quantities must be finite and nonnegative. Other
problems retain the configured opacity law and zero photon heating.

The external photon source changes the radiation energy equation to

    dEr/dt = transport - eta Q + eta H.

It adds no direct gas-energy or momentum source. In particular, prescribed
photons continue to be injected where absorption is exactly zero. The
callback's opacity supplies absorption for constant-material problems;
scattering, when configured, also contributes to face fluxes and optical
depth diagnostics. The radiation flux cone still uses physical `c`.

Photon heating is included in the known driving of the source-aware midpoint
and final implicit solves. The accepted final source uses midpoint quadrature
over its physical interval. Predictor, coarse halo forecast, and shadow
evaluations do not count as physical injections. Only completed physical leaf
updates contribute to the source ledger; a rejected interval restores its
previous value along with the evolved states. Regridding preserves that
cumulative ledger.

`Runtime::radiationSourceEnergy()` and the CSV column
`radiation_source_energy_erg` record the actual prescribed energy added to
`Er`, namely the discrete counterpart of `integral eta H dV dt`.
`rsla_source_energy_erg` divides it by `eta`, matching the source in the
weighted combined-energy budget. Corrected radiation and combined-energy
columns subtract the appropriate prescribed source separately from boundary
inflow/outflow. Radiation energy alone can still change through exchange with
gas; the combined weighted budget is the conservation check. Silo reports the
local fields `radiationOpacity` (Planck absorption),
`radiationRosselandAbsorption`, `radiationScatteringOpacity`,
`radiationFluxOpacity` (their sum), and `photonHeating`, with
`photonHeating` showing physical `H` before reduced-speed scaling.

The [`photon-source` verification problem](../problem/radiation_tests/PhotonSource/README.md)
uses uniform, fixed `H=1 erg/(cm^3 s)`. Its transparent solution is analytic;
positive-opacity cases test the accepted source ledger across coupling, AMR
subcycling, regridding, and rollback. These checks validate the source
infrastructure, not a stellar equilibrium configuration.

## Pressure balance and thick-cell transport

For static gas, the flux equation is

    dF/dt = -c chat div P - chat chi F.

The two terms balance in diffusion. A pressure-gradient kick followed by
independent exponential damping gives
`(F_old-c chat dt div P) exp(-chat chi dt)`, which incorrectly vanishes as
the optical timestep grows. In contrast, explicit pressure forcing with
backward-Euler damping gives

    F_new = (F_old-c chat dt div P*) / (1+chat chi dt)
          -> -c/chi div P*.

Thus an implicit spatial pressure solve is not required. The implemented
higher-order source solve accepts the known transport increment as constant
forcing throughout its stages. An initial numerical-flux probe supplies the
transport drive for the coupled half-step source solve. The resulting midpoint
cell averages are reconstructed for the accepted numerical face fluxes, which
then drive the final full-interval source solve. This balances pressure transport
against stiff relaxation without a global implicit solve. The timestep still
obeys the explicit radiation transport CFL condition.

The runtime retains both initial and midpoint halo exchanges, each using two
ghost layers for reconstruction. Coupled AMR intervals remain capped at twice
the finest leaf's CFL limit because every leaf supplies the source-aware
midpoint donor state; active substeps also obey their own CFL bounds. The
self-contained patch calculation used for shadows and diffusion checks receives
four initial ghost layers to compute its enlarged midpoint stencil locally.

The runtime cache removes repeated reads within these stages. For coupled
transport it retains only the first layer of initial gas and radiation ghosts
after spatial interpolation and physical boundary treatment; those are the
original states needed again by the corrector. This compact cache does not
reduce the reconstruction stencil. Gravity without radiation retains the full
raw initial halo, including prolongation donors, because its source-rate
midpoint must be formed before spatial interpolation. Midpoint and material
halo buffers persist in each worker workspace and share the fetched gas values
between hydro and radiation calculations.

Each locality reserves a pool of initial-halo buffers for its owned blocks plus
its worker count. Additional buffers acquired for stolen probes return to that
pool for later reuse. A probe may be stolen; the associated corrector executes
on the locality holding its cache. Cache identity includes the input bank,
physical time, reference interval, source-rate field, and coarse donor times.
Replacement probes and rollback invalidate retained data. An opening gas kick
requires a fresh initial gas halo while the unchanged radiation halo remains
reusable.

For a uniform unforced coupled step, the request count changes from
`5 G + 3 R` to `2 G + 2 R` scalar values per donor cell, where `G` and `R` are the
stored gas and radiation component counts. An opening kick requires `3 G + 2 R`.
The 1D and 3D serial fixtures verify these counts, second-order temporal
convergence, conservation, regridding, and rollback. The counts do not measure
whole-application network or wall-time improvements. See
[halo communication accounting](halo-communication-accounting.md) for the
measurement scope.

Time integration alone does not remove excess HLL diffusion in an opaque
cell. With face extinction `chi_f=(chi_L+chi_R)/2`, `tau=chi_f dx`, and
`b=1/(1+tau^2)`, the radiation energy face flux is corrected to

    J_E = b J_thin + (1-b) J_material
          - chat tau/[3(1+tau^2)] (Er_R-Er_L).

`J_material` contains the radiation flux on the moving material source
equilibrium and the ALE mesh advection, upwinded by their net energy
transport speed. Neither is attenuated with opacity. The pressure flux
likewise approaches the source-equilibrium M1 tensor. Consequently static
opaque cells approach `J_E=-chat/(3 chi_f) grad Er`, while the physical
radiation flux approaches `F=-c/(3 chi) grad Er`. The extinction average
represents two half-cell optical resistances in series.

This interpolation is an application-specific asymptotic matching choice,
not a verbatim implementation of a published flux. The general need to
modify numerical diffusion and balance the flux source is established in
[Bloch et al. (2021)](https://doi.org/10.1051/0004-6361/202038579).
At zero opacity the original flux is recovered exactly; its transparent
limit correction is quadratic in cell width for resolved smooth states.
At fixed cell width and `chi=chi_tilde/epsilon`, `b=O(epsilon^2)` suppresses
the HLL energy viscosity below the `O(epsilon)` physical diffusive flux.
The gradient uses midpoint cell averages, rather than the much smaller jump
between two reconstructed face states. The correction applies to
both high- and low-order radiation fluxes before the shared realizability
limiter. If the corrected low-order update violates the radiation cone,
the limiter can fall back to the original dissipative flux. The thick-cell
coefficient is therefore not promised in cells where that fallback is active.
The same qualification applies whenever the high-order shared-face limiter
changes the corrected flux, even if it does not require the original fallback.

The moving source equilibrium satisfies both
`aT^4=Er-v dot F/c^2` and `F=(aT^4 I+P) v`. With `beta=|v|/c`, the
Lorentz-isotropic M1 manifold at fixed inertial energy is

    F = 4 Er v/(3+beta^2)
    aT^4 = 3(1-beta^2) Er/(3+beta^2).

This face construction requires subluminal material speed. It is algebraically
the Lorentz transform of isotropic thermal radiation; the retained evolution
equations are still a nonrelativistic approximation. Its subluminal domain
does not establish physical accuracy at large speeds. At stellar velocities
the flux reduces to `F=(4/3) Er v + O((v/c)^3 c Er)`.

The associated M1 tensor satisfies both retained source balances exactly.
Its static limit is `P=Er I/3`.
This also removes a transient source-free pressure predictor's anisotropy
from the opaque numerical pressure flux. The full method still requires
source-aware midpoint states and pressure forcing in the implicit stages.
The measured scope and reproduction commands are recorded in
[radiation-validation.md](radiation-validation.md).

## Conservation

For a local exchange step, density and gravitational binding energy are
unchanged. With `w = c/chat`, the paired increments satisfy

    delta Egas = -w * delta Er
    delta (rho v) = -w * delta F / c^2.

Consequently, the source conserves `Egas + w Er` and `rho v + w F/c^2` to
floating-point rounding. At full light speed these are physical energy and
momentum. At reduced light speed they are the invariants of the approximation;
the physical sums can change. The diagnostic output reports both, with
distinct names rather than interpreting the weighted balance as physical
conservation. Self gravity adds `integral rho phi/2 dV` and the existing
potential boundary transport to the energy budget. See
[conservation.md](conservation.md).

Momentum exchange also changes kinetic energy. The thermal increment is
`delta u = delta Egas - delta K`, with the stable fixed-density expression

    delta K = delta p dot (p_old + delta p/2) / rho.

The entropy auxiliary must receive this nonadiabatic increment even when the
usual total-energy synchronization threshold is not met. It must not be used
to overwrite independently conserved gas energy. Gravity's provisional work
and mass-flux energy corrections remain separate from the local exchange.

An imposed acceleration must also do work on momentum transferred from
radiation during the step. Its kinetic kicks therefore include
`dt * acceleration dot delta_m_source_mid`, where `delta_m_source_mid` is
the source-only midpoint impulse, with the predicted radiation transport
increment removed. On a globally kicked interval this correction is applied
once in the closing kick, after the numerical midpoint prediction has been computed. The same correction applies to the optional naive
self-gravity work. Conservative Mullen self-gravity already uses the
source-aware midpoint mass flux and does not receive a duplicate correction.
The dual-energy auxiliary receives the additional raw-drive work before the
final coupled solve accounts for the full kinetic-energy change.

## Nonfatal reduced-speed diagnostics

Silo includes three dimensionless cell fields whenever gas and radiation are
both enabled:

| Field | Definition |
| --- | --- |
| `cellOpticalDepth` | `tau_cell = rho * kappa_R,total * dx` |
| `radiationTrappingParameter` | `|v| * tau_L / c` |
| `rslaCriterion` | `(|v| + a_bound) * max(1, tau_L) / chat` |

Here `v` is the inertial gas velocity and

    tau_L = rho * kappa_R,total * L
    a_bound = sqrt((gamma * Pgas + 4 Er/9) / rho).

The acoustic expression includes a radiation-pressure contribution and is a
conservative timescale proxy, not a computed characteristic speed of the full
coupled system. `radiation.diagnosticLength` sets the fixed physical length
`L` in cm. Its default, zero, selects the domain side length. Unlike `dx`, this
length stays fixed under mesh refinement, so refining a cell does not by
itself improve the reported RSLA criterion. The default domain scale is a
conservative choice for a uniform medium, but local `rho * kappa_R,total * L` is not
a rigorous column optical-depth bound in a heterogeneous flow. Choose a
physical scale appropriate to the phenomenon and inspect the field, not only
the maximum.

The trapping parameter compares physical radiation diffusion and material
advection timescales. Values near or above one indicate a regime where
advection and radiation work can matter; they do not constitute a numerical
failure. A reduced-speed calculation should have `rslaCriterion` well below
one on the relevant scales. A value near or above one flags a possible
timescale conflict and motivates a comparison with a larger `chat`. It is a
heuristic diagnostic, not a measured error or a universal validity threshold.
In particular, at `chat=c`, a large value can describe physical radiation
trapping and is not a failure of the physical-light-speed solver.

The run continues for all finite values of these diagnostics. No automatic
parameter changes, clipping, or timestep rejection are attached to the
criterion. `conservation.csv` reports the respective global maxima as
`maximum_cell_optical_depth`, `maximum_radiation_trapping_parameter`, and
`maximum_rsla_criterion` on every completed timestep, including when Silo
output is disabled.

## Accuracy checks are separate from balance checks

A conservative source update alone does not establish correct optically thick
transport. In the static diffusion limit the physical flux must recover

    F = -c / (3 rho kappa_R,total) * grad Er.

Reduced-speed transport changes the energy evolution timescale, while this
steady physical-flux relation should retain physical `c`. Test that limit
independently of the paired exchange identities, along with transparent
transport, LTE equilibrium, moving-medium work, and temporal convergence.
Passing a global conservation test alone cannot establish the accuracy of a
diffusion solution or justify a reduced light speed.

## Helmholtz material thermodynamics

With `hydro.eos=helmholtz`, hydro and the local exchange solver use the gas-only
EOS, excluding equilibrium photons. The nonlinear source evaluates the EOS
temperature and heat capacity with the advected composition. Emission derivatives
use `1/(rho*cv)` rather than the ideal-gas proportionality. See
[Helmholtz hydro](helmholtz-hydro.md) for inversion, entropy transport and the
explicit floor source that must be subtracted from conservation budgets.
