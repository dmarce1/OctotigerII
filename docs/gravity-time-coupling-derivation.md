# Gravity time refinement: finite-volume coupling derivation

Derivation and implementation design, 2026-09-25. This document describes the
finite-volume adaptation implemented here. Identities below are exact where
stated; temporal-accuracy claims have explicit smoothness and stage-accuracy
assumptions. Smooth fixed-mesh regressions support second-order convergence
for the tested case; they do not establish accuracy outside those assumptions.

Baseline correction: the gravity-conservation update from
`/home/dmarce1/Downloads/OctotigerII-energy-conservation-source.tar.gz` has now
been merged into this working checkout. It includes the symmetric adaptive
traversal, a direct potential-reciprocity regression, and the independent
`gravity.energyTreatment` and `gravity.conserveRegridEnergy` options. Global
timesteps and conservative regridding balance gas plus binding energy to
roundoff in the tested cases. The multirate scheme below should preserve those
properties; it must not reintroduce the older asymmetric gravity walk.

The result is a separation of responsibilities: hierarchical pair interactions
determine momentum impulses, while the numerical continuity equation determines
gravitational energy work. Moving mass between Eulerian timestep groups requires
no additional momentum source. It does require accounting for both sides of
every mass flux and, if binding energy is split into groups, transfers between
those binding-energy components.

References:

- [Pelupessy, Janes & Portegies Zwart (2012), section 2.1](https://arxiv.org/html/1205.5668):
  recursive slow/fast interaction splitting.
- [GADGET-4 code paper, sections 4.1-4.2](https://wwwmpa.mpa-garching.mpg.de/gadget4/gadget4-code-paper.pdf):
  conventional and hierarchical gravity schedules.
- [Mullen, Hanawa & Gammie (2021), equations 52-59](https://arxiv.org/html/2012.01340):
  gravitational work using the Riemann mass flux and time-centered gravity.

The finite-volume derivation and synchronization corrections below are
an adaptation for this repository, not an algorithm claimed by those papers.

**1. Fix the discrete system during a synchronization interval.**

Let cell i have fixed volume V_i and center x_i. Use extensive variables

    m_i = V_i rho_i,   P_i = V_i (rho v)_i,   E_i = V_i e_i.

E is gas total energy, excluding gravitational binding energy. Only active leaf
cells enter these vectors; covered coarse cells are not additional masses.
Regridding is a separate synchronized operation.

Orient each interior face from cell i to cell j. On an AMR interface, treat each
fine subface as a separate face. Define the accepted mass transfer

    q_f = integral_[T,T+H] A_f F_rho,f(t) dt.                    (1)

Then that face contributes -q_f to m_i and +q_f to m_j. In a subcycled step q_f
is the sum of accepted substep flux integrals. Reflux must replace the coarse
estimate with exactly that sum. When species define density, use their actual
total mass flux, as the existing gravity-energy coupling already does.

For a fixed linear gravity discretization, write

    phi = K m,                W = (1/2) m^T K m.               (2)

The mass form absorbs unequal cell volumes. Reciprocity is K = K^T; in density
variables it is volume-weighted reciprocity. The force operator has a separate
requirement, action/reaction symmetry. A reciprocal potential does not by itself
prove that an approximate acceleration operator conserves momentum.

**2. Split forces according to pair cadence.**

Let b_i be a cell's fixed timestep-bin index for the interval, with larger b
meaning a smaller timestep, H / 2^b. Spatial level and timestep bin are distinct.
Define nested sets A_l = {i: b_i >= l} and diagonal masks D_l. Empty finest sets
have zero mask. The interaction shell at cadence l is

    K_l = D_l K D_l - D_(l+1) K D_(l+1),
    W_l = (1/2) m^T K_l m.                                    (3)

Thus sum_l K_l = K and sum_l W_l = W. Each pair occurs once, in shell
l = min(b_i,b_j): its interaction has the slower member's cadence. These are
slow-slow and slow-fast pairs; fast-fast pairs go into the next nested set.

For isolated direct gravity, define the force on i from j by

    F_ij = G m_i m_j (x_j-x_i) / |x_j-x_i|^3,
    F_ji = -F_ij.                                              (4)

A shell kick of duration tau uses one common mass state and applies

    delta P_i = tau sum_[j in shell(i)] F_ij.                  (5)

The two endpoints of each interaction receive opposite extensive impulses.
Consequently sum_i delta P_i = 0 even though the masses differ from those used
in an earlier kick. There is no dm/dt correction to (5): advected momentum is
already in the conservative hydro momentum flux. Adding another transfer
impulse would count it twice.

The particle HOLD map realizes this impulse using an opening and closing
half-kick. The finite-volume adaptation instead retains physical gas states,
adds recorded provisional source impulses during hydro substeps, and replaces
their accumulated value at shell closure. Its accepted impulse uses the same
paired endpoint quadrature,

    I_ij = (h_l/2) [F_ij(t_a) + F_ij(t_b)].                    (6)

The closing mass state includes accepted interface mass fluxes and reflux.
The same I_ij is then applied with opposite signs. The accepted shell impulse
agrees with HOLD's endpoint quadrature; the intermediate states and hydro map
are not the particle HOLD composition. A frozen opening acceleration may be
used for the provisional impulse only if that impulse is recorded and replaced.
Keeping it multiplied by later target masses as the final accepted impulse
would break the pairing.

For an FMM, the corresponding requirement is mutual, balanced force exchange
for every selected interaction. Independent source/target masks alone do not
establish it. Potential reciprocity and force action/reaction are separate
properties; the corrected archive's potential-reciprocity fix should not by
itself be read as a proof about momentum conservation for newly masked forces.
The current backend additionally keeps an auxiliary force local through degree
p+1. Its gradient has target degree p, matching the source degree p; retaining
only the scalar order-p local would leave target force degree p-1. See the
[scalar-kernel and mutual-force proof](parallel-fmm.md#scalar-reciprocity-and-mutual-force)
for the resulting symmetric force construction and its image-kernel treatment.
Periodic conservative pair kernels admit the same argument; reflecting images
can exchange momentum with the boundary and require a different boundary budget.

For comparison, conventional time refinement integrates the total force on i
at i's own cadence. Its quadrature for i<-j generally differs from that for
j<-i. Equation (4) then need not cancel after time integration, even with an
exact spatial force solver.

**3. Changing group mass prevents a literal particle-Hamiltonian substitution.**

The hydro transport operator changes m_i and couples cells across group
boundaries. Separate slow and fast hydro updates therefore do not commute in
the way particle drifts of disjoint particle sets do. The particle argument
does not prove a symplectic Eulerian hydro integrator.

The energy consequence can be seen using one fast set F. Let phi_F denote the
potential of its masses evaluated everywhere, and chi_F its fixed spatial mask.
Its binding energy and mass derivative are

    W_FF = (1/2) m^T D_F K D_F m,
    dW_FF/dm_i = chi_F(i) phi_F(i).                            (7)

Its physical force is chi_F g_F, whereas a spatial gradient of the derivative
in (7) contains an additional group-boundary term:

    -grad(chi_F phi_F) = chi_F g_F - phi_F grad(chi_F).         (8)

This term describes reassignment of binding energy when mass crosses the mask;
it is not an extra gravitational force on the gas.

In the continuum notation for this bookkeeping identity, with mass current J,
the physical fast-fast gas work satisfies

    dW_FF/dt + Qdot_FF = -integral_boundary(F) phi_F J.n dA.   (9)

The normal is outward from F. The complementary slow-containing binding
energy has the opposite interface transfer. For a mass q moving from S into F,
the transfer is approximately +q phi_F,face into the W_FF budget, with the
opposite sign in the complementary budget. Its time-integrated value must be
shared by the two budgets. Independently sampling it at different timesteps
does not cancel exactly.

For a concrete discrete example, take two cells, K_12 = K_21 = -1, omit self
terms, and let their masses be a and b. Moving q from the slow cell to the fast
cell changes total binding energy by

    delta W = -(a-q)(b+q) + ab = (b-a)q + q^2.                (10)

Both density changes are essential. For a=2, b=3, q=0.1, the change is +0.11;
using only the fast-cell term with the old potential gives -0.2.

This rules out treating the existing energy formula on independently masked
potentials as independent physical gas-heating terms without transfer accounting.

**4. A synchronization energy identity avoids the group-boundary ambiguity.**

Use the full potentials at accepted synchronization endpoints:

    phi_bar_i = [phi_i(T) + phi_i(T+H)] / 2.                  (11)

For an interior face i->j choose one shared face potential. Linear interpolation
between centers gives

    phi_bar_f = (d_j phi_bar_i + d_i phi_bar_j)/(d_i+d_j),     (12)

where d_i and d_j are their normal distances to the face. The conservation
identity only requires a shared face value; consistency determines how it
should be interpolated. Define gravitational gas-energy increments

    Q_i,f = -q_f (phi_bar_f - phi_bar_i),
    Q_j,f = +q_f (phi_bar_f - phi_bar_j).                     (13)

For equal widths each cell receives half of -q_f(phi_bar_j-phi_bar_i).
For unequal widths use the same fine subface and q_f on both sides, and divide
by each cell's volume only when converting back to energy density.

Adding (13) over cells cancels every interior face potential:

    sum_i Q_i = -phi_bar^T delta m                            (14)

for a closed domain. For reciprocal, fixed K, polarization of (2) gives

    W(T+H)-W(T) = (1/2)(m_new+m_old)^T K (m_new-m_old)
                 = phi_bar^T delta m.                        (15)

Therefore

    sum_i Q_i + delta W = 0.                                 (16)

This is independent of the number of hydro substeps, the gravity kick schedule,
and the amount of mass crossing between timestep groups. Every accepted mass
change must be represented by the face ledger. A final instantaneous flux,
or a coarse flux before reflux, is not sufficient.

For outward boundary mass q_f, (13) becomes

    Q_i,f = -q_f(phi_bar_boundary - phi_bar_i).

The extra potential-energy boundary ledger is B_phi = sum_boundary q_f phi_bar_boundary.
Together with the accepted outward hydro energy flux B_gas,

    delta(sum_i E_i + W) + B_gas + B_phi = 0.                 (17)

This describes the discrete boundary convention already used in this repository;
it does not evolve the gravitational field of matter after it leaves the domain.

For arbitrary computed endpoint potentials, the exact algebra instead gives

    R = (1/2)[m_old^T phi_new - m_new^T phi_old],
    delta(sum_i E_i + W) + B_gas + B_phi = R.                 (18)

R vanishes for a fixed reciprocal linear operator. Otherwise it measures the
endpoint defect separately from the work discretization, as the existing tests
already do. The merged `gravity.conserveRegridEnergy=on` path separately
transfers gas plus binding energy through regridding; time refinement should
continue to close its work ledgers before invoking that path.

**5. Physical gas states and provisional sources must be reconciled explicitly.**

The globally synchronized reference path uses kinetic-work kicks. At fixed
mass such a kick preserves internal energy by temporarily adding

    delta E_kin = delta P . (P_before + delta P/2) / m.       (19)

Record actual sequential kick increments; separate kicks have cross terms and
must not each be evaluated from the same original momentum. With simultaneous
self and external impulses, the self-gravity share is

    delta E_self = delta P_self . (P_before+delta P_total/2)/m.

The refined Mullen path instead evolves physical gas momentum and energy
throughout. Its provisional momentum source and its provisional mass-flux
energy source are separate; it does not impose an additional kinetic-energy
increment when applying a momentum correction. The naive option retains
kinetic-work energy updates. In either construction, let S_i be the accumulated
self-gravity energy sources actually added to cell i
during the interval. This includes temporary kick work and any provisional
mass-flux work corrections, with their actual signs. Then the endpoint operation

    E_i <- E_i + Q_i - S_i                                  (20)

leaves exactly the source (13) in the finite-volume update. External work is not
subtracted. S_i records source additions, not a passively advected material
quantity; the hydro energy-flux ledger already records transport of the evolved
gas energy. All refluxed hydro energy fluxes must be included in that ledger.

Equation (20) changes no momentum or density. The auxiliary entropy can then be
synchronized from admissible thermodynamics. Clipping conservative energy or
repairing positivity would be an additional energy source and cannot be hidden
inside the conservation claim.

**6. Close the energy components on their own timesteps.**

The synchronization identity also gives a genuinely multirate work ledger. For
each shell l and each of its intervals I=[a,b], define

    phi_l(t) = K_l m(t),
    phi_bar_l,I = [phi_l(a)+phi_l(b)]/2,
    q_f,I = sum of accepted mass transfers on f during I.     (21)

All cells in A_l have completed their mass updates at a and b. Cells outside
A_l have zero entries in phi_l, so their asynchronous density states are not
needed for this shell solve. Fluxes crossing the edge of A_l are still needed.
Each face must have a single time-integrated flux history shared by its two
neighbors; use that history to obtain q_f,I, not each neighbor's separate
provisional flux. Bin alignment and flux ownership must make these intervals
well defined.

Apply the same linear face interpolation as (12) to phi_bar_l,I and form the
accounting increment

    Q_i,l,I = -sum_[f outward from i] q_f,I
                  (phi_bar_l,I,f - phi_bar_l,I,i).            (22)

This sum includes cells just outside A_l when a face touches A_l. Dropping that
side of a crossing flux would destroy the proof. Each coarse/fine subface is
again a shared transaction, including when it crosses a timestep-group boundary.

For symmetric K_l, exactly the same polarization argument gives

    sum_i Q_i,l,I = -[W_l(b)-W_l(a)].                         (23)

Sum first over all intervals of each shell, then over shells. Every shell's
intermediate binding energies telescope, and (3) gives

    sum_[i,l,I] Q_i,l,I = -[W(T+H)-W(T)].                     (24)

Thus fine-fine energy work can close on fine timesteps while gas plus total
binding energy is conserved at the common synchronization endpoints. It is not
necessary to average every energy interaction over the root timestep.

The masked potential in (22) is the mass derivative of the binding-energy
component, not the physical potential whose gradient should drive a shell
momentum kick. It includes the group-transfer bookkeeping identified in (8)-(9).
Individual Q_i,l,I must not be used as independent physical heating rates.
They are components of the final energy source. A conservative implementation
uses the following forecast-and-correction construction. Let N_l = D_l K D_l
denote the nested operator, so K_l = N_l - N_(l+1). The current runtime nests
occupied spatial levels, with dyadic substeps selected within each enclosing
interval; masks remain fixed while the corresponding interval is open.

1. At opening time a, solve the nested field N_l m(a) and shell field K_l m(a).
   Keep each open ancestor shell at its own opening field. When level l advances,
   its source predictor is the sum of those ancestor shell fields and its
   current nested field. The nested field forecasts the descendants that have
   not yet advanced. Faster levels subsequently use their own nested field and
   every still-open ancestor shell. No density-rate gravity solve is required.

2. Use physical, time-interpolated gas states for the hydro update. A zero-duration
   Riemann-flux probe supplies the full discrete numerical-flux divergence L(U)
   and the discrete mass-flux energy-work rate from the assembled potential.
   Add rho g to the momentum derivative and the work rate to the energy
   derivative, forming R(U)=L(U)+S(U). Predict an explicit midpoint stage

       U_half = U_start + (h/2) R(U_start).

   Form midpoint halo donor states before spatial prolongation, so donor slopes
   contain the same time evolution as their centers. Reconstruct the midpoint
   cells and evaluate the shared Riemann solver with the Hancock time predictor
   disabled; this stage must not receive a second half-step flux-divergence
   prediction. The conservative update starts from the original physical gas
   state. Reconstruction and flux limiting remain shared with the global
   transport path. The predictor uses the numerical-flux operator whose
   midpoint flux is integrated, rather than assuming that a divergence of
   reconstructed physical fluxes is its time derivative at fixed mesh spacing.
   It also retains the discrete energy source used by the corrector instead
   of substituting p.g for the Mullen work rate.
   Initialize these combined derivatives on every level at the root opening,
   then refresh an active level before it advances. An inactive donor's
   derivative may be old by O(H); its O(H) error is multiplied by half a substep,
   giving the O(H^2) midpoint-state error allowed below.

   The coarse forecast used by fine halos and conventional inactive-source
   interpolation is a separate value,

       U_forecast(a+H) = U(a) + H R_canonical(U(a)).

   Here R_canonical uses the initial fine subface numerical flux on both sides
   of each AMR interface, including a probe-time reflux correction to the full
   hydro derivative. Store this forecast in bank 3 and interpolate it with the
   old state during the interval. The unrefluxed coarse transport endpoint has
   a different boundary flux and can differ by O(H) at fixed spatial resolution;
   it is not the donor endpoint. It remains in its normal update bank until
   accepted reflux and shell closure. Bank 2 retains the root checkpoint.

   Both midpoint and coarse forecast states include the convex kinetic-energy
   remainder

       C = |p_pred - rho_pred p_old/rho_old|^2 / (2 rho_pred).

   For smooth positive-density states C=O(h^2). It helps retain a physical cold
   accelerating forecast without changing its first-order time derivative.
   This addition is predictor-only: it is not inserted into the accepted
   conservative gas-energy source or counted as accepted gravitational work.

3. On each accepted hydro substep, provisionally integrate every open shell's
   momentum source using its frozen opening acceleration and the substep's old
   and provisional new densities. Record these impulses. At shell closure,
   replace their sum by

       delta P_i,l = (h_l/2) [m_i(a) g_i,l(a)
                              + m_i(b) g_i,l(b)],
       P_i <- P_i + delta P_i,l - Itilde_i,l.                (25)

   The closing density includes the accepted reflux correction. This replacement
   establishes the paired impulse (6); the provisional impulses themselves need
   not be exactly paired. At a fast substep endpoint an ancestor shell can still
   be open, so its final momentum balance is obtained when that ancestor closes.

4. For energy, record the work added with the actual substep mass flux. On its
   own level, frame l records a nested-field forecast Qtilde_i,N_l. On faster
   levels the same open frame records only its shell provisional work
   Qtilde_i,l. Open ancestor shells each keep their own contributions. At common
   opening times the nested field plus ancestor shells reconstructs the full
   field; at intermediate times it approximates that sum. This is why the
   individual masked shell works are bookkeeping components rather than
   independent physical heating rates.

5. Every open frame accumulates an all-face mass-flux register over its interval.
   A stored interval-average flux multiplied by h_l represents q_f,I. At closure,
   the work calculation uses the fine subface register on both sides of each
   coarse/fine interface, exactly as mass reflux does. Provisional coarse work
   uses the coarse block's own flux and is subsequently replaced. A final flux
   snapshot cannot replace these interval registers.

6. After descendants have closed and the current level has refluxed, solve the
   closing shell field and compute Q_i,l,I from (22) for every block. Work on
   cells outside A_l is added to a deferred energy register. It is not dropped
   and is not immediately inserted into a slower level's provisional endpoint.
   Let D_i,l be the descendants' deferred work received by a cell on level l.
   The energy correction when frame l closes is

       level l:  E_i <- E_i + Q_i,l,I + D_i,l - Qtilde_i,N_l,
       faster:   E_i <- E_i + Q_i,l,I         - Qtilde_i,l.   (26)

   The first line replaces the nested forecast by the exact current-shell work
   plus the descendants' transfers across the mask. The second line replaces
   only this open ancestor shell's provisional work; the fast cell's own shell
   has already closed. Consume the deferred register exactly once when its
   recipient level closes. It does not create a second hydro flux or an extra
   momentum impulse. Reconcile any separately added temporary kinetic work as
   in (19)-(20); the refined Mullen update adds only its recorded flux work.

This construction is local: it requires no redistribution of a measured global
energy error. At the root endpoint no deferred work or shell ledger remains
open, and the net energy source is precisely the sum in (24). Linear signed
density-rate solves could improve provisional field prediction, but are not a
requirement for the second-order smooth-time argument below.

For a nonreciprocal solver, replace (23) by its endpoint defect:

    R_l,I = [m(a)^T phi_l(b) - m(b)^T phi_l(a)]/2.            (27)

The final energy-budget residual is sum_[l,I] R_l,I, plus the boundary terms
handled as in (17). It is not generally the single root residual (18). Shared
shell operators and diagnostic potentials must satisfy the decomposition (3);
unrelated approximations in different masked solves cannot be presumed to do so.

**7. Conservation is proved above; second-order multirate accuracy has conditions.**

For a smooth scalar mass-transfer rate f(t) and a smooth face potential difference
psi(t), expansion about the interval midpoint gives

    psi_bar_endpoints * integral f(t) dt - integral f(t)psi(t) dt
      = H^3/12 [f_mid psi''_mid - f'_mid psi'_mid] + O(H^5).  (28)

Thus endpoint work has second-order accuracy over repeated synchronization
intervals on a fixed smooth problem. This result concerns time quadrature;
spatial reconstruction and the gravity operator have their own errors.

The provisional field need not have an O(H^2) pointwise error. For a smooth
shell field frozen at its opening time a,

    g_l(t) - g_l(a) = O(H),   phi_l(t) - phi_l(a) = O(H)

through an enclosing interval of length H. Integrating the corresponding source
error gives an O(H^2) error in provisional physical gas states. The nested
forecast on a slower level supplies the first-order source contribution of
unopened descendants; deferred descendant transfers later replace this
forecast through (26). This keeps those states first-order correct in time
instead of omitting a source for an entire coarse interval.

The finite-volume explicit midpoint construction then obtains midpoint states
with O(H^2) temporal error: errors in the opening derivative are multiplied by
half a local step, while the physical gas state supplied to that predictor
already has O(H^2) error. The derivative uses the same numerical Riemann-flux
operator as the corrector. Applying it to donor states before prolongation and
reconstruction includes its spatial variation in their slopes. For a fixed
spatial discretization with bounded derivatives, the resulting hydro flux error
is O(H^2). Its integral over the enclosing interval is O(H^3).
There is therefore no requirement that a provisional frozen-field source
itself integrate to O(h_l^3): its source error can be O(H^2), provided (25)-(26)
replace it at the correct closure and its effect on the intervening hydro
fluxes has the O(H^3) bound above.

With O(H^2) physical-state errors at the force endpoints, accepted endpoint
quadrature has O(h H^2+h^3) error on a shell interval of length h. The work
quadrature (28) has the corresponding bound. Summing intervals covering H
gives O(H^3) local error for fixed, finite refinement ratios and bounded source
and flux derivatives. Linear extrapolation using K_l mdot could improve the
provisional field, but is unnecessary for this second-order argument.

This reasoning requires the following implementation conditions:

- Accepted hydro and species mass fluxes use the same face history, including
  reflux, and the conservative gas-energy flux ledger records the evolved gas
  energy rather than an uncorrected proxy.
- All derivatives used by the midpoint predictor, including halo donors, use
  the same discrete numerical-flux and energy-work operators as the corrector
  and have error at most O(H). Root initialization followed by active-level
  refresh allows inactive derivatives to be O(H) old. Substituting kinetic
  power for flux work would change that operator rather than merely age its
  evaluation. Midpoint donor states must precede spatial prolongation; copying
  only the donor center's derivative omits its slope contribution at fixed mesh.
- Opening and closing shell solves use common-time densities for their active
  sets. Closing endpoints include accepted reflux. Inactive-source interpolation
  uses the dedicated canonical-initial-RHS forecast, not the unrefluxed coarse
  transport endpoint. Fine hydro halos use the same forecast.
- Provisional momentum, nested work, shell work, and deferred transfers are
  recorded separately and replaced exactly once. All close before regridding.
- Step selection resolves each shell's density and force evolution. A slow
  cell label alone is not a bound on rapid incoming mass flux.

The usual stability assumptions are also needed for a global second-order
conclusion. These estimates hold for smooth solutions on a fixed mesh with
bounded derivatives and fixed step ratios as H tends to zero. They do not give
a uniform bound for unbounded time-scale separation, shocks, limiter switches,
or arbitrary changes to an open shell's membership. In particular, exact energy
closure alone is not a temporal-accuracy proof. The implemented scheduler and
its explicit midpoint predictor have smooth fixed-mesh temporal regressions
with observed orders approximately 2.04–2.09 for both modes. Their scope and
results are in the [validation report](validation/gravity-time-integration.txt).
This evidence supports the conditional smooth-case argument; it is not a
general stability theorem. No symplectic property is claimed for this
finite-volume map.

There is an accuracy limitation in using (11) over the whole root interval:
it can underresolve rapidly changing fast-fast work even when momentum is
subcycled accurately. Equation (28) does not give an error bound independent of
fast timescales. Full-potential endpoint work is a conservative reference
coupling, not sufficient evidence that the desired multirate efficiency has
been achieved.

Equations (21)-(26) address that limitation by integrating fast-fast work at its
own cadence. The same smooth analysis applies with H replaced by h_l. Terms
assigned to a slow shell must still vary slowly enough for that shell's step;
fast mass transfer across its boundary can require shortening that step. This
requirement is physical, not removed by an exact global energy identity.

**8. Implementation consequences and checks.**

The derivation establishes the endpoint momentum pairing, mass-transfer
bookkeeping, and local energy-work replacements whose sum exactly balances
binding energy for reciprocal operators. The accepted impulses follow HOLD's
pair cadence, while the physical gas states and flux-based work form an Eulerian
adaptation. A literal substitution of hydro updates for Gadget's independent
particle drifts is neither assumed nor justified.

The runtime design provides common-time density views, nested and shell force
fields, all-face interval mass-flux registers, explicit midpoint hydro prediction,
recorded provisional impulses/work, and deferred energy transfers. Each open
shell retains its opening density and fields until its refluxed endpoint solve
and closure. A fourth state bank separates canonical coarse forecasts from
actual unrefluxed transport endpoints and the third-bank checkpoint. The
conservative work helper computes endpoint work and the
physical-boundary ledger without changing the hydro state. Its provisional
mode reads the current block's own flux; its canonical mode reads the shared
fine subface fluxes used by reflux.

Nested endpoint fields can be reused when the source set and physical time
match, avoiding duplicate partial solves between aligned frames. The partial
source hierarchy's upward pass is still unpruned, and the source/work ledgers
still have central coordination. These remaining costs preclude a speedup or
scalability claim from the cadence argument alone.

`gravity.timeIntegration=hierarchical` applies (25) on each shell's cadence.
`conventional` uses the full source field for each active level's own endpoint
impulse, including time-interpolated inactive sources. Conventional impulses do
not satisfy the pair-cancellation argument because opposite targets can use
different time quadratures. Both modes can use the same shell energy ledger;
its conservation identity depends on accepted mass fluxes and reciprocal
potential operators, not on the momentum schedule. `gravity.energyTreatment`
and `gravity.conserveRegridEnergy` remain independent choices. The energy
conservation claims in this document concern the Mullen option.

Masks and level membership stay fixed while a shell ledger is open. Dyadic
subdivision at aligned boundaries can shorten a level's next interval without
changing its spatial mask. A change that would move a cell or interaction
between masks must close/restart or reject/subdivide the enclosing interval.
Regridding occurs only after every impulse, work register, and deferred transfer
has closed.

Earlier standalone arithmetic checks of the energy identities and smooth
endpoint quadrature, with eight unequal masses, several timestep bins, changing
cell masses, and unequal face weights, found:

- Reciprocal-kernel energy identity error: 7.98e-16 in dimensionless units.
- Deliberately nonreciprocal kernel: energy defect 1.85434894e-4, agreeing with
  independently computed R within 2.33e-15.
- Net paired momentum impulse: maximum component 1.01e-16.
- Endpoint-work local error decreased by a factor of eight per timestep halving,
  consistent with (28).
- Four shells closing at different cadences, with mass flux across their masks:
  energy-budget error at most 1.43e-15. The summed magnitudes of component work
  corrections decreased with local orders 2.962, 2.981, and 2.991 under successive
  timestep halvings, approaching the expected cubic local scaling.

These earlier checks tested algebra and quadrature, including a linearly
predicted provisional-field construction. They are not validation of the
frozen-field physical-state scheduler described here, nor hydro convergence,
stability, FMM conservation, or HPX execution tests for that scheduler.
The implemented scheduler's separate numerical evidence is recorded in the
[gravity integration validation report](validation/gravity-time-integration.txt).
