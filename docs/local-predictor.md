# Local source-aware prediction

**Historical, rejected experiment:** this implementation fails the required
fixed-mesh second-order gravity time convergence. Production now uses cached
numerical midpoint with the original accuracy gates restored. The
joint-refinement results below do not satisfy that requirement. See the
[time-order recovery analysis](time-order-recovery.md) for the defect and
`validation/cached-midpoint/` for current checks. The local-predictor test
targets below are no longer registered and should not be used with the restored
owning-patch API.

The coupled transport path reconstructs at the beginning of the interval and
predicts each reconstructed face locally. It uses two initial ghost layers and
does not exchange midpoint hydro or radiation states. This is a second-order
MUSCL–Hancock construction with a local source solve, not a high-order ADER
space–time discretization.

## Stencil and update

For cell `i`, PLM forms initial left/right states using `i-1,i,i+1`. The first
ghost cell therefore needs the second ghost cell. Prediction uses the physical
fluxes of that cell's reconstructed states:

    D_i = sum_d [ f_d(U_i,d,-) - f_d(U_i,d,+) ] / dx.

Each face and the cell center evolve over `h=dt/2` with their own initial state
and the cell's constant transport forcing:

    dU/dt = D_i + S(U).

The gas and radiation are paired in each local forced implicit source solve.
Gravity supplies its opening acceleration. No neighbor's time derivative or
initial numerical Riemann flux is required. An inadmissible prediction reduces
the initial face offsets and recomputes their physical flux divergence together.
The last fallback retains local source evolution with zero spatial slope.

One numerical Riemann solve per face then supplies the accepted transport
increment. The existing positivity/realizability limiter, composition flux,
boundary ledgers, and AMR reflux use that shared accepted flux. The final local
radiation source solve starts from the original cell average and is driven by
the accepted full-interval transport and other accepted force work. Predictor
increments are never added as extra accepted transport.

For subcycled neighbors, the ordinary initial halo still interpolates the
coarse time endpoints. This does not require a midpoint exchange. Gravity's
separate coarse forecast now uses its local initial physical derivative;
accepted mass-flux gravitational work and endpoint force corrections remain.

## Stiff sources

Radiation pressure forcing remains inside the source solve. For a simplified
relaxation component `F' = D_F - chat*chi*F`, its equilibrium is
`F = D_F/(chat*chi)`. Projecting to the homogeneous source equilibrium `F=0`
would erase diffusion. The existing thick-cell flux correction is also retained
because implicit temporal damping alone does not remove excess numerical
diffusion from the Riemann flux.

The source integrator is L-stable. Unresolved stiff transients can lose temporal
order; this implementation does not claim uniform second order through an
unprepared relaxation layer. Its rational damping leaves an algebraic startup
remainder. The homogeneous stiff test demonstrates stability and conservative
exchange, not exponential accuracy of that remainder. Diffusion tests separately
measure the accepted physical decay rate.

For nonstiff smooth solutions, the accuracy contract is joint space–time
refinement at fixed CFL. The local derivative differs from the fixed-mesh
numerical-flux operator. Accordingly, refining only `dt` at fixed `dx` can reveal
first-order time error even while joint refinement converges at second order.
Source-only temporal order and gravity endpoint quadrature order are separate
questions. Historical numerical-midpoint convergence results do not establish
the new method's order.

## Evidence and reproduction

The [radiation and communication validation records](validation/local-predictor/README.md)
and [gravity validation records](validation/local-gravity-predictor.txt) contain
the measured scope and raw results. The rejected predictor sources and tests
are not included in the production tree, and their old test targets are no
longer registered. The logs document the experiment; they do not imply that it
can be rebuilt from this checkout.

The local radiation tests verify two layers against four-layer storage whose
outer layers contain NaNs, exact agreement across a block split, conservative
stiff relaxation, sharp opacity transitions, and fallback behavior. Independent
small-amplitude scattering-wave errors at 16, 32, 64 cells were approximately
`0.028024, 0.006735, 0.001546`, giving refinement ratios `4.16, 4.35`.
The accelerated smooth Euler wave gave density orders `1.86, 2.01` and the
expected integrated mass, momentum, and force-work balances.

The existing coupled thick Fourier test measured decay-rate ratios to its
analytic discrete expectation of `1.00607, 1.00067, 1.00009` at cell optical
depths `10, 100, 1000`, for both full and reduced light speed. These tests have
different purposes: convergence, diffusion balance, conservation, and stencil
dependence must each pass.

[Halo accounting](halo-communication-accounting.md) defines the measured
communication quantities and the production before/after fixture. Fewer halo
requests do not establish an application speedup. The local predictor performs
more source solves per cell; the next performance comparison must measure full
steps at matched error, mesh, physics, placement, and physical elapsed time.
Gravity communication, AMR reflux, and non-halo traffic remain outside that
fixture's bandwidth claim.
