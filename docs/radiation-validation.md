# Radiation coupling and diffusion validation

These checks were run on 2026-09-26–27 against the opacity-aware face flux,
coupled source-aware midpoint integrator, and application runtime. They
establish the specific limits and budgets listed below. They do not establish
a rotating-star equilibrium or arbitrary dynamic-diffusion accuracy.

## Source and application checks

All 16 `radiationCouplingChecks` pass in the 3D serial build. Independent
references include a scalar thermal-equilibrium root, a linear momentum
relaxation solution, temporal refinement, and the forced opaque flux balance.
The tests also cover moving LTE, a cold absorbing streaming beam, full M1
anisotropy, high-kinetic-energy dual-energy heating, and failure atomicity.
A gas-dominated cell with `rho=1 g/cm^3`, `T=1 K`, and initially zero radiation
recovers the very small LTE radiation energy, approximately
`7.56e-15 erg/cm^3`, instead of losing it to a gas-scaled solver tolerance.

The full-light-speed transverse-diffusion regression boosts comoving M1
moments into the laboratory frame and independently balances the meridional
radiation force. Its 32 cases cover axis-aligned/oblique velocities and source
optical intervals through `c chi dt=1e12`. At the highest stiffness, maximum
relative gas thermal and radiation-energy drifts are `8.88e-16` and `3.33e-15`;
the flux change divided by `c aT^4` is at most `1.10e-16`. Resolving a tiny
transverse flux inside much larger Cartesian advective components has a
cancellation floor: at comoving flux factor `1e-12`, its relative change can
reach `7.78e-6`. This is a different normalization from the thermal error and
is not a claim of machine-relative-precision diffusion flux in that limit.

A captured dilute atmosphere cell from the opaque-star startup exposed a
nonlinear initial-guess failure when transport brightened it by five orders
of magnitude. An admissible transported starting guess, selected only when it
reduces the stage residual, fixes that failure without changing the equations,
tolerances, atmosphere, or timestep. The exact cell is retained as a
conservation regression. The complete star startup and its limitations are
recorded in [opaque-star validation](validation/opaque-star/README.md).

Production periodic-grid tests at full and reduced light speed give relative
combined-energy and momentum budget residuals below `2e-16`. At reduced
speed these are the weighted invariants defined below. Mixed-level
subcycling, refinement and coarsening, and positive/negative mesh rotation
also meet those budgets in the tested cases. Zero opacity reproduces the
uncoupled update bit for bit. A deliberately infeasible numerical step
restores the published state; the next valid step matches a clean run.

The coupled polytrope test with self gravity checks globally synchronized,
hierarchical, and conventional time stepping, full and reduced light speed,
regridding, and a rotating case. The largest scaled combined-energy residual
in that serial test is `3.23e-16`, including gravitational binding energy and
the boundary ledger. This verifies a conserved budget during evolution; the
initial polytrope is not a gas/radiation equilibrium and this is not evidence
that such an equilibrium remains stationary.

Temporal refinement on a fixed spatial grid gives these errors relative to
a finer-timestep reference. Halving the timestep from left to right leaves
spatial discretization unchanged:

| Mesh | Coarse timestep | Half timestep | Quarter timestep |
| --- | ---: | ---: | ---: |
| Uniform | 0.00788585 | 0.00136863 | 0.000312280 |
| Mixed levels 1 and 2 | 0.000786152 | 0.000182610 | 0.0000441965 |

The finest error ratios are 4.38 and 4.13, consistent with second-order
temporal convergence in these smooth, finite-stiffness problems. They are
not a uniform-in-stiffness order proof. Two separate three-level tests use
the runtime's reported CFL interval, including the new predictor cap, for
both smooth material and a sharp stellar surface. Both verify admissibility,
expected level substeps, and mass/combined-energy budgets.

A separate imposed-acceleration test checks work against an independently
implemented long-double RK4 solution of the homogeneous gas/radiation ODEs.
Its periodic box is wide enough that transport during the comparison is
negligible; the entropy auxiliary supplies the gas temperature. With 8, 16,
and 32 steps, relative work errors are `6.06e-11`, `1.53e-11`, and `3.83e-12`
at full light speed, and `2.82e-10`, `7.12e-11`, and `1.81e-11` at `eta=0.2`.
The reference changes by less than `1e-13` relative when its own resolution
doubles from 512 to 1024 steps. This specifically checks acceleration work
on the radiation-induced midpoint momentum, including auxiliary-selected
thermodynamics.

Five conservation/diagnostic tests pass in both serial and HPX builds,
including reading the new fields back from Silo. Three HPX serialization
tests pass with profiling disabled. All ten application integration tests
pass in the serial build. Under HPX, nine pass together on two independent
TCP localities on this host; the longer temporal-refinement test passes
separately on one locality. The distributed gravity cases have the same
reported energy residuals as the serial cases. These checks establish
correctness for the tested partitioning and schedules, not cluster scaling
or whole-application performance. A radiation-only build with hydro disabled
also succeeds.

The supplied `radiation-matter` example
completes 11 steps to `t=0.001 s`; its final combined-energy drift is zero
at output precision and its largest intermediate scaled drift is at roundoff.
An intentionally inadequate reduced-speed example with `eta=0.01`, opacity
`1000 cm^2/g`, and diagnostic length `1e9 cm` also completes: its maximum
`rslaCriterion` is about 1848 and is recorded without stopping execution.
Completion demonstrates the requested nonfatal behavior, not physical
accuracy at that reduced speed.

## Governing reference

For stationary material, opaque gray radiation has

    F = -c/(3 chi) grad Er,       chi = rho*kappa.

The energy transport speed uses `chat=eta*c`; the stored physical flux still
uses `c`. In LTE, ideal-gas internal energy density is `u=rho*c_v*T` and
`Er=a*T^4`. At fixed density, linearizing the conserved weighted energy
`u+Er/eta` gives a temperature perturbation with diffusion coefficient

    D_eff = chat/(3 chi) * 4 Er/(4 Er + eta*u).

For a periodic Fourier mode on a uniform grid, the centered discrete
Laplacian has eigenvalue `-4 sin^2(k*dx/2)/dx^2`. The reference decay rate
is therefore `D_eff * 4 sin^2(k*dx/2)/dx^2`. This analytical reference is
independent of the numerical source root solver.

## Measured checks

`radiationDiffusionChecks` exercises the production M1 reconstruction,
face correction, shared realizability limiter, and, for the evolution test,
`advanceCoupledPatch` and `coupleForced`.

The spatial Fourier check holds an eight-cell grid fixed while increasing
cell optical depth. It measures the energy face-flux divergence before a
time update. Initial `F=0` deliberately exposes HLL numerical diffusion.
The same corrected flux survives the production limiter unchanged at a
nonzero radiation CFL. Both `eta=1` and `eta=0.1` give these ratios to the
expected discrete diffusion coefficient:

| Cell optical depth | Measured / reference |
| ---: | ---: |
| 10 | 1.01521 |
| 100 | 1.00244 |
| 1,000 | 1.00025 |
| 10,000 | 1.00003 |

The full planar gas/radiation evolution check uses 16 periodic cells across
1 cm, `rho=1e-17 g/cm^3`, `T=100 K`, `gamma=5/3`, mean molecular weight 1,
and a temperature mode of relative amplitude `1e-5`. Gas and radiation start
in LTE, with the centered diffusion flux. Eight steps use
`dt=0.3 dx/chat`. The measured logarithmic decay of the weighted total-energy
mode is compared with the analytical LTE heat-capacity factor above:

| Cell optical depth | Decay / reference, eta=1 | Decay / reference, eta=0.1 |
| ---: | ---: | ---: |
| 10 | 1.0068 | 1.0068 |
| 100 | 1.00066 | 1.00066 |
| 1,000 | 1.00008 | 1.00009 |

Periodic weighted total energy meets a relative `5e-13` tolerance. This is
a short-time decay measurement, not an integration over many diffusion
times. Hydrodynamics is active; its mechanical response is negligible on
the tested light-crossing intervals. The planar evolution test runs in the
1D executable; 2D and 3D explicitly skip it.

Additional tests require:

- Exactly unchanged numerical flux at zero extinction.
- Quadratic thin-limit flux perturbation with cell width for a smooth state.
- Static opaque pressure `Er/3` despite a transient input reduced flux of 0.7.
- Unattenuated material and ALE advection for a uniform moving source
  equilibrium, including oblique velocities, full and reduced light speed,
  and cell optical depths through `1e8`.
- A realizable finite predictor across a sharp opaque-core/transparent-
  atmosphere interface with energy contrast `1e14` and cell optical depths
  100 and `1e-10`. The test also confirms that omitting the finite predictor
  limiter interval produces an infeasible state from the instantaneous flux.
- Positive gas thermal states, realizable radiation, and periodic weighted
  energy conservation in nonlinear transition-regime runs described below.

The six spatial/flux/predictor tests passed in both 1D and 3D; all eight tests
passed in 1D. All 13 `radiationChecks` tests passed in both dimensions
after the stricter face-cone acceptance, reconstruction, and interpolation
changes. These runs use the explicitly HPX-free serial build.

The initial zero-step flux probe applies the shared face limiter over the
actual half-step predictor interval without advancing its geometry or
publishing a transport update. The face limiter requires its convex states
to lie inside the radiation cone. A tolerance based on the brighter neighbor
would otherwise allow an overshoot exceeding the faint cell's own roundoff
budget. Final summed transport updates retain their scale-aware roundoff
repair. Reconstruction removes only tolerated face-state cone overshoot
before a subsequent rotation of its components; physically invalid trial
states remain invalid for slope rejection. AMR interpolation checks the exact
cone when choosing its common slope fraction; a regression checks both
conserved child averages and subsequent rotation of coarse/fine corner
states. No independent child or source-state clipping is introduced.
These restrictions can activate the original low-order flux near a
sharp interface; a diffusion-coefficient claim is not made there.

## Transition-regime stability check

The nonlinear planar test uses 16 cells with density varying by 40%,
temperature by 20%, and sinusoidal velocity of amplitude `0.001 chat`.
It starts on the moving energy/momentum-source equilibrium and runs 24 steps at
`dt=0.3 dx/chat`, using mass opacity so extinction evolves with density.
Nominal cell depths are 0.03, 0.3, 1, and 3, each at `eta=1` and `eta=0.1`.
The assertions check every accepted gas/radiation state and weighted total
energy, with the same `5e-13` conservation tolerance. This probes the flux
interpolation around optical depth one; it is not a general nonlinear
stability proof or an accuracy reference for that regime.

## Dynamic diffusion, force balance, and conservation limits

The opaque energy flux retains material advection, and the implicit source
stages include the accepted pressure-transport increment. These are necessary
for moving, trapped radiation. The tests above establish a static diffusion
coefficient and a uniform moving face-flux identity; they do not measure
advection, compression work, or wave speeds in a spatially varying opaque
flow. In particular, L-stability of the local source solve does not by itself
prove second-order accuracy uniformly as the source relaxation time vanishes.
An advected opaque radiation pulse, radiation-modified acoustic modes, and a
compression/work comparison are still needed before claiming general dynamic
diffusion accuracy. Those comparisons must vary timestep, spatial resolution,
and reduced light speed independently.

At full light speed, `|v| tau_L/c >= 1` can describe physical radiation
trapping. It does not invalidate a simulation merely because the diagnostic
is large. At reduced light speed, even a correctly measured stationary
physical flux `F=-c grad Er/(3 chi)` is insufficient: the energy evolution,
thermal capacity, and radiation inertia follow the modified equations.
The `eta` dependence in the measured LTE diffusion coefficient above is a
concrete example. A successful `eta=0.1` test against its modified reference
is not evidence of agreement with `eta=1` for the same physical transient.
The static-diffusion timescale criterion and its limitations are discussed
in [Skinner and Ostriker (2013), Section 3.2](https://arxiv.org/abs/1306.0010).

Both choices of light speed still use the documented equal-opacity,
nonrelativistic mixed-frame source approximation, including all equal-opacity
first-order terms in the energy and momentum sources. The algebraic moving
M1 equilibrium used by the flux is the Lorentz transform of isotropic thermal
radiation and satisfies both retained source balances. The subluminal domain
of that construction does not establish relativistic evolution accuracy.
The importance of terms beyond the retained approximation
depends on both velocity and optical-depth scaling, especially when large
work and emission/absorption terms nearly cancel. No extension to large
velocities or arbitrary opacity models follows from the present tests.

Pressure-gradient forcing balances stiff radiation momentum relaxation
locally. This is different from preserving an entire discrete stellar
equilibrium. The code does not construct a shared equilibrium stencil for
gas pressure, radiation pressure, self gravity, and centrifugal support.
Neither a constant luminosity across abrupt opacity changes nor a rotating
hydrostatic star has been verified here. The arithmetic face extinction has
the correct two-half-cell resistance interpretation; that fact alone does
not prove force balance or accuracy at an opacity discontinuity. A global
claim of an asymptotic-preserving, well-balanced moving-flow scheme would
therefore exceed the evidence recorded here.

The paired local source increments preserve
`Egas+Er/eta` and `rho v+F/(eta c^2)`. With `coupleForced`, their change is
exactly the corresponding sum of the supplied transport/gravity increments,
up to rounding; these quantities are not held fixed in a cell experiencing
transport. Periodic global balances cancel common face fluxes. External
boundary fluxes, applied forces, and gravitational energy require their
respective budget terms. At `eta=1` the source invariants are physical total
energy and momentum; at `eta<1` they are modified invariants and should not
be reported as physical conservation. Conservation is necessary evidence,
but does not bound diffusion, radiation work, or force-balance error.

## Diagnostic interpretation checked against the implementation

Silo cell fields and conservation-CSV maxima call the same diagnostic
function. `tau_cell=rho*kappa*dx` uses the same local extinction as the face
correction. The other indicators use `tau_L=rho*kappa*L`, with a fixed
physical length, so refinement alone cannot improve the indicated RSLA
timescale ordering. `tau_L` is a local estimate, not a ray-integrated optical
depth; a maximum of these local indicators is not a general bound on a
heterogeneous star's transport error. Their velocities are inertial, so
rotating the mesh basis does not remove bulk material advection.

The term `4 Er/(9 rho)` in the acoustic estimate is motivated by isotropic,
trapped radiation. The resulting `a_bound` is a timescale proxy, not a
computed characteristic speed or a universal bound for anisotropic or
nonequilibrium radiation. A small reported RSLA criterion must still be
checked by increasing `chat` on the physical problem of interest. A large
full-light-speed value is recorded without rejection, as intended.

## Reproduction

Use a configured build containing both hydro and radiation. For example,
with dependencies installed or their prefixes supplied:

```bash
cmake -S . -B /tmp/octo-radiation-validation \
  -DOCTOTIGERII_WITH_HPX=OFF -DOCTOTIGERII_BUILD_TESTS=ON
cmake --build /tmp/octo-radiation-validation \
  --target radiationDiffusionChecks-1d radiationDiffusionChecks-3d \
    radiationChecks-3d radiationCouplingChecks-3d \
    radiationIntegrationChecks-3d radiationDepthChecks-3d \
    radiationConservationChecks-3d octoII-3d -j2
ctest --test-dir /tmp/octo-radiation-validation \
  -R '(1d|3d)\.radiation(Diffusion|Coupling|Integration|Depth|Conservation)?Checks' \
  --output-on-failure
/tmp/octo-radiation-validation/octoII-3d \
  --config examples/radiation-matter.ini \
  --output.directory=/tmp/octo-radiation-example
/tmp/octo-radiation-validation/octoII-3d \
  --config examples/radiation-matter.ini \
  --radiation.lightSpeedRatio=0.01 --radiation.opacity=1000 \
  --radiation.diagnosticLength=1e9 --runtime.stopTime=1e-5 \
  --output.directory=/tmp/octo-radiation-diagnostic-example
```

With an HPX-enabled build at `/tmp/octo-radiation-hpx`, the distributed and
single-locality temporal checks are:

```bash
env OCTOTIGERII_TEST_LOCALITIES=2 \
  GTEST_FILTER=-RadiationIntegration.CoupledTransportConvergesAtSecondOrderInTime \
  python3 tests/distributed.py \
  /tmp/octo-radiation-hpx/3d/tests/radiationIntegrationChecks-3d
/tmp/octo-radiation-hpx/3d/tests/radiationIntegrationChecks-3d \
  --hpx:threads=2 --hpx:bind=none \
  --gtest_filter=RadiationIntegration.CoupledTransportConvergesAtSecondOrderInTime
```

The distributed launcher needs local TCP sockets. The tests avoid GoogleTest
thread-local scoped traces across HPX waits, since a resumed HPX task can run
on a different operating-system thread.

The initial checks reported here were compiled directly with GCC against
the serial build's generated dimension headers and GoogleTest, using the
current radiation transport/source translation units. The registered CMake
targets above reproduce the same test source through the normal build.

The asymptotic flux construction is documented in
[radiation-coupling.md](radiation-coupling.md). If the shared face limiter
changes a corrected flux, conservation and admissibility take priority;
these checks do not establish the physical diffusion coefficient in such
cells. No diffusive timestep acceleration is claimed: the method retains
the radiation transport CFL restriction.
