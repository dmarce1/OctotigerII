# Analytic verification

Every problem implements `problemReference(Config const&)` alongside its
initial conditions. Analytic comparison is **automatic by default** when a
reference is available. For gravity builds, the default reference is instead the direct sum of the
current discrete mass distribution, including for evolving collapse. This
checks gravity fields only. A non-gravity problem without a continuum reference
reports `unavailable`, with a reason. An unavailable reference is never reported as a passed accuracy test.

## Available references

| Problem | Dimensions | Exact solution and restrictions |
| --- | --- | --- |
| `sod` | 1, 2, 3 | Ideal-gas Riemann solution for the configured `hydro.gamma`; valid until the first wave reaches an outflow/reflecting x boundary; analytic data on both x faces removes that limit. Periodic x faces disable this reference. |
| `gravity-sphere` | 3 | Direct discrete reference by default; optional continuum isolated uniform sphere, including interior and exterior potential and acceleration; potential vanishes at infinity. |
| `gravity-gaussian` | 3 | Direct discrete reference by default; optional continuum isolated spherical Gaussian truncated at radius half the domain width, including the contribution of exterior shells to the interior potential. |
| `streaming` | 1, 2, 3 | Periodic translation along the diagonal of the active dimensions, with transport speed `lightSpeedRatio*c` and physical flux magnitude `c*E`. |
| `polytrope` | 3 | Direct gravity by default; `verification.gravityReference=continuum` compares against the Lane–Emden equilibrium profile (tenuous atmosphere approximation). |
| `collapse` | 3 | Direct gravity reference at the current time; no analytic hydrodynamic solution. |
| `kelvin-helmholtz`, `radiation-pulse` | As supported by their manifests | No full analytic reference currently implemented; status is explicitly unavailable. |

The Gaussian benchmark now has a **spherical cutoff**, rather than the old
implicit cubical truncation. Its center density remains `1e4 g/cm^3`, its
standard deviation is `0.15` times the domain width, and its cutoff radius is
`0.5` times that width. Initial conditions and the reference use the same cutoff.
This is a deliberate change to this benchmark's density distribution.

## Running a comparison

Use the executable for the chosen problem and dimension. For example, starting
from HOME, with the checkout at `~/workspace/OctotigerII`:

```bash
cd "$HOME/workspace/OctotigerII"
cmake -S . -B release \
  -DCMAKE_BUILD_TYPE=Release
cmake --build release -j
./release/octoII-1d --problem.name=sod \
  --problem.name=sod --config=bin/hydro_tests/Sod/inputs \
  --verification.analytic=on \
  --mesh.cells=32 --mesh.level=2 \
  --output.directory=output/sod-analytic
ctest --test-dir release --output-on-failure
```

Supply the normal HPX dependency paths for your installation, or configure with
`-DOCTOTIGERII_WITH_HPX=OFF` for the serial backend. Gravity builds use
`--problem.name=gravity-sphere` or `--problem.name=gravity-gaussian` with `octoII-3d`.
Their input files are `bin/gravity_tests/Sphere/inputs` and
`bin/gravity_tests/Gaussian/inputs`. Gravity tests require `runtime.stopTime=0`.

| Option | Default | Meaning |
| --- | --- | --- |
| `verification.analytic` | `auto` | `auto`: compare when available; `on`: require a valid reference; `off`: disable. |
| `verification.relativeL1Tolerance` | `-1` | A nonnegative value requires every nonzero-reference field's relative L1 error to be at most this value. Failure exits nonzero. `-1` disables the gate. |
| `verification.absoluteTolerance` | `1e-12` | Maximum absolute Linf for identically zero reference fields when the gate is enabled, in each field's CGS units. |

These options support both dotted INI keys and `[verification]` sections.
Tolerances are opt-in: the resolution, problem, and desired accuracy determine
appropriate values. Reporting runs automatically without an accuracy gate.
A gate requires an available reference and cannot be combined with `off`.

## Output and norms

The final console summary lists each field's **relative** L1, L2, and Linf error.
With a continuum reference, each written Silo frame additionally contains `<field>Exact`
and `<field>Error`, where the signed error is numerical minus exact. For Sod,
these include density, pressure, total gas energy, velocity, and momentum;
for gravity, potential and acceleration; for streaming, radiation energy and
flux. Momentum, velocity, acceleration, and radiation flux retain separate
scalar components over the active dimensions. For example,
`radiationFluxXExact` and `radiationFluxXError` are directly selectable in
Pseudocolor, with the latter storing the signed x-flux error. Embedded VisIt
vector expressions combine these scalars in x, y, z order, exposing names
such as `radiationFlux`, `radiationFluxExact`, and `radiationFluxError` for
Vector plots. The expressions add no duplicate field arrays; inactive
directions in 1D/2D use zone-centered zero expressions. The console and JSON
error norms continue to report each component separately.

`analytic-errors.json` in the output directory contains the latest written
frame's time, resolution, reference status, relative `L1`, `L2`, and `Linf`
norms, and `absoluteL1`, `absoluteL2`, `absoluteLinf`, and reference norms in
each field's stated CGS units. Schema version 3 retains dimensionless `L1`,
`L2`, and `Linf` and adds reference kind, population/source counts, selected
global target indices, RNG/seed and sampling uncertainty estimates. `cells`
is the number of compared targets; `totalCells` is the population size. The final frame is always written when
output is enabled. With `output.enabled=off`, the console comparison still runs
and any accuracy gate still applies; no output files are created.

References are evaluated at **cell centers**, at the actual synchronized
snapshot time, including nonzero times for Sod and streaming. Let `e_i` be the
numerical value minus its exact point value and `V_i` the cell measure. Then

- `L1 = sum(V_i*abs(e_i)) / sum(V_i*abs(f_i))`;
- `L2 = sqrt(sum(V_i*e_i^2) / sum(V_i*f_i^2))`;
- `Linf = max(abs(e_i)) / max(abs(f_i))`.

Here `f_i` is the analytic field value. Each denominator uses the same norm as
the numerator, and the volume measure cancels where appropriate. If the
analytic field is identically zero, a relative norm is undefined; JSON uses
`null` and the console reports the absolute norms with CGS units instead.
The numerical hydro fields are finite-volume values; these comparisons use
exact point samples, **not exact volume averages**. Discontinuities therefore
have sampling error, and Linf need not decrease monotonically as they move
relative to cell centers. The convergence tests check L1 reduction.

## Direct gravity reference (default)

`verification.gravityReference=direct` applies to every gravity-enabled build.
The reference uses the current snapshot density, cell mass `rho*dx^3`, isolated
`-G*m/r` potential and its acceleration, with the same self-interaction exclusion
as the solver. Every selected target is summed against **all nonzero sources**;
only the targets are sampled. Vacuum target cells remain in the population.
No continuum density approximation or FMM operator enters the reference sum.
CGS contributions are accumulated in `long double` to reduce summation roundoff.

| Option | Default | Meaning |
| --- | --- | --- |
| `gravity.multipoleOrder` | `5` | Integer expansion order from 1 through 10. |
| `verification.gravityReference` | `direct` | `direct` or `continuum`. The latter retains the smooth problem-provided reference. |
| `verification.directSamples` | `0` | `0`: automatic; positive: requested target count, clamped to the number of cells and overriding the automatic pair budget. Set it at least as large as the population to force a full comparison. |
| `verification.directMaxPairs` | `20000000` | Automatic work budget. Compare every target if `Ntargets*Nsources` fits; otherwise sample `min(1024, max(1, budget/Nsources))` targets. At least one complete target sum is always done. |
| `randomSeed` | `5489` | Global nonnegative 64-bit signed integer seed (0 through 9223372036854775807). |

For example, append these to the usual sphere run:

```bash
--gravity.multipoleOrder=10 --verification.directSamples=2048 --randomSeed=5489
```

The generator is `std::mt19937_64`. Floyd sampling without replacement uses an
explicit rejection/remainder mapping, avoiding library-dependent distribution
algorithms. Targets are identified by global x-fast index `(z*N+y)*N+x` and
sources are summed in that same canonical order. For a fixed mesh, target count
and seed, changing order, opening angle, block order or HPX scheduling preserves
the sample. Explicit sample counts also keep the same targets as density evolves.
Automatic counts may change if the number of nonzero sources changes.

Console and JSON direct-reference norms use only the selected targets with
equal cell weights. On a mixed-level AMR mesh these are cell-sampled norms,
not volume-weighted continuum errors. Source masses use actual cell volumes. For a sampled comparison, the
console reports approximate 95% half-widths for relative L1 and L2, and warns
when either exceeds the corresponding measured error. A one-target sample or
zero sample reference norm produces an unavailable-uncertainty warning.

The uncertainty estimate uses a ratio delta method, including the covariance
of error and reference magnitudes and the finite-population correction. With
`k` sampled targets out of `N`, for L1 let `a_i=abs(error_i)`, `b_i=abs(ref_i)`,
`R=mean(a)/mean(b)`, and `z_i=(a_i-R*b_i)/mean(b)`. The estimated half-width is
`1.96*sqrt((1-k/N)*sum(z_i^2)/(k*(k-1)))`. L2 applies the same estimator to
squared errors and reference values, then divides the half-width for their
ratio by twice the relative L2 norm. A full census has zero sampling uncertainty.
These are approximate normal intervals, not guarantees: a small sample can miss
rare large errors, and zero observed errors cannot exclude unobserved errors.
Linf is explicitly the **sample maximum**, with no full-population confidence
bound. Accuracy gates likewise apply to the sampled norms.

Direct comparisons populate the console/JSON report; they do not create dense
`Exact`/`Error` Silo arrays, since unsampled cells have no direct reference.
The usual numerical gravity Silo fields remain available. Select
`verification.gravityReference=continuum` to obtain the original dense analytic
Silo fields and to measure density-discretization plus FMM error against the
smooth sphere/Gaussian solution. Direct mode isolates the FMM error for the
actual cell-center point masses.

## Automated checks

`analytic.references` independently checks:

- Sod star-state values and shock jump conditions at a second gamma;
- spherical potentials and forces against shell integration, including the
  Gaussian center, small radii, cutoff, and exterior;
- streaming wraparound with reduced light speed;
- zero-error synthetic fields, deliberately injected errors, availability,
  disabled mode, and invalid Sod reference times.

`analytic.convergence` runs each supported compiled problem at two resolutions,
checks every nonzero reference field's relative L1 decrease, and requires a
bounded final error. Gravity uses 16 and 32 cells per axis and opening angle
0.35; 3D transport uses 16 and 32, other transport builds use 32 and 64.
Silo tests read back exact/error arrays, their units, and multivar registration.
Application tests parse the JSON report and check nonzero exit on failed gates.
HPX serialization tests cover the verification options and global seed.
`gravity.directReference` checks point-mass fields and self exclusion, vacuum,
explicit/automatic sampling, reproducibility and block-order independence,
plus the uncertainty estimator and warnings. The `gravity` check compares all
orders 1 through 10 against an independent full direct sum. The existing mesh
convergence test explicitly uses continuum gravity.

## Adding a problem

Implement `problemReference` in that problem's `problem.cpp`. Return a named
`verification::Reference` with an evaluator returning typed CGS `ExactState`
values for all enabled physics. Give any time restriction in `validUntil` and
explain it in `reason`. Return an empty evaluator and a useful reason when no
analytic result is available. The common reporting and Silo paths then work
without changes. Add independent reference identities and a problem-appropriate
resolution test to `tests/analyticChecks.cpp` and register it in CMake.

Sod's reference matches the rarefaction and shock curves using a bracketed
pressure solve, then evaluates the self-similar regions. It is independent of
the evolution code's approximate HLLC solver. For the derivation, see
[Ketcheson, LeVeque, and del Razo, Euler equations](https://www.clawpack.org/riemann_book/html/Euler.html).
Spherical gravity follows the Newtonian shell theorem, including outer shells
in the potential and enclosed mass in the force.

Streaming also supports paired analytic faces per axis; see [boundary conditions](boundaries.md).
