# Recovering fixed-mesh second-order time accuracy

Investigation, 2026-09-28. The rejected local predictor has been replaced in
production by cached numerical midpoint. The original gravity and
rotating-gravity temporal convergence assertions and test names are restored.
The previous first-order results are failures of those requirements; historical
passing logs used weakened assertions. Current validation is recorded in
`validation/cached-midpoint/`.

## Mechanism

Write the smooth fixed-mesh problem as `U' = K_h(U)`, including the numerical
transport and discretized source terms. Write reconstructed face states as
`R(U)`. A midpoint flux needs

    R(U(t + dt/2)) = R(U) + (dt/2) DR(U)[K_h(U)] + O(dt^2).

The local predictor instead evolves each initial face with its cell's physical
flux divergence and a local source. This generally differs from the expression
above even in smooth nonstiff flow. The difference includes both the numerical
versus physical flux divergence and the time derivative of the reconstruction's
slopes and primitive-variable conversion. An O(dt) face error produces O(dt^2)
error per step and O(dt) accumulated time error at fixed mesh spacing.

Changing only the common center derivative while freezing the slopes is
insufficient. This issue exists without gravity or radiation; gravity-specific
force/work staging and AMR interpolation impose additional consistency conditions.

### Independent check

[time-order-analysis.py](validation/local-predictor/time-order-analysis.py) uses
periodic scalar advection with unlimited centered PLM slopes and an upwind flux.
For a Fourier mode it compares each method against the exact exponential of the
same spatial operator, eliminating spatial-reference error. At 32 fixed cells,
mode 4, final time 0.25, and 128/256/512/1024 steps:

| Method | Observed orders |
| --- | --- |
| Local Hancock | 0.9840, 0.9925, 0.9964 |
| Numerical center derivative, frozen slopes | 0.9723, 0.9869, 0.9937 |
| Numerical midpoint, then reconstruct | 2.0012, 2.0007, 2.0003 |
| AB2 with midpoint startup | 1.9996, 1.9993, 1.9995 |

These establish the mechanism and the scalar remedies. They are not new
production self-gravity, radiation, or AMR convergence results.

## Recommended recovery: cache and share the numerical midpoint stages

1. Obtain the initial two-layer gas/radiation halo once per block and stage.
   Reuse it for gas transport, radiation transport, and radiation material data.
2. Evaluate the initial shared numerical fluxes. Form a consistent cell-centered
   numerical RHS, including the gravity work prescription. Retain the original
   AMR treatment of coarse donor derivatives and initial flux reflux.
3. Form source-aware midpoint **cell averages**. Radiation exchange remains a
   local forced implicit solve driven by the numerical transport RHS.
4. Exchange two layers of midpoint cell averages, prolonging the midpoint donor
   states themselves. Reconstruct these averages with the extra Hancock time
   prediction disabled. Apply physical boundaries at their stage times.
5. Use one shared set of accepted midpoint fluxes for conservative transport,
   composition, reflux, mass-flux gravity work, and the full local source solve.
   Retain the existing gravity endpoint closure and consistent coarse forecasts.

For a smooth source `S`, the forced half-step gives
`Uhalf = U + dt/2*(L_h(U)+S(U)) + O(dt^2)`. Evaluating numerical transport at
`Uhalf` and integrating the full local source problem with that constant
transport forcing supplies both transport and source cross terms at second
order. This Taylor argument assumes resolved, smooth source dynamics; it does
not establish uniform order through an unresolved stiff initial layer.

Stage caches must be tied to block, field bank, donor times, and generation,
and remain available across the predictor/corrector scheduling boundary.
Worker scratch alone is insufficient when tasks move between workers. Gravity
and radiation must use the same source stages and work discretization as the
accepted update. In particular, analytic `p dot g` is not automatically the
fixed-mesh derivative of the discrete mass-flux energy work.

### Expected communication, not a new measurement

The recorded homogeneous coupled fixture used five gas plus three radiation
halo reads per block before the local-predictor conversion. Sharing and caching
the two stages targets two gas plus two radiation reads. With the same halo
plans, this is a calculated reduction of 55.56% in 1D (288 to 128 requested
scalar values per step) and 53.85% in 3D (186,368 to 86,016).
The 3D old baseline itself was calculated, not measured. These estimates exclude
gravity traffic, species, reflux, AMR endpoint multiplicity, cache migration,
and other communication. They do not establish achieved bandwidth or speedup.

## Can we retain one exchange and two layers?

For the current stateless one-step method and unchanged numerical spatial
operator, generally no. The reconstructed numerical RHS has stencil radius
two. Its second time derivative contains `DK_h[K_h]`, which generically reaches
radius four. Two initial ghost layers do not contain that information. A single
exchange of four initial layers permits redundant local stage computation on
a uniform mesh, but may increase bytes: a full cubic halo around a four-cell
3D block has 448 cells at width two and 1,664 at width four. Actual sparse
stencils and AMR communication would need separate accounting.

Using history changes the answer. AB2 integrates the numerical face flux using
`1.5*F_n - 0.5*F_previous` for equal timesteps, requiring one fresh two-layer
halo per step. It is a candidate for a later design, with separate requirements
for variable step coefficients, initialization, regridding, subcycling,
conservative shared histories, positivity, and stiff source coupling. Its
oscillatory stability also needs scrutiny; see the primary
[MITgcm time-stepping documentation](https://mitgcm.org/public/r2_manual/latest/online_documents/node35.html).
It is a different integration method, not a correction to the local predictor.

## Required evidence before accepting a production repair

- Original fixed-mesh gravity and rotating-gravity convergence thresholds.
- Conservative gravity work and AMR interface behavior, including subcycling.
- Nonstiff coupled temporal order and the existing stiff/diffusion limits.
- Actual halo counts for the repaired path, followed by comparisons at matched
  physical accuracy. The current local predictor's lower counts cannot be
  transferred to a different algorithm.

No production solver changes or production test runs were made during this
investigation. The scalar diagnostic was run; original accuracy gates were
restored.
