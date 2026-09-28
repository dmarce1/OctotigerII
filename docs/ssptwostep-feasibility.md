# Feasibility of the supplied SSP two-step Runge–Kutta draft

2026-09-28. This is a numerical and architecture audit of the user-supplied
draft, not a production implementation. No runtime changes or production
builds were made for this audit.

## Finding

The third-order transport method is a credible candidate for implementation.
The printed formulas reproduce the supplied coefficient table and satisfy the
nonlinear third-order conditions in the numerical audit. Independent scalar
calculations recover third order with fixed and varying timesteps.

It requires **two fresh spatial RHS evaluations and normally two halo exchanges
per step**, each using two ghost layers. The previous RHS is retained locally.
This saves a stage relative to SSPRK3; it has the same number of fresh stages
as a conventional two-stage RK2 method. It does not meet a one-exchange target
by itself. The existing gravity and radiation source integration is second
order; third-order transport alone does not make the complete solver third
order.

## Reproducible coefficient and convergence audit

Run [ssptwostep-audit.py](validation/local-predictor/ssptwostep-audit.py) with
Python 3.11 or later. It uses the supplied formulas, standard-library arithmetic,
and no simulation dependencies.

At equal timesteps, the fitted coefficients are approximately

    c1 = 0.666107436829     c2 = 0.333892563171
    c3 = 0.003680245579     c4 = 0.586930483271
    c5 = 0.409389271151
    v  = 1.368738564140     v5 = 0.651913215516

The intermediate state's time is

    t_stage = t_n + h * (c1*v - alpha*c2).

At alpha=1 this is `t_n + 0.577834373477*h`. Boundaries, rotating geometry,
material properties, and time-dependent forces must use that stage time.

Let `A=c1*v-alpha*c2`, `B=c2*alpha^2/2`, and `C=-c2*alpha^3/6`.
Besides the unit sums of the convex weights, the nonlinear third-order
conditions are

    c4*(A+v) + c5*(v5-alpha) = 1
    c4*(B+v*A) + c5*(alpha^2/2-alpha*v5) = 1/2
    c4*(C+v*B) + c5*(-alpha^3/6+v5*alpha^2/2) = 1/6
    c4*(C+v*A^2/2) + c5*(-alpha^3/6+v5*alpha^2/2) = 1/6.

The last two are distinct nonlinear order conditions. On alpha=0.7 through 2
at spacing 0.001, the largest residual was 5.22e-15 and the smallest coefficient
was 0.00143193. This is a sampled numerical check, not proof of positivity over
an arbitrary continuous domain.

| Independent calculation | Last measured temporal order |
| --- | ---: |
| Nonlinear `y'=y^2`, constant steps | 2.9869 |
| Nonlinear `y'=y^2`, alternating step weights 0.8/1.2 | 2.9856 |
| Fixed 32-cell PLM/upwind advection Fourier mode | 2.99999 |

The nonlinear tests use exact startup to isolate the multistep formula.
The advection test uses SSPRK3 startup and the exact exponential of its spatial
operator as the reference. These are not production Euler, self-gravity,
radiation, shock, or AMR results.

## SSP and performance qualifications

The dimensional SSP condition involves the **forward-Euler timestep limit**,
not simply a mesh spacing. A sufficient common-bound form is

    h <= h_FE / max(v,v5).

More precisely, `v*h` must satisfy the forward-Euler bound for the current and
intermediate states, and `v5*h` must satisfy the historical state's bound.
The spatial RHS must itself have the required forward-Euler property. Reusing
a RHS whose positivity limiter depended on a different interval requires care.

For alpha=1, `C_SSP=0.730599711442`. Against SSPRK3 with coefficient 1, assuming
RHS cost dominates and each method operates at its SSP limit, the ideal speed
ratio is `3*C_SSP/2 = 1.09590`: about 9.6% faster, or 8.75% less RHS work per
physical interval. At the same admissible timestep, two versus three RHS calls
gives one-third less RHS work. Neither calculation is an application benchmark.
Matched-error comparisons may give different results.

The supplied table matches the formulas, but every table row has negative c3;
those rows verify arithmetic and do not demonstrate SSP. The fitted alpha=1
case above has positive coefficients.

The empirical fit is not valid globally. For example, at alpha=10 it gives
`c5=-0.00368736` and `v5=-0.254645`. Its domain must be bounded and checked,
or coefficient selection must enforce positivity directly. Step doubling gives
alpha=0.5, the draft's problematic limit; use controlled step growth or a
one-step restart. Startup, failed-step rollback, and regridding also need
explicit history handling.

### Simpler coefficient evaluation

An independent rearrangement of the order conditions permits choosing `v`
directly, instead of evaluating the long radicals:

    c1 = alpha*(alpha+2*v)/(alpha+v)^2
    c4 = (alpha+1)^2*(alpha+v)^2 / (2*alpha^2*v*(alpha+2*v)^2)

Use the draft's c3, `c2=1-c1`, `c5=1-c3-c4`, and
`v5=(1-c4*(A+v)+alpha*c5)/c5`. Validate the denominators, positivity, and order
residuals during coefficient selection. This also avoids a removable singularity
in the printed c1 expression at `alpha=3*v`.
At alpha=1, an exact nonnegative choice is
`v=(1+sqrt(3))/2`, `c1=2/3`, `c3=0`, `c4=9*sqrt(3)-15`, and
`v5=(5+2*sqrt(3))/13`. Its SSP coefficient is `sqrt(3)-1`, approximately
0.732051. Zero convex weights are allowed in the SSP argument.

## Conservation and coupling

Define `beta=c4*c2+c5`, `w_n=c4*c1*v`, `w_stage=c4*v`, and `w_old=c5*v5`.
Expanding the two stages gives

    U_new-U_n = beta*(U_old-U_n)
                + h*(w_n*K_n + w_stage*K_stage + w_old*K_old).

If `I_old` denotes the previous accepted time-integrated numerical face flux,
then the new interval's effective integrated flux is

    I_new = h*(w_n*F_n + w_stage*F_stage + w_old*F_old) - beta*I_old.

Boundary ledgers, species transport, and AMR reflux must include this history
term. Counting only the two fresh flux evaluations is incorrect. Source
increments have corresponding history contributions. Shared interfaces require
consistent history and timestep coefficients.

Closed uniform-grid linear conservation follows from consistent conservative
RHS evaluations and state combinations. Gravitational potential energy depends
quadratically on density, so that argument does not establish gas-plus-gravity
energy conservation. The discrete mass-flux work and force stages need a
compatible derivation. Existing second-order shell closure cannot be assumed
to deliver third order.

The draft is explicit. Applying it directly to stiff radiation reinstates the
relaxation timestep limit. Substituting implicit solves for its forward-Euler
pieces does not automatically preserve the order conditions. A compatible
implicit-explicit or other source treatment must be derived and checked.

AMR also needs histories consistent with reflux and gravity corrections,
donor states at the intermediate stage time, and history migration/checkpoint
support. Third-order subcycling requires more than the current linear temporal
interpolation. Restarting after history invalidation is a reasonable first
prototype policy; its frequency belongs in performance accounting.

## Bounded next experiment

Implement an optional uniform-grid hydro transport prototype with a global
timestep, validated coefficients, correct stage times, SSPRK3 startup, and
conservative flux history. Require smooth nonlinear temporal convergence,
conservation, positivity/shock checks, and measured halo counts. Compare with
RK2 and SSPRK3 at matched error before committing to the coupled AMR extension.
The production second-order requirements remain unchanged.

This belongs to the established SSP two-step Runge–Kutta family; related
[Ketcheson, Gottlieb, and Macdonald work](https://arxiv.org/abs/1106.3626)
provides relevant order and SSP theory. Variable-step coefficient selection and
timestep restrictions also have established
[SSP multistep literature](https://arxiv.org/abs/1504.04107).
These references establish context, not novelty or validation of this draft's
particular variable-step formula.
