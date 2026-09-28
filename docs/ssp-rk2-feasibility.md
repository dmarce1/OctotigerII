# SSP-RK2 with shared stage data: coupling investigation

2026-09-28. The user authorized an implementation only if the original
accuracy could be retained. **No SSP-RK2 production implementation was made.**
The straightforward source-coupling candidates below do not preserve the
original requirements. This is not an impossibility result for all IMEX or
SSP methods.

## What works for transport and gravity

Heun's two spatial evaluations give second-order fixed-mesh transport with
two initial-state layers at each stage. Each stage can share its gas halo with
radiation material calculations. For source-free Heun, only the original
interior state must survive the first stage; its initial ghosts are unnecessary
for the final convex average.

For gravity, form a first-order endpoint predictor from the canonical numerical
RHS, including the initial flux reflux and discrete gravity energy derivative.
Use the average of initial and endpoint numerical fluxes for accepted transport.
The existing endpoint force quadrature can retain second order without an
additional gravity solve, and the energy work must use the same accepted mass
flux. The accepted transport base remains the original state: adding a gravity
source inside the predictor and subsequently applying the full SSP state average
plus the existing gravity correction would count part of that source twice.
This consistency argument does not establish SSP for the complete gravity map.

## Source insertion after each Euler substep fails temporal order

Consider source-only decay `q'=-q`. Even with an EXACT source solve for each
full Euler substep, followed by Heun's convex average, the step becomes

    q_new = (1 + exp(-2*h))/2 * q_old.

Its expansion is `1-h+h^2+O(h^3)` instead of `1-h+h^2/2+O(h^3)`.
It is first order globally, and in the infinitely stiff limit retains half
the initial disturbance instead of relaxing it away. The independent check
at 20/40/80/160 steps gives orders 1.0086, 1.0044, 1.0022.

## Standard IMEX-SSP2(2,2,2) fails the required forced balance

The published second-order method couples explicit Heun with an L-stable
implicit tableau, using `gamma=1-1/sqrt(2)`. Its general order does not imply
preservation of a particular stiff transport-source balance.

Use the scalar radiation-momentum model

    q' = d - lambda*q,

where `d` is a constant pressure-gradient drive. Begin at the exact equilibrium
`q=d/lambda`. Treat `d` explicitly and relaxation implicitly. The stages are

    q1 = q/(1+gamma*h*lambda)
    q2 = [q+h*d-(1-2*gamma)*h*lambda*q1]/(1+gamma*h*lambda)
    qnew = q+h*d-(h*lambda/2)*(q1+q2).

The true solution is constant. Normalize it to 1:

| lambda*h | Exact / forced SDIRK | IMEX-SSP2(2,2,2) |
| --- | ---: | ---: |
| 1 | 1 | 0.963711 |
| 10 | 1 | -2.92966 |
| 100 | 1 | -65.1187 |
| 10,000 | 1 | -7065.24 |

At fixed `h,d` with `lambda` tending to infinity, this IMEX update approaches
`-h*d/sqrt(2)` instead of the equilibrium `d/lambda` tending to zero. This is
a failure for prepared equilibrium data, not merely an unresolved initial
transient. A negative flux component is not itself inadmissible, but the wrong
sign and O(h) magnitude violate this force balance.

The existing forced SDIRK includes pressure driving in both implicit stages
and preserves this constant equilibrium. The production test
`RadiationCoupling.ForcedStiffDiffusionBalancePreservesFlux` expresses the
corresponding coupled requirement; it was not weakened or rerun for this audit.

The scalar argument does not include the full code's opacity-dependent spatial
flux correction. It establishes failure of the candidate forced relaxation
update; it is not a measured whole-code diffusion-rate error.

## AMR and alternative coupling work

A full coarse-step Euler predictor applied to all fine donors can exceed their
forward-Euler CFL, even when the old half-step predictor was admissible. A
global donor bound would reduce the available subcycling benefit. Alternatively,
the coarse endpoint flux could be evaluated after fine subcycling completes,
using completed fine endpoint data. That can retain second-order consistency,
but requires reorganizing the current gravity/AMR schedule and verifying its
conservation identities.

A full forced predictor followed by a final forced solve driven by averaged
transport can satisfy resolved second-order Taylor conditions. It is not the
simple SSP convex-combination update. Averaging an unrelaxed initial transport
flux into the accepted flux can also spoil an unresolved stiff transient's
effect on conserved variables. It needs a separate stiff-limit and positivity
analysis before implementation.

Other IMEX formulations remain possible. The investigated SSP-RK2 replacements
have not met the conditions for replacing the original source-aware numerical
midpoint method. Sharing data between that method's existing stages remains a
separate optimization that avoids changing these numerical balances.

## Reproduction and sources

Run [ssp-rk2-source-audit.py](validation/local-predictor/ssp-rk2-source-audit.py).
It uses only Python's standard library and computes the scalar results above.
No production build or convergence run was performed for this investigation.

- [Pareschi and Russo, IMEX Runge–Kutta schemes](https://arxiv.org/abs/1009.2757):
  established SSP explicit / L-stable implicit methods and their relaxation
  limits. Their general AP results do not establish preservation of the
  constant forced balance tested here.
- [Boscarino and Russo review, appendix](https://link.springer.com/article/10.1007/s40324-024-00351-x):
  the SSP2-IMEX(2,2,2) tableau used in the calculation.
