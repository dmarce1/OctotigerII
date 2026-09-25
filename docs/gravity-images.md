# Periodic and reflecting gravity

The 3D gravity solver follows `mesh.boundary` on every face. Existing all-outflow
inputs still select isolated gravity. Periodic faces must be paired. Analytic
gravity faces are not implemented. These are Newtonian 1/r interactions in three
physical dimensions even when only one or two directions are periodic.

| Faces on an axis | Gravitational source extension |
| --- | --- |
| outflow or inflow / outflow or inflow | No source images on that axis; free-space gravity. |
| periodic / periodic | Infinite image lattice with the domain length as its period. |
| reflecting / outflow or inflow | Original sources plus their same-sign mirror across the lower face. |
| outflow or inflow / reflecting | Original sources plus their same-sign mirror across the upper face. |
| reflecting / reflecting | Even extension with period twice the domain length; includes the entire infinite sequence of mirrors. |

Reflections in different directions compose, including edges and corners. Every
source copy carries the same mass. A reflected source moment changes sign for
an odd exponent in each reflected coordinate. The physical particle/cell self
term is omitted; its mirror and periodic images are included. A reflecting
**gravity** face here means even potential and zero normal gravitational field at
the plane. Hydro uses its existing reflecting rule separately.

When all three directions repeat (whether by periodic faces or paired mirrors),
the solver uses the usual mean-subtracted Poisson equation, with a uniform
background density `-M/V` in the extended cell. This is required for a compatible
periodic/closed-Neumann problem with positive mass. For the original uniform
cubic domain, this is also minus its volume-averaged density. With one or two
repeating directions, no density background is subtracted: the line/sheet zero
mode remains. Thus the potential may grow logarithmically or linearly in an open
direction. Outflow refers to a source-free exterior, not zero potential gradient.

## Examples

Build a 3D gravity problem as usual, then choose the faces at runtime. For a slab:

```sh
./release/octoII-3d --problem.name=gravity-sphere --mesh.periodic=off \
  --mesh.boundary.xLower=periodic --mesh.boundary.xUpper=periodic \
  --mesh.boundary.yLower=periodic --mesh.boundary.yUpper=periodic
```

For one periodic axis and a lower reflecting z face:

```sh
./release/octoII-3d --problem.name=gravity-sphere --mesh.periodic=off \
  --mesh.boundary.xLower=periodic --mesh.boundary.xUpper=periodic \
  --mesh.boundary.zLower=reflecting
```

`--mesh.periodic=on` selects three periodic axes. The six explicit face options
also work in INI files. Their existing precedence and serialization are retained.

## Interactions and acceptance

The Newtonian walk uses the nearest image only along repeating coordinates.
Reflected distances are measured to the actual reflected source center. Its
cached diagonal translation therefore depends on that image separation; source
moments have the reflection parity applied before the transform.

A separate walk supplies the Ewald correction `psi_periodic(r) + 1/|r|`, excluding
exactly that nearest image. Its acceptance distance is the nearest **remaining**
image distance. For a nearest-image displacement r and periods P_d in cell widths:

```
nextDistance² = min_d (|r|² + P_d² - 2 P_d |r_d|), over repeating axes d.
```

Both walks accept when `openingAngle * dimensionlessDistance > 1` (equivalently,
the physical distance times the opening angle exceeds the cell width). A multipole pair
that crosses a nearest-image selection plane is opened until that choice is
uniform across the source/target supports. This ensures independent acceptance
levels never subtract different images. Leaf point interactions resolve ties
with an antisymmetric displacement convention.

The scalar correction M2L is a dense `(p+1)²` by `(p+1)²` contraction, O(p⁴).
The auxiliary force local adds one target degree, giving a `(p+2)²` by `(p+1)²`
contraction. Both are cached by source order, local kind, separation, and lattice
periods measured in cell widths. The auxiliary gradient has degree p at both
ends of a pair, preserving action/reaction; see the
[expansion derivation](parallel-fmm.md#scalar-reciprocity-and-mutual-force).
Newtonian image
interactions retain the cached diagonal translations. Ewald initialization is
more expensive than subsequent applications. Caches are local to each process;
cache size grows with the distinct geometries encountered. The same source-only
exchange and disjoint-target HPX scheduling used by isolated gravity are retained.
No complete density field is gathered onto a single locality for the solve.

The existing `gravity.multipoleOrder=1..10` range and opening-angle range remain.
`Statistics::ewaldPairs` and `reflectedPairs` count directed image interactions;
for image solves, `directPairs` and `multipolePairs` also count directed work.
The legacy isolated pair counts are unchanged. APEX samples include these image
counters when profiling is enabled.

## Ewald kernels and traces

The real-space kernel uses erfc(alpha*r)/r. Its radial derivative recurrence is
`D[k+1] = ((2*k+1)*D[k] + 2^(k+1)*alpha^(2*k+1)*exp(-alpha²*r²)/sqrt(pi))/r²`.
The central term is evaluated directly as erf(alpha*r)/r using a nonsingular
Gaussian integral, including at r=0. No singular quantities are subtracted.

The reciprocal part uses the ordinary discrete Fourier series in 3P, the analytic
slab erfc expressions in 2P, and the Gaussian integral representation of the
incomplete-Bessel coefficients in 1P. The latter is evaluated by 64-point
Gauss-Legendre quadrature; Gaussian/Hermite recurrences generate all derivatives.
Sums use order-dependent exponential cutoffs. Kernel tests vary both alpha and
quadrature resolution through degree 20. This quadrature is for coefficient
setup, not an approximation of the infinite lattice by a fixed number of images.

Three-periodic gravity has `Laplacian(psi_correction) = -4*pi/V`. Harmonic-only
moment storage would lose its trace. Image source vectors therefore carry one
additional raw z²/2 moment, and local vectors carry an additional Laplacian.
M2M, dense M2L, and L2L transport these terms explicitly. Ordinary harmonic
reduction still applies to higher derivatives because that Laplacian is constant.

The 3P potential uses the mean-zero periodic Green function (with physical self
removed). The 2P gauge is the usual slab expression. The 1P finite-part gauge
uses `log(alpha*minimumPeriod)` in the zero mode; changing alpha leaves the
potential unchanged. These conventions fix otherwise arbitrary constants in
1P/2P. The GADGET-4 slab formula has an additional position-independent constant,
which has no effect on forces.

## Verification

The default sampled direct reference now sums the same physical/image sources
and evaluates the Ewald pair correction without multipole approximation. Its
kernel is shared with the solver, so separate independent kernel tests are also
provided. Isolated sphere/Gaussian continuum references are unavailable for image
boundaries; select `verification.gravityReference=direct` (the default).

`ewaldChecks` checks alpha/quadrature independence, potential/gradient/Hessian
consistency, the 1P force against a long direct sum, the 2P force and potential
against an independent slab Fourier solution, the 3P cubic lattice constant,
periodic seams, self limits, moment traces, image geometry, reflections,
convergence with expansion order, and uniform lattice force cancellation.
`gravityParallelChecks` exercises periodic/mixed/reflected problems through
partitioned storage, checks repeat solves and publication, and runs under the
existing two-/three-locality launcher.

```sh
ctest --test-dir <build-directory> --output-on-failure -R 'Ewald|GravityImages|GravityBoundaries|PartitionedGravity|gravity.distributed'
```

## References

* [GADGET-4 paper](https://wwwmpa.mpa-garching.mpg.de/gadget4/gadget4-code-paper.pdf),
  section 2.8 and equation (39): slab Ewald correction and nearest-image split.
* [GADGET-4 Ewald source](https://wwwmpa.mpa-garching.mpg.de/gadget4/doxygen/ewald_8cc_source.html):
  slab derivatives and regular central-image treatment.
* [Tornberg (2016), The Ewald sums for singly, doubly and triply periodic electrostatic systems](https://arxiv.org/abs/1404.3534):
  the mixed Fourier representation, with the gravitational sign and retained
  transverse zero mode for nonneutral mass distributions used here.
