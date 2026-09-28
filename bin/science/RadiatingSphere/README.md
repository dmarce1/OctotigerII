# Spherical radiative equilibrium control

This problem evolves gas, full M1 radiation, and isolated self-gravity. It is
the nonrotating control for the rotating-star work, not a rotating equilibrium
or a demonstrated long-term stable star. Start from
[`examples/radiating-sphere.ini`](../../../examples/radiating-sphere.ini).

The initializer integrates the spherical hydrostatic, gravitational, and M1
momentum equations together. In the opaque core the structural relation is

```
pg + a T^4/3 = K rho^(1 + 1/n),  pg = rho Rgas T,  E = a T^4.
```

The default `n=3.5` describes the structure; the evolving gas still has
`gamma=5/3`. The radiation pressure used in mechanical balance is the full M1
tensor, so the structural sum above is not identified with its radial stress.

This is a deliberately constructed opacity and heating test. Let `rho_ref(r)`
be the initial reference density, `alpha` the structural length scale, and
`rho_cut` the selected matching density. The prescribed opacity is

```
kappa_star(r) = kappa0 max(0, 1 - rho_cut/rho_ref(r))^2,
kappa0 rho_c alpha = radiatingStar.opticalDepthScale.
```

It is held fixed as a function of position throughout evolution. It is **not**
recomputed from the evolving density. The outer gas envelope is transparent
and follows a gas polytrope matched in pressure and pressure derivative at the
transition. It has a finite zero-density surface; radiation continues through
it and into a vacuum M1 solution. This avoids extrapolating LTE to a cold
surface carrying finite luminosity, which would violate `|F| <= c E`.

The fixed photon heater is `j_star(r) = div F_star`. It is positive in the
opaque core, includes a distributed outer tail, and is zero in the transparent
envelope. It is not a point source or a compact central heater. The source
adds radiation energy with no imposed gas force or momentum, and does not
restore evolving fields to their initial values. At reduced light speed the
code injects `(chat/c) j_star` into radiation, consistently with its equations.
The source-corrected energy budget is recorded in the diagnostic CSV.

The default reference has `r_transition/alpha = 3.5708999`,
`R_surface/alpha = 6.3941644`, central radial optical depth about `286.41`, and
transition flux factor about `0.71534`. A constant numerical gas floor outside
the material surface is explicitly an approximation, not part of the exact
atmosphere. Its mass and influence must decrease with the floor. Transport
uses outflow boundaries with isolated gravity; boundary placement is also a
convergence parameter.

`radiation.enabled` should remain off: this problem's manifest already enables
radiation and its initializer supplies explicit M1 moments. That separate
option adds LTE radiation to problems which otherwise initialize gas only.
The size is derived from central density, gas pressure fraction, composition,
and structural index; specifying `star.radius` simultaneously is rejected.

Validation is in `radiatingStarAtmosphereChecks` and `radiatingSphereChecks`:
independent differential force/Poisson checks, positive heating versus escaping
luminosity, reference integration convergence, and mesh comparison of actual
coupled evolution with energy and mass ledgers. Short evolution checks do not
establish thermal, convective, or global stability. A rotating result also
requires the two-dimensional balances in
[`docs/radiating-star-design.md`](../../../docs/radiating-star-design.md).
