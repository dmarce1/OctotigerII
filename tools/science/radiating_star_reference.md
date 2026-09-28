# Rotating gas–radiation reference generator

**Deferred research, 2026-09-27.** The active benchmark now uses an opaque,
source-free SCF star with physical light speed. This generator retains the
earlier leaking/heated construction for later investigation. The
[research record](../../docs/research/radiating-star-2026-09/README.md) preserves
the experiments, reference arrays, and limitations; no stable evolving rotating
result is claimed for this table-based model.

`generate_radiating_star_reference.py` constructs an axisymmetric reference with
self-consistent gravity, gas pressure, full M1 radiation pressure, and a
prescribed photon energy source. It is an offline initializer, not an extra
force or relaxation term in the evolution solver.

The current family fixes the structural core index to 3.5, central gas-pressure
fraction to 0.8, opacity coefficient to 100, and opacity cutoff density to 0.015.
It uses units `4*pi*G = rho_c = Rgas*T_c = 1`; radiation flux is expressed as
`F/c`. The opacity is a prescribed spatial function of the spherical reference
profile, `kappa(r) = 100 max(0, 1 - 0.015/rho_spherical(r))^2`. The equatorial
angular velocity is the requested spin fraction of `sqrt(M_spherical/R_spherical^3)`.
The required differential rotation elsewhere follows from radial force balance.

The opaque core obeys the gas barotrope induced by
`P_g + a*T^4/3 = 1.25*rho^(9/7)`. At the opacity cutoff, a transparent gas
polytrope matches gas pressure and its density derivative. Radiation then
satisfies source-free M1 transport independently of gas temperature; this avoids
a material wall or an isothermal gas atmosphere of infinite radius. Opacity is
zero before the gas surface. The photon heating is distributed through the
opaque region, rather than compactly supported at the center.

## Reproduction

Python with NumPy and SciPy is required. From the repository root:

```sh
python3 tools/science/generate_radiating_star_reference.py \
  --nr 193 --nm 65 --extent 4 --spin 0.1 \
  --output /tmp/radiating-star-reference.npz \
  --include /tmp/radiating-star-reference.inc
```

The `.npz` file retains the nonlinear unknowns and metadata for independent
checks. `--resume` accepts a state on the identical radial and angular grid.
`--static-only` is for diagnosing the meridional zeroth-order calculation and
cannot export a shipped include. The normal path also solves the toroidal
radiation momentum equation and the retained moving-gas thermal/force terms.
It uses the physical light speed for `rho_c=1 g/cm^3`, mean molecular weight
0.6, and central gas fraction 0.8: `c/sqrt(Rgas*T_c) = 520.3917149`.

## Discretization and checks

The equations are solved for corrections to the known spherical M1 reference:
`F_h(U) - F_h(U_spherical) = 0`. This consistent analytic-background
formulation preserves the spherical solution without adding any balancing
force to a running simulation. Gravity uses the spherical potential plus
isolated multipoles of the density contrast. Radial pressure derivatives are
second-order outward differences. The nonuniform radial mesh places the
opacity cutoff exactly on a node and concentrates cells around its narrow M1
critical layer. Central/equatorial derivatives use the even quadratic limit.

A small nonlinear residual is not sufficient validation: centered collocation
experiments produced tiny algebraic residuals but nonconvergent continuum
errors. The accepted outward-difference calculation was independently
differentiated with piecewise cubic interpolants. The spherical background and
the two-dimensional correction both require separate interpolation branches
at the cutoff. Inside the opaque region, temperature and radiation energy must
satisfy the retained thermal balance exactly after interpolation.

For the static meridional spin-0.1 sequence, doubling 97x33 to 193x65 nodes
gave these dimensionless RMS continuum residuals over the opaque region:

| Equation | 97x33 | 193x65 |
| --- | ---: | ---: |
| Radial radiation momentum | 1.313e-5 | 3.430e-6 |
| Polar radiation momentum | 5.464e-6 | 1.355e-6 |
| Vertical gas force balance | 4.592e-7 | 1.162e-7 |

The fine model has positive angular velocity squared and positive opaque-core
heating. A separate dense 241x129 near-cut sample gave minimum heating
`6.834e-6`, including points within `1e-8` of the cutoff. The retained
moving-equation solve has maximum algebraic residual `8.19e-10`; this is close
to the Poisson stencil's floating-point floor and much smaller than the
continuum discretization error. These are reference-construction results, not
a claim that the evolving Cartesian/AMR calculation preserves equilibrium to
machine precision. Evolution stability, drift, conservation, and resolution
checks remain separate acceptance tests. Escaping photon angular momentum is
retained; no compensating gas torque is applied.

## Include format

The include is one C++ raw string literal. Its contents are:

1. `OCTOII_RADIATING_STAR 1`.
2. Seven metadata numbers: index, central gas fraction, opacity coefficient,
   cutoff density, spin fraction, outer-radius/cutoff-radius, physical
   dimensionless light speed.
3. Core radial count, exterior radial count, angular count, field count (6).
4. The mu axis, core radius axis, and exterior radius axis.
5. Core then exterior numeric blocks, with radius as the outer index.

Each node contains six fields, each followed by four numbers:
`value, d/dr, d/dmu, d²/(dr dmu)`. The fields are corrections to gas enthalpy,
`log(E/E_spherical)`, radial flux factor, the polar factor `q`, potential, and
the toroidal flux variable `W`. Here
`F_theta/c = -sqrt(1-mu²) E q` and
`F_phi/c = sqrt(1-mu²) W`. Enthalpy and potential use the physical
dimensionless units above, not the nonlinear solver's central-enthalpy scaling.
The cutoff node appears in both pieces with exactly equal values and separate
one-sided derivatives. No boundary condition sets `q` or `W` to zero on the
polar axis.
