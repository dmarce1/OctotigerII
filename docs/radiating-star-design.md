# Opaque rotating gas–radiation star

Implementation and measured startup results are recorded in
[opaque-star validation](validation/opaque-star/README.md).

The first benchmark is a mechanically balanced, uniformly rotating star with
an optically thick interior, physical light speed, and **no photon heater**.
Radiation energy is confined at the computational boundary. Slow thermal
evolution and adjustment of the last few stellar cells are allowed. This is
not a claim of exact radiative equilibrium or a calibrated stellar-evolution
model.

The earlier leaking/heated equilibrium investigation is preserved in the
[research record](research/radiating-star-2026-09/README.md), including equations,
reference solutions, failed approaches, and the unfinished table initializer.

## Structure and thermodynamics

Use a monatomic ideal gas with `Rgas=kB/(mu m_u)` and a structural polytrope:

\[
 p_g=\rho\mathcal RT,\quad u_g=\frac32p_g,\quad E_0=aT^4,
 \qquad P=p_g+\frac{E_0}{3}=K\rho^{\gamma_s},
 \quad\gamma_s=1+\frac1n=\frac97.
\]

The structural index is `n=3.5`. The evolved gas still has adiabatic index 5/3;
the gas and radiation energies evolve separately. At specified central density
and gas-pressure fraction `beta_c=pg_c/P_c`,

\[
 T_c^3=\frac{3\rho_c\mathcal R(1-\beta_c)}{a\beta_c},\qquad
 K=\frac{\rho_c\mathcal RT_c}{\beta_c\rho_c^{\gamma_s}}.
\]

The local temperature is the unique positive root of the pressure equation.
The isolated multipole self-consistent-field calculation solves

\[
 \nabla^2\Phi=4\pi G\rho,\qquad
 H+\Phi-\frac12\Omega^2R^2=C,\qquad
 H=\int\frac{dP}{\rho}=(n+1)K\rho^{1/n}.
\]

`H` is a barotropic integral, not thermodynamic enthalpy. Entropy varies with
radius. This deliberately stably stratified family has `Gamma1>=4/3>9/7` and
outward-increasing entropy. Uniform rotation gives outward-increasing specific
angular momentum and no vertical rotational shear. These are adiabatic
stability properties; measured stability under the evolved radiation/hydro
equations remains a numerical test.

Defaults are `rho_c=1 g/cm³`, `beta_c=0.8`, `mu=0.6`, and
`Omega=0.2 sqrt(G M_spherical/R_spherical³)`. The spin normalization uses the
spherical reference, not an assumption that the final rotating radius is
unchanged. The numerical atmosphere has density `1e-12 rho_c`.

## Opacity, initial radiation, and boundaries

Set constant gray **absorption** opacity `kappa=0.34 cm²/g`, zero prescribed
photon power, and `chat=c`. This opacity is a controlled gray approximation;
sharing a numerical value with electron-scattering opacity does not turn the
implemented absorption/emission coupling into scattering.

For the thick interior, initialize the comoving diffusion flux from the SCF:

\[
 \Psi=\Phi-\frac12\Omega^2R^2,\quad b=\frac{dp_r}{dP},\qquad
 \mathbf F_0=-\frac{c}{3\kappa\rho}\nabla E_0
             =\frac c\kappa b\nabla\Psi.
\]

Limit this flux to a realizable value in the last thin surface layer. Do not
inject the divergence of this flux as heating: the structure helper's `heating`
return is useful only for estimating an initial thermal-adjustment time here.

Radiation moments are stored in the inertial frame. For azimuthal velocity
`v=Omega R e_phi`, the meridional diffusion flux satisfies `v·F0=0`. Define

\[
 f_0=\frac{|\mathbf F_0|}{cE_0},\quad
 A=\frac{E_0(1-f_0^2)}{1+\sqrt{4-3f_0^2}},\quad
 \beta_v^2=\frac{v^2}{c^2},\quad\gamma_v^2=\frac1{1-\beta_v^2}.
\]

Initialize the Lorentz-boosted M1 moments

\[
 E=\gamma_v^2(E_0+\beta_v^2 A),\qquad
 \mathbf F=\gamma_v\mathbf F_0+\gamma_v^2(E_0+A)\mathbf v.
\]

These moments give `Q=0` and `v·G=0` for the retained equal-opacity coupling,
with meridional force `G=rho kappa F0/(gamma_v c)`. Thus the initializer avoids
an artificial thermal transient from inconsistent moving LTE moments. The SCF
uses isotropic pressure and Newtonian mechanics: finite-flux M1 anisotropy,
the small `1/gamma_v` correction, spatial discretization, and the surface limiter
still prevent exact discrete stationarity. For the default star, `v_max/c`
is about `4.13e-4`.

`radiation.closedBoundary=true` sets the accepted shared radiation **energy**
face flux exactly to zero at domain walls, after the optical-depth correction
and through any realizability fallback. Interior diffusion and gas/radiation
energy exchange remain active. Hydro uses outflow boundaries and gravity stays
isolated. The retained radiation pressure flux can exert wall forces; momentum
and angular-momentum boundary ledgers remain necessary. Energy confinement does
not imply a torque-free wall.

The first configuration uses an inertial grid (`frame.omega=0`) while the star
rotates. The closed-energy boundary initially requires that choice: a moving
reflecting wall can perform work, so its treatment would require a separate
boundary derivation. This is a numerical insulating box, not a model of a real
stellar photosphere or a constraint applied at the deforming stellar surface.

## Scales and limitations

The 256-radial-cell SCF reference gives the following physical scales. The
[calculation and full output](research/radiating-star-2026-09/opaque_star_scales.txt)
are saved with the research record.

| Quantity | Default reference |
| --- | ---: |
| Mass | 7.8971e34 g (about 39.7 solar masses) |
| Equatorial / polar radius | 1.4491e12 / 1.4189e12 cm |
| Central temperature | 2.3950e7 K |
| Dynamical time `sqrt(Re³/GM)` | 2.4028e4 s |
| Rotation period | 7.3465e5 s |
| Central radial optical depth | 6.66e10 |
| Diffusion estimate `3 integral(kappa rho r dr)/c` | 9.09e11 s, or 3.78e7 dynamical times |

The diffusion integral is a scale estimate, not a computed global thermal
eigenmode. The formal initial local ratio `(ug+E0)/|div F0|` falls from millions
of dynamical times in the bulk to about two at `0.95 Re`. Only `2.76e-7` of the
mass lies outside that radius. In a 128-cell-wide box of width `2.4 Re`, this is
about three surface cells. Their cell optical depths at `0.95, 0.98, 0.99 Re`
are approximately `878, 31.7, 2.69`; the interior remains very thick.
This uniform-grid illustration is not the example's adaptive surface
resolution: the example derives its box from the nominal spherical radius and
uses density-based refinement, concentrating the finest cells in the interior.

The full light speed retains the radiation CFL cost even though the local
matter exchange is implicit. Brief smoke tests cannot establish many-dynamical-
time stability. Validation should measure density/radius drift, meridional
motion, thermal redistribution, and total gas+radiation+gravity energy with
accepted boundary fluxes. Compare mesh/time resolutions before interpreting a
small change as physical evolution. The exterior floor, Cartesian force error,
and surface-flux limiter are explicit numerical approximations to vary.

The reusable EOS/SCF tests check analytic thermodynamic derivatives, the
spherical Lane–Emden limit, rotational deformation, virial convergence, and the
diffusion gradient. The radiation source tests independently check boosted M1
thermal balance; boundary tests check energy confinement with the actual
transport, coupling, AMR, and conservation ledgers. These tests support the
construction without replacing the long evolution experiment.
