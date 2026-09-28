# Rotating radiating white dwarf

`radiatingWhiteDwarf` is a carbon/oxygen white-dwarf reference with cold,
relativistically degenerate electrons, a small thermal-ion component,
uniform rotation, isolated Newtonian self-gravity, and full-speed M1 radiation.
Use [`examples/radiating-white-dwarf.ini`](../../../examples/radiating-white-dwarf.ini).

## Governing model

For mass density `rho` and electron molecular weight `mu_e=2`, the electron
number density is `n_e=rho/(mu_e m_u)` and the Fermi momentum parameter is
`x=hbar (3 pi² n_e)^(1/3)/(m_e c)`. The hydro EOS is

    P = P_deg(rho) + (gamma_i-1) u_i,
    E = u_deg(rho) + u_i + rho |v|²/2,
    c_s² = dP_deg/drho + gamma_i P_i/rho.

`P_deg` is the zero-temperature Fermi-electron pressure; its energy satisfies
`u_deg=rho H_deg-P_deg`, with `H_deg=(m_e c²/(mu_e m_u)) [sqrt(1+x²)-1]`.
This is the ideal degenerate-electron WD approximation used in
[Mathew and Nandy (2017)](https://arxiv.org/abs/1401.0819).
The thermal ions have `gamma_i=5/3` and default mean mass per ion 14 atomic
mass units. The dual-energy auxiliary tracks *thermal-ion* energy, since it
is small compared with electron degeneracy energy. The finite-volume solver
uses the EOS derivative in its wave speeds, not `sqrt(gamma P/rho)`.

The reference uses `P_i=epsilon P_deg` with `epsilon=1e-4`; hence its
hydrostatic barotrope is `P_ref=(1+epsilon)P_deg`. A spherical enthalpy
integration establishes the mass, radius, and breakup frequency. An
axisymmetric self-consistent-field solve then enforces

    H_ref + Phi - Omega² R²/2 = constant,
    nabla² Phi = 4 pi G rho.

The default `Omega` is 0.2 times the spherical reference breakup rate.
The initialized radiation field has LTE comoving energy `a T⁴`, the interior
diffusion flux limited to a realizable M1 state, and an azimuthal Lorentz
boost. It evolves with full-speed M1 and gray absorption. No photon heater
is imposed. The numerical radiation boundary is insulating and the gas
boundary is outflow, as in the existing `radiatingStar` benchmark.

## Reference values

For central density `1e9 g/cm³`, `mu_e=2`, and the defaults above, the
reference constructor reports:

| Quantity | Value |
| --- | ---: |
| Mass | `2.7354e33 g = 1.3756 M_sun` |
| Equatorial radius | `2.5116e8 cm = 2512 km` |
| Polar radius | `2.4558e8 cm = 2456 km` |
| Central sound speed | `8.0808e8 cm/s = 0.02695 c` |
| Central ion temperature | `8.186e6 K` |
| Angular velocity | `0.69396 s^-1` |
| Rotation period | `9.054 s` |
| Equatorial dynamical time `sqrt(R³/GM)` | `0.2946 s` |

With 256 radial cells and 12 multipoles, the SCF density and Bernoulli
residuals are `7.4e-10` and `3.7e-10`; the discretized scalar virial residual
is `5.7e-5`. Those are reference-solver diagnostics, not evidence that a
Cartesian evolution preserves the star for many dynamical times.

At the center, LTE radiation pressure is only about `2.3e-14` of degeneracy
pressure. The M1 field is dynamically weak for a cool white dwarf. A physical
opacity would require composition-, density-, and temperature-dependent
absorption/scattering and conduction; the fixed gray coefficient here is a
controlled transport approximation. Electron captures, nuclear burning,
general relativity, and a resolved photosphere are not included. The example's
`0.03 s` stop time is a startup run of about one tenth of a dynamical time.
Electron captures in particular limit high-density WD stability, as discussed
by [Chamel et al. (2021)](https://arxiv.org/abs/2110.11038).
