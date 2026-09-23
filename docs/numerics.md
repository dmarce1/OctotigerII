# Numerical conventions {#numerical_conventions}

## Stored quantities and ordering

Every stored scalar is a Boost.Units quantity in CGS. There is no configurable
code-unit scale. Columns hold one scalar field over contiguous interior cells;
`units::State` reconstructs the heterogeneous tuple used by a kernel.

| System | Column order | Units |
| --- | --- | --- |
| Hydro | ρ, ndim momentum components, E_gas | g cm⁻³; g cm⁻² s⁻¹; erg cm⁻³ |
| Radiation | Eᵣ, ndim flux components | erg cm⁻³; erg cm⁻² s⁻¹ |
| Gravity | Φ, gₓ, gᵧ, g_z | cm² s⁻²; cm s⁻² |

The x coordinate is contiguous. All coordinates and traversal bounds have exactly
`ndim` entries, selected at compile time. There are no inactive coordinate axes.
CGS integrated totals use unit transverse measure in 1D/2D. Persistent interiors
have no ghosts; transport uses two temporary ghost layers in each of the `ndim`
directions. `PatchView` joins direct interior views to compact ghost storage.
See [problem builds](../BUILDING.md) for dimension-dependent state sizes.

## Hydro and reconstruction

The ideal-gas closure is P = (γ − 1)(E_gas − |ρv|²/(2ρ)). Hydro reconstructs
(ρ, vₓ, vᵧ, v_z, P). See @ref ref_vanleer1979 "van Leer (1979)" for the
piecewise-linear reconstruction lineage. The current unsplit predictor uses
all active directional flux divergences at the half step, then corrects with
shared face fluxes. This is MUSCL–Hancock, not a multi-stage Runge–Kutta update.

The main hydro Riemann solver is HLLC, following @ref ref_toro1994 "Toro et al. (1994)".
Degenerate or inadmissible star states use HLL; see
@ref ref_harten1983 "Harten et al. (1983)". A common face blend toward a
local Lax–Friedrichs flux enforces admissibility. The blending approach is
related to @ref ref_hu2013 "Hu et al. (2013)"; our coefficient is found by bisection.

## Radiation: physical storage and calculation variables

Stored radiation state is **(E, F)**. `radiativeFlux(axis)` returns the actual
physical flux, with dimensions energy density times velocity. Snapshots,
communication buffers, and Silo all use this convention.

The M1 adapter forms Q = F/c temporarily. `M1` operates on (E, Q), and
`RadiationSystem` maps its fluxes back to the physical stored convention.
The stored-state transport flux in direction n is:

```text
energy component:       (ĉ/c) F_n
flux-vector component: c ĉ P_in
```

Thus the latter has units erg cm⁻¹ s⁻²: it transports a stored radiation flux.
Multiplication by Δt/Δx returns a change in (E, F) with the correct units.
The reconstruction variables remain (E, F/c). See
@ref ref_levermore1984 "Levermore (1984)" for the closure and
@ref ref_skinner2013 "Skinner and Ostriker (2013)" for the transport formulation.

Physical c always sets |F| ≤ cE and the conversion to Q; ĉ only changes transport
speeds. With ĉ/c = 1/4, a streaming state still stores F = cE. Tests check that
state, its transport flux, and the independently translated analytic profile.

## Gravity

`gravity::solve` treats each finest-level cell as a point mass, excludes self
interaction, and uses isolated boundaries. It compares with that discrete
problem, not a volume-integrated continuum Green function. Physical G and CGS
units are applied at the typed solver interface.

The far field uses the diagonal factorization of
@ref ref_greengard1997 "Greengard and Rokhlin (1997)". Multipoles are
harmonic-equivalent coefficients; retaining only z exponents 0 and 1 does not
mean discarding the other raw moments. Local coefficients are potential
derivatives in normalized coordinates. `fmm.hpp` documents shift signs and scales.

The run loop uses symmetric gravity kicks around transport, following the
composition principle of @ref ref_strang1968 "Strang (1968)". A kick adjusts total
gas energy by the kinetic-energy change, preserving its internal energy.
This alone does not establish global gas-plus-gravity energy conservation.

## Fixtures and constants

The shock tube follows @ref ref_sod1978 "Sod (1978)". The shear layer and Gaussian
fixtures are project regression problems; they do not claim to reproduce a
published benchmark specification exactly.

Constants follow the 2022 CODATA adjustment, published by
@ref ref_mohr2025 "Mohr et al. (2025)". Decimal literals are converted to CGS once
in `units/constants.hpp`. c, k_B, and h have exact SI defining values; their
CGS equivalents are rounded only by the floating-point representation.
