# Bibliography {#bibliography}

Code comments use author–year citations linked to the records below. One author:
Levermore (1984); two: Skinner and Ostriker (2013); three or more: Harten et al.
(1983). Use the publication year, not the download date. If the same displayed
citation identifies more than one paper, assign a/b/c in title order and use the
suffix in both the displayed year and its stable anchor everywhere. No current
entries require that disambiguation.

This follows the author–year convention familiar from AASTeX. It does not require
a LaTeX or BibTeX installation. Each record gives the full author list, title,
publication, volume and pages or article number, and a DOI or primary-source URL.
The scope notes distinguish implemented methods from related design ideas.

## Biddiscombe et al. (2017) {#ref_biddiscombe2017}

Biddiscombe, J., Bikineev, A., Heller, T., & Kaiser, H. 2017,
“Zero Copy Serialization Using RMA in the HPX Distributed Task-Based Runtime,”
in *Proceedings of the 14th International Conference on Applied Computing*,
ed. P. Isaías & H. Weghorn (IADIS), 151–158. ISBN 978-989-8533-69-2.
[Publisher record](https://www.iadisportal.org/digital-library/zero-copy-serialization-using-rma-in-the-hpx-distributed-task-based-runtime).

Scope: HPX reference-backed serialization and buffer lifetimes. Our TCP tests
verify eligible pointer chunks and transferred values, not hardware RMA or a
copy-free network stack. Author order follows the publisher's record.

## Greengard and Rokhlin (1997) {#ref_greengard1997}

Greengard, L., & Rokhlin, V. 1997,
“A new version of the Fast Multipole Method for the Laplace equation in three dimensions,”
*Acta Numerica*, **6**, 229–269.
[doi:10.1017/S0962492900002725](https://doi.org/10.1017/S0962492900002725).

Scope: the diagonal plane-wave translation factorization, especially Section 7.
OctotigerII uses its own harmonic-equivalent Cartesian coefficient representation,
rotated quadrature, and finite-order implementation; it does not reproduce every
algorithmic choice in the paper.

## Harten et al. (1983) {#ref_harten1983}

Harten, A., Lax, P. D., & van Leer, B. 1983,
“On Upstream Differencing and Godunov-Type Schemes for Hyperbolic Conservation Laws,”
*SIAM Review*, **25**, 35–61.
[doi:10.1137/1025002](https://doi.org/10.1137/1025002).

Scope: the two-wave HLL numerical flux used for M1 transport and as a hydro fallback.

## Hu et al. (2013) {#ref_hu2013}

Hu, X. Y., Adams, N. A., & Shu, C.-W. 2013,
“Positivity-preserving method for high-order conservative schemes solving compressible Euler equations,”
*Journal of Computational Physics*, **242**, 169–180.
[doi:10.1016/j.jcp.2013.01.024](https://doi.org/10.1016/j.jcp.2013.01.024).

Scope: the high-order/low-order flux blending idea. The implementation uses
bisection of a common face coefficient against its own admissibility tests;
the radiation extension tests the M1 cone. Those details are application choices,
not a claim to implement the paper's Euler limiter verbatim.

## Kaiser et al. (2014) {#ref_kaiser2014}

Kaiser, H., Heller, T., Adelstein-Lelbach, B., Serio, A., & Fey, D. 2014,
“HPX: A Task Based Programming Model in a Global Address Space,”
in *Proceedings of the 8th International Conference on Partitioned Global Address Spaces*
(PGAS '14; New York: ACM), article 6, 1–11.
[doi:10.1145/2676870.2676883](https://doi.org/10.1145/2676870.2676883).
[Author manuscript](https://stellar-group.org/pubs/pgas14.pdf).

Scope: HPX components, actions, futures, and distributed execution. OctotigerII's
field partition component, two-bank publication policy, and locality work queues
are application designs. Author order follows the paper's title page.

## Levermore (1984) {#ref_levermore1984}

Levermore, C. D. 1984, “Relating Eddington factors to flux limiters,”
*Journal of Quantitative Spectroscopy and Radiative Transfer*, **31**, 149–160.
[doi:10.1016/0022-4073(84)90112-2](https://doi.org/10.1016/0022-4073(84)90112-2).

Scope: the M1 Eddington closure for a Lorentz-boosted isotropic distribution.

## Mohr et al. (2025) {#ref_mohr2025}

Mohr, P. J., Newell, D. B., Taylor, B. N., & Tiesinga, E. 2025,
“CODATA recommended values of the fundamental physical constants: 2022,”
*Reviews of Modern Physics*, **97**, 025002.
[doi:10.1103/RevModPhys.97.025002](https://doi.org/10.1103/RevModPhys.97.025002).
[NIST constants table](https://physics.nist.gov/cuu/Constants/Table/allascii.txt).

Scope: G, c, atomic mass unit, Boltzmann constant, and Planck constant, converted
to CGS. The radiation constant is derived as 8π⁵k_B⁴/(15h³c³). The reference year
is 2025; “2022” identifies the adjustment.

## Skinner and Ostriker (2013) {#ref_skinner2013}

Skinner, M. A., & Ostriker, E. C. 2013,
“A Two-moment Radiation Hydrodynamics Module in Athena Using a Time-explicit Godunov Method,”
*The Astrophysical Journal Supplement Series*, **206**, 21 (29 pp.).
[doi:10.1088/0067-0049/206/2/21](https://doi.org/10.1088/0067-0049/206/2/21).
[arXiv:1306.0010](https://arxiv.org/abs/1306.0010).

Scope: M1 closure, characteristic speeds, explicit transport, and the reduced
speed of light convention. Persistent radiation fields here are (E, F); the
calculation adapter converts to (E, F/c) using physical c. This project currently
implements uncoupled transport, not the paper's full implicit coupling or
radiation subcycling machinery.

## Sod (1978) {#ref_sod1978}

Sod, G. A. 1978,
“A survey of several finite difference methods for systems of nonlinear hyperbolic conservation laws,”
*Journal of Computational Physics*, **27**, 1–31.
[doi:10.1016/0021-9991(78)90023-2](https://doi.org/10.1016/0021-9991(78)90023-2).

Scope: the standard shock-tube initial discontinuity. OctotigerII's fixture uses
CGS values with the same density and pressure ratios.

## Strang (1968) {#ref_strang1968}

Strang, G. 1968, “On the Construction and Comparison of Difference Schemes,”
*SIAM Journal on Numerical Analysis*, **5**, 506–517.
[doi:10.1137/0705041](https://doi.org/10.1137/0705041).

Scope: symmetric source/transport/source composition in the gravity-coupled
run loop. This citation does not assert discrete gravitational energy conservation.

## Toro et al. (1994) {#ref_toro1994}

Toro, E. F., Spruce, M., & Speares, W. 1994,
“Restoration of the contact surface in the HLL-Riemann solver,”
*Shock Waves*, **4**, 25–34.
[doi:10.1007/BF01414629](https://doi.org/10.1007/BF01414629).

Scope: the contact-resolving HLLC hydro flux. Degenerate or inadmissible star
states fall back to HLL in this implementation.

## van Leer (1979) {#ref_vanleer1979}

van Leer, B. 1979,
“Towards the ultimate conservative difference scheme. V. A second-order sequel to Godunov's method,”
*Journal of Computational Physics*, **32**, 101–136.
[doi:10.1016/0021-9991(79)90145-1](https://doi.org/10.1016/0021-9991(79)90145-1).

Scope: piecewise-linear conservative reconstruction and limiting. The present
unsplit MUSCL–Hancock predictor combines all active directional divergences;
it is not the directionally split Lagrangian/remap algorithm of this paper.
