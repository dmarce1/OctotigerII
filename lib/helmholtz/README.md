# Standalone Helmholtz EOS library

This is a C++ translation of **F. X. (Frank) Timmes's** Helmholtz evaluator and
of the direct electron–positron routines used to construct its table. See
[NOTICE.md](NOTICE.md) for scientific credit, source hashes, and references.
Original Fortran is retained as comments in the C++ source.

The library is `octotigerII_helmholtz`, with CMake alias
`OctotigerII::helmholtz`. It has separate translation units for the tabulated
EOS (`helmholtz.cpp`), direct Fermi–Dirac integration (`direct.cpp`), and table
I/O/generation (`table.cpp`). It requires C++17 and the standard library.
There is no Fortran runtime dependency. Fortran is only needed for optional
reference comparisons. OctoII's main CMake build also includes these targets.

The library exposes the forward `(rho,T,Abar,Zbar)` EOS. OctoII now has an
optional [hydro adapter](../../docs/helmholtz-hydro.md) with temperature inversion,
species-dependent composition, entropy dual energy, and gas-only radiation
coupling. The library remains independently buildable.

## Build and test

From the repository root:

```sh
cmake -S lib/helmholtz -B build-helmholtz -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build-helmholtz -j 4
ctest --test-dir build-helmholtz --output-on-failure
```

The production library does not depend on HPX or introduce any threads, locks,
lazy initialization, or shared mutable caches. Construct an `Eos` explicitly
before parallel evaluation; evaluations use stack-local work and immutable
shared tables. The offline generator is serial.

## Default table

`Eos eos;` loads Timmes's precomputed `data/helm_table.dat`, copied unchanged
from the user-supplied archive. CMake verifies its SHA-256 against NOTICE.md.
The default path is absolute and tied to this source checkout; when relocating
an application, ship the table and pass its path explicitly to `Eos`.
No table generation or Fortran compiler is required for normal use. The supplied
table also has a negative electron thermal derivative at a tested cold-degenerate
state; see VALIDATION.md for the measured value and scope.

## Experimental table generation

```sh
cmake --build build-helmholtz --target helmholtz_table
```

This computes `build-helmholtz/helm_table.experimental.dat` from direct Fermi–Dirac integrals,
without reading the supplied table. The default is 541 density nodes by 201
temperature nodes, covering `log10(rho*Ye) = [-12,15]` and
`log10(T) = [3,13]`, at 20 intervals per decade. The four data blocks match
Timmes's ASCII format. A `.meta.json` sidecar records the grid, method,
derivative step, source hash and sign diagnostics. Experimental build products are ignored by Git and never replace the bundled
Timmes table. The generator executable is excluded from normal builds.

**The full default table is a compatibility/research artifact, not a validated
production table.** Our driver uses `F_TT=-s_T`; the original direct routine
can lose accuracy in that entropy derivative at cold, strongly degenerate
states. Its separate energy derivative `e_T` is better conditioned in the
tested example. Our full generation found 3,622 of 108,741 nodes with
nonpositive heat capacity inferred from `T*s_T`. The tool reports this
in both its output and metadata. Successful process completion means generation
finished, not that the entire domain is scientifically validated. Positive
sign checks alone would not establish accuracy either.

A smaller table in the domain exercised by the basic thermodynamic tests:

```sh
cmake --build build-helmholtz --target helmholtz_build_table
build-helmholtz/helmholtz_build_table build-helmholtz/helm-warm.dat 21 21 4 5 7 8
```

Arguments are:

```text
OUTPUT [ND NT LOG_D_MIN LOG_D_MAX LOG_T_MIN LOG_T_MAX [DERIVATIVE_STEP]]
```

`D = rho*Ye`, not the physical density unless `Ye=1`. Custom tables must be
loaded with their matching `Grid`; the legacy ASCII format does not embed grid
coordinates. The reader checks finite entries, exact entry count, and bounds,
but cannot infer grid limits from this format. Keep the metadata with the data.

### Table construction method

Timmes's EOS web page identifies the direct Timmes EOS as the source of the
Helmholtz free-energy table. His supplied archive has no generator, and the
separately distributed direct EOS also has no table-writing driver. This
implementation therefore translates his numerical kernels and supplies an
explicit OctoII driver; it does not claim bytewise reproduction of his original
unpublished table-generation process.

The direct translation includes `xneroot`, `etages`, `dfermi`, `fdfunc1`,
`fdfunc2`, `dqleg020`, and `dqlag020`. It uses the same charge-neutrality Newton
solve, Aparicio integration intervals, and 20-point Gauss–Legendre/Laguerre
quadratures. The root solve has a finite-value failure path and avoids division
by zero at zero chemical potential. Direct calls use full ionization, no
ionization-potential contribution, and `Ye=1`, so the density is the table's
`D`. The source's constants and positron threshold are retained.

For the specific electron–positron free energy, the driver constructs

```text
F = e - T*s
F_D = P/D^2                  F_T = -s
F_DD = P_D/D^2 - 2P/D^3      F_TT = -s_T
F_DT = P_T/D^2
F_DDT = P_DT/D^2 - 2P_T/D^3
F_DTT = P_TT/D^2
```

The distributed direct code supplies second thermodynamic derivatives.
`F_DDTT` and the auxiliary pressure derivative table's `P_DDT` are completed
using fourth-order centered temperature differences of `F_DDT` and `P_DD`,
respectively. The relative step defaults to `1e-3`; tests also halve it. This
stencil is new driver code, explicitly distinguished from the translated
Fortran. All four blocks are written: free energy, pressure density derivative,
chemical potential, and total electron-plus-positron number density.

## Use the library

```cpp
#include <octotigerII/helmholtz/helmholtz.hpp>

using namespace octotigerII::helmholtz;
Eos eos; // Timmes's bundled precomputed table.
auto q = eos.evaluate(6.0e4, 2.718e7, 12.0, 6.0);
// rho [g/cm^3], T [K], mean nuclear mass Abar, mean nuclear charge Zbar.
// q.ptot [erg/cm^3], q.etot [erg/g], q.stot [erg/g/K].
// q.pgas, q.egas, q.sgas exclude photons.
```

Link `OctotigerII::helmholtz` from CMake. `Result` preserves Timmes's 126 output
names without `_row`, including component contributions, derivatives,
adiabatic indices, specific heats, sound speeds and consistency diagnostics.
Derivative suffixes `t`, `d`, `a`, and `z` mean differentiation with respect to
`T`, `rho`, `Abar`, and `Zbar`, holding the others fixed.

The tabulated `pele/eele/sele` already include positrons; the separate `pos`
outputs are zero, as in Timmes's evaluator. `xne` denotes total electrons plus
positrons, while `xnem` is their charge-neutrality difference. `cs` and `cs_gas`
are the source's **relativistic** sound speeds. A Newtonian hydro solver instead
needs `sqrt(Gamma1*P/rho)` for the selected components. Radiation contributions
must be accounted for once when coupling to separately evolved radiation.

Inputs must be finite with positive density, temperature, and `Abar`, and
`0 < Zbar <= Abar`. Out-of-table states throw; there is no extrapolation or
silent input clamping. Timmes's internal `Ye` floor, pressure-derivative floor,
and Coulomb suppression when total pressure or energy would become nonpositive
are preserved. This is a fully ionized EOS; partial ionization is not exposed.
Demo printing and the unused ion-table reader were not translated into runtime
features.

## Verification and limitations

[VALIDATION.md](VALIDATION.md) records measured results and their scope. The
ordinary CTest runs analytic limits, table construction, interpolation and
thermodynamic identities, finite-difference derivatives, stencil-step checks,
pair-region resolution convergence, and invalid input/table checks.

For the independent Fortran comparison, extract the two original archives into
separate source directories, then run:

```sh
python3 lib/helmholtz/tests/compare_fortran.py \
  --helmholtz-source /path/to/helmholtz \
  --timmes-source /path/to/timmes \
  --probe build-helmholtz/helmholtz_probe
```

This compiles the original routines with `gfortran`, replacing only their demo
programs with input/output drivers in a temporary directory. It changes neither
reference numerical routines nor the supplied table. Source archives and
`gfortran` are not required for normal library builds. No test downloads data.

Before production hydro use, establish a validated density/temperature/composition
envelope and check interpolation accuracy throughout that envelope. The direct
routine's cold-degenerate thermal issue concerns experimental table generation;
it is not a prerequisite for using Timmes's supplied table. The hydro adapter
adds inversion and general-EOS dual energy, sound speed, and conservation tests;
stellar equilibrium and merger validation remain separate. No merger accuracy
or performance claim is made by this isolated translation.
