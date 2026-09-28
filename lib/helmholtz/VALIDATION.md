# Helmholtz translation validation

Environment: GCC/GFortran 13.3.0, Linux x86-64, C++17 standalone Release build.
AddressSanitizer, LeakSanitizer and UndefinedBehaviorSanitizer checks also passed
in a separate Debug build (outside the sandbox, which prevents leak tracing).
Source identities are recorded in [NOTICE.md](NOTICE.md). Commands are in
[README.md](README.md). No fast-math flags were used.

## Selected runtime data

The default C++ evaluator now loads the unchanged bundled Timmes table. CMake
verifies its source-record SHA-256. The default constructor was tested from the
build directory; normal builds exclude the experimental generator executable.

At rho=2e9 g/cm^3, T=1e4 K, Abar=12, Zbar=6 (D=1e9), the supplied table gives
`deept=-9356.262817670478` erg/g/K while total `cv=134805215.25472507` erg/g/K.
This state is included in the reference comparison. Thus supplied data do not
by themselves eliminate the cold-degenerate electron thermal-derivative issue.
Translation parity is distinct from physical accuracy.

## Implementation comparisons

- **Tabulated evaluator:** all 126 returned quantities at 421 states matched
  the original Fortran numerically exactly (53,046 comparisons, no NaNs).
  The fixed-seed sample covers the supplied table's density/temperature range,
  six compositions, table endpoints, and off-node states. Both evaluators use
  the exact user-supplied table in this comparison. This establishes translation
  agreement, not independent physical validation of that table.
- **Direct quadrature:** 18 electron/positron quantities at 25 states compared
  with the original direct Fortran EOS. The grid uses
  `D = [1e-6, 1, 1e3, 1e6, 1e9]` g/cm^3 and
  `T = [1e4, 1e6, 1e8, 1e9, 1e10]` K. The largest difference normalized by
  the documented dimensional/cancellation scales was `1.71e-14`.
  These are not uniformly tiny *relative errors in small derivatives*.
  Some derivative outputs are differences of much larger terms, and raw
  relative differences can be large near zero. Entropy uses `e/T` as its
  cancellation scale; derivative scales are documented in the comparison script.

## Analytic and generated-table checks

- Dilute, nonrelativistic electron pressure and energy recover `n k T` and
  `3 k T/(2 m_u)` to the `5e-4` test tolerance at the tested state.
- Cold nonrelativistic electron pressure agrees with
  `(3 pi^2)^(2/3) hbar^2 n_e^(5/3)/(5 m_e)` to the 1% test tolerance, which
  allows the finite-temperature and mildly relativistic corrections at that
  test state.
- A newly generated 21x21 table over `logD=[4,5]`, `logT=[7,8]` agrees with
  direct quadrature at the tested node/interior combinations within `2e-6`
  for electron pressure, specific energy, and entropy, and within `2e-5` for
  the density derivative of pressure.
- The three normalized Maxwell residuals are below `2e-12` in those tests.
- The pressure temperature derivative agrees with a centered perturbation
  within `1e-6`. Halving the mixed-derivative stencil step changes the tested
  pressure/energy by less than `1e-10` and electron heat capacity by less than
  `1e-8`.
- Invalid inputs, out-of-domain states, mismatched table sizes and truncated
  tables are rejected.

At `D=6131.90350185317 g/cm^3`, `T=1.653076520037604e9 K`, the electron
compressibility error against direct quadrature decreases with table refinement:

| Intervals per decade | Relative error in dP/dD |
|---|---:|
| 10 | 2.54816e-3 |
| 20 | 4.64812e-4 |
| 40 | 2.19332e-5 |

This is a targeted convergence check across the pair-production regime, not a
uniform accuracy bound on the full domain. Thermodynamic consistency alone does
not establish accuracy.

## Full table generation and unresolved thermal accuracy

The generated full-table SHA-256 on this build is
`304d4df4065ada67708c5e31c8a1f2dcd8d633906e20c93664974b868ee0ddbd`.

The standalone target generated all four table blocks on the full default
541x201 grid, from quadrature alone. The metadata reports **3,622 of 108,741
nodes with nonpositive electron heat capacity**, and zero nodes with negative
electron compressibility. These numerical sign failures prevent treating the
full generated grid as production validated.

The generator currently sets `F_TT=-s_T`, so its sign diagnostic measures
`Cv=T*s_T`. It does **not** measure the direct routine's separately calculated
`Cv=e_T`. At `D=1e9 g/cm^3`, `T=1e4 K`, the original Fortran gives
`T*s_T=-18712.5256 erg g^-1 K^-1` but `e_T=+137.707389 erg g^-1 K^-1`.
Our initial description of this as the original routine returning a negative
heat capacity was too broad: its entropy derivative is inconsistent with its
energy derivative at that state. The new table driver selected the less stable
expression. No silent floor is imposed on the generated derivatives.

Agreement with either original expression to a cancellation-scaled tolerance
is insufficient evidence of thermal accuracy. The controlled experiment below
separates precision, quadrature order and choice of derivative expression.

The supplied table and a newly generated table are not assumed identical.
The exact original table-building driver was absent from both archives, and
the two missing higher derivatives use the documented OctoII stencil. Preserve
this distinction in future scientific comparisons.

The pure C++ library also builds through the top-level OctoII CMake project with
HPX disabled and tests disabled. No application dynamics were changed or tested
as part of this standalone EOS step.

## Long-double investigation

Run the isolated experiment (C++17 compiler; `mpmath` for the higher-order rule):

```sh
python3 lib/helmholtz/tests/precision_experiment.py
```

It generates temporary C++ variants and CSV measurements under
`build-helmholtz/precision-experiment`. It changes no production EOS sources
or tables. All internal floating-point arithmetic, constants, quadrature
coefficients and literals are promoted. Timmes's full decimal quadrature
coefficients are recovered from the retained Fortran comments. The 40-point
rule uses coefficients computed at 60 decimal digits, with the same integration
intervals and integrands. Newton stopping tolerances are held fixed.

On this platform `double` has 53 mantissa bits and `long double` has 64, rather
than the 113 bits of IEEE binary128. At `D=1e9 g/cm^3`, `T=1e4 K`, with `Ye=1`
normalization, the C++ experiment finds:

| Arithmetic and quadrature | `T*s_T` | `e_T` |
|---|---:|---:|
| double, 20 points | -18712.8313 | 137.438953 |
| long double, 20 points | -1945.60483 | 137.851961 |
| long double, 40 points | 134.655190 | 137.844097 |

Heat capacities are in erg g^-1 K^-1. The leading relativistic Sommerfeld
low-temperature value is **137.843797** at this state; `kT/E_F=1.85e-7`, making
it a strong independent asymptotic check. In the 40-point long-double test,
`e_T` differs from this limit by about `2.18e-6` relatively, while `T*s_T`
still differs by about 2.31%. Increasing precision alone does not repair the
current entropy-derivative route.

The experiment also covers 15 other density/temperature combinations. The
coldest, densest tests still fail badly even with long double and 40 points;
positivity at one state is not a full-domain validation. A better candidate is
`F_TT=-e_T/T` with extended-precision integration, followed by checks of precision,
quadrature convergence and the cold-degenerate limit over the full intended
domain. This candidate has not been promoted to the production table builder.
