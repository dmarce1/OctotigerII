# Timmes EOS attribution and source record

The original Fortran and the scientific algorithms translated here are the work
of **F. X. (Frank) Timmes** and the authors cited below. The Helmholtz interpolation
method is credited to **F. X. Timmes and F. Douglas Swesty**. Original Fortran
statements and explanatory comments are preserved in the C++ implementation.
The translation, immutable C++ interface, build integration, table-building
driver, validation code, and derivative-completion stencil are OctotigerII work.
They are not presented as original EOS physics or as software endorsed by Timmes.

References:

- Timmes, F. X., & Arnett, D. (1999), *The Accuracy, Consistency, and Speed of
  Five Equations of State for Stellar Hydrodynamics*, ApJS **125**, 277.
  DOI: [10.1086/313271](https://doi.org/10.1086/313271).
- Timmes, F. X., & Swesty, F. D. (2000), *The Accuracy, Consistency, and Speed of
  an Electron–Positron Equation of State Based on Table Interpolation of the
  Helmholtz Free Energy*, ApJS **126**, 501.
  DOI: [10.1086/313304](https://doi.org/10.1086/313304).
- Aparicio, J. M. (1998), *A Simple and Accurate Method for the Calculation of
  Generalized Fermi Functions*, ApJS **117**, 627. This is the integration method cited
  by Timmes's `dfermi` routine (which gives page 632 in its comment).
  DOI: [10.1086/313121](https://doi.org/10.1086/313121).
- The Coulomb terms retain Timmes's attribution to Yakovlev & Shalybkov (1989).

Sources:

- User-supplied `helmholtz.tbz`, originally distributed through
  [Timmes's stellar EOS page](https://cococubed.com/code_pages/eos.shtml).
- `timmes.tbz` retrieved from
  [Timmes's direct EOS distribution](https://cococubed.com/codes/eos/timmes.tbz)
  during this implementation. The attached Helmholtz archive has an evaluator
  and a table, but no table-generation driver. Timmes's page states that the
  Helmholtz table is calculated from the direct Timmes EOS.

The bundled `data/helm_table.dat` is an unchanged copy of the supplied table,
credited to Timmes. CMake checks its digest before building.

SHA-256 identifiers:

| Source | SHA-256 |
|---|---|
| Supplied `helmholtz.tbz` | `fc55ca3b188598ed19f9dbf63bacf1033676a22d2065f420b67375a493b91eeb` |
| `helmholtz.f90` | `e23a6603e8963301e27c3e98e54fc589c63f0e07873a95635bbd664d5f884501` |
| Supplied `helm_table.dat` | `c9a57c26c6fd2b2b378b9d5295ca1214022f6fec6289d038b47bf8c8938881a1` |
| Downloaded `timmes.tbz` | `e20c7d27e66c240486a3397a649c49673e33100284f0905cd6fa9893dbad30a9` |
| `eosfxt.f90` | `4c3e45924c00cff751885377c936cdd44793838c3b38db6ba4c5fec1c58710d3` |

The upstream web page requests citation of the relevant references when these
codes or modified versions are used in publications, and appropriate offers of
coauthorship. Neither inspected archive contains a separate explicit license
file. This notice records attribution and provenance; it does not assign a new
license to Timmes's original code. The project's own licensing must not be
mistaken for an upstream license grant.
