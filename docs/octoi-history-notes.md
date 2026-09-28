# OctoI historical EOS and radiation references

Inspected 2026-09-28 in the full local repository at
`/home/dmarce1/octotiger/src/octotiger`, across locally available refs.
Repository: https://github.com/STEllAR-GROUP/octotiger.
No historical implementation was restored or modified.

## Helmholtz

- `23b8238895e9ed2822f70cba63bcced95bcb9281` (2016-10-26): adds
  `src/helmholtz.cpp`, a double-precision table reader/evaluator (1081 x 401).
- `ccedddad9d62c27ce45e905b01ace3f9a225814b` (2016-11-08): removes it.
- `77ec6243bd846a5f9cdea868b206b32debbce777` (2022-11-01), Helmholtz merge:
  `src/stellar_eos/helmholtz.cpp`, reader/evaluator crediting Timmes.
- `37378fac6193d727f868d45a7df0e596d610207b` (2022-11-14): adds
  `helmholtz.table.dat` (433481 lines).

No table-generation driver or Fermi-Dirac integration routine was found in the
searched history. This does not rule out another repository or uncommitted work.

## VET / short-characteristic-style intensity sweeps: bookmarks only

All listed changes are authored by Dominic Marcello. Classification is based on
source inspection: octant sweeps interpolate neighboring upstream intensities,
apply exponential attenuation, accumulate angular moments, and normalize the
pressure tensor by energy. This identifies the relevant implementation without
asserting its numerical correctness or validation status.

- [587a1ecc](https://github.com/STEllAR-GROUP/octotiger/commit/587a1ecc5492e46305816f718aca045c8dbce274)
  (2016-06-02), Working on radiation module: intensity implementation.
- [1acb0242](https://github.com/STEllAR-GROUP/octotiger/commit/1acb024222f4388db2a3f0cda70d6509a5a06573)
  (2016-06-03), AMR for eddington tensor: AMR intensity exchange and tensor
  construction in `src/rad_grid.cpp`.
- [e3fd15ca](https://github.com/STEllAR-GROUP/octotiger/commit/e3fd15ca001ff2954959c64584a9c98811fd3c64)
  (2016-06-04), Added triangle decomp for sphere points.
- [1390008a](https://github.com/STEllAR-GROUP/octotiger/commit/1390008ade83c79d9bea0897f8cf3abd03363bb9)
  (2017-02-07), Adding radiation: `compute_intensity`, `accumulate_intensity`,
  sphere quadrature, and HPX/AMR boundary exchanges.
- `3deb2f01` (2017-02-09), Radiation test working in uni-grid: still contains
  the intensity sweep and accumulation. Commit title is historical, not a new
  verification of the test.
- [c4bbe7d7](https://github.com/STEllAR-GROUP/octotiger/commit/c4bbe7d7fe1312f1e55787c02193df37f119aed6)
  (2017-02-15), UPdating radiation transfer w/o intensity based eddington:
  deletes sphere-point files.
- [da7d22c7](https://github.com/STEllAR-GROUP/octotiger/commit/da7d22c78895b432ba6b513e032fee5d45e4cf31)
  (2017-02-15), Correcting github misuse: disables/removes intensity machinery.
