# Dimension executables and runtime problems

Validated with GCC 13.3, C++20, the serial backend, and Debug builds.

- One configuration builds octoII-1d, octoII-2d, and octoII-3d.
- Build and installation directories contain the relative octoII -> octoII-3d link.
- Application checks passed in all dimensions for four module configurations:
  all ON; HYDRO OFF; only HYDRO ON; and all OFF.
- Application checks cover all nine problem names, dimension restrictions,
  unavailable modules, unknown problems, rejected runtime dimension settings,
  layered INI/CLI problem selection, simulation startup/steps, Silo output,
  and the selected problem/dimension in verification JSON.
- Compilation databases confirm disabled numerical solver sources are omitted.
- The initial full regression run passed 408 of 414 checks. Five tests required
  migration from implicit build-selected physics to explicit runtime problem
  selection; those five passed after correction. The new 1D streaming-profile
  check also passed. The final focused CTest run passed all 12 selected checks.
- Additional hydro-disabled checks passed: 37 option tests, three Silo tests,
  ten AMR tests and the corrected adaptive-gravity test. Two density/hydro AMR
  tests correctly skipped when the test's selected problem was radiation-only.

Remaining validation limits:

- The long RayleighTaylor.PerturbationGrowsUnderGravity regression exceeded its
  existing 300-second timeout in Debug. It did not complete; its physics result
  is not claimed as validated. The other four Rayleigh-Taylor checks passed.
- HPX was unavailable in this environment. Distributed execution and HPX
  serialization were not run. Config now serializes the selected problem name,
  and the existing HPX test registrations remain available per dimension.

Use `ctest --test-dir release --output-on-failure` for the installed toolchain.
`python3 tests/run_matrix.py --module-matrix` exercises all eight module
combinations, with all three dimensions in each configuration.
