# GoogleTest validation — 2026-09-23

All **1172 CTest executions passed across 14 supported problem/dimension
builds**: 1144 GoogleTest case executions and 28 existing whole-process
example/application checks. Counts include the same shared tests exercised in
multiple compiled configurations; they are not counts of unique test definitions.
No cases were skipped in this serial matrix.

Environment: GCC 13.3.0, CMake 3.28.3, GoogleTest 1.14.0, C++20, Release,
`OCTOTIGERII_WITH_HPX=OFF`. GoogleTest assertions are active with `NDEBUG`.
The per-build JUnit reports are in [testing-results/serial](testing-results/serial/).

| Problem | Dimension | CTest executions | Result |
| --- | ---: | ---: | --- |
| collapse | 3 | 96 | Passed |
| gravity-gaussian | 3 | 80 | Passed |
| gravity-sphere | 3 | 80 | Passed |
| kelvin-helmholtz | 2 | 81 | Passed |
| kelvin-helmholtz | 3 | 81 | Passed |
| radiation-pulse | 1 | 84 | Passed |
| radiation-pulse | 2 | 84 | Passed |
| radiation-pulse | 3 | 84 | Passed |
| sod | 1 | 82 | Passed |
| sod | 2 | 82 | Passed |
| sod | 3 | 82 | Passed |
| streaming | 1 | 86 | Passed |
| streaming | 2 | 85 | Passed |
| streaming | 3 | 85 | Passed |
| **Total** | | **1172** | **Passed** |

## Findings

The new negative-input tests caught a production validation defect:
`HydroSystem` accepted positive infinity for gamma and either positivity floor.
The constructor now requires finite gamma greater than one and finite positive
floors. All hydro configurations passed the regression after this change.
No numerical accuracy or convergence tolerance in the migrated solver checks
was relaxed. A new time-metadata test uses an eight-epsilon relative tolerance
for subtracting successive floating-point times rather than bitwise equality.

## Additional build-path checks

- Installed GoogleTest discovery and CTest's individual case discovery passed.
- The offline FetchContent source override was configured using an unpacked
  GoogleTest tree; its math test executable built and passed all three cases.
- Disabling tests configured successfully without looking up GoogleTest.
- A missing GoogleTest with automatic fetching disabled produced the intended
  actionable configuration error.
- Python syntax, build-script shell syntax, and patch whitespace checks passed.

## Scope

HPX was not installed in this validation environment. HPX initialization,
serialization, zero-copy archive checks and the registered 2-/3-locality TCP
tests were adapted or retained but were **not compiled or executed here**.
Serial partitioned-gravity checks cover field decomposition and publication,
not actual HPX thread scheduling, network behavior, or cluster performance.
APEX/PAPI-enabled execution and a fresh online GoogleTest download were not
exercised. The existing profiling configuration was preserved.

For commands to reproduce the matrix and execute the HPX tests, see
[testing.md](testing.md).
