# Source reorganization validation — 2026-09-26

This change separates existing responsibilities without changing the numerical
method. See the [implemented plan](../source-reorganization-plan.md),
[current inventory](../source-inventory.md), and
[highlighted change map](../source-inventory-proposed.md).

## Scope and source equivalence

- Removed the unused `octotigerII/subgrid/exchange.hpp`. No maintained source,
  test, or build definition included it; installation rules do not export these
  headers. `fluxPacket.hpp` remains the single `FieldFluxPacket` definition.
- Extracted AMR selection into `regridSelection.hpp/.cpp`. Selection and
  geometry helper bodies were copied unchanged; startup, runtime regridding,
  and AMR tests now include the selection header explicitly.
- Split the 1,992-line runtime into its lifecycle facade and private execution,
  stage, transport, gravity, and regrid files. The public `Runtime` header is
  unchanged.
- Compared the split to the complete source snapshot taken immediately before
  this reorganization. All 42 extracted method bodies match after whitespace
  normalization. An independent token comparison also matched all 382 brace
  bodies inside the implementation, ignoring comments and whitespace; only
  their enclosing namespace/class outlines changed.
- Shared HPX action argument types now have one named namespace across all
  translation units. The five action identifiers and `LocalExecutor` component
  identity are unchanged. Action declarations live in the private header;
  component and action registrations remain in exactly one implementation file.
- Reviewed exception draining, rollback order, compile guards, constructor
  initialization, and final field-bank publication. HPX mutexes remain selected
  by `OCTOTIGERII_WITH_HPX`; the global gravity fallback still calls the public
  methods without first taking their lock.
- The user's separate `examples/polytrope.ini` edit matches the snapshot byte
  for byte. Existing rotating-frame work is retained.

The snapshot is local at `/tmp/octotigerII-before-source-reorganization`.
It includes the uncommitted source, inputs, Git status, and pre-refactor diff;
it is a comparison aid, not an installed dependency.

## Builds and regression tests

All 1D, 2D, and 3D applications and test targets built in both the serial and
HPX release configurations.

| Configuration | Result |
|---|---|
| Serial, full CTest suite | 625 passed, 3 expected skips, 0 failed; 563.44 s |
| HPX, full single-locality CTest suite | 634 passed, 3 expected skips, 0 failed; 769.16 s |
| HPX, two-locality rotating gravity | 4 passed, 0 failed; 152.03 s |
| HPX 3D, profiling disabled | Clean build and rotating-star smoke run passed |

Each suite's three skips are the dimension-inapplicable rotating boundary and
rotating Silo cases in 1D. The rotating-gravity temporal convergence test ran
and passed in both suites (429.46 s serial, 365.06 s HPX); it was not filtered
out of these full suites.

Commands used for the full builds and suites:

```sh
cmake --build /tmp/octotigerII-conservation-serial --parallel 3
ctest --test-dir /tmp/octotigerII-conservation-serial -j3 --output-on-failure
cmake --build release --parallel 2
ctest --test-dir release -LE distributed -j2 --output-on-failure
```

The two-locality run covers reciprocal rotation work, all three gravity modes
with both signs of rotation, refinement/coarsening energy preservation, and
independent naive work. The long convergence case was excluded only from this
distributed run, having passed in both full suites:

```sh
OCTOTIGERII_TEST_LOCALITIES=2 \
GTEST_FILTER=-RotatingGravityIntegration.TemporalConvergence \
  python3 tests/distributed.py release/3d/tests/rotatingGravityChecks-3d
```

The launcher required execution outside the sandbox because local TCP sockets
are blocked inside it. Both localities exited successfully. Its log is
`/tmp/octo-reorg-gravity-distributed.log`.

The separate HPX build at `/tmp/octotigerII-reorg-hpx-no-profile` used
`OCTOTIGERII_WITH_HPX=ON`, `OCTOTIGERII_WITH_PROFILING=OFF`, and tests disabled.
Its compile commands retain `OCTOTIGERII_WITH_HPX=1` in all six runtime
translation units, without `OCTOTIGERII_PROFILE_HPX`. A rotating-star run
through 0.05 s on 64 level-2 blocks, with two HPX worker threads, completed
with zero corrected mass and combined-energy drift; the scaled energy-budget
residual was −4.07×10⁻¹⁷.

Local build and test logs are `/tmp/octo-reorg-serial-{build,tests}.log`,
`/tmp/octo-reorg-hpx-{build,tests}.log`, and
`/tmp/octo-reorg-hpx-no-profile-{configure,build,star}.log`.

## Rotating-star comparison

The serial corotating star was rerun through 0.25 s and compared with the saved
pre-existing validation CSV in `rotating-frame/`. All 255 diagnostic entries
(five rows, 51 columns) were compared. Steps, times, corrected mass, peak
density, kinetic and grid potential energy, grid angular momentum, and scaled
combined-energy drift match exactly. Grid gas/thermal energies differ by at
most 4.63×10⁻¹⁶ relative; corrected momentum differences are at most
3.06×10⁻¹⁷ of their diagnostic norms. Nearly cancelled momentum residuals are
compared against those norms rather than divided by their near-zero totals.

Both runs take four steps, starting with 512 level-3 blocks and then using
56 level-2 plus 64 level-3 blocks. Final corrected mass drift is zero; scaled
combined-energy drift is −9.27×10⁻¹⁷. The energy-budget residual changes by at
most 8.18×10⁻¹⁸ between the saved run and this run. These are roundoff-scale
differences, within the existing validation tolerances.

```sh
/tmp/octotigerII-conservation-serial/octoII-3d \
  --config=examples/rotating-star.ini --runtime.stopTime=0.25 \
  --verification.analytic=off --output.enabled=off \
  --output.directory=/tmp/octo-reorg-star-corotating
```

The stationary-grid HPX star comparison checks all 204 diagnostic entries
(four rows, 51 columns) against the corresponding saved validation CSV.
Its three steps, all times, hierarchy, corrected mass, peak density, kinetic
energy, and grid angular momentum match exactly. Grid and corrected energies
differ by at most 4.34×10⁻¹⁶ relative. Corrected momentum differences are at
most 2.11×10⁻¹⁷ of their diagnostic norms. Scaled combined-energy drift differs
by at most 9.28×10⁻¹⁷, and the energy-budget residual by at most 8.75×10⁻¹⁷.

Final scaled mass drift is 1.47×10⁻¹⁶, combined-energy drift is −4.64×10⁻¹⁷,
and the energy-budget residual is 2.20×10⁻¹⁷. The same initial 512 level-3
blocks become 56 level-2 plus 64 level-3 blocks. This comparison also stays
within the existing tolerances.

```sh
./release/octoII-3d \
  --config=examples/rotating-star.ini --runtime.stopTime=0.25 \
  --frame.omega=0 --verification.analytic=off --output.enabled=off \
  --output.directory=/tmp/octo-reorg-star-inertial \
  --hpx:localities=1 --hpx:threads=2 --hpx:bind=none
```

The candidate CSVs, complete per-column comparison reports, and application
logs are retained locally under `/tmp/octo-reorg-star-{corotating,inertial}*`.
Reference CSVs are in [the rotating-frame validation directory](rotating-frame/).
These short comparisons check preservation of the existing behavior; they do
not extend the original claim to long-term stellar equilibrium accuracy.

Both inventories were checked against the maintained filesystem set: 187
unique file rows, no omissions, no broken links. `git diff --check` passed.
