# AMR validation

See `VALIDATION.md` at the repository root for results and scope. The logs here
come from Release builds with GCC 13.3.0 and HPX 1.11.0 where enabled.

Build any supported problem/dimension with tests enabled and run:

```bash
ctest --test-dir BUILD_DIRECTORY --output-on-failure -j 2
BUILD_DIRECTORY/tests/amrChecks
```

For HPX, the serial test command needs `--hpx:threads=2 --hpx:bind=none`.
Distributed AMR checks use:

```bash
OCTOTIGERII_TEST_LOCALITIES=2 python3 tests/distributed.py BUILD_DIRECTORY/tests/amrChecks
OCTOTIGERII_TEST_LOCALITIES=3 python3 tests/distributed.py BUILD_DIRECTORY/tests/amrChecks
```

The serial matrix checked Sod 1D, Rayleigh–Taylor 2D, streaming 3D, and gravity
sphere 3D. The coupled collapse 3D build ran the AMR suite. Distributed runs
checked Rayleigh–Taylor, streaming, and gravity sphere, all in 3D.

The RT smoke run used the new example input:

```bash
BUILD_DIRECTORY/octotigerII-rayleigh-taylor-2d \
  --config=examples/rayleigh-taylor-amr.ini \
  --runtime.stopTime=0.2 --output.every=4
```

The coupled gravity smoke used:

```bash
BUILD_DIRECTORY/octotigerII-collapse-3d \
  --config=bin/science/Collapse/inputs \
  --amr.enabled=on --amr.minLevel=1 --amr.maxLevel=3 \
  --amr.maxCellMass=1e28 --amr.shadowTolerance=0.1 --amr.regridEvery=2 \
  --runtime.stopTime=0.5 --output.enabled=off \
  --verification.analytic=on --verification.gravityReference=direct \
  --verification.directSamples=32 --gravity.multipoleOrder=4 --gravity.openingAngle=0.3
```

The longer mixed-level collapse run uses the same command with:

```bash
  --amr.maxCellMass=2e27 --amr.shadowTolerance=0 --amr.regridEvery=1 \
  --amr.bufferCells=0 --timestep.cfl=0.002 --runtime.stopTime=2 \
  --output.enabled=on --output.directory=output/collapse-amr --output.every=1
```

Replace the corresponding arguments in the first command; duplicate option
names in one command are rejected. This run isolates the mass criterion and
exercises coarsening and gravity-tree rebuilding during time integration.

The Silo tests construct a genuinely mixed-level mesh and read back each
leaf's coordinates and refinement level. The RT application log exercises
repeated shadow updates and regrids, and writes a VisIt time series.
