# Helmholtz hydro coupling validation

2026-09-28, GCC 13.3, Linux x86-64, Release, no fast-math. The default remains
ideal gas. This validates the opt-in Helmholtz adapter, using Timmes's unchanged
precomputed table (digest in `lib/helmholtz/NOTICE.md`).

Measured checks:

- Nine Helmholtz tests passed in serial 1D. In both 2D and 3D, eight passed and
  the recorded 1D source-state regression was skipped. They include energy,
  pressure and entropy round trips; gas/total selection; the independent
  constant-entropy pressure derivative for sound speed; floor/entropy recovery;
  material composition and tracer exclusion; stiff LTE/relaxation; and coupled
  species/energy conservation through AMR subcycling and coarsening.
- All nine checks passed in the 1D HPX build with profiling disabled,
  with two worker threads on one locality and with two TCP localities. The
  distributed run exercises serialization of weighted species field handles.
- The suite includes an accepted-patch floor ledger check and a regression
  for a stiff, composition-dependent radiation
  drive. Tightening the temperature inversion bracket from 2e-11 to 4e-15
  relative width resolved that source convergence failure.
- Existing 1D ideal-gas hydro and dual-energy suites passed. Radiation coupling:
  14 passed, two dimension-dependent tests skipped. Species: 7 passed;
  storage: 2; units: 4; options: 43. These guard the existing default closure.
- The hydro-disabled 1D executable built and completed a radiation streaming
  smoke run, checking the optional-module build guards.
- `examples/helmholtz-advection.ini` completed its 0.01 s interval in the serial
  1D executable and wrote conservation diagnostics. With radiation and opacity
  1 cm^2/g it completed 35 steps to the same time; the reported combined-energy
  drift and floor corrections were both zero at the output precision.

Reference comparison of the standalone C++ evaluator was previously exact for
all 126 outputs at 421 states (53,046 comparisons), including the documented
cold-degenerate derivative failure. That is translation parity, not independent
physical validation of every table state.

Reproduction (dependencies already installed):

```sh
cmake -S . -B /tmp/octoii-helm-coupling -DOCTOTIGERII_WITH_HPX=OFF \
  -DOCTOTIGERII_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/octoii-helm-coupling --target \
  helmholtzHydroChecks-1d helmholtzHydroChecks-2d helmholtzHydroChecks-3d -j 4
/tmp/octoii-helm-coupling/1d/tests/helmholtzHydroChecks-1d
/tmp/octoii-helm-coupling/2d/tests/helmholtzHydroChecks-2d
/tmp/octoii-helm-coupling/3d/tests/helmholtzHydroChecks-3d
```

For HPX, configure `OCTOTIGERII_WITH_HPX=ON`, set `HPX_DIR`, and use
`OCTOTIGERII_WITH_PROFILING=OFF` to exercise the independent synchronization
guard. Run the executable with `--hpx:threads=2 --hpx:bind=none`. For two
localities, run `python3 tests/distributed.py /path/to/helmholtzHydroChecks-1d`.
Local TCP sockets require execution outside the restrictive sandbox.

The EOS inversion path is not performance tuned. These checks do not establish
shock convergence, production merger accuracy, or equilibrium of a stellar
profile built with a different EOS. The library's cold-degenerate table
limitation remains. Read `docs/helmholtz-hydro.md` for the physical closure,
composition treatment, temperature floor and its conservation accounting.
