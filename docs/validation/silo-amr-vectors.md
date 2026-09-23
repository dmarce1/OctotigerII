# Adaptive Silo metadata and vector output (2026-09-23)

The reported RT disappearance was reproduced with the 3D AMR example:

```sh
release/rayleigh-taylor/3d/octotigerII-rayleigh-taylor-3d \
  --config=examples/rayleigh-taylor-amr.ini \
  --runtime.stopTime=0.014 --output.every=4
```

| Frame | Step | Time (s) | Leaf blocks | Written volume (cm³) | Written mass (g) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 0 | 0.000000e+00 | 64 | 1.000000e+00 | 1.500000e+00 |
| 1 | 4 | 8.825578e-03 | 64 | 1.000000e+00 | 1.500000e+00 |
| 2 | 8 | 1.323793e-02 | 512 | 1.000000e+00 | 1.500000e+00 |
| 3 | 9 | 1.400000e-02 | 512 | 1.000000e+00 | 1.500000e+00 |

Every frame covers [-0.5, 0.5] cm on all three axes. At frame 2, the first
64 Morton-ordered blocks cover only [-0.5, 0]³: one eighth of the volume and
one quarter of a central slice. The original files omit VisIt's
`MetadataIsTimeVarying` flag, allowing its domain selection to remain at
64 blocks after refinement. The simulation and stored mesh retain the full
domain. VisIt's Silo reader checks this flag in
`avtSiloFileFormat::CheckForTimeVaryingMetadata` and separately checks
`ConnectivityIsTimeVarying` in `GetTimeVaryingInformation`.

Both flags are now integer 1 at the root of every AMR frame, including the
first. The new regression test failed on the original writer because the
metadata flag was absent, then passed with the corrected writer. It writes
coarse, mixed, fine, and coarsened frames, checks both flags and block counts,
and verifies exactly-once coverage on a finest-level block lattice using
coordinates read from disk.

Native vectors use `DBPutQuadvar`; scalar fields use `DBPutQuadvar1`.
Multivars explicitly name their mesh and scalar/vector rank. The tests read
every component, CGS units, centering, active dimension count, exact values,
and signed errors. Test inputs give components distinct values to expose
permutations and sign mistakes. No separate X/Y/Z scalar quadvars or
multivars remain for vector quantities. JSON component error norms retain
their existing format.

All three Silo tests passed in each of these Release configurations:

| Problem | Dimensions | Backend | Passed |
| --- | ---: | --- | ---: |
| Sod | 1 | Serial | 3 |
| Sod | 2 | Serial | 3 |
| Rayleigh–Taylor | 3 | Serial | 3 |
| Streaming radiation | 3 | Serial | 3 |
| Gravity sphere | 3 | Serial | 3 |
| Rayleigh–Taylor | 3 | HPX 1.11, two threads | 3 |

The RT run above was repeated after the fix. Its console diagnostics matched
the original run, and readback confirmed complete volume, mass, bounds, and
the new metadata flags across all four frames. VisIt itself was unavailable
in the validation environment; these are Silo-level checks and inspection of
the VisIt 3.4.2 reader source.

Build validation rejects dimensions 1 and 2 for RT, gravity sphere, gravity
Gaussian, and collapse, while accepting 3D. A Sod 1D invocation with nonzero
`hydro.acceleration.x` is also rejected.

References:

- [VisIt 3.4.2 Silo reader](https://github.com/visit-dav/visit/blob/v3.4.2/src/databases/Silo/avtSiloFileFormat.C)
- [VisIt report of the same AMR time-series metadata issue](https://github.com/visit-dav/visit/issues/4845)
- [Silo quad variables](https://silo.readthedocs.io/latest/objects.html#dbputquadvar)
- [Silo multivars](https://silo.readthedocs.io/latest/parallel.html#dbputmultivar)
