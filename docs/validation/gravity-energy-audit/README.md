# Gravity energy audit — 2026-09-24

Based on the uploaded `OctotigerII-source.tar(10).gz` and `conservation.csv`.
All pre-existing working-tree changes and Git history were preserved.

## Finding

The supplied CSV ends at step 135, t=7.396116 s, with scaled mass drift
-3.605464e-13 and combined gas/gravity energy drift -1.913069e-09. The largest
energy drift magnitude is 5.714435e-09 at step 10. The initial mass, timestep,
and first three energy-drift values were reproduced using the shipped polytrope
example with the original adaptive interaction walk. The later 135-step history
was not rerun; the comparisons here stop at t=0.25 s (four timesteps).

The reproduced energy defect agrees with

    R = (1/2) sum[V (rho_old*phi_new - rho_new*phi_old)]

to approximately 1e-17 of the energy scale. The flux work and physical-boundary
ledger therefore close independently of the potential-operator error. In this
short reproduction, regrid energy changes are much smaller than the original
1e-9 drift, although the old remap is a real, separate conservation gap.

The adaptive walk previously split the target first for equally sized nodes.
Reversing source/target could therefore accept different interaction groups,
violating reciprocal potential coefficients. Descending both nodes together
restores the pair symmetry; no multipole-order increase was made. Independent
bilinear tests at p=2 and p=5 use unrelated mass distributions on a mixed-level
mesh to check this without any hydrodynamic energy correction.

## Changes

- `gravity.energyTreatment=mullen|naive`, default `mullen`. The naive path matches
  the old unbracketed kick evolution exactly in a three-step comparison, including
  every hydro component. Both modes retain the endpoint momentum kicks and record
  the same physical-boundary potential-flux convention.
- `gravity.conserveRegridEnergy=on|off`, default `on`, independent of the first
  option. Restriction/prolongation remap E+rho*phi/2 as a signed scalar; after the
  new field solve, recover E with the new density and potential. Gas and entropy
  remain separate from the remap register. There is no global redistribution.
- Symmetric adaptive gravity-tree traversal. This changes the approximate fields
  slightly but leaves the momentum-update algorithm and multipole order intact.
- Compensated grid-integral summation. The approximately 3.6e-13 mass offset in
  the original short reproduction disappears at printed precision with this
  summation; zero in a CSV is not a claim of exact real-arithmetic conservation.
- CSV columns for cumulative potential reciprocity defect, actual regrid energy
  change, and the remaining scaled budget residual. Original drift columns retain
  their physical meaning and include all defects.

## Comparisons

GCC 13.3.0, serial Release build, p=5, openingAngle=0.5, shipped
`examples/polytrope.ini`, analytic verification disabled, t=0.25 s.
Both the initial and evolved meshes are adaptive. Scales are those defined in
`docs/conservation.md`, not division by the signed net binding-plus-gas energy.

| Case | Final scaled mass drift | Final scaled energy drift |
| --- | ---: | ---: |
| Original asymmetric walk; Mullen; conservative remap | 0.000000e+00 | -2.968104e-09 |
| Symmetric walk; Mullen; conservative remap | 0.000000e+00 | -6.252829e-17 |
| Symmetric walk; naive; conservative remap | 0.000000e+00 | 1.092302e-03 |
| Symmetric walk; Mullen; gas-only remap | -2.061443e-16 | 1.394381e-13 |
| Symmetric walk; naive; gas-only remap | 0.000000e+00 | 1.092302e-03 |

The forced refinement/coarsening regression separately exercises a substantial
change in binding energy, all four control combinations, species enabled/disabled,
and local recovery of the conservatively transferred combined energy in each
cell. It checks total-energy preservation at 2e-14 of the energy scale with the
remap enabled, and verifies that gas-only transfer does not preserve total energy.
A repeated gravity solve must not apply the remap correction again.

The symmetric walk incurs more interactions in this four-step example:
3,628,848 multipole pairs and 27,261,312 direct pairs, versus 3,577,760 and
24,142,464 previously (all solves, including startup/regrid). These are work
counts, not a controlled wall-time performance benchmark.

## Validation

**510/510 serial tests passed across 1D, 2D, and 3D**, including isolated and
periodic adaptive energy budgets, direct-field comparisons on mixed meshes and
image boundaries, forced refinement/coarsening, dual energy, Silo, options,
and the naive compatibility check. The complete CTest log and five comparison
CSVs accompany this note. Build completed without compiler warnings.

HPX is not installed here. Distributed execution and the updated HPX option
serialization regression remain untested in this environment. No new standard
threading/synchronization facilities were introduced in HPX paths.

Reproduce:

```bash
cmake -S . -B build -DOCTOTIGERII_WITH_HPX=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
ctest --test-dir build --output-on-failure -j4
for treatment in mullen naive; do
  for remap in on off; do
    ./build/octoII-3d --problem.name=polytrope --config=examples/polytrope.ini \
      --runtime.stopTime=0.25 --verification.analytic=off --output.enabled=off \
      --gravity.energyTreatment="$treatment" --gravity.conserveRegridEnergy="$remap" \
      --output.directory="output/energy-${treatment}-${remap}"
  done
done
```
