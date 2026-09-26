# A rigidly rotating grid with inertial conserved quantities

`frame.omega` is a constant angular velocity in radians per second about the
z axis through the coordinate origin. It defaults to zero. For nonzero values,
all physical faces must use `outflow` (free) boundaries. Reflecting, periodic,
inflow, and analytic boundary configurations are rejected. Two- and
three-dimensional transport support rotation; self-gravity remains 3D.

Logical coordinates xi are fixed to the mesh. Physical coordinates are
`x = R_z(omega*t) xi`, with zero initial phase, and mesh velocity is
`w = omega cross x`. Cell volumes and face areas do not change. Hydro momentum,
gas energy, and radiation moments retain their **inertial Cartesian** meaning.
There are no Coriolis or centrifugal source terms in these stored equations.

## Moving-face transport

For a conserved inertial state U, the face flux is

```
F_ALE = F(U) dot n - (w dot n) U.
```

The hydro fluxes are therefore

```
mass:      rho (v_n - w_n)
momentum:  rho v (v_n - w_n) + p n
energy:    E (v_n - w_n) + p v_n.
```

The pressure work uses inertial velocity. Hydro Riemann states are temporarily
expressed in the face basis and boosted by its normal speed; the resulting
flux is returned to inertial components. Limiting uses the same moving-face
fluxes. Species and entropy use the relative mass flux.

Radiation stores physical E_r and F_r. With reduced light speed c_hat, its
moving-face flux is

```
energy:  (c_hat/c) F_r dot n - w_n E_r
flux:    c*c_hat P_r dot n - w_n F_r.
```

Only components are rotated. Radiation moments are not velocity-boosted, and
the physical realizability bound remains `|F_r| <= c E_r`. HLL characteristic
speeds are shifted by the face speed.

Predictors use beginning-of-stage geometry and final fluxes use midpoint
geometry. Each substep accumulates inertial vector fluxes directly in the
coarse/fine registers. CFL uses relative characteristic speeds, including the
variation of mesh speed across a face. A further `|omega| dt <= 0.1` bound
resolves changing normals even when the gas corotates.

The gas outflow diode tests velocity relative to the moving face. A clamped
nonzero-speed ghost state preserves its internal energy. For rotating radiation
boundaries, the exterior Riemann state is vacuum; this avoids inventing an
unrealizable flux when mesh speed exceeds the reduced light speed. Boundary
budgets integrate the actual moving-face fluxes. At zero rotation the original
boundary and transport paths are retained.

## Conservative gravitational rotation work

The gravity solver continues to operate in logical coordinates. Its potential
is a scalar and its acceleration is expressed in grid components. Accelerations
are rotated into inertial components at the actual source stage; a field held
fixed in grid components must not be held fixed in inertial components.

On a fixed mesh the reciprocal scalar operator gives

```
phi = K m,  K^T = K,
W = (1/2) m^T K m,
Delta W = ((phi_a + phi_b)/2)^T Delta m.
```

Rigid rotation leaves K time independent in logical coordinates. The existing
face-work calculation uses the accepted **relative** mass transfers, including
canonical fine-face transfers at refinement boundaries. Locally this accounts
for `rho (v-w) dot g`. The additional rotation work should approximate
`rho w dot g` while producing zero net energy in an isolated system.

Let F_ij be the force on cell i from j. The existing mutual force construction
satisfies F_ji = -F_ij. We use the power exchange

```
B_ij = ((w_i + w_j)/2) dot F_ij,  B_ji = -B_ij.
```

Its sum is zero regardless of approximate gravitational torque. For central
Newtonian forces it equals `w_i dot F_ij`, because

```
(w_i-w_j) dot F_ij = omega dot ((x_i-x_j) cross F_ij) = 0.
```

For approximate noncentral forces it changes the local energy work by half the
pair's spurious torque power. This is a modification at the force approximation's
error order, not a timestep error. It does not change accelerations or establish
angular-momentum conservation. No measured global energy residual is redistributed.

To evaluate the exchange, write each acceleration component as `g_d = A_d m`,
where the mass-coordinate operator satisfies `A_d^T = -A_d`. In addition to the
ordinary field, evaluate signed density fields `(x/L)rho` and `(y/L)rho`, using
one fixed positive reference length L and the same source/target selections:

```
g_xweighted = A [(x/L)m]
g_yweighted = A [(y/L)m]
B_i = omega*m_i/2 * [x_i*g_y - y_i*g_x
                    + L*(g_xweighted_y - g_yweighted_x)].
```

The coordinate-weighted matrix in brackets is skew-symmetric, which proves
`sum B_i = 0`. Signed sources must not be clipped or skipped merely because a
node's monopole vanishes. The first implementation performs two extra field
solves per required field evaluation. Zero rotation skips these solves.

## Rungs and time refinement

A rung denotes a timestep group. Interactions assigned to a rung include
interactions within that group and with faster groups. Each pair belongs to
its slower endpoint's group. Existing internal `shell` identifiers denote these
interaction sets, not spatial shells.

Each symmetric rung operator retains force reciprocity, so the rotation-work
identity also holds separately on each rung. Accepted work uses matched
beginning/end masses and fields on both sides of an interaction. Provisional
substep work is recorded and replaced at interval closure, after reflux. Both
`gravity.timeIntegration=hierarchical` and `conventional` use this conservative
energy accounting, even though conventional momentum updates have a different
cadence. The global-step source/transport/source API also includes matched
rotation work at its endpoints. Existing regrid energy accounting remains active.

## Diagnostics and validation

Silo files keep mesh coordinates in the rotating grid frame, so the mesh stays
aligned with the grid across output times. The cell field arrays are written as
before. Conservation CSV files include `angular_momentum_z_g_cm2_s_grid` and
`maximum_density_g_cm3`, in addition to existing mass and energy budgets. The
angular momentum column is a grid integral; it is not corrected for material
that leaves the domain and is not promised to be conserved to roundoff.

See [the original rotating-star benchmark](rotating-star.md) for stationary-grid
and corotating examples, and the [validation record](validation/rotating-frame/README.md)
for measured errors and reproducible commands. Check total energy independently of stellar structure,
angular momentum, and convergence. Roundoff energy closure by itself does not
establish an accurate stellar evolution.

Korobkin et al. (2021), [Conservation of Angular Momentum in the Fast Multipole
Method](https://arxiv.org/abs/2107.07166), describe centralizing implicit pair
forces. That would be a separate improvement to the force operator. Their
low-order formulas are not substituted for the current arbitrary-order harmonic
FMM, and no angular-momentum-conserving force change is included here.
