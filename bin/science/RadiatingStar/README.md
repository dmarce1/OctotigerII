# Opaque rotating gas–radiation star

`radiatingStar` initializes a self-consistent, uniformly rotating stellar
interior with gas pressure, radiation pressure, and isolated self-gravity.
It evolves with the physical speed of light, constant gray absorption opacity,
and **zero imposed photon source**. Radiation energy cannot leave the numerical
box; gas retains outflow boundaries. The grid is inertial.

The structural total pressure is `P=K rho^(1+1/n)`, with `n=3.5`. Temperature
comes from `P=rho Rgas T+aT^4/3`, and the central gas-pressure fraction is 0.8.
The gas evolution still uses gamma=5/3. Entropy increases outward, and the
uniform angular velocity introduces no vertical shear. The SCF solve includes
the centrifugal potential and the nonspherical gravitational potential.

Default parameters give approximately 39.7 solar masses, equatorial radius
`1.4491e12 cm`, polar radius `1.4189e12 cm`, central temperature `2.395e7 K`,
and dynamical time `2.403e4 s`. The spin is 0.2 times the spherical model's
breakup frequency, giving an 8.50-day rotation period. `star.radius` is derived
from the EOS; it must not be supplied independently. The default box width is
2.4 times the nominal spherical radius.

The opacity `0.34 cm²/g` is a constant gray **absorption** coefficient in this
test. Its value is numerically familiar from scattering models, but this
benchmark does not implement scattering opacity. The central optical depth is
`6.66e10`; a diffusion-time estimate exceeds `3.7e7` dynamical times. Only the
low-density outer few cells have a short thermal adjustment time: material
outside 0.95 equatorial radii contains about `2.8e-7` of the mass.

Initialization uses LTE comoving radiation energy and the interior diffusion
flux. The comoving flux is limited only where the formal surface diffusion
formula becomes unrealizable. Its exact perpendicular M1 boost adds the
azimuthal radiation advection and preserves local thermal balance for the
retained mixed-frame source equations. There is no compensating heater,
acceleration, or torque. The insulating numerical boundary sets accepted
radiation energy transport to zero while retaining radiation momentum flux.

This is a mechanical equilibrium with a long bulk thermal time, not an exact
stationary radiative stellar atmosphere. Surface adjustment and internal
redistribution are expected. The initial SCF reference is useful for measuring
drift; longer evolution and mesh refinement are needed to demonstrate stability.

Use `examples/radiating-star.ini`. It requests only a short startup run; the
full-light-speed transport timestep makes longer runs more expensive. The
`radiatingStarChecks` tests verify the initial mixed-frame moments, force
balance, configuration, and a short coupled evolution with no photon heater or
escaping radiation energy. Existing SCF tests independently check the spherical
Lane–Emden limit, rotational deformation, virial convergence, and EOS stability.
