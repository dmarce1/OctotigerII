# Prescribed photon-source verification

`photon-source` starts homogeneous gas at rest with density
`1e-10 g/cm^3`, temperature `1e4 K`, and initially thermal radiation in a
periodic box. The prescribed isotropic photon emissivity is fixed at
`H = 1 erg/(cm^3 s)`. It is independent of the evolving temperature and
radiation field. It supplies no direct momentum or gas heating.

With `eta = chat/c`, radiation receives `dEr/dt = eta H` in addition to
gas/radiation exchange. Set `radiation.opacity=0` for the exact transparent
solution `Er(t)=Er(0)+eta H t`, with unchanged gas and zero radiation flux.
For positive opacity, absorption and emission redistribute this supplied
energy. The combined source budget is

    Egas(t) + Er(t)/eta = Egas(0) + Er(0)/eta + H t

per unit volume. The defaults use opacity `1 cm^2/g`; the same problem can
check full and reduced light speed by changing `radiation.lightSpeedRatio`.

`conservation.csv` records the actual injected radiation energy and its
RSLA-weighted counterpart separately from boundary transport. Corrected
combined energy subtracts the supplied photon energy. The automated checks
also exercise AMR subcycling, regridding, and restoration after a rejected
interval so predictor and forecast stages cannot be counted as injections.

For example, from a build's dimension directory:

```sh
./octoII-3d --problem.name=photon-source --radiation.opacity=0 --radiation.lightSpeedRatio=0.2 --output.directory=output-photon-source
```

The `photonHeatingChecks-3d` test executable also checks the production
source-aware patch update against a time-linear prescribed source, for which
midpoint injection is exact, and verifies that source energy remains separate
from boundary transport in the CSV budget. All source units and reduced-speed
weights are the same in the serial and HPX implementations.
