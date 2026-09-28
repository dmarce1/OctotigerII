# Original rotating-star equilibrium

`equilibrium.inc` contains the positive cylindrical-radius/vertical-coordinate
quadrant of the `rotating_star.bin` equilibrium distributed with Octo-Tiger's
published rotating-star benchmark. The other quadrants are mirror copies.
Hexadecimal floating-point literals preserve every stored double exactly.
The table is compiled into the executable: execution requires no data-file path,
network access, or particular working directory.

Source repository revision: `3e2f949cf8fb9dc54ce101d05546713ec4ab9034`.

- [Original binary](https://github.com/STEllAR-GROUP/octotiger/blob/3e2f949cf8fb9dc54ce101d05546713ec4ab9034/test_problems/rotating_star/rotating_star_paper/rs_3_w/rotating_star.bin)
- [SCF generator](https://github.com/STEllAR-GROUP/octotiger/blob/3e2f949cf8fb9dc54ce101d05546713ec4ab9034/tools/gen_rotating_star_init/make_bin.cpp)
- [Original initializer](https://github.com/STEllAR-GROUP/octotiger/blob/3e2f949cf8fb9dc54ce101d05546713ec4ab9034/src/test_problems/rotating_star/rotating_star.cpp)
- [Published benchmark, section 3.7](https://doi.org/10.1093/mnras/stab937)

Original binary SHA256:

```
9958252bae60348e08d4d043bf33954f75c43286133a99f49977e58d3b72da08
```

The binary contains two little-endian 32-bit counts (100, 100), one double
angular velocity (0.5155532816213834), and 200 by 200 pairs of density and
internal-energy doubles. Positive samples are at `(i+0.5)/100, (k+0.5)/100`.
The generator uses G=1, n=1.5, and nominal equatorial/polar edges 0.9/0.6;
the selected surface grid centers are 0.905/0.605. Its pressure constant is
approximately 0.1681244461473046. There are 51 SCF iterations.

To reproduce the include file from an already downloaded original binary:

```bash
python3 problem/science/RotatingStar/import_equilibrium.py /path/to/rotating_star.bin
```

The original data generator and initializer are Copyright (c) 2019 AUTHORS,
distributed under the Boost Software License, Version 1.0. This extracted data
and adaptation retain that license; see the repository's `LICENSE` file.
