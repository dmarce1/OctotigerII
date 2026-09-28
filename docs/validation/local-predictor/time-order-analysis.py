"""Exact Fourier-mode diagnostic for fixed-mesh time order; no dependencies.

Scalar periodic advection u_t + u_x = 0, centered unlimited PLM slopes,
positive-speed upwind numerical flux. The reference is the exact solution
of that spatial discretization, exp(T*L), not the continuum PDE solution.
This isolates time error. It does not validate the production gravity solver.
Run: python3 docs/validation/local-predictor/time-order-analysis.py
"""

import cmath
import math

cells, mode, stop = 32, 4, 0.25
dx = 1 / cells
theta = 2 * math.pi * mode / cells
jump = 1 - cmath.exp(-1j * theta)
reconstruction = 1 + 0.5j * math.sin(theta)
spatial_operator = -jump * reconstruction / dx
hancock_quadratic = 0.5 * jump * (1j * math.sin(theta)) / dx**2
reference = cmath.exp(stop * spatial_operator)

for method in (
    "Local Hancock",
    "Numerical center derivative with frozen slopes",
    "Numerical midpoint followed by reconstruction",
    "AB2 with midpoint startup",
):
    errors = []
    for steps in (128, 256, 512, 1024):
        dt = stop / steps
        z = dt * spatial_operator
        if method == "Local Hancock":
            result = (1 + z + dt**2 * hancock_quadratic) ** steps
        elif method == "Numerical center derivative with frozen slopes":
            result = (1 + z - 0.5 * dt**2 * jump * spatial_operator / dx) ** steps
        elif method == "Numerical midpoint followed by reconstruction":
            result = (1 + z + 0.5 * z**2) ** steps
        else:
            previous, current = 1 + 0j, 1 + z + 0.5 * z**2
            for _ in range(1, steps):
                previous, current = current, current + z * (1.5 * current - 0.5 * previous)
            result = current
        errors.append(abs(result - reference))
    orders = [math.log2(a / b) for a, b in zip(errors, errors[1:])]
    print(method)
    print("  errors:", " ".join(f"{value:.9g}" for value in errors))
    print("  orders:", " ".join(f"{value:.6f}" for value in orders))
