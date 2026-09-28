"""Independent scalar checks of proposed SSP-RK2 source couplings.

No production solver is called or changed. The IMEX tableau is the standard
SSP2(2,2,2), gamma=1-1/sqrt(2), explicit Heun weights (1/2,1/2), implicit
rows (gamma,0) and (1-2*gamma,gamma), implicit weights (1/2,1/2).
"""

import math

gamma = 1 - 1 / math.sqrt(2)


def imex_forced_relaxation(value, drive, rate, step):
    denominator = 1 + gamma * step * rate
    first = value / denominator
    second = (value + step * drive - (1 - 2 * gamma) * step * rate * first) / denominator
    return value + step * drive - step * rate * (first + second) / 2


def forced_sdirk(value, drive, rate, step):
    # The constant drive participates in BOTH implicit source stages.
    first = (value + gamma * step * drive) / (1 + gamma * step * rate)
    base = value + (1 - gamma) / gamma * (first - value)
    return (base + gamma * step * drive) / (1 + gamma * step * rate)


print("Exact steady forced relaxation: q'=d-rate*q, d=rate, q(0)=1.")
print("rate*dt          IMEX SSP2(2,2,2)       forced SDIRK")
for stiffness in (1e-3, 1, 10, 100, 1e4, 1e8):
    imex = imex_forced_relaxation(1, stiffness, stiffness, 1)
    original = forced_sdirk(1, stiffness, stiffness, 1)
    print(f"{stiffness:<16g}{imex:<24.12g}{original:.12g}")

print("\nSource solve after each Heun Euler substep, then convex average:")
print("q'=-q, exact source substeps, q(0)=1, final time 1")
errors = []
for steps in (20, 40, 80, 160):
    step = 1 / steps
    result = ((1 + math.exp(-2 * step)) / 2) ** steps
    errors.append(abs(result - math.exp(-1)))
print("errors:", " ".join(f"{error:.9g}" for error in errors))
print("orders:", " ".join(f"{math.log2(a / b):.6f}" for a, b in zip(errors, errors[1:])))
