"""Independent checks of the offline moving-star equations.

Run with ``python3 tests/radiatingStarRotationChecks.py``. The manufactured
field is differentiated as a Cartesian tensor, so this checks the spherical
geometry terms without copying their implementation into the expectation.
"""
from pathlib import Path
from types import SimpleNamespace
import sys
import unittest

import numpy as np
from scipy.sparse import csr_matrix, diags, eye, kron

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "science"))
from radiating_star_rotation import moving_residual


def cartesian_state(position):
    x, y, z = position
    radius2 = x*x + y*y
    radial = .12 + .03*z*z
    vertical = .17 + .02*radius2
    azimuthal = .11 + .01*z*z
    energy = 1.3 + .13*radius2 + .09*z*z + .03*radius2*z*z
    flux = np.array([x*radial-y*azimuthal, y*radial+x*azimuthal, z*vertical])
    beta = np.array([-y, x, 0]) * (.08 + .02*z*z)
    rho = 1 + .1*radius2 + .05*z*z
    temperature = .9 + .04*radius2 + .02*z*z
    return energy, flux, beta, rho, temperature


def pressure_tensor(position):
    energy, flux, _, _, _ = cartesian_state(position)
    flux2 = np.sum(flux*flux)
    factor2 = flux2/(energy*energy)
    eddington = (3+4*factor2)/(5+2*np.sqrt(4-3*factor2))
    return energy*((1-eddington)/2*np.eye(3)
                   + (3*eddington-1)/2*np.outer(flux, flux)/flux2)


def derivative(nodes):
    result = np.zeros((len(nodes), len(nodes)))
    for i in range(len(nodes)):
        start = min(max(i-2, 0), len(nodes)-5)
        sites = np.arange(start, start+5)
        offsets = nodes[sites]-nodes[i]
        moments = np.array([offsets**power for power in range(5)])
        result[i, sites] = np.linalg.solve(moments, [0, 1, 0, 0, 0])
    return csr_matrix(result)


def manufactured_model(count):
    radii = np.linspace(.7, 1.3, count)
    cosines = np.linspace(.2, .7, count)
    r, mu = np.meshgrid(radii, cosines, indexing="ij")
    symmetry = 1-mu*mu
    cylindrical = r*np.sqrt(symmetry)
    z = r*mu
    radial = .12 + .03*z*z
    vertical = .17 + .02*cylindrical*cylindrical
    azimuthal = .11 + .01*z*z
    energy = 1.3+.13*cylindrical**2+.09*z*z+.03*cylindrical**2*z*z
    flux = r*(symmetry*radial+mu*mu*vertical)
    polar = r*mu*(vertical-radial)
    W = r*azimuthal
    beta = cylindrical*(.08+.02*z*z)
    rho = 1+.1*cylindrical**2+.05*z*z
    temperature = .9+.04*cylindrical**2+.02*z*z
    zero = np.zeros_like(r)
    dr = kron(derivative(radii), eye(count), format="csr")
    dm = kron(eye(count), derivative(cosines), format="csr")
    inverse = diags(1/r.ravel())
    size = count*count
    empty = np.array([], dtype=int)
    model = SimpleNamespace(
        nr=count, nm=count, N=size, r=r, MU=mu, S=symmetry,
        x=np.array([]), omega2=0, hscale=1,
        kappa=np.full_like(r, .4), forceScale=np.ones_like(r),
        radScale=np.ones_like(r), mask=np.ones_like(r, dtype=bool),
        rb=dr, me=dm, mo=dm, ir=inverse,
        gas=diags(mu.ravel())@dr+inverse@diags(symmetry.ravel())@dm,
        lap=csr_matrix((size, size)), eqint=csr_matrix((0, size)),
        center=empty, equator=empty, outer=empty, phib=np.array([]),
        baseline=np.zeros(5*size),
        decode=lambda unused: (zero, zero, rho, temperature, energy, flux, polar))
    return model, W, beta


class MovingReferenceEquations(unittest.TestCase):
    def test_full_m1_meridional_divergence_matches_cartesian_tensor(self):
        errors = []
        for count in (9, 17, 33):
            model, W, beta = manufactured_model(count)
            actual = moving_residual(model, None, W, beta).reshape(5, count, count)
            error = 0.
            for i in range(2, count-2):
                for j in range(2, count-2):
                    mu = model.MU[i, j]
                    sine = np.sqrt(1-mu*mu)
                    point = np.array([model.r[i, j]*sine, 0., model.r[i, j]*mu])
                    energy, flux, velocity, rho, temperature = cartesian_state(point)
                    pressure = pressure_tensor(point)
                    divergence = np.zeros(3)
                    for dimension in range(3):
                        displaced = point.astype(complex)
                        displaced[dimension] += 1e-28j
                        divergence += pressure_tensor(displaced).imag[:, dimension]/1e-28
                    emission = .75*temperature**4
                    force = .4*rho*(flux-pressure@velocity-emission*velocity)
                    source_energy = .4*rho*(energy-emission-velocity@flux)
                    expected_radial = np.array([sine, 0, mu])@(divergence+force)
                    expected_polar = -np.array([mu, 0, -sine])@(divergence+force)/sine
                    error = max(error, abs(actual[2, i, j]-expected_radial),
                                abs(actual[3, i, j]-expected_polar))
                    self.assertAlmostEqual(actual[0, i, j], -force[2]/rho, places=13)
                    expected_thermal = (source_energy-velocity@force)/(.4*rho*energy)
                    self.assertAlmostEqual(actual[1, i, j], expected_thermal, places=13)
            errors.append(error)
        print("Cartesian-tensor residual errors:", errors)
        self.assertGreater(errors[0]/errors[1], 8)
        self.assertGreater(errors[1]/errors[2], 8)
        self.assertLess(errors[2], 5e-9)

    def test_retained_thermal_balance_is_the_comoving_energy(self):
        for point in ([.8, 0, .4], [.6, .4, .9], [1.1, -.3, .2]):
            energy, flux, beta, _, _ = cartesian_state(np.array(point))
            pressure = pressure_tensor(np.array(point))
            emission = (energy-2*beta@flux+beta@pressure@beta)/(1-beta@beta)
            source_energy = energy-emission-beta@flux
            force = flux-pressure@beta-emission*beta
            self.assertGreater(emission, 0)
            self.assertAlmostEqual(source_energy-beta@force, 0, places=14)


if __name__ == "__main__":
    unittest.main()
