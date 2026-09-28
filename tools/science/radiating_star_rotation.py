"""Moving-gas terms for the offline axisymmetric radiation-star reference.

Radiation flux is normalized by c times the pressure unit. The toroidal
unknown W has F_phi/c = sqrt(1-mu**2) W. No compensating gas torque is added.
The linear toroidal solve is followed by the same retained mixed-frame
thermal and meridional force terms used by OctoII's matter coupler.
"""
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve


def light_speed(rho_c=1.0, mean_molecular_weight=0.6, beta_c=0.8):
    """Physical c / sqrt(Rgas Tc), using OctoII's CGS constants."""
    gas_constant = 1.380649e-16 / (mean_molecular_weight * 1.66053906892e-24)
    temperature = (3 * rho_c * gas_constant * (1-beta_c)
                   / (7.565733250280004e-15 * beta_c)) ** (1/3)
    return 2.99792458e10 / np.sqrt(gas_constant * temperature)


def angular_velocity_squared(model, state, W=None, c=None):
    """Recover centrifugal support with exact spherical balance subtracted.

    Negative support in material is an invalid reference, never clipped.
    The virtual barotropic potential in vacuum need not describe rotation.
    """
    h, potential, rho, _, energy, flux, polar = model.decode(state)
    h0, potential0, _, _, _, flux0, _ = model.decode(model.initial.ravel())
    delta_s = (h + potential - h0 - potential0).ravel()
    inverse_r = model.ir.diagonal()
    mu = model.MU.ravel()
    opacity = model.kappa.ravel()
    force = opacity * ((flux-flux0).ravel() - mu*polar.ravel())
    numerator = model.re @ delta_s - mu*inverse_r*(model.me @ delta_s) - force
    omega2 = (inverse_r*numerator).reshape(model.nr, model.nm)
    omega2[0] = omega2[1]
    if W is not None:
        if c is None:
            raise ValueError("A physical light speed is required for moving-gas corrections")
        # Picard iteration of the small M1 off-diagonal force correction.
        for _ in range(4):
            active = rho > 1e-12
            if np.any(omega2[active] < 0):
                raise ValueError("Reference requires negative angular velocity squared")
            beta = model.r * np.sqrt(np.where(active, omega2, 0)) / c
            phi_flux = np.sqrt(model.S)*W
            f2 = (flux*flux + model.S*polar*polar + phi_flux*phi_flux)/(energy*energy)
            closure = 3/(2+np.sqrt(4-3*f2))
            factor = model.S*beta*W*closure/energy
            correction = opacity*(factor*(flux-model.MU*polar)).ravel()
            omega2 = (inverse_r*(numerator+correction)).reshape(model.nr, model.nm)
            omega2[0] = omega2[1]
    if np.any(omega2[rho > 1e-12] < 0):
        raise ValueError("Reference requires negative angular velocity squared")
    return np.where(rho > 1e-12, omega2, 0)


def toroidal_flux(model, state, omega2, c, previous=None):
    """Outgoing azimuthal M1 momentum balance, with no momentum heater.

    ar=b f_r, aq=b q, where b=3/(2+sqrt(4-3 f**2)). The equation is
      (ar W)_r + [3 ar W + (1-mu²)(aq W)_mu - 4 mu aq W]/r
        + chi W = chi (P_phi_phi+B) r Omega/c.
    At the center W=0; polar regularity permits finite W on the axis.
    """
    _, _, rho, temperature, energy, flux, polar = [a.ravel() for a in model.decode(state)]
    mu = model.MU.ravel()
    symmetry = model.S.ravel()
    W0 = np.zeros(model.N) if previous is None else previous.ravel()
    fr, q = flux/energy, polar/energy
    phi_flux = np.sqrt(symmetry)*W0
    f2 = fr*fr + symmetry*q*q + (phi_flux/energy)**2
    if np.any(f2 >= 1):
        raise ValueError("Toroidal correction exceeds the radiation flux cone")
    closure = 3/(2+np.sqrt(4-3*f2))
    ar, aq = closure*fr, closure*q
    p_phi = energy*(1/3-closure*f2/3) + closure*phi_flux*phi_flux/energy
    blackbody = 0.75*temperature**4
    chi = model.kappa.ravel()*rho
    matrix = (model.rb @ diags(ar)
              + model.ir @ (diags(3*ar) + model.sm @ model.mo @ diags(aq)
                            - diags(4*mu*aq)) + diags(chi)).tolil()
    rhs = chi*(p_phi+blackbody)*model.r.ravel()*np.sqrt(omega2.ravel())/c
    matrix[model.center] = 0
    matrix[model.center, model.center] = 1
    rhs[model.center] = 0
    W = spsolve(matrix.tocsr(), rhs).reshape(model.nr, model.nm)
    if not np.all(np.isfinite(W)):
        raise ValueError("Nonfinite azimuthal radiation solution")
    return W


def moving_residual(model, state, W, beta):
    """Five meridional equations with frozen toroidal flux and velocity.

    Used for a Picard correction, so the existing sparse meridional Jacobian
    pattern remains valid. beta is the actual v_phi/c, including sin(theta).
    The static spherical residual is subtracted only in reference generation.
    """
    h, potential, rho, temperature, energy, flux, polar = [a.ravel() for a in model.decode(state)]
    beta = beta.ravel()
    phi_flux = np.sqrt(model.S.ravel())*W.ravel()
    opacity = model.kappa.ravel()
    mu, symmetry = model.MU.ravel(), model.S.ravel()
    ir = model.ir.diagonal()
    fs, rs = model.forceScale.ravel(), model.radScale.ravel()
    fr, q, fphi = flux/energy, polar/energy, phi_flux/energy
    f2 = fr*fr + symmetry*q*q + fphi*fphi
    if np.any(f2 >= 1):
        raise ValueError("Moving reference exceeds the radiation flux cone")
    b = 3/(2+np.sqrt(4-3*f2))
    A = energy*(1/3-b*f2/3)
    P = A+energy*b*fr*fr
    C, D = energy*b*fr*q, energy*b*q*q
    Pphi = A+energy*b*fphi*fphi
    B = 0.75*temperature**4
    force_factor = 1-beta*b*fphi
    kap_force = opacity*force_factor
    S = h+potential
    gas = (model.gas @ S - kap_force*(mu*flux+symmetry*polar))/fs
    eq = model.equator
    gas[eq] = (S[eq]-model.hscale-potential[0] - 0.5*model.omega2*model.x**2
               - model.eqint @ (kap_force*flux))/model.hscale
    gas[model.center] = h[model.center]/model.hscale-1
    radial = (model.rb @ P + ir*(2*P-2*A-symmetry*D
               - energy*b*fphi*fphi + symmetry*(model.mo @ C)-2*mu*C)
               + kap_force*rho*flux)/rs
    polar_equation = (model.rb @ C + ir*(3*C+model.me @ A
                      + symmetry*(model.me @ D)-3*mu*D
                      + mu*energy*b*(W.ravel()/energy)**2)
                      + kap_force*rho*polar)/rs
    radial[model.center] = fr[model.center]
    polar_equation[eq] = q[eq]
    polar_equation[model.center] = q[model.center]
    div_flux = model.rb @ flux + ir*(2*flux+symmetry*(model.mo @ polar)-2*mu*polar)
    # Exact internal-energy balance for the retained equal-opacity equations.
    thermal = np.where(model.mask.ravel(),
                       (energy-B-2*beta*phi_flux+beta*beta*(Pphi+B))/np.maximum(energy,1e-8),
                       div_flux/rs)
    gravity = model.lap @ potential-rho
    gravity[model.center[1:]] = (potential[model.center[1:]]-potential[0])/model.hscale
    gravity[model.outer] = (potential[model.outer]-model.phib)/model.hscale
    return np.array([gas,thermal,radial,polar_equation,gravity]).ravel()-model.baseline
