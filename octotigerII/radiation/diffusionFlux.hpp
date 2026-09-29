/** @file
 * @brief Opacity-aware matching of M1 transport to its equilibrium diffusion flux.
 */
#pragma once
#include "octotigerII/radiation/radiationTransport.hpp"

namespace octotigerII::radiation {

using Extinction = units::Quantity<-1, 0, 0>;
using MaterialVelocity = std::array<units::Velocity, ndim>;

/// Lorentz-isotropic M1 moments at fixed inertial energy. With
/// aT^4=E-v.F/c^2 these satisfy both equal-opacity mixed-frame source balances.
/// The velocity and returned moments use the same orthonormal basis. The
/// manifold is subluminal; the evolution equations remain nonrelativistic.
inline RadiationSystem::State materialEquilibriumMoments(units::EnergyDensity energy, MaterialVelocity const& velocity) {
	Real beta2 = 0;
	for (auto v : velocity) {
		if (!units::finite(v)) throw std::invalid_argument("Nonfinite diffusion material velocity");
		Real const beta = v / constants::c;
		beta2 += beta * beta;
	}
	if (!(beta2 < Real(1))) throw std::invalid_argument("Mixed-frame equilibrium requires subluminal material speed");
	// For an isotropic rest-frame radiation field, f=4*beta/(3+beta^2),
	// P_parallel/E=(1+3*beta^2)/(3+beta^2), and aT^4=3(1-beta^2)E/(3+beta^2).
	Real const fluxFactor = Real(4) / (Real(3) + beta2);
	RadiationSystem::State state{};
	state.energy() = energy;
	for (int d = 0; d < ndim; ++d) { state.radiativeFlux(d) = energy * velocity[d] * fluxFactor; }
	M1::checkState(RadiationSystem::toCalculationState(state));
	return state;
}

/// Match the original high- or low-order numerical flux to diffusion, retaining
/// the actual cell-center separation in the diffusive gradient. All states,
/// material velocities, and the face normal must use the same basis.
///
/// With tau=mean(chi_L,chi_R)*dx and b=1/(1+tau^2), the energy flux is
/// b*J_thin+(1-b)*(J_material-ALE-chat*delta(E)/(3*tau)). The material/ALE
/// term is upwinded with the M1 equilibrium tensor. It is never attenuated by
/// opacity. The flux-moment transport tends to c*chat*P_eq-w*F_eq, preventing
/// a source-free Hancock predictor's transient F from imposing a spurious
/// anisotropic pressure in opaque cells.
///
/// At fixed opacity, 1-b=O(dx^2), retaining the thin scheme's consistency.
/// At fixed dx and chi~1/epsilon, b=O(epsilon^2): numerical light-speed
/// viscosity is subleading to physical diffusion, whose coefficient is
/// chat/(3*chi), or c/(3*chi) for physical F. This is a spatial asymptotic
/// matching construction, not a transcription of the first-order HLL scheme
/// of Bloch et al. (2020), https://arxiv.org/abs/2011.13926. Their section 6.2
/// explains why attenuating the entire flux fails in moving matter.
///
/// Full diffusion consistency additionally requires source-aware midpoint
/// states and a forced implicit source update containing the pressure-flux
/// divergence. A standalone source/transport Strang split is insufficient.
/// A conservative shared-face realizability limiter must follow this function;
/// where that limiter changes the result, the asymptotic argument does not
/// establish the local diffusion coefficient.
inline RadiationSystem::Flux diffusionCorrectedFlux(RadiationSystem const& system,
	RadiationSystem::Flux const& thinFlux, RadiationSystem::State const& centerLeft,
	RadiationSystem::State const& centerRight, RadiationSystem::State const& faceLeft,
	RadiationSystem::State const& faceRight, Extinction chiLeft, Extinction chiRight,
	MaterialVelocity const& velocityLeft, MaterialVelocity const& velocityRight,
	int normal, units::Length cellWidth, units::Velocity faceSpeed = {}) {
	if (!(chiLeft >= Extinction{} && chiRight >= Extinction{} && cellWidth > units::Length{}) ||
		!units::finite(chiLeft) || !units::finite(chiRight) || !units::finite(cellWidth) || !units::finite(faceSpeed))
		throw std::invalid_argument("Invalid diffusion face extinction or geometry");
	if (normal < 0 || normal >= ndim) throw std::invalid_argument("Invalid diffusion face normal");
	Extinction const chi = Real(0.5) * chiLeft + Real(0.5) * chiRight;
	if (chi == Extinction{}) return thinFlux;
	Real const tau = chi * cellWidth;
	if (!std::isfinite(tau)) throw std::invalid_argument("Nonfinite optical depth at radiation face");
	// Evaluate both weights without squaring a potentially huge optical depth.
	Real const inverse = tau > 1 ? Real(1) / tau : tau;
	Real const small = inverse * inverse / (Real(1) + inverse * inverse);
	Real const thinWeight = tau > 1 ? small : Real(1) - small;
	Real const thickWeight = tau > 1 ? Real(1) - small : small;
	MaterialVelocity velocity{};
	for (int d = 0; d < ndim; ++d) { velocity[d] = Real(0.5) * (velocityLeft[d] + velocityRight[d]); }
	auto const meanEnergy = Real(0.5) * (faceLeft.energy() + faceRight.energy());
	auto const equilibrium = materialEquilibriumMoments(meanEnergy, velocity);
	auto const equilibriumFlux = system.physicalFlux(equilibrium, normal, faceSpeed);
	units::Velocity advectiveSpeed{};
	if (meanEnergy > units::EnergyDensity{}) advectiveSpeed = equilibriumFlux.energy() / meanEnergy;
	auto const upwindEnergy = advectiveSpeed >= units::Velocity{} ? faceLeft.energy() : faceRight.energy();
	auto const materialFlux = advectiveSpeed * upwindEnergy;
	// thickWeight/tau = tau/(1+tau^2), with no 0/0 or overflowing tau^2.
	Real const diffusionWeight = tau > 1 ? inverse / (Real(1) + inverse * inverse) : tau / (Real(1) + tau * tau);
	RadiationSystem::Flux corrected = thinWeight * thinFlux + thickWeight * equilibriumFlux;
	corrected.energy() = thinWeight * thinFlux.energy() + thickWeight * materialFlux
		- (system.reducedLightSpeed() * (diffusionWeight / Real(3))) * (centerRight.energy() - centerLeft.energy());
	return corrected;
}

} // namespace octotigerII::radiation
