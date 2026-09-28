/** @file
 * @brief Nonfatal optical-depth and reduced-light-speed diagnostic estimates.
 */
#pragma once

#include <algorithm>
#include "octotigerII/hydro/hydroSystem.hpp"
#include "octotigerII/radiation/radiationTransport.hpp"
#include "octotigerII/radiation/opacity.hpp"

namespace octotigerII::radiation {

struct CouplingDiagnostics {
	Real cellOpticalDepth = 0;
	Real trappingParameter = 0;
	Real rslaCriterion = 0;
};

/// Local estimates, not error bounds. The physical scale stays fixed under AMR.
/// All velocities are inertial; rotation of the mesh does not change these values.
inline CouplingDiagnostics couplingDiagnostics(hydro::ConservedState const& gas,
	RadiationSystem::State const& radiation, units::Length cellWidth, Config const& config,
	units::Quantity<2, -1, 0> opacity) {
	auto const inverseMeanFreePath = gas.density() * opacity;
	auto const length = config.radiation.diagnosticLength > 0 ?
		units::Length::from_value(config.radiation.diagnosticLength) : config.mesh.upper - config.mesh.lower;
	Real const cellDepth = units::value(inverseMeanFreePath * cellWidth);
	Real const scaleDepth = units::value(inverseMeanFreePath * length);
	units::Velocity speed{};
	for (int d = 0; d < ndim; ++d) speed = units::hypot(speed, gas.momentum(d) / gas.density());
	auto const gasSoundSpeed = hydro::HydroSystem(config).adiabaticSoundSpeed(gas);
	auto const soundSpeed = units::sqrt(gasSoundSpeed*gasSoundSpeed + (Real(4) / 9) * radiation.energy() / gas.density());
	return {cellDepth, Real(speed / constants::c) * scaleDepth,
		Real((speed + soundSpeed) / (config.radiation.lightSpeedRatio * constants::c)) * std::max(Real(1), scaleDepth)};
}

inline CouplingDiagnostics couplingDiagnostics(hydro::ConservedState const& gas,
	RadiationSystem::State const& radiation, units::Length cellWidth, Config const& config) {
	return couplingDiagnostics(gas, radiation, cellWidth, config,
		opacityLaw(config, Opacity::from_value(config.radiation.opacity))
			.evaluate(gas, hydro::HydroSystem(config)).fluxExtinction);
}

} // namespace octotigerII::radiation
