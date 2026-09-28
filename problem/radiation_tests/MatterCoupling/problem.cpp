/** @file
 * @brief Smooth moving gas/radiation relaxation fixture in physical CGS units.
 */
#include <cmath>
#include "octotigerII/problems.hpp"
#include "octotigerII/verification/analytic.hpp"

namespace octotigerII::radiation_matter {

ProblemBoundary problemBoundary(Config const&) { return {}; }
verification::Reference problemReference(Config const&) { return {}; }
void validateProblem(Config const&) {}

void problemDefaults(Config& c) {
	c.mesh.cells = 8;
	c.mesh.level = 0;
	c.mesh.lower = units::Length::from_value(-5e7);
	c.mesh.upper = units::Length::from_value(5e7);
	c.mesh.boundary = finiteVolume::BoundaryConditions::periodic();
	c.hydro.gamma = Real(5) / 3;
	c.radiation.enabled = true;
	c.radiation.opacity = 1;
	c.radiation.initialEnergyRatio = 0.5;
	c.amr.shadowTolerance = 0;
	c.runtime.stopTime = units::Time::from_value(0.001);
}

void initializeProblem(Snapshot& data, Config const& c, bool) {
	hydro::HydroSystem const system(c);
	auto const length = c.mesh.upper - c.mesh.lower;
	data.layout.forEachInterior([&](auto const& cell, std::size_t i) {
		auto const x = data.layout.cellCenter(data.lower, data.cellWidth, cell);
		Real const phase = 2 * piR * Real((x[0] - c.mesh.lower) / length);
		hydro::PrimitiveState gas;
		gas.density() = units::Density::from_value(1e-7 * (1 + 0.1 * std::sin(phase)));
		auto const temperature = units::Temperature::from_value(1e5 * (1 + 0.03 * std::cos(phase)));
		gas.pressure() = gas.density() * constants::k_B * temperature / (c.hydro.meanMolecularWeight * constants::m_u);
		gas.velocity(0) = units::Velocity::from_value(1e6);
		if constexpr (ndim >= 2) gas.velocity(1) = units::Velocity::from_value(-5e5);
		data.hydro.values()[i] = system.conservedState(gas);
		// initialSnapshot supplies the inertial radiation moments after gas initialization.
	});
}
} // namespace octotigerII::radiation_matter
