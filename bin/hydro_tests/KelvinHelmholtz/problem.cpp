/** @file
 * @brief kelvin-helmholtz defaults, validation, and physical CGS initial conditions.
 * @ingroup runtime
 */
#include <cmath>
#include <stdexcept>
#include "octotigerII/problems.hpp"
#include "octotigerII/subgrid/subgrid.hpp"
#include "octotigerII/verification/analytic.hpp"

namespace octotigerII {
verification::Reference problemReference([[maybe_unused]] Config const& c) {
	return {};
}

void problemDefaults(Config& c) {
	c.mesh.periodic = true;
	c.runtime.stopTime = units::Time::from_value(0.1);
}

void validateProblem(Config const& c) {
	if (!c.mesh.periodic) throw std::invalid_argument("Kelvin–Helmholtz requires periodic boundaries");
}

/// Initialize this problem in the executable's compile-time dimension.

void initializeProblem(Snapshot& data, Config const& c) {
	using std::sin;
	using std::tanh;

	hydro::HydroSystem gas(c.hydro.gamma);
	auto const length = c.mesh.upper - c.mesh.lower;
	data.layout.forEachInterior([&](mesh::Coordinates const& cell, std::size_t i) {
		auto const point = data.layout.cellCenter(data.lower, data.cellWidth, cell);
		Real const x = (point[0] - c.mesh.lower) / length, y = (point[1] - c.mesh.lower) / length;
		Real const layer = 0.5 * (tanh((y - 0.25) / 0.025) - tanh((y - 0.75) / 0.025));
		hydro::PrimitiveState primitive{};
		primitive.density() = units::Density::from_value(1 + layer);
		primitive.pressure() = units::Pressure::from_value(2.5);
		primitive.velocity(0) = units::Velocity::from_value(layer - 0.5);
		primitive.velocity(1) = units::Velocity::from_value(0.01 * sin(4 * piR * x));
		data.hydro.values()[i] = gas.conservedState(primitive);
	});
}
}	 // namespace octotigerII
