/** @file
 * @brief collapse defaults, validation, and physical CGS initial conditions.
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
	c.mesh.cells = 4;
	c.mesh.lower = units::Length::from_value(-1e9);
	c.mesh.upper = units::Length::from_value(1e9);
	c.runtime.stopTime = units::Time::from_value(1);
}

void validateProblem(Config const&) {}

/// Initialize this problem in the executable's compile-time dimension.

void initializeProblem(Snapshot& data, Config const& c) {
	using std::exp;

	hydro::HydroSystem gas(c.hydro.gamma);
	auto const length = c.mesh.upper - c.mesh.lower;
	data.layout.forEachInterior([&](mesh::Coordinates const& cell, std::size_t i) {
		auto const point = data.layout.cellCenter(data.lower, data.cellWidth, cell);
		Real radius2 = 0;
		for (int axis = 0; axis < ndim; ++axis) {
			Real const r = (point[axis] - (c.mesh.lower + c.mesh.upper) / 2.0) / length;
			radius2 += r * r;
		}
		hydro::PrimitiveState primitive{};
		primitive.density() = units::Density::from_value(1e4 * (0.01 + exp(-radius2 / (2 * 0.15 * 0.15))));
		primitive.pressure() = units::Pressure::from_value(1e12);
		data.hydro.values()[i] = gas.conservedState(primitive);
	});
}
}	 // namespace octotigerII
