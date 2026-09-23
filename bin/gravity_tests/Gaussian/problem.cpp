/** @file
 * @brief gravity-gaussian defaults, validation, and physical CGS initial conditions.
 * @ingroup runtime
 */
#include <cmath>
#include <stdexcept>
#include "octotigerII/problems.hpp"
#include "octotigerII/subgrid/subgrid.hpp"
#include "octotigerII/verification/analytic.hpp"

namespace octotigerII {
verification::Reference problemReference([[maybe_unused]] Config const& c) {
	return verification::sphereReference(c, true);
}

void problemDefaults(Config& c) {
	c.mesh.cells = 4;
	c.mesh.lower = units::Length::from_value(-1e9);
	c.mesh.upper = units::Length::from_value(1e9);
	c.runtime.stopTime = units::Time::from_value(0);
}

void validateProblem(Config const& c) {
	if (c.runtime.stopTime != units::Time{}) throw std::invalid_argument("Static gravity problems require runtime.stopTime=0");
}

/// Initialize this problem in the executable's compile-time dimension.

void initializeProblem(Snapshot& data, Config const& c) {
	using std::exp;


	auto const length = c.mesh.upper - c.mesh.lower;
	data.layout.forEachInterior([&](mesh::Coordinates const& cell, std::size_t i) {
		auto const point = data.layout.cellCenter(data.lower, data.cellWidth, cell);
		Real radius2 = 0;
		for (int axis = 0; axis < ndim; ++axis) {
			Real const r = (point[axis] - (c.mesh.lower + c.mesh.upper) / 2.0) / length;
			radius2 += r * r;
		}
		data.density.values()[i] = units::Density::from_value(radius2 < 0.5 * 0.5 ? 1e4 * exp(-radius2 / (2 * 0.15 * 0.15)) : 0);
	});
}
}	 // namespace octotigerII
