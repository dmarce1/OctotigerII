/** @file
 * @brief streaming defaults, validation, and physical CGS initial conditions.
 * @ingroup runtime
 */
#include <cmath>
#include <stdexcept>
#include "octotigerII/problems.hpp"
#include "octotigerII/subgrid/subgrid.hpp"
#include "octotigerII/verification/analytic.hpp"

namespace octotigerII {

ProblemBoundary problemBoundary(Config const& c) {
	return [c](mesh::PhysicalCoordinates const& position, units::Time time) { return verification::streamingState(c, position, time); };
}


verification::Reference problemReference([[maybe_unused]] Config const& c) {
	return verification::streamingReference(c);
}

void problemDefaults(Config& c) {
	c.mesh.lower = units::Length::from_value(-3e10);
	c.mesh.upper = units::Length::from_value(3e10);
	c.runtime.stopTime = units::Time::from_value(0.4);
	c.mesh.boundary = physics::BoundaryConditions::periodic();
}

void validateProblem(Config const&) {}

/// Initialize this problem in the executable's compile-time dimension.

void initializeProblem(Snapshot& data, Config const& c) {
	using std::exp;
	using std::round;
	using std::sqrt;


	auto const length = c.mesh.upper - c.mesh.lower;
	data.layout.forEachInterior([&](mesh::Coordinates const& cell, std::size_t i) {
		auto const point = data.layout.cellCenter(data.lower, data.cellWidth, cell);
		Real distance2 = 0;
		for (int axis = 0; axis < ndim; ++axis) {
			auto distance = point[axis] - (c.mesh.lower + 0.25 * length);
			if (c.mesh.boundary.periodic(axis)) distance -= round(distance / length) * length;
			distance2 += distance * distance / (0.08 * length * 0.08 * length);
		}
		auto& state = data.radiation.values()[i];
		state.energy() = units::EnergyDensity::from_value(1e-6 + exp(-0.5 * distance2));
		for (int axis = 0; axis < ndim; ++axis)
			state.radiativeFlux(axis) = constants::c * state.energy() / sqrt(Real(ndim));
	});
}
}	 // namespace octotigerII
