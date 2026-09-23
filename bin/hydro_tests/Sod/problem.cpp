/** @file
 * @brief sod defaults, validation, and physical CGS initial conditions.
 * @ingroup runtime
 */
#include <cmath>
#include <stdexcept>
#include "octotigerII/problems.hpp"
#include "octotigerII/subgrid/subgrid.hpp"
#include "octotigerII/verification/analytic.hpp"

namespace octotigerII {

ProblemBoundary problemBoundary(Config const& c) {
	return verification::sodReference(c).evaluate;
}


verification::Reference problemReference([[maybe_unused]] Config const& c) {
	auto reference = verification::sodReference(c);
	if (c.mesh.boundary.periodic(0)) {
		reference.evaluate = {};
		reference.reason = "The isolated Sod Riemann reference does not apply to periodic x boundaries";
	} else if (c.mesh.boundary.lower[0] == physics::BoundaryCondition::Analytic && c.mesh.boundary.upper[0] == physics::BoundaryCondition::Analytic) {
		reference.validUntil = units::Time::from_value(std::numeric_limits<Real>::infinity());
	}
	return reference;
}

void problemDefaults(Config&) {}

void validateProblem(Config const&) {}

/// Initialize this problem in the executable's compile-time dimension.
/// Shock-tube states follow @ref ref_sod1978 "Sod (1978)".
void initializeProblem(Snapshot& data, Config const& c) {
	hydro::HydroSystem gas(c.hydro.gamma);

	data.layout.forEachInterior([&](mesh::Coordinates const& cell, std::size_t i) {
		auto const point = data.layout.cellCenter(data.lower, data.cellWidth, cell);
		hydro::PrimitiveState primitive{};
		bool const left = point[0] < (c.mesh.lower + c.mesh.upper) / 2.0;
		primitive.density() = units::Density::from_value(left ? 1 : 0.125);
		primitive.pressure() = units::Pressure::from_value(left ? 1 : 0.1);
		data.hydro.values()[i] = gas.conservedState(primitive);
	});
}
}	 // namespace octotigerII
