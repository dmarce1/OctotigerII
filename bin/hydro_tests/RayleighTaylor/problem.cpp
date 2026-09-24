/** @file
 * @brief Rayleigh–Taylor layers in a uniform downward gravitational field.
 * @ingroup runtime
 */
#include <cmath>
#include <stdexcept>
#include "octotigerII/problems.hpp"
#include "octotigerII/subgrid/subgrid.hpp"
#include "octotigerII/verification/analytic.hpp"

namespace octotigerII::rayleigh_taylor {

ProblemBoundary problemBoundary(Config const&) {
	return {};
}

verification::Reference problemReference(Config const&) {
	verification::Reference result;
	result.reason = "The nonlinear Rayleigh-Taylor problem has no implemented exact solution";
	return result;
}

void problemDefaults(Config& c) {
	c.mesh.lower = units::Length::from_value(-0.5);
	c.mesh.upper = units::Length::from_value(0.5);
	c.mesh.boundary = physics::BoundaryConditions::periodic();
	c.mesh.boundary.lower[ndim - 1] = c.mesh.boundary.upper[ndim - 1] = physics::BoundaryCondition::Reflecting;
	c.hydro.acceleration[ndim - 1] = units::Acceleration::from_value(-0.1);
	c.runtime.stopTime = units::Time::from_value(10);
}

void validateProblem(Config const& c) {
	auto const& rt = c.rayleighTaylor;
	if (!units::finite(rt.densityLower) || !units::finite(rt.densityUpper) || !(rt.densityLower >= units::Density::from_value(1e-14)) ||
		!(rt.densityUpper > rt.densityLower) || !units::finite(rt.interfacePressure) || !units::finite(rt.perturbation) || rt.perturbation < units::Velocity{})
		throw std::invalid_argument("Rayleigh-Taylor requires finite upper > lower > 0 densities, pressure, and a nonnegative perturbation");
	for (int axis = 0; axis < ndim - 1; ++axis)
		if (c.hydro.acceleration[axis] != units::Acceleration{}) throw std::invalid_argument("Rayleigh-Taylor gravity must act along the last active axis");
	auto const acceleration = c.hydro.acceleration[ndim - 1];
	if (acceleration > units::Acceleration{}) throw std::invalid_argument("Rayleigh-Taylor gravity must point downward (or be zero for a control run)");
	auto const topPressure = rt.interfacePressure + rt.densityUpper * acceleration * (c.mesh.upper - c.mesh.lower) / 2.0;
	if (!units::finite(topPressure) || !(topPressure >= units::Pressure::from_value(1e-14)))
		throw std::invalid_argument("Rayleigh-Taylor hydrostatic pressure must remain positive through the upper wall");
}

void initializeProblem(Snapshot& data, Config const& c, [[maybe_unused]] bool refinementProbe) {
	using std::cos;

	hydro::HydroSystem gas(c.hydro.gamma);
	auto const length = c.mesh.upper - c.mesh.lower;
	auto const center = (c.mesh.lower + c.mesh.upper) / 2.0;
	auto const& rt = c.rayleighTaylor;
	data.layout.forEachInterior([&](mesh::Coordinates const& cell, std::size_t i) {
		auto const point = data.layout.cellCenter(data.lower, data.cellWidth, cell);
		auto const height = point[ndim - 1] - center;
		hydro::PrimitiveState primitive{};
		primitive.density() = height < units::Length{} ? rt.densityLower : rt.densityUpper;
		primitive.pressure() = rt.interfacePressure + primitive.density() * c.hydro.acceleration[ndim - 1] * height;
		// A single transverse Fourier mode with a smooth envelope vanishing at both walls.
		Real mode = cos(piR * Real(height / length));
		mode *= mode;
		for (int axis = 0; axis < ndim - 1; ++axis)
			mode *= cos(2 * piR * Real((point[axis] - center) / length));
		primitive.velocity(ndim - 1) = rt.perturbation * mode;
		data.hydro.values()[i] = gas.conservedState(primitive);
	});
}

}	 // namespace octotigerII
