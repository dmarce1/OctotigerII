/** @file
 * @brief Isolated, spherical Lane-Emden polytrope in hydrostatic equilibrium.
 */
#include <cmath>
#include <stdexcept>
#include "octotigerII/problems.hpp"
#include "octotigerII/problems/laneEmden.hpp"
#include "octotigerII/verification/analytic.hpp"

namespace octotigerII::polytrope {
namespace {
	verification::ExactState state(Config const& c, problems::Polytrope const& star, mesh::PhysicalCoordinates const& position) {
		using std::pow;
		units::Length radius{};
		for (int d = 0; d < ndim; ++d) radius = units::hypot(radius, position[d] - c.star.center[d]);
		auto const profile = star(radius);
		verification::ExactState result;
		result.hydro.density() = std::max(profile.density, c.star.atmosphereFraction * c.star.centralDensity);
		result.hydro.pressure() = std::max(profile.pressure, star.centralPressure() * pow(c.star.atmosphereFraction, 1 + 1 / c.star.polytropicIndex));
		result.gravity.potential() = profile.potential;
		if (radius > units::Length{})
			for (int d = 0; d < ndim; ++d)
				result.gravity.acceleration(d) = -constants::G * profile.enclosedMass / (radius * radius) * Real((position[d] - c.star.center[d]) / radius);
		return result;
	}
}

ProblemBoundary problemBoundary(Config const& c) {
	problems::Polytrope const star(c.star.polytropicIndex, c.star.radius, c.star.centralDensity);
	return [c, star](auto const& x, units::Time) { return state(c, star, x); };
}

verification::Reference problemReference(Config const& c) {
	verification::Reference result;
	result.name = "Lane-Emden hydrostatic equilibrium (tenuous atmosphere approximation)";
	result.evaluate = polytrope::problemBoundary(c);
	return result;
}

void problemDefaults(Config& c) {
	c.mesh.cells = 4;
	c.mesh.lower = -2.0 * c.star.radius;
	c.mesh.upper = 2.0 * c.star.radius;
	c.hydro.gamma = 1 + 1 / c.star.polytropicIndex;
	c.amr.refineDensity = 0.01 * c.star.centralDensity;
	c.amr.shadowTolerance = 0;
	c.runtime.stopTime = units::Time::from_value(0.1);
}

void validateProblem(Config const& c) {
	using std::isfinite;
	using std::pow;
	for (auto x : c.star.center)
		if (!units::finite(x)) throw std::invalid_argument("Star center must be finite");
	if (!(c.star.radius > units::Length{}) || !units::finite(c.star.radius) ||
		!(c.star.centralDensity > units::Density{}) || !units::finite(c.star.centralDensity) ||
		!isfinite(c.star.polytropicIndex) || !(c.star.polytropicIndex > 0 && c.star.polytropicIndex < 5) ||
		!isfinite(c.star.atmosphereFraction) || !(c.star.atmosphereFraction > 0 && c.star.atmosphereFraction < 1))
		throw std::invalid_argument("Invalid star radius, central density, polytropic index, or atmosphere fraction");
	problems::Polytrope const star(c.star.polytropicIndex, c.star.radius, c.star.centralDensity);
	// Resolve the core scale a=R/xi1, which can be much smaller than R for
	// centrally concentrated indices. Intermediate levels may widen it.
	if (c.amr.enabled && c.amr.maxLevel >= 0 && c.amr.maxLevel <= 16 && c.mesh.cells > 0 &&
		star.scaleLength() < (c.mesh.upper - c.mesh.lower) / Real(c.mesh.cells * (1 << c.amr.maxLevel)))
		throw std::invalid_argument("Star core R/xi1 is unresolved at amr.maxLevel; increase maxLevel or mesh.cells, or reduce the domain");
	if (c.star.atmosphereFraction * c.star.centralDensity < units::Density::from_value(1e-14) ||
		star.centralPressure() * pow(c.star.atmosphereFraction, 1 + 1 / c.star.polytropicIndex) < units::Pressure::from_value(1e-14))
		throw std::invalid_argument("Polytrope atmosphere falls below the hydro density or pressure floor");
}

void initializeProblem(Snapshot& data, Config const& c, [[maybe_unused]] bool refinementProbe) {
	// Preserve the core on unresolved startup grids so density tagging sees it,
	// including for highly concentrated polytropic indices.
	problems::Polytrope star(c.star.polytropicIndex, c.star.radius, c.star.centralDensity);
	if (refinementProbe) {
		auto const scale = initialFeatureWidth(star.scaleLength(), data.cellWidth);
		star = problems::Polytrope(c.star.polytropicIndex, c.star.radius * Real(scale / star.scaleLength()), c.star.centralDensity);
	}
	hydro::HydroSystem const gas(c.hydro.gamma);
	data.layout.forEachInterior([&](mesh::Coordinates const& cell, std::size_t i) {
		auto const x = data.layout.cellCenter(data.lower, data.cellWidth, cell);
		data.hydro.values()[i] = gas.conservedState(state(c, star, x).hydro);
	});
}
} // namespace octotigerII
