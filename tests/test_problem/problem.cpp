/** @file
 * @brief Problem-independent initial state for the all-physics test build.
 * @ingroup runtime
 */
#include "octotigerII/problems.hpp"
#include "octotigerII/subgrid/subgrid.hpp"
#include "octotigerII/verification/analytic.hpp"

namespace octotigerII {

ProblemBoundary problemBoundary(Config const&) {
	return {};
}

verification::Reference problemReference(Config const&) {
	return {};
}

void problemDefaults(Config&) {}

void validateProblem(Config const&) {}

void initializeProblem(Snapshot& data, Config const& c, [[maybe_unused]] bool refinementProbe) {
	if constexpr (build::hydro) {
		hydro::HydroSystem gas(c.hydro.gamma);
		hydro::PrimitiveState primitive{};
		primitive.density() = units::Density::from_value(1e-10);
		primitive.pressure() = units::Pressure::from_value(1e2);
		auto const state = gas.conservedState(primitive);
		data.hydro.values().assign(data.layout.interiorCellCount(), state);
	}
	if constexpr (build::radiation) {
		radiation::RadiationSystem::State state{};
		state.energy() = units::EnergyDensity::from_value(1e-8);
		data.radiation.values().assign(data.layout.interiorCellCount(), state);
	}
}

} 	// namespace octotigerII
