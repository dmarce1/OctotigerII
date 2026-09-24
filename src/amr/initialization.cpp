#include "octotigerII/amr/hierarchy.hpp"
#include <algorithm>
#include <limits>

namespace octotigerII::amr {

InitialMesh initializeMesh(Config const& config, refinement::Criteria const& criteria, Initializer const& initialize) {
	config.validate();
	// mesh.level is the required startup base resolution. minLevel controls
	// subsequent coarsening and may be lower than this startup base.
	auto selection = config;
	selection.amr.minLevel = config.mesh.level;
	std::vector<Snapshot> candidate{initialize(config, {0, {}})};
	InitialMesh result;
	// No coarsening is allowed, so every change adds leaves and the maximum
	// level bounds the process. New analytic detail can expose neighboring
	// tags late; the number of passes need not equal the maximum tree depth.
	for (;;) {
		Hierarchy const hierarchy(config, candidate);
		auto dt = units::Time::from_value(std::numeric_limits<Real>::infinity());
		for (auto const& block : candidate) {
			if (build::hydro && config.hydroEnabled())
				dt = std::min(dt, hydro::Solver(hydro::HydroSystem(config.hydro.gamma)).stableTimestep(block.hydro, config.timestep.cfl));
			if (build::radiation && config.radiationEnabled())
				dt = std::min(dt, radiation::Solver(radiation::RadiationSystem(config.radiation.lightSpeedRatio * constants::c))
					.stableTimestep(block.radiation, config.timestep.cfl));
		}
		// Recompute the lookahead after every new initialization. A coarse-grid
		// timestep must not pad the final fine mesh by many fine cell widths.
		dt = std::min(dt, config.runtime.stopTime);
		if (!units::finite(dt)) dt = {};
		auto const selected = selectMesh(selection, candidate, hierarchy, Real(config.amr.regridEvery) * dt, criteria, false, nullptr, true);
		if (!selected.changed) {
			result.leaves = selected.leaves;
			result.signalSpeed = selected.signalSpeed;
			result.timestep = dt;
			return result;
		}
		candidate.clear();
		for (auto const& leaf : selected.leaves)
			candidate.push_back(initialize(config, leaf));
	}
}

} // namespace octotigerII::amr
