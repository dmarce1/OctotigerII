#pragma once
#include "octotigerII/config.hpp"
#include "octotigerII/gravity/solver.hpp"
#include "octotigerII/runtime.hpp"
#include <functional>

namespace octotigerII {
struct Diagnostics {
	Real time = 0, mass = 0, gasEnergy = 0, radiationEnergy = 0;
	Real minimumDensity = 0, minimumPressure = 0, minimumRadiationEnergy = 0,
		 maximumReducedFlux = 0;
	std::array<Real, 3> momentum{};
};
Diagnostics diagnose(std::vector<Snapshot> const& snapshots, Config const& config);
struct RunResult {
	Diagnostics initial, final;
	int steps = 0;
	gravity::Statistics gravityWork;
	std::vector<Snapshot> snapshots;
};
using Observer = std::function<void(std::vector<Snapshot> const&, int, Diagnostics const&)>;
RunResult run(Config const& config, Observer const& observer = {});
} // namespace octotigerII
