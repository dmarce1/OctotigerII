#include "octotigerII/simulation.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace octotigerII {
namespace {
std::size_t globalIndex(mesh::BlockLocation const& location, mesh::Coordinates cell, int cells,
						int n) {
	for (int d = 0; d < 3; ++d)
		cell[d] += location.coordinates[d] * cells;
	return (static_cast<std::size_t>(cell[2]) * n + cell[1]) * n + cell[0];
}
gravity::Statistics solveGravity(Runtime& runtime, std::vector<Snapshot> const& snapshots,
								 Config const& c) {
	int const n = c.cells * (1 << c.level);
	std::vector<Real> density(static_cast<std::size_t>(n) * n * n);
	for (auto const& block : snapshots)
		block.layout.forEachInterior([&](mesh::Coordinates const& cell, std::size_t i) {
			density[globalIndex(block.location, cell, c.cells, n)] =
				block.hydroEnabled ? block.hydro.values()[i].density() : block.density.values()[i];
		});
	auto solution =
		gravity::solve(density, n, (c.upper - c.lower) / n, c.multipoleOrder, c.openingAngle);
	std::vector<std::vector<gravity::State>> fields(snapshots.size());
	for (std::size_t b = 0; b < snapshots.size(); ++b)
		snapshots[b].layout.forEachInterior([&](mesh::Coordinates const& cell, std::size_t) {
			fields[b].push_back(
				solution.fields[globalIndex(snapshots[b].location, cell, c.cells, n)]);
		});
	runtime.setGravity(fields);
	return solution.statistics;
}
} // namespace
Diagnostics diagnose(std::vector<Snapshot> const& snapshots, Config const& c) {
	if (snapshots.empty())
		throw std::invalid_argument("Empty snapshot directory");
	Diagnostics d;
	d.time = snapshots.front().time;
	d.minimumDensity = d.minimumPressure = d.minimumRadiationEnergy =
		std::numeric_limits<Real>::infinity();
	hydro::HydroSystem gas(c.gamma);
	radiation::RadiationSystem rad(c.lightSpeedRatio * physicalLightSpeed);
	for (auto const& block : snapshots) {
		if (!transportExchange::sameTime(d.time, block.time))
			throw std::runtime_error("Unsynchronized snapshots");
		Real const volume = block.layout.cellMeasure(block.cellWidth);
		block.layout.forEachInterior([&](mesh::Coordinates const&, std::size_t i) {
			if (block.hydroEnabled) {
				auto const& u = block.hydro.values()[i];
				if (!gas.admissible(u))
					throw std::runtime_error("Inadmissible gas state");
				d.mass += volume * u.density();
				d.gasEnergy += volume * u.totalEnergy();
				d.minimumDensity = std::min(d.minimumDensity, u.density());
				d.minimumPressure =
					std::min(d.minimumPressure, gas.reconstructionVariables(u).pressure());
				for (int axis = 0; axis < 3; ++axis)
					d.momentum[axis] += volume * u.momentum(axis);
			} else if (block.gravityEnabled)
				d.mass += volume * block.density.values()[i];
			if (block.radiationEnabled) {
				auto const& u = block.radiation.values()[i];
				if (!rad.admissible(u))
					throw std::runtime_error("Inadmissible radiation state");
				d.radiationEnergy += volume * u[0];
				d.minimumRadiationEnergy = std::min(d.minimumRadiationEnergy, u[0]);
				Real const f = u[0] > 0 ? std::hypot(u[1], u[2], u[3]) / u[0] : 0;
				d.maximumReducedFlux = std::max(d.maximumReducedFlux, f);
			}
			if (block.gravityEnabled)
				for (int f = 0; f < 4; ++f)
					if (!std::isfinite(block.gravity.values()[i][f]))
						throw std::runtime_error("Nonfinite gravity field");
		});
	}
	if (!c.hydroEnabled())
		d.minimumDensity = d.minimumPressure = 0;
	if (!c.radiationEnabled())
		d.minimumRadiationEnergy = 0;
	if (!std::isfinite(d.mass) || !std::isfinite(d.gasEnergy) || !std::isfinite(d.radiationEnergy))
		throw std::runtime_error("Nonfinite global totals");
	return d;
}
RunResult run(Config const& c, Observer const& observer) {
	c.validate();
	Runtime runtime(c);
	RunResult result;
	auto snapshots = runtime.snapshots();
	if (c.gravityEnabled()) {
		result.gravityWork = solveGravity(runtime, snapshots, c);
		snapshots = runtime.snapshots();
	}
	result.initial = diagnose(snapshots, c);
	result.final = result.initial;
	if (observer)
		observer(snapshots, 0, result.initial);
	while (result.final.time < c.stopTime) {
		if (result.steps >= c.maxSteps)
			throw std::runtime_error("Maximum steps reached before requested stop time");
		Real const dt = std::min(runtime.stableTimestep(), c.stopTime - result.final.time);
		if (!(dt > 0) || !std::isfinite(dt) || result.final.time + dt == result.final.time)
			throw std::runtime_error("Timestep cannot advance physical time");
		if (c.gravityEnabled()) {
			runtime.kickGravity(dt / 2);
			snapshots = runtime.snapshots();
		}
		runtime.advance(snapshots, dt);
		snapshots = runtime.snapshots();
		if (c.gravityEnabled()) {
			auto const work = solveGravity(runtime, snapshots, c);
			result.gravityWork.multipolePairs += work.multipolePairs;
			result.gravityWork.directPairs += work.directPairs;
			runtime.kickGravity(dt / 2);
			snapshots = runtime.snapshots();
		}
		++result.steps;
		result.final = diagnose(snapshots, c);
		if (observer)
			observer(snapshots, result.steps, result.final);
	}
	result.snapshots = std::move(snapshots);
	return result;
}
} // namespace octotigerII
