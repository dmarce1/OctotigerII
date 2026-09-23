#include "octotigerII/simulation.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include "octotigerII/verification/analytic.hpp"


namespace octotigerII {


namespace {
#if OCTOTIGERII_GRAVITY
	std::size_t globalIndex(mesh::BlockLocation const& location, mesh::Coordinates cell, int cells, int n) {
		for (int d = 0; d < ndim; ++d)
			cell[d] += location.coordinates[d] * cells;
		return mesh::linearIndex(cell, mesh::filledCoordinates(n));
	}

	gravity::Statistics solveGravity(Runtime& runtime, std::vector<Snapshot> const& snapshots, Config const& c) {
		int const n = c.mesh.cells * (1 << c.mesh.level);
		std::vector<units::Density> density(static_cast<std::size_t>(n) * n * n);
		for (auto const& block : snapshots)
			block.layout.forEachInterior([&](mesh::Coordinates const& cell, std::size_t i) {
				density[globalIndex(block.location, cell, c.mesh.cells, n)] =
					block.hydroEnabled ? block.hydro.values()[i].density() : block.density.values()[i];
			});
		auto solution = gravity::solve(density, n, (c.mesh.upper - c.mesh.lower) / Real(n), c.gravity.multipoleOrder, c.gravity.openingAngle);
		std::vector<std::vector<gravity::State>> fields(snapshots.size());
		for (std::size_t b = 0; b < snapshots.size(); ++b)
			snapshots[b].layout.forEachInterior([&](mesh::Coordinates const& cell, std::size_t) {
				fields[b].push_back(solution.fields[globalIndex(snapshots[b].location, cell, c.mesh.cells, n)]);
			});
		runtime.setGravity(fields);
		return solution.statistics;
	}

#else
	[[maybe_unused]] gravity::Statistics solveGravity(Runtime&, std::vector<Snapshot> const&, Config const&) {
		throw std::logic_error("Gravity is not part of this executable");
	}
#endif

}	 // namespace


Diagnostics diagnose(std::vector<Snapshot> const& snapshots, Config const& c) {
	if (snapshots.empty()) throw std::invalid_argument("Empty snapshot directory");
	Diagnostics d;
	d.time = snapshots.front().time;
	d.minimumDensity = units::Density::from_value(std::numeric_limits<Real>::infinity());
	d.minimumPressure = d.minimumRadiationEnergy = units::Pressure::from_value(std::numeric_limits<Real>::infinity());
	for (auto const& block : snapshots) {
		if (!units::finite(block.time) || d.time != block.time) throw std::runtime_error("Unsynchronized snapshots");
		auto const volume = block.layout.cellMeasure(block.cellWidth);
		block.layout.forEachInterior([&](mesh::Coordinates const&, std::size_t i) {
			if constexpr (build::hydro) {
				hydro::HydroSystem gas(c.hydro.gamma);
				auto const& u = block.hydro.values()[i];
				if (!gas.admissible(u)) throw std::runtime_error("Inadmissible gas state");
				d.mass += volume * u.density();
				d.gasEnergy += volume * u.totalEnergy();
				d.minimumDensity = std::min(d.minimumDensity, u.density());
				d.minimumPressure = std::min(d.minimumPressure, gas.reconstructionVariables(u).pressure());
				for (int axis = 0; axis < ndim; ++axis)
					d.momentum[axis] += volume * u.momentum(axis);
			} else if constexpr (build::gravity)
				d.mass += volume * block.density.values()[i];
			if constexpr (build::radiation) {
				radiation::RadiationSystem rad(c.radiation.lightSpeedRatio * constants::c);
				auto const& u = block.radiation.values()[i];
				if (!rad.admissible(u)) throw std::runtime_error("Inadmissible radiation state");
				d.radiationEnergy += volume * u.energy();
				d.minimumRadiationEnergy = std::min(d.minimumRadiationEnergy, u.energy());
				units::EnergyFlux magnitude{};
				for (int axis = 0; axis < ndim; ++axis)
					magnitude = units::hypot(magnitude, u.radiativeFlux(axis));
				Real const f = u.energy() > units::EnergyDensity{} ? Real(magnitude / (constants::c * u.energy())) : 0;
				d.maximumReducedFlux = std::max(d.maximumReducedFlux, f);
			}
			if constexpr (build::gravity)
				if (!finite(block.gravity.values()[i])) throw std::runtime_error("Nonfinite gravity field");
		});
	}
	if (!c.hydroEnabled()) {
		d.minimumDensity = {};
		d.minimumPressure = {};
	}
	if (!c.radiationEnabled()) d.minimumRadiationEnergy = {};
	if (!units::finite(d.mass) || !units::finite(d.gasEnergy) || !units::finite(d.radiationEnergy)) throw std::runtime_error("Nonfinite global totals");
	return d;
}

RunResult run(Config const& c, Observer const& observer) {
	c.validate();
	verification::reference(c, c.runtime.stopTime);
	Runtime runtime(c);
	RunResult result;
	auto snapshots = runtime.snapshots();
	if constexpr (build::gravity) {
		result.gravityWork = solveGravity(runtime, snapshots, c);
		snapshots = runtime.snapshots();
	}
	result.initial = diagnose(snapshots, c);
	result.final = result.initial;
	if (observer) observer(snapshots, 0, result.initial);
	while (result.final.time < c.runtime.stopTime) {
		if (result.steps >= c.runtime.maxSteps) throw std::runtime_error("Maximum steps reached before requested stop time");
		auto const dt = std::min(runtime.stableTimestep(), c.runtime.stopTime - result.final.time);
		if (!(dt > units::Time{}) || !units::finite(dt) || result.final.time + dt == result.final.time)
			throw std::runtime_error("Timestep cannot advance physical time");
		if constexpr (build::gravity) {
			runtime.kickGravity(dt / 2.0);
			snapshots = runtime.snapshots();
		}
		runtime.advance(dt);
		snapshots = runtime.snapshots();
		if constexpr (build::gravity) {
			auto const work = solveGravity(runtime, snapshots, c);
			result.gravityWork.multipolePairs += work.multipolePairs;
			result.gravityWork.directPairs += work.directPairs;
			runtime.kickGravity(dt / 2.0);
			snapshots = runtime.snapshots();
		}
		++result.steps;
		result.final = diagnose(snapshots, c);
		if (observer) observer(snapshots, result.steps, result.final);
	}
	result.snapshots = std::move(snapshots);
	return result;
}
}	 // namespace octotigerII
