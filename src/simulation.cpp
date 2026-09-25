#include "octotigerII/simulation.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include "octotigerII/profiling.hpp"
#include "octotigerII/verification/analytic.hpp"

namespace octotigerII {

Diagnostics diagnose(std::vector<Snapshot> const& snapshots, Config const& c) {
	profiling::Region profile("diagnostics");
	if (snapshots.empty()) throw std::invalid_argument("Empty snapshot directory");
	Diagnostics d, compensation;
	auto add = [](auto& sum, auto& error, auto value) {
		auto const y = value - error;
		auto const next = sum + y;
		error = (next - sum) - y;
		sum = next;
	};
	d.time = snapshots.front().time;
	d.minimumDensity = units::Density::from_value(std::numeric_limits<Real>::infinity());
	d.minimumPressure = d.minimumRadiationEnergy = units::Pressure::from_value(std::numeric_limits<Real>::infinity());
	for (auto const& block : snapshots) {
		if (!units::finite(block.time) || d.time != block.time) throw std::runtime_error("Unsynchronized snapshots");
		auto const volume = block.layout.cellMeasure(block.cellWidth);
		block.layout.forEachInterior([&](mesh::Coordinates const&, std::size_t i) {
			if (build::hydro && c.hydroEnabled()) {
				hydro::HydroSystem gas(c.hydro);
				auto const& u = block.hydro.values()[i];
				if (!gas.admissible(u)) throw std::runtime_error("Inadmissible gas state");
				add(d.mass, compensation.mass, volume * u.density());
				add(d.norm.mass, compensation.norm.mass, volume * units::abs(u.density()));
				add(d.norm.gasEnergy, compensation.norm.gasEnergy, volume * units::abs(u.totalEnergy()));
				add(d.gasEnergy, compensation.gasEnergy, volume * u.totalEnergy());
				units::EnergyDensity kinetic{};
				for (int axis = 0; axis < ndim; ++axis)
					kinetic += 0.5 * u.momentum(axis) * u.momentum(axis) / u.density();
				add(d.kineticEnergy, compensation.kineticEnergy, volume * kinetic);
				add(d.thermalEnergy, compensation.thermalEnergy, volume * gas.internalEnergy(u));
				add(d.gasGravityNorm, compensation.gasGravityNorm, volume * units::abs(u.totalEnergy()));
				if (build::gravity && c.gravityEnabled()) {
					auto const potential = 0.5 * volume * u.density() * block.gravity.values()[i].potential();
					add(d.potentialEnergy, compensation.potentialEnergy, potential);
					add(d.gasGravityNorm, compensation.gasGravityNorm, units::abs(potential));
				}
				d.minimumDensity = std::min(d.minimumDensity, u.density());
				d.minimumPressure = std::min(d.minimumPressure, gas.reconstructionVariables(u).pressure());
				for (int axis = 0; axis < ndim; ++axis) {
					add(d.momentum[axis], compensation.momentum[axis], volume * u.momentum(axis));
					add(d.norm.momentum[axis], compensation.norm.momentum[axis], volume * units::abs(u.momentum(axis)));
				}
			} else if (build::gravity && c.gravityEnabled()) {
				add(d.mass, compensation.mass, volume * block.density.values()[i]);
				add(d.norm.mass, compensation.norm.mass, volume * units::abs(block.density.values()[i]));
			}
			if (build::radiation && c.radiationEnabled()) {
				radiation::RadiationSystem rad(c.radiation.lightSpeedRatio * constants::c);
				auto const& u = block.radiation.values()[i];
				if (!rad.admissible(u)) throw std::runtime_error("Inadmissible radiation state");
				add(d.radiationEnergy, compensation.radiationEnergy, volume * u.energy());
				add(d.norm.radiationEnergy, compensation.norm.radiationEnergy, volume * units::abs(u.energy()));
				for (int axis = 0; axis < ndim; ++axis) {
					add(d.radiationFlux[axis], compensation.radiationFlux[axis], volume * u.radiativeFlux(axis));
					add(d.norm.radiationFlux[axis], compensation.norm.radiationFlux[axis], volume * units::abs(u.radiativeFlux(axis)));
				}
				d.minimumRadiationEnergy = std::min(d.minimumRadiationEnergy, u.energy());
				units::EnergyFlux magnitude{};
				for (int axis = 0; axis < ndim; ++axis)
					magnitude = units::hypot(magnitude, u.radiativeFlux(axis));
				Real const f = u.energy() > units::EnergyDensity{} ? Real(magnitude / (constants::c * u.energy())) : 0;
				d.maximumReducedFlux = std::max(d.maximumReducedFlux, f);
			}
			if (build::gravity && c.gravityEnabled())
				if (!finite(block.gravity.values()[i])) throw std::runtime_error("Nonfinite gravity field");
		});
	}
	if (!c.hydroEnabled()) {
		d.minimumDensity = {};
		d.minimumPressure = {};
	}
	if (!c.radiationEnabled()) d.minimumRadiationEnergy = {};
	if (!units::finite(d.mass) || !units::finite(d.gasEnergy) || !units::finite(d.radiationEnergy)) throw std::runtime_error("Nonfinite global totals");
	d.gasGravityEnergy = d.gasEnergy + d.potentialEnergy;
	if (!units::finite(d.gasGravityEnergy)) throw std::runtime_error("Nonfinite gravitational energy integral");
	return d;
}

namespace {
units::Energy reciprocityDefect(std::vector<Snapshot> const& old, std::vector<Snapshot> const& now) {
	if (old.size() != now.size()) throw std::logic_error("Gravity energy endpoints must use the same mesh");
	units::Energy result{}, error{};
	for (std::size_t b = 0; b < old.size(); ++b) {
		if (old[b].location != now[b].location) throw std::logic_error("Gravity energy endpoint block mismatch");
		auto const volume = old[b].layout.cellMeasure(old[b].cellWidth);
		for (std::size_t i = 0; i < old[b].hydro.values().size(); ++i) {
			auto const rho = old[b].hydro.values()[i].density();
			auto const phi = old[b].gravity.values()[i].potential();
			auto const delta = 0.5 * volume * (rho * (now[b].gravity.values()[i].potential() - phi)
				- phi * (now[b].hydro.values()[i].density() - rho));
			auto const y = delta - error;
			auto const next = result + y;
			error = (next - result) - y;
			result = next;
		}
	}
	return result;
}
}

RunResult run(Config const& c, Observer const& observer) {
	profiling::Elapsed profile("simulation.wall_ns");
	c.validate();
	verification::reference(c, c.runtime.stopTime);
	Runtime runtime(c);
	bool const kick = c.hydroEnabled() && (c.gravityEnabled() || c.hasExternalAcceleration());
	RunResult result;
	auto countGravity = [&](gravity::Statistics const& work) {
		result.gravityWork.multipolePairs += work.multipolePairs;
		result.gravityWork.directPairs += work.directPairs;
		result.gravityWork.workerTasks += work.workerTasks;
		result.gravityWork.ewaldPairs += work.ewaldPairs;
		result.gravityWork.reflectedPairs += work.reflectedPairs;
		result.gravityWork.localityCells = work.localityCells;
	};
	[[maybe_unused]] auto solveGravity = [&] { countGravity(runtime.solveGravity()); };
	if (build::gravity && c.gravityEnabled()) solveGravity();
	auto snapshots = runtime.snapshots();
	result.initial = diagnose(snapshots, c);
	result.final = result.initial;
	if (observer) observer(snapshots, 0, result.initial);
	units::Energy solverDefect{}, regridChange{};
	while (result.final.time < c.runtime.stopTime) {
		if (result.steps >= c.runtime.maxSteps) throw std::runtime_error("Maximum steps reached before requested stop time");
		auto dt = std::min(runtime.stableTimestep(), c.runtime.stopTime - result.final.time);
		if (runtime.regrid(dt)) {
			if (build::gravity && c.gravityEnabled()) solveGravity();
			snapshots = runtime.snapshots();
			if (c.hydroEnabled() && c.gravityEnabled())
				regridChange += diagnose(snapshots, c).gasGravityEnergy - result.final.gasGravityEnergy;
			dt = std::min(runtime.stableTimestep(), c.runtime.stopTime - result.final.time);
		}
		if (!(dt > units::Time{}) || !units::finite(dt) || result.final.time + dt == result.final.time)
			throw std::runtime_error("Timestep cannot advance physical time");
		if (c.hydroEnabled() && c.gravityEnabled() && c.amr.enabled && c.timestep.refinement) {
			countGravity(runtime.advanceGravity(dt));
		} else {
			if (c.hydroEnabled() && c.gravityEnabled()) runtime.beginGravityEnergy();
			if (kick) runtime.kickGravity(dt / 2.0);
			runtime.advance(dt);
			if (build::gravity && c.gravityEnabled()) solveGravity();
			if (kick) runtime.kickGravity(dt / 2.0);
			if (c.hydroEnabled() && c.gravityEnabled()) runtime.finishGravityEnergy(dt);
		}
		auto nextSnapshots = runtime.snapshots();
		if (c.hydroEnabled() && c.gravityEnabled()) solverDefect += reciprocityDefect(snapshots, nextSnapshots);
		snapshots = std::move(nextSnapshots);
		++result.steps;
		result.final = diagnose(snapshots, c);
		result.final.boundary = runtime.boundaryTransport();
		result.final.gravityReciprocityDefect = solverDefect;
		result.final.gravityRegridEnergyChange = regridChange;
		if (observer) observer(snapshots, result.steps, result.final);
	}
	result.snapshots = std::move(snapshots);
	return result;
}
}	 // namespace octotigerII
