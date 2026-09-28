// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#include "octotigerII/amr/hierarchy.hpp"
#include <bit>
#include "octotigerII/composition/transport.hpp"
#include "octotigerII/amr/interpolation.hpp"
#include "octotigerII/problems.hpp"
#include "octotigerII/verification/analytic.hpp"
#include "octotigerII/radiation/coupledPatch.hpp"

namespace octotigerII::amr {
namespace {
mesh::PhysicalCoordinates logicalCenter(Config const& config, mesh::BlockLocation const& cell) {
	mesh::PhysicalCoordinates x{};
	auto const width = (config.mesh.upper - config.mesh.lower) / Real(std::uint64_t(1) << cell.level);
	for (int d = 0; d < ndim; ++d) x[d] = config.mesh.lower + (cell.coordinates[d] + 0.5) * width;
	return x;
}
template <typename System>
units::Velocity meshSignalSpeed(System const& system, typename System::State const& state, int axis,
	Config const& config, mesh::BlockLocation const& cell, units::Time time) {
	physics::RotatingFrame const frame(config.frame.omega);
	auto const meshSpeed = frame.normalSpeed(logicalCenter(config, cell), axis);
	auto result = system.maximumSignalSpeed(frame.toGridState(state, time), axis, meshSpeed);
	if (frame.active() && axis < 2) {
		auto const width = (config.mesh.upper - config.mesh.lower) / Real(std::uint64_t(1) << cell.level);
		result += 0.5 * units::abs(frame.omega()) * width;
	}
	return result;
}
}


Hierarchy::Values& Hierarchy::Values::operator+=(Values const& other) {
	hydro += other.hydro;
	radiation += other.radiation;
	gravity += other.gravity;
	density += other.density;
	gasGravityEnergy += other.gasGravityEnergy;
	if (species.empty()) species.resize(other.species.size());
	if (species.size() != other.species.size()) throw std::logic_error("Shadow species size mismatch");
	for (std::size_t s = 0; s < species.size(); ++s) species[s] += other.species[s];
	return *this;
}

Hierarchy::Values Hierarchy::Values::operator*(Real weight) const {
	Values result{hydro * weight, radiation * weight, gravity * weight, density * weight, gasGravityEnergy * weight, species};
	for (auto& q : result.species) q *= weight;
	return result;
}

Hierarchy::Hierarchy(Config const& config, std::vector<Snapshot> const& leaves)
  : config_(config)
  , radiationMaterial_(problemRadiationMaterial(config))
  , blockBits_(std::countr_zero(unsigned(config.mesh.cells))) {
	if (leaves.empty()) throw std::invalid_argument("Cannot build an empty shadow hierarchy");
	time_ = leaves.front().time;
	std::unordered_map<mesh::BlockLocation, bool, mesh::BlockLocationHash> shadowBlocks;
	for (auto const& block : leaves) {
		auto location = block.location;
		for (;;) {
			shadowBlocks.emplace(location, true);
			if (location.isRoot()) break;
			location = location.parent();
		}
		if (block.time != time_) throw std::invalid_argument("Shadow hierarchy requires synchronized leaves");
		block.layout.forEachInterior([&](auto const& coordinate, std::size_t i) {
			Values value;
			for (auto const& s : block.species) value.species.push_back(s.values()[i]);
			if (build::hydro && config_.hydroEnabled()) value.hydro = block.hydro.values()[i];
			if (build::radiation && config_.radiationEnabled()) value.radiation = block.radiation.values()[i];
			if (build::gravity && config_.gravityEnabled()) {
				value.gravity = block.gravity.values()[i];
				if (!config_.hydroEnabled()) value.density = block.density.values()[i];
			}
			if (build::hydro && config_.hydroEnabled()) value.density = value.hydro.density();
			if (config_.hydroEnabled() && config_.gravityEnabled())
				value.gasGravityEnergy = value.hydro.totalEnergy() + 0.5 * value.density * value.gravity.potential();
			auto scale = [](auto& maximum, auto const& state) {
				state.forEach([&](auto f, auto v) { maximum.template get<f>() = std::max(maximum.template get<f>(), units::abs(v)); });
			};
			scale(scale_.hydro, value.hydro);
			scale(scale_.radiation, value.radiation);
			mesh::BlockLocation cell{cellLevel(block.location.level), {}};
			for (int d = 0; d < ndim; ++d)
				cell.coordinates[d] = block.location.coordinates[d] * config.mesh.cells + coordinate[d];
			for (;;) {
				cells_[cell] += value;
				if (cell.isRoot()) break;
				cell = cell.parent();
				value = value * (Real(1) / (1 << ndim));
			}
		});
	}
	for (auto const& [block, present] : shadowBlocks) {
		(void) present;
		shadowBlocks_.push_back(block);
	}
}

void Hierarchy::advance(units::Time dt) {
	if (!config_.timestep.refinement || config_.gravityEnabled() || config_.hasExternalAcceleration()) {
		advanceOnce(dt);
		return;
	}
	// The independent shadow hierarchy includes fine covered cells. A coarse
	// synchronization interval must never be applied to all of them in one step.
	auto remaining = dt;
	auto step = dt;
	while (remaining > units::Time{}) {
		std::array<units::Velocity, ndim> speed{};
		auto width = config_.mesh.upper - config_.mesh.lower;
		for (auto const& [cell, value] : cells_) {
			width = std::min(width, (config_.mesh.upper - config_.mesh.lower) / Real(std::uint64_t(1) << cell.level));
			for (int d = 0; d < ndim; ++d) {
				if (build::hydro && config_.hydroEnabled()) speed[d] = std::max(speed[d], meshSignalSpeed(hydro::HydroSystem(config_), value.hydro, d, config_, cell, time_));
				if (build::radiation && config_.radiationEnabled()) speed[d] = std::max(speed[d], meshSignalSpeed(radiation::RadiationSystem(config_.radiation.lightSpeedRatio * constants::c), value.radiation, d, config_, cell, time_));
			}
		}
		units::InverseTime rate{};
		for (auto v : speed) rate += v / width;
		if (!(rate > units::InverseTime{}) || !units::finite(rate)) throw std::runtime_error("Invalid shadow CFL rate");
		auto const limit = std::min(config_.timestep.cfl / rate, physics::RotatingFrame(config_.frame.omega).maximumTimestep());
		while (step > limit || step > remaining) step /= 2;
		if (!(step > units::Time{}) || time_ + step == time_) throw std::runtime_error("Shadow timestep cannot advance time");
		advanceOnce(step);
		remaining -= step;
		if (remaining <= 32 * epsilonR * dt) remaining = {};
	}
}

void Hierarchy::advanceGravity(units::Time dt) {
	if (!(dt > units::Time{}) || !units::finite(dt))
		throw std::invalid_argument("Shadow gravity timestep must be positive and finite");
	auto const end = time_ + dt;
	if (!units::finite(end) || end == time_) throw std::runtime_error("Shadow timestep cannot advance time");
	Real elapsed = 0, fraction = 1;
	while (elapsed < 1) {
		std::array<units::Velocity, ndim> speed{};
		auto width = config_.mesh.upper - config_.mesh.lower;
		units::Acceleration maximumAcceleration{};
		auto const external = physics::RotatingFrame(config_.frame.omega).toGrid(config_.hydro.acceleration, time_);
		for (auto const& [cell, value] : cells_) {
			width = std::min(width, (config_.mesh.upper - config_.mesh.lower) / Real(std::uint64_t(1) << cell.level));
			units::Acceleration acceleration{};
			for (int d = 0; d < ndim; ++d) {
				if (build::hydro && config_.hydroEnabled()) {
					speed[d] = std::max(speed[d], meshSignalSpeed(hydro::HydroSystem(config_), value.hydro, d, config_, cell, time_));
					acceleration += units::abs(external[d]
						+ (config_.gravityEnabled() ? value.gravity.acceleration(d) : units::Acceleration{}));
				}
				if (build::radiation && config_.radiationEnabled())
					speed[d] = std::max(speed[d], meshSignalSpeed(radiation::RadiationSystem(config_.radiation.lightSpeedRatio * constants::c), value.radiation, d, config_, cell, time_));
			}
			maximumAcceleration = std::max(maximumAcceleration, acceleration);
		}
		units::InverseTime rate{};
		for (auto v : speed) rate += v / width;
		if (!(rate > units::InverseTime{}) || !units::finite(rate)) throw std::runtime_error("Invalid shadow CFL rate");
		auto limit = std::min(config_.timestep.cfl / rate, physics::RotatingFrame(config_.frame.omega).maximumTimestep());
		if (maximumAcceleration > units::Acceleration{}) {
			// Include the speed gained during the source step, using the same
			// displacement and acceleration constraints as the physical cells.
			auto const a = 0.5 * maximumAcceleration / width;
			limit = std::min(limit, 2.0 * config_.timestep.cfl /
				(rate + units::sqrt(rate * rate + 4.0 * a * config_.timestep.cfl)));
			limit = std::min(limit, 0.2 * units::sqrt(width / maximumAcceleration));
		}
		while (fraction * dt > limit || fraction > 1 - elapsed) fraction /= 2;
		auto const step = fraction * dt;
		if (!(step > units::Time{}) || time_ + step == time_) throw std::runtime_error("Shadow timestep cannot advance time");
		// These are independent error-estimation states, not physical leaves.
		// Their field stays fixed until refreshGravity; no extra gravity solve
		// or physical energy/momentum budget is introduced by this forecast.
		kick(step / 2.0);
		advanceOnce(step);
		kick(step / 2.0);
		elapsed += fraction;
	}
	time_ = end;
}

void Hierarchy::advanceOnce(units::Time dt) {
	// A half-sized patch on every block (including ancestors) covers exactly
	// its coarse shadow cells. Different levels advance their own solution.
	// All stencils read the immutable old hierarchy; publication follows all work.
	std::unordered_map<mesh::BlockLocation, Values, mesh::BlockLocationHash> next;
	int const n = config_.mesh.cells / 2;
	bool const coupled = build::hydro && build::radiation && config_.hydroEnabled() && config_.radiationEnabled()
		&& (config_.radiation.opacity > 0 || problemHasRadiationMaterial(config_));
	int const ghosts = coupled ? 4 : 2;
	for (auto const& block : shadowBlocks_) {
		int const level = cellLevel(block.level) - 1;
		auto const width = (config_.mesh.upper - config_.mesh.lower) / Real(1 << level);
		mesh::PhysicalCoordinates lower{};
		for (int d = 0; d < ndim; ++d)
			lower[d] = config_.mesh.lower + Real(block.coordinates[d] * n) * width;
		mesh::MeshLayout layout(n, ghosts);
		mesh::PatchData<hydro::ConservedState> gas(layout, width, lower);
		mesh::PatchData<radiation::RadiationSystem::State> radiation(layout, width, lower);
		mesh::forEachCoordinate(layout.extents(), [&](auto const& cell) {
			mesh::BlockLocation global{level, {}};
			for (int d = 0; d < ndim; ++d)
				global.coordinates[d] = block.coordinates[d] * n + cell[d] - ghosts;
			auto const value = average(global);
			if (build::hydro && config_.hydroEnabled()) gas.atStorage(cell) = value.hydro;
			if (build::radiation && config_.radiationEnabled()) radiation.atStorage(cell) = value.radiation;
		});
		auto destination = [&](auto const& cell) -> Values& {
			mesh::BlockLocation global{level, {}};
			for (int d = 0; d < ndim; ++d)
				global.coordinates[d] = block.coordinates[d] * n + cell[d];
			return next.try_emplace(global, average(global)).first->second;
		};
		hydro::Solver::Workspace gasWork;
		if constexpr (build::hydro && build::radiation) if (coupled) {
			radiation::CoupledPatchWorkspace work;
			auto midpointBoundary = [&](hydro::Fields& midGas, radiation::Fields& midRad, units::Time at) {
				auto const oldGas = midGas;
				auto const oldRad = midRad;
				physics::RotatingFrame const frame(config_.frame.omega);
				hydro::HydroSystem const gasSystem(config_);
				radiation::RadiationSystem const radSystem(config_.radiation.lightSpeedRatio * constants::c,
					{config_.radiation.closedBoundary, config_.mesh.lower, config_.mesh.upper});
				mesh::forEachCoordinate(layout.extents(), [&](auto const& cell) {
					mesh::BlockLocation global{level, {}};
					for (int d = 0; d < ndim; ++d) global.coordinates[d] = block.coordinates[d] * n + cell[d] - ghosts;
					auto const mapped = config_.mesh.boundary.map(global.coordinates, 1 << level);
					bool exterior = false;
					for (int d = 0; d < ndim; ++d) exterior = exterior || (!config_.mesh.boundary.periodic(d) &&
						(global.coordinates[d] < 0 || global.coordinates[d] >= (1 << level)));
					if (!exterior) return;
					auto const position = logicalCenter(config_, global);
					if (mapped.analytic) {
						auto const prescribed = problemBoundary(config_)(position, at);
						midGas.atStorage(cell) = gasSystem.conservedState(prescribed.hydro);
						midRad.atStorage(cell) = prescribed.radiation;
						return;
					}
					mesh::Coordinates source{};
					for (int d = 0; d < ndim; ++d) source[d] =
						(config_.mesh.boundary.periodic(d) ? global.coordinates[d] : mapped.source[d]) - block.coordinates[d] * n + ghosts;
					midGas.atStorage(cell) = physics::transformBoundary(oldGas.atStorage(source), mapped.reflectionMask,
						mapped.outflowLowerMask, mapped.outflowUpperMask, gasSystem, frame, position, at);
					midRad.atStorage(cell) = physics::transformBoundary(oldRad.atStorage(source), mapped.reflectionMask,
						mapped.outflowLowerMask, mapped.outflowUpperMask, radSystem, frame, position, at);
				});
			};
			radiation::advanceCoupledPatch(gas, radiation, hydro::HydroSystem(config_),
				radiation::RadiationSystem(config_.radiation.lightSpeedRatio * constants::c,
					{config_.radiation.closedBoundary, config_.mesh.lower, config_.mesh.upper}),
				radiation::Opacity::from_value(config_.radiation.opacity), dt, work,
				[&](auto const& cell, auto const& g, auto const& r) {
					auto& target = destination(cell);
					target.hydro = g;
					target.radiation = r;
					target.density = g.density();
				}, physics::RotatingFrame(config_.frame.omega), time_, midpointBoundary, radiationMaterial_);
			gasWork = std::move(work.hydro);
		}
		if constexpr (build::hydro) if (config_.hydroEnabled()) {
			if (!coupled) {
				hydro::Solver(hydro::HydroSystem(config_), physics::RotatingFrame(config_.frame.omega), time_).advanceInto(gas, dt, gasWork, [&](auto const& cell, auto const& value) {
					auto& target = destination(cell);
					target.hydro = value;
					hydro::HydroSystem(config_).synchronize(target.hydro);
					target.density = value.density();
				});
			}
			if (config_.massFractions.enabled) {
				auto read = [&](std::size_t s, auto const& cell) {
					mesh::BlockLocation global{level, {}};
					for (int d = 0; d < ndim; ++d) global.coordinates[d] = block.coordinates[d] * n + cell[d];
					return average(global).species.at(s);
				};
				auto flux = composition::fluxes(config_.massFractions, layout, read,
					[&](int d, auto const& face) { return gasWork.fluxes[d][layout.faceIndex(d, face)].template get<0>(); });
				layout.forEachInterior([&](auto const& cell, std::size_t) {
					auto& target = destination(cell);
					for (std::size_t s = 0; s < flux.size(); ++s)
						target.species[s] = composition::update(read(s, cell), flux[s], layout, cell, dt / width);
					target.density = target.hydro.density() = composition::totalDensity(config_.massFractions, target.species);
					if constexpr (build::hydro) hydro::HydroSystem(config_).setComposition(target.hydro,target.species,config_.massFractions);
				});
			}
		}
		if constexpr (build::radiation) if (config_.radiationEnabled() && !coupled) {
			radiation::Solver::Workspace work;
			radiation::Solver(radiation::RadiationSystem(config_.radiation.lightSpeedRatio * constants::c,
				{config_.radiation.closedBoundary, config_.mesh.lower, config_.mesh.upper}), physics::RotatingFrame(config_.frame.omega), time_)
				.advanceInto(radiation, dt, work, [&](auto const& cell, auto const& value) { destination(cell).radiation = value; });
		}
	}
	for (auto const& [cell, value] : next)
		cells_.at(cell) = value;
	time_ += dt;
}

void Hierarchy::kick(units::Time dt) {
	if (build::hydro && config_.hydroEnabled())
		for (auto& [cell, value] : cells_) {
			(void) cell;
			auto const gravity = physics::RotatingFrame(config_.frame.omega).toInertialState(value.gravity, time_);
			for (int d = 0; d < ndim; ++d) {
				auto const acceleration = config_.hydro.acceleration[d] + (config_.gravityEnabled() ? gravity.acceleration(d) : units::Acceleration{});
				auto const impulse = dt * value.hydro.density() * acceleration;
				value.hydro.totalEnergy() += impulse * (value.hydro.momentum(d) + 0.5 * impulse) / value.hydro.density();
				value.hydro.momentum(d) += impulse;
			}
			hydro::HydroSystem(config_).synchronize(value.hydro);
		}
}

void Hierarchy::refreshLeaves(std::vector<Snapshot> const& leaves) {
	// Active cells supply boundary data for independent coarse evolution. Covered
	// cells are deliberately preserved until the next error estimate and regrid.
	for (auto const& block : leaves)
		block.layout.forEachInterior([&](auto const& coordinate, std::size_t i) {
			mesh::BlockLocation cell{cellLevel(block.location.level), {}};
			for (int d = 0; d < ndim; ++d)
				cell.coordinates[d] = block.location.coordinates[d] * config_.mesh.cells + coordinate[d];
			auto& value = cells_.at(cell);
			value.species.clear();
			for (auto const& s : block.species) value.species.push_back(s.values()[i]);
			if (build::hydro && config_.hydroEnabled()) {
				value.hydro = block.hydro.values()[i];
				value.density = value.hydro.density();
			}
			if (build::radiation && config_.radiationEnabled()) value.radiation = block.radiation.values()[i];
			auto maximum = [](auto& scale, auto const& state) {
				state.forEach([&](auto f, auto q) { scale.template get<f>() = std::max(scale.template get<f>(), units::abs(q)); });
			};
			maximum(scale_.hydro, value.hydro);
			maximum(scale_.radiation, value.radiation);
		});
}

void Hierarchy::refreshGravity(std::vector<Snapshot> const& leaves) {
	Hierarchy current(config_, leaves);
	for (auto& [cell, value] : cells_)
		value.gravity = current.cells_.at(cell).gravity;
}

Hierarchy::Values Hierarchy::average(mesh::BlockLocation cell) const {
	auto const position = logicalCenter(config_, cell);
	physics::RotatingFrame const frame(config_.frame.omega);
	auto const mapped = config_.mesh.boundary.map(cell.coordinates, 1 << cell.level);
	if (mapped.analytic) {
		mesh::PhysicalCoordinates position{};
		auto const h = (config_.mesh.upper - config_.mesh.lower) / Real(1 << cell.level);
		for (int d = 0; d < ndim; ++d)
			position[d] = config_.mesh.lower + (cell.coordinates[d] + 0.5) * h;
		auto const sample = problemBoundary(config_)(position, time_);
		Values result;
		if (build::hydro && config_.hydroEnabled()) result.hydro = hydro::HydroSystem(config_).conservedState(sample.hydro);
		if (build::radiation && config_.radiationEnabled()) result.radiation = sample.radiation;
		if (build::hydro && config_.hydroEnabled()) result.density = result.hydro.density();
		if (config_.massFractions.enabled) result.species = composition::initialDensities(config_.massFractions, result.density);
		return result;
	}
	cell.coordinates = mapped.source;
	auto found = cells_.find(cell);
	while (found == cells_.end()) {
		if (cell.isRoot()) throw std::logic_error("Incomplete shadow hierarchy");
		cell = cell.parent();
		found = cells_.find(cell);
	}
	auto result = found->second;
	if (build::hydro && config_.hydroEnabled())
		result.hydro = physics::transformBoundary(
			result.hydro, mapped.reflectionMask, mapped.outflowLowerMask, mapped.outflowUpperMask, hydro::HydroSystem(config_), frame, position, time_);
	if (build::radiation && config_.radiationEnabled())
		result.radiation = physics::transformBoundary(result.radiation, mapped.reflectionMask, mapped.outflowLowerMask, mapped.outflowUpperMask,
			radiation::RadiationSystem(config_.radiation.lightSpeedRatio * constants::c,
				{config_.radiation.closedBoundary, config_.mesh.lower, config_.mesh.upper}), frame, position, time_);
	return result;
}

Hierarchy::Values Hierarchy::reconstructFrom(mesh::BlockLocation coarse, mesh::BlockLocation fine) const {
	auto const center = average(coarse);
	std::array<hydro::ConservedState, ndim> gasSlopes{};
	std::array<radiation::RadiationSystem::State, ndim> radiationSlopes{};
	std::array<Real, ndim> offset{};
	std::array<units::EnergyDensity, ndim> energySlopes{};
	Real const ratio = Real(1 << (fine.level - coarse.level));
	for (int d = 0; d < ndim; ++d) {
		auto left = coarse, right = coarse;
		--left.coordinates[d];
		++right.coordinates[d];
		auto const a = average(left), b = average(right);
		energySlopes[d] = physics::limitedSlope(center.gasGravityEnergy - a.gasGravityEnergy,
			b.gasGravityEnergy - center.gasGravityEnergy, physics::Limiter::Minmod);
		gasSlopes[d] = slope(a.hydro, center.hydro, b.hydro);
		radiationSlopes[d] = slope(a.radiation, center.radiation, b.radiation);
		offset[d] = (fine.coordinates[d] + 0.5) / ratio - (coarse.coordinates[d] + 0.5);
	}
	auto result = center;
	for (int d = 0; d < ndim; ++d) result.gasGravityEnergy += offset[d] * energySlopes[d];
	if (build::hydro && config_.hydroEnabled()) result.hydro = interpolate(center.hydro, gasSlopes, offset, hydro::HydroSystem(config_));
	if (build::radiation && config_.radiationEnabled())
		result.radiation = interpolate(center.radiation, radiationSlopes, offset, radiation::RadiationSystem(config_.radiation.lightSpeedRatio * constants::c));
	if (build::hydro && config_.hydroEnabled()) result.density = result.hydro.density();
	if (config_.massFractions.enabled) {
		// Inject concentrations and prolong hydro density. Summed child partial
		// densities preserve each parent species mass and the hydro mass together.
		for (std::size_t s = 0; s < result.species.size(); ++s)
			result.species[s] = center.species[s] * Real(result.density / center.density);
		result.density = result.hydro.density() = composition::totalDensity(config_.massFractions, result.species);
		if constexpr (build::hydro) hydro::HydroSystem(config_).setComposition(result.hydro,result.species,config_.massFractions);
	}
	// Prescribed gravity-only densities use conservative, positive injection.
	return result;
}

Hierarchy::Values Hierarchy::reconstruct(mesh::BlockLocation cell) const {
	if (auto found = cells_.find(cell); found != cells_.end()) return found->second;
	auto ancestor = cell;
	do {
		ancestor = ancestor.parent();
	} while (!cells_.contains(ancestor));
	return reconstructFrom(ancestor, cell);
}

Hierarchy::Values Hierarchy::shadow(mesh::BlockLocation cell) const {
	return cell.isRoot() ? average(cell) : reconstructFrom(cell.parent(), cell);
}

refinement::CellView Hierarchy::sample(mesh::BlockLocation block, mesh::Coordinates const& coordinate) const {
	mesh::BlockLocation cell{cellLevel(block.level), {}};
	for (int d = 0; d < ndim; ++d)
		cell.coordinates[d] = block.coordinates[d] * config_.mesh.cells + coordinate[d];
	auto const value = average(cell), coarse = shadow(cell);
	refinement::CellView result;
	result.width = (config_.mesh.upper - config_.mesh.lower) / Real(1 << cell.level);
	result.volume = mesh::MeshLayout(config_.mesh.cells).cellMeasure(result.width);
	result.time = time_;
	result.level = block.level;
	result.mass = value.density * result.volume;
	result.hydro.value = value.hydro;
	result.hydro.shadow = coarse.hydro;
	result.hydro.scale = scale_.hydro;
	result.radiation.value = value.radiation;
	result.radiation.shadow = coarse.radiation;
	result.radiation.scale = scale_.radiation;
	result.hasHydro = config_.hydroEnabled();
	result.hasRadiation = config_.radiationEnabled();
	for (int d = 0; d < ndim; ++d) {
		result.center[d] = config_.mesh.lower + (cell.coordinates[d] + 0.5) * result.width;
		auto left = cell, right = cell;
		--left.coordinates[d];
		++right.coordinates[d];
		auto const a = average(left), b = average(right);
		result.hydro.gradient[d] = (b.hydro - a.hydro) / (2.0 * result.width);
		result.radiation.gradient[d] = (b.radiation - a.radiation) / (2.0 * result.width);
		if (build::hydro && config_.hydroEnabled()) result.signalSpeed[d] = meshSignalSpeed(hydro::HydroSystem(config_), value.hydro, d, config_, cell, time_);
		if (build::hydro && config_.hydroEnabled())
			result.acceleration[d] = physics::RotatingFrame(config_.frame.omega).toGrid(config_.hydro.acceleration, time_)[d] + (config_.gravityEnabled() ? value.gravity.acceleration(d) : units::Acceleration{});
		if (build::radiation && config_.radiationEnabled()) {
			auto speed = config_.radiation.lightSpeedRatio * constants::c
				+ units::abs(physics::RotatingFrame(config_.frame.omega).normalSpeed(logicalCenter(config_, cell), d));
			if (d < 2) speed += 0.5 * units::abs(config_.frame.omega) * result.width;
			result.signalSpeed[d] = std::max(result.signalSpeed[d], speed);
		}
	}
	return result;
}

std::vector<units::EnergyDensity> Hierarchy::transferGasGravityEnergy(mesh::BlockLocation location) const {
	mesh::MeshLayout const layout(config_.mesh.cells);
	std::vector<units::EnergyDensity> result(layout.cellCount());
	layout.forEachInterior([&](auto const& coordinate, std::size_t i) {
		mesh::BlockLocation cell{cellLevel(location.level), {}};
		for (int d = 0; d < ndim; ++d)
			cell.coordinates[d] = location.coordinates[d] * config_.mesh.cells + coordinate[d];
		result[i] = reconstruct(cell).gasGravityEnergy;
	});
	return result;
}

Snapshot Hierarchy::transfer(mesh::BlockLocation location) const {
	Snapshot result;
	result.location = location;
	result.layout = mesh::MeshLayout(config_.mesh.cells);
	result.cellWidth = (config_.mesh.upper - config_.mesh.lower) / Real(config_.mesh.cells * (1 << location.level));
	result.time = time_;
	for (int d = 0; d < ndim; ++d)
		result.lower[d] = config_.mesh.lower + Real(location.coordinates[d] * config_.mesh.cells) * result.cellWidth;
	result.hydroEnabled = config_.hydroEnabled();
	result.radiationEnabled = config_.radiationEnabled();
	result.gravityEnabled = config_.gravityEnabled();
	if (build::hydro && config_.hydroEnabled()) result.hydro = hydro::Fields(result.layout, result.cellWidth, result.lower);
	if (build::radiation && config_.radiationEnabled()) result.radiation = radiation::Fields(result.layout, result.cellWidth, result.lower);
	if (build::gravity && config_.gravityEnabled()) {
		result.gravity = gravity::Fields(result.layout, result.cellWidth, result.lower);
		if (!config_.hydroEnabled()) result.density = mesh::PatchData<units::Density>(result.layout, result.cellWidth, result.lower);
	}
	if (config_.massFractions.enabled)
		for (std::size_t s = 0; s < config_.massFractions.species.size(); ++s) result.species.emplace_back(result.layout, result.cellWidth, result.lower);
	result.layout.forEachInterior([&](auto const& coordinate, std::size_t i) {
		mesh::BlockLocation cell{cellLevel(location.level), {}};
		for (int d = 0; d < ndim; ++d)
			cell.coordinates[d] = location.coordinates[d] * config_.mesh.cells + coordinate[d];
		auto const value = reconstruct(cell);
		for (std::size_t s = 0; s < result.species.size(); ++s) result.species[s].values()[i] = value.species.at(s);
		if (build::hydro && config_.hydroEnabled()) result.hydro.values()[i] = value.hydro;
		if (build::radiation && config_.radiationEnabled()) result.radiation.values()[i] = value.radiation;
		if (build::gravity && config_.gravityEnabled()) {
			result.gravity.values()[i] = value.gravity;
			if (!config_.hydroEnabled()) result.density.values()[i] = value.density;
		}
	});
	return result;
}


}	 // namespace octotigerII::amr
