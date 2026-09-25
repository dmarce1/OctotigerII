// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#include "octotigerII/amr/hierarchy.hpp"
#include <bit>
#include "octotigerII/composition/transport.hpp"
#include <set>
#include "octotigerII/amr/interpolation.hpp"
#include "octotigerII/problems.hpp"
#include "octotigerII/verification/analytic.hpp"

namespace octotigerII::amr {

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
				if (build::hydro && config_.hydroEnabled()) speed[d] = std::max(speed[d], hydro::HydroSystem(config_.hydro).maximumSignalSpeed(value.hydro, d));
				if (build::radiation && config_.radiationEnabled()) speed[d] = std::max(speed[d], radiation::RadiationSystem(config_.radiation.lightSpeedRatio * constants::c).maximumSignalSpeed(value.radiation, d));
			}
		}
		units::InverseTime rate{};
		for (auto v : speed) rate += v / width;
		if (!(rate > units::InverseTime{}) || !units::finite(rate)) throw std::runtime_error("Invalid shadow CFL rate");
		auto const limit = config_.timestep.cfl / rate;
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
		for (auto const& [cell, value] : cells_) {
			width = std::min(width, (config_.mesh.upper - config_.mesh.lower) / Real(std::uint64_t(1) << cell.level));
			units::Acceleration acceleration{};
			for (int d = 0; d < ndim; ++d) {
				if (build::hydro && config_.hydroEnabled()) {
					speed[d] = std::max(speed[d], hydro::HydroSystem(config_.hydro).maximumSignalSpeed(value.hydro, d));
					acceleration += units::abs(config_.hydro.acceleration[d]
						+ (config_.gravityEnabled() ? value.gravity.acceleration(d) : units::Acceleration{}));
				}
				if (build::radiation && config_.radiationEnabled())
					speed[d] = std::max(speed[d], radiation::RadiationSystem(config_.radiation.lightSpeedRatio * constants::c).maximumSignalSpeed(value.radiation, d));
			}
			maximumAcceleration = std::max(maximumAcceleration, acceleration);
		}
		units::InverseTime rate{};
		for (auto v : speed) rate += v / width;
		if (!(rate > units::InverseTime{}) || !units::finite(rate)) throw std::runtime_error("Invalid shadow CFL rate");
		auto limit = config_.timestep.cfl / rate;
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
	for (auto const& block : shadowBlocks_) {
		int const level = cellLevel(block.level) - 1;
		auto const width = (config_.mesh.upper - config_.mesh.lower) / Real(1 << level);
		mesh::PhysicalCoordinates lower{};
		for (int d = 0; d < ndim; ++d)
			lower[d] = config_.mesh.lower + Real(block.coordinates[d] * n) * width;
		mesh::MeshLayout layout(n, 2);
		mesh::PatchData<hydro::ConservedState> gas(layout, width, lower);
		mesh::PatchData<radiation::RadiationSystem::State> radiation(layout, width, lower);
		mesh::forEachCoordinate(layout.extents(), [&](auto const& cell) {
			mesh::BlockLocation global{level, {}};
			for (int d = 0; d < ndim; ++d)
				global.coordinates[d] = block.coordinates[d] * n + cell[d] - 2;
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
		if constexpr (build::hydro) if (config_.hydroEnabled()) {
			hydro::Solver::Workspace work;
			hydro::Solver(hydro::HydroSystem(config_.hydro)).advanceInto(gas, dt, work, [&](auto const& cell, auto const& value) {
				auto& target = destination(cell);
				target.hydro = value;
				hydro::HydroSystem(config_.hydro).synchronize(target.hydro);
				target.density = value.density();
			});
			if (config_.massFractions.enabled) {
				auto read = [&](std::size_t s, auto const& cell) {
					mesh::BlockLocation global{level, {}};
					for (int d = 0; d < ndim; ++d) global.coordinates[d] = block.coordinates[d] * n + cell[d];
					return average(global).species.at(s);
				};
				auto flux = composition::fluxes(config_.massFractions, layout, read,
					[&](int d, auto const& face) { return work.fluxes[d][layout.faceIndex(d, face)].template get<0>(); });
				layout.forEachInterior([&](auto const& cell, std::size_t) {
					auto& target = destination(cell);
					for (std::size_t s = 0; s < flux.size(); ++s)
						target.species[s] = composition::update(read(s, cell), flux[s], layout, cell, dt / width);
					target.density = target.hydro.density() = composition::totalDensity(config_.massFractions, target.species);
				});
			}
		}
		if constexpr (build::radiation) if (config_.radiationEnabled()) {
			radiation::Solver::Workspace work;
			radiation::Solver(radiation::RadiationSystem(config_.radiation.lightSpeedRatio * constants::c))
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
			for (int d = 0; d < ndim; ++d) {
				auto const acceleration = config_.hydro.acceleration[d] + (config_.gravityEnabled() ? value.gravity.acceleration(d) : units::Acceleration{});
				auto const impulse = dt * value.hydro.density() * acceleration;
				value.hydro.totalEnergy() += impulse * (value.hydro.momentum(d) + 0.5 * impulse) / value.hydro.density();
				value.hydro.momentum(d) += impulse;
			}
			hydro::HydroSystem(config_.hydro).synchronize(value.hydro);
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
	auto const mapped = config_.mesh.boundary.map(cell.coordinates, 1 << cell.level);
	if (mapped.analytic) {
		mesh::PhysicalCoordinates position{};
		auto const h = (config_.mesh.upper - config_.mesh.lower) / Real(1 << cell.level);
		for (int d = 0; d < ndim; ++d)
			position[d] = config_.mesh.lower + (cell.coordinates[d] + 0.5) * h;
		auto const sample = problemBoundary(config_)(position, time_);
		Values result;
		if (build::hydro && config_.hydroEnabled()) result.hydro = hydro::HydroSystem(config_.hydro).conservedState(sample.hydro);
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
			result.hydro, mapped.reflectionMask, mapped.outflowLowerMask, mapped.outflowUpperMask, hydro::HydroSystem(config_.hydro));
	if (build::radiation && config_.radiationEnabled())
		result.radiation = physics::transformBoundary(result.radiation, mapped.reflectionMask, mapped.outflowLowerMask, mapped.outflowUpperMask,
			radiation::RadiationSystem(config_.radiation.lightSpeedRatio * constants::c));
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
	if (build::hydro && config_.hydroEnabled()) result.hydro = interpolate(center.hydro, gasSlopes, offset, hydro::HydroSystem(config_.hydro));
	if (build::radiation && config_.radiationEnabled())
		result.radiation = interpolate(center.radiation, radiationSlopes, offset, radiation::RadiationSystem(config_.radiation.lightSpeedRatio * constants::c));
	if (build::hydro && config_.hydroEnabled()) result.density = result.hydro.density();
	if (config_.massFractions.enabled) {
		// Inject concentrations and prolong hydro density. Summed child partial
		// densities preserve each parent species mass and the hydro mass together.
		for (std::size_t s = 0; s < result.species.size(); ++s)
			result.species[s] = center.species[s] * Real(result.density / center.density);
		result.density = result.hydro.density() = composition::totalDensity(config_.massFractions, result.species);
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
		if (build::hydro && config_.hydroEnabled()) result.signalSpeed[d] = hydro::HydroSystem(config_.hydro).maximumSignalSpeed(value.hydro, d);
		if (build::hydro && config_.hydroEnabled())
			result.acceleration[d] = config_.hydro.acceleration[d] + (config_.gravityEnabled() ? value.gravity.acceleration(d) : units::Acceleration{});
		if (build::radiation && config_.radiationEnabled()) result.signalSpeed[d] = std::max(result.signalSpeed[d], config_.radiation.lightSpeedRatio * constants::c);
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

namespace {
	using Location = mesh::BlockLocation;
	class LocationLess {
	public:
		bool operator()(Location const& a, Location const& b) const {
			auto const x = mesh::mortonKey(a), y = mesh::mortonKey(b);
			return x != y ? x < y : a.level < b.level;
		}
	};
	class Box {
	public:
		std::array<Real, ndim> lower{}, upper{};
		int level = 0;
	};
	Box bounds(Location const& location) {
		Box box;
		box.level = location.level;
		for (int d = 0; d < ndim; ++d) {
			box.lower[d] = Real(location.coordinates[d]) / Real(1 << location.level);
			box.upper[d] = Real(location.coordinates[d] + 1) / Real(1 << location.level);
		}
		return box;
	}
	bool overlaps(Box const& a, Box const& b, physics::BoundaryConditions const& bc, bool touching = false) {
		for (int d = 0; d < ndim; ++d) {
			bool hit = false;
			for (int shift = bc.periodic(d) ? -1 : 0; shift <= (bc.periodic(d) ? 1 : 0); ++shift) {
				if (touching ? a.lower[d] <= b.upper[d] + shift && b.lower[d] + shift <= a.upper[d] :
							   a.lower[d] < b.upper[d] + shift && b.lower[d] + shift < a.upper[d])
					hit = true;
			}
			if (!hit) return false;
		}
		return true;
	}
}	 // namespace

RegridResult selectMesh(Config const& c, std::vector<Snapshot> const& snapshots, Hierarchy const& hierarchy, units::Time horizon,
	refinement::Criteria const& criteria, bool allowCoarsening, Hierarchy const* restricted, bool oneLevel) {
	std::vector<Box> seeds;
	std::array<units::Velocity, ndim> speed{};
	std::unordered_map<Location, Real, mesh::BlockLocationHash> scores;
	std::set<Location, LocationLess> leaves;
	auto const length = c.mesh.upper - c.mesh.lower;
	for (auto const& block : snapshots) {
		leaves.insert(block.location);
		Real maximum = 0;
		Box seed;
		seed.lower.fill(1);
		seed.upper.fill(0);
		seed.level = std::min(c.amr.maxLevel, block.location.level + 1);
		bool flagged = false;
		block.layout.forEachInterior([&](auto const& cell, std::size_t) {
			auto const view = hierarchy.sample(block.location, cell);
			auto const value = refinement::score(view, criteria);
			maximum = std::max(maximum, value);
			for (int d = 0; d < ndim; ++d)
				speed[d] = std::max(speed[d], view.signalSpeed[d] + units::abs(view.acceleration[d]) * horizon);
			if (value <= 1) return;
			flagged = true;
			for (int d = 0; d < ndim; ++d) {
				Real const center = Real((view.center[d] - c.mesh.lower) / length);
				Real const radius = (Real(c.amr.bufferCells) + 0.5) * Real(view.width / length);
				seed.lower[d] = std::min(seed.lower[d], center - radius);
				seed.upper[d] = std::max(seed.upper[d], center + radius);
			}
		});
		scores[block.location] = maximum;
		if (flagged) seeds.push_back(seed);
	}
	for (auto& seed : seeds)
		for (int d = 0; d < ndim; ++d) {
			Real const travel = Real(c.amr.signalBuffer * speed[d] * horizon / length);
			seed.lower[d] -= travel;
			seed.upper[d] += travel;
		}
	auto desired = [&](Location const& leaf) {
		int level = c.amr.minLevel < 0 ? c.mesh.level : c.amr.minLevel;
		auto const box = bounds(leaf);
		for (auto const& seed : seeds)
			if (overlaps(box, seed, c.mesh.boundary)) level = std::max(level, seed.level);
		return level;
	};
	RegridResult result;
	result.signalSpeed = speed;
	int const minimum = c.amr.minLevel < 0 ? c.mesh.level : c.amr.minLevel;
	if (allowCoarsening) {
		std::set<Location, LocationLess> parents;
		for (auto const& leaf : leaves)
			if (leaf.level > minimum) parents.insert(leaf.parent());
		for (auto const& parent : parents) {
			bool safe = desired(parent) <= parent.level;
			for (int slot = 0; slot < (1 << ndim); ++slot)
				safe = safe && leaves.contains(parent.child(slot)) && scores.at(parent.child(slot)) < c.amr.coarsenFactor;
			if (!safe) continue;
			// Test the proposed coarse cells too: summing children can violate
			// the mass limit even when every child is below the threshold.
			mesh::forEachCoordinate(mesh::filledCoordinates(c.mesh.cells), [&](auto const& cell) {
				if (refinement::score((restricted ? *restricted : hierarchy).sample(parent, cell), criteria) >= c.amr.coarsenFactor) safe = false;
			});
			if (!safe) continue;
			for (int slot = 0; slot < (1 << ndim); ++slot)
				leaves.erase(parent.child(slot));
			leaves.insert(parent);
			++result.coarsened;
		}
	}
	bool changed = true;
	while (changed) {
		changed = false;
		std::vector<Location> split;
		for (auto const& leaf : leaves)
			if (leaf.level < desired(leaf)) split.push_back(leaf);
		for (auto const& leaf : split) {
			leaves.erase(leaf);
			for (int slot = 0; slot < (1 << ndim); ++slot)
				leaves.insert(leaf.child(slot));
			++result.refined;
			changed = true;
		}
		if (oneLevel) break;
	}
	// Full 2:1 balance includes edges, corners, and periodic neighbors, because
	// multidimensional reconstruction uses all of those ghost dependencies.
	changed = true;
	while (changed) {
		changed = false;
		std::set<Location, LocationLess> split;
		for (auto const& fine : leaves)
			mesh::forEachCoordinate(mesh::filledCoordinates(3), [&](auto const& direction) {
				Location neighbor = fine;
				for (int d = 0; d < ndim; ++d) {
					neighbor.coordinates[d] += direction[d] - 1;
					if ((neighbor.coordinates[d] < 0 || neighbor.coordinates[d] >= (1 << fine.level)) && !c.mesh.boundary.periodic(d)) return;
				}
				neighbor.coordinates = c.mesh.boundary.map(neighbor.coordinates, 1 << fine.level).source;
				while (!leaves.contains(neighbor) && !neighbor.isRoot())
					neighbor = neighbor.parent();
				if (leaves.contains(neighbor) && neighbor.level + 1 < fine.level) split.insert(neighbor);
			});
		for (auto const& leaf : split) {
			leaves.erase(leaf);
			for (int slot = 0; slot < (1 << ndim); ++slot)
				leaves.insert(leaf.child(slot));
			++result.refined;
			changed = true;
		}
	}
	result.leaves.assign(leaves.begin(), leaves.end());
	result.changed = result.leaves.size() != snapshots.size();
	if (!result.changed)
		for (auto const& block : snapshots)
			if (!leaves.contains(block.location)) result.changed = true;
	return result;
}

}	 // namespace octotigerII::amr
