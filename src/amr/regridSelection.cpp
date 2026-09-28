// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#include "octotigerII/amr/regridSelection.hpp"
#include <set>
#include "octotigerII/amr/hierarchy.hpp"

namespace octotigerII::amr {
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
	bool overlaps(Box const& a, Box const& b, finiteVolume::BoundaryConditions const& bc, bool touching = false) {
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
