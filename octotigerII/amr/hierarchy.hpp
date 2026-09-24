// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#pragma once
#include <unordered_map>
#include "octotigerII/refinement/criteria.hpp"
#include "octotigerII/subgrid/subgrid.hpp"

namespace octotigerII::amr {

/// Cell-centered shadows on every dyadic level, including levels within a block.
/// Initialize ancestors by volume restriction, then evolve covered coarse
/// states independently. Shadows never count as additional physical cells.
class Hierarchy {
public:
	class Values {
	public:
		hydro::ConservedState hydro;
		radiation::RadiationSystem::State radiation;
		gravity::State gravity;
		units::Density density{};
		Values& operator+=(Values const& other);
		Values operator*(Real weight) const;
	};

	Hierarchy(Config const&, std::vector<Snapshot> const& leaves);
	Values average(mesh::BlockLocation cell) const;
	Values reconstruct(mesh::BlockLocation cell) const;
	Values shadow(mesh::BlockLocation cell) const;
	refinement::CellView sample(mesh::BlockLocation block, mesh::Coordinates const& cell) const;
	Snapshot transfer(mesh::BlockLocation block) const;
	/// Advance covered coarse states independently, without averaging leaves down.
	void advance(units::Time dt);
	void kick(units::Time dt);
	void refreshLeaves(std::vector<Snapshot> const&);
	void refreshGravity(std::vector<Snapshot> const&);
	std::size_t size() const {
		return cells_.size();
	}
	int cellLevel(int blockLevel) const {
		return blockLevel + blockBits_;
	}

private:
	Config config_;
	units::Time time_{};
	int blockBits_ = 0;
	Values scale_;
	std::vector<mesh::BlockLocation> shadowBlocks_;
	std::unordered_map<mesh::BlockLocation, Values, mesh::BlockLocationHash> cells_;
	Values reconstructFrom(mesh::BlockLocation coarse, mesh::BlockLocation fine) const;
};

class RegridResult {
public:
	std::vector<mesh::BlockLocation> leaves;
	std::size_t refined = 0, coarsened = 0;
	std::array<units::Velocity, ndim> signalSpeed{};
	bool changed = false;
};

RegridResult selectMesh(Config const&, std::vector<Snapshot> const&, Hierarchy const&, units::Time horizon, refinement::Criteria const&,
	bool allowCoarsening = true, Hierarchy const* restricted = nullptr, bool oneLevel = false);

class InitialMesh {
public:
	std::vector<mesh::BlockLocation> leaves;
	std::array<units::Velocity, ndim> signalSpeed{};
	units::Time timestep{};
};

/// Start at the root and sample the initial condition again after each split.
/// The optional initializer supports independent tests of the startup algorithm.
using Initializer = std::function<Snapshot(Config const&, mesh::BlockLocation)>;
InitialMesh initializeMesh(Config const&, refinement::Criteria const&, Initializer const& = [](Config const& c, mesh::BlockLocation l) { return initialSnapshot(c, l, true); });

}	 // namespace octotigerII::amr
