// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#pragma once
#include "octotigerII/refinement/criteria.hpp"
#include "octotigerII/subgrid/subgrid.hpp"

namespace octotigerII::amr {

class Hierarchy;

class RegridResult {
public:
	std::vector<mesh::BlockLocation> leaves;
	std::size_t refined = 0, coarsened = 0;
	std::array<units::Velocity, ndim> signalSpeed{};
	bool changed = false;
};

RegridResult selectMesh(Config const&, std::vector<Snapshot> const&, Hierarchy const&, units::Time horizon, refinement::Criteria const&,
	bool allowCoarsening = true, Hierarchy const* restricted = nullptr, bool oneLevel = false);

}	 // namespace octotigerII::amr
