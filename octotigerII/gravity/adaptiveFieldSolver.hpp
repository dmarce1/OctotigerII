// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#pragma once
#include "octotigerII/gravity/fieldSolver.hpp"

namespace octotigerII::gravity {
/// Sparse cell octree independent of hydro block ancestry and field ownership.
/// Geometric centers and dyadic scales keep image translations reusable.
class AdaptiveFieldSolver {
public:
	AdaptiveFieldSolver(Config const&, std::vector<Subgrid> const&, FieldDirectory const&, std::vector<storage::Locality> const&);
	~AdaptiveFieldSolver();
	Statistics solve(unsigned bank);
	Statistics solve(FieldSolveRequest const& request);

private:
	class Impl;
	std::unique_ptr<Impl> impl_;
};
}	 // namespace octotigerII::gravity
