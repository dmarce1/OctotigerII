// Distributed under the Boost Software License, Version 1.0.
#pragma once
#include "octotigerII/gravity/gravityFields.hpp"
#include <cstdint>
#include <vector>

namespace octotigerII::gravity {
inline constexpr Real gravitationalConstant = 6.67430e-8; // cm^3 g^-1 s^-2

struct Statistics {
	std::uint64_t multipolePairs = 0;
	std::uint64_t directPairs = 0;
};

struct Solution {
	std::vector<State> fields;
	Statistics statistics;
};

// Isolated Newtonian gravity on a uniform cubic cell octree. Each leaf source
// is a cell-centered mass, with its self term omitted. x is contiguous.
// The only far-field operator is the diagonal plane-wave M2L.
Solution solve(std::vector<Real> const& density, int cellsPerAxis, Real cellWidth, int order = 5,
			   Real openingAngle = 0.5);
} // namespace octotigerII::gravity
