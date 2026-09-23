/** @file
 * @brief Isolated Newtonian gravity on a uniform cubic cell hierarchy.
 * @ingroup numerics
 */
// Distributed under the Boost Software License, Version 1.0.
#pragma once
#include <cstdint>
#include <vector>
#include "octotigerII/gravity/gravityFields.hpp"
#include "octotigerII/units/constants.hpp"


namespace octotigerII::gravity {


/// Number of accepted multipole pairs and direct cell pairs in the gravity solve.
/// @ingroup numerics
class Statistics {
public:

	std::uint64_t multipolePairs = 0;
	std::uint64_t directPairs = 0;
};


/// Gravity fields in the same x-contiguous order as the input density, plus work counts.
/// @ingroup numerics
class Solution {
public:

	std::vector<State> fields;
	Statistics statistics;
};


// Isolated Newtonian gravity on a uniform cubic cell octree. Each leaf source
// is a cell-centered mass, with its self term omitted. x is contiguous.
// The only far-field operator is the diagonal plane-wave M2L.
/// Solve isolated gravity for cell-centered point masses, excluding self terms.
/// Input and output are x-contiguous; cellsPerAxis must be a power of two.
/// Far interactions use @ref ref_greengard1997 "Greengard and Rokhlin (1997)".
Solution solve(std::vector<units::Density> const& density, int cellsPerAxis, units::Length cellWidth, int order = 5, Real openingAngle = 0.5);
}	 // namespace octotigerII::gravity
