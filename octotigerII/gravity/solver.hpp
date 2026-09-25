/** @file
 * @brief Newtonian gravity with open, periodic and reflecting image boundaries.
 * @ingroup numerics
 */
// Distributed under the Boost Software License, Version 1.0.
#pragma once
#include <cstdint>
#include <vector>
#include "octotigerII/gravity/gravityFields.hpp"
#include "octotigerII/physics/boundary.hpp"
#include "octotigerII/units/constants.hpp"

namespace octotigerII::gravity {

/// Number of accepted multipole pairs and direct cell pairs in the gravity solve.
/// @ingroup numerics
class Statistics {
public:
	std::uint64_t multipolePairs = 0;
	std::uint64_t directPairs = 0;
	/// Bounded worker tasks dispatched by the field solver, including exchanges.
	std::uint64_t workerTasks = 0;
	std::uint64_t ewaldPairs = 0, reflectedPairs = 0;
	/// Leaves evaluated on each locality; empty for the standalone serial reference.
	std::vector<std::uint64_t> localityCells;

	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & multipolePairs & directPairs & workerTasks & localityCells & ewaldPairs & reflectedPairs;
	}
};

/// Gravity fields in the same x-contiguous order as the input density, plus work counts.
/// @ingroup numerics
class Solution {
public:
	std::vector<State> fields;
	Statistics statistics;
};

// Newtonian gravity on a uniform cubic cell octree. Each leaf source
// is a cell-centered mass, with its self term omitted. x is contiguous.
// Newtonian images use diagonal plane waves; Ewald uses dense harmonic M2Ls.
// Scalar potentials retain order-p source/local expansions. Accelerations use
// a separate order-(p+1) local so the gradient has degree p in both endpoints;
// source-target reversal then produces equal and opposite integrated forces.
/// Solve image-boundary gravity for cell-centered masses, excluding only the physical self.
/// Input and output are x-contiguous; cellsPerAxis must be a power of two.
/// Far interactions use @ref ref_greengard1997 "Greengard and Rokhlin (1997)".
Solution solve(std::vector<units::Density> const& density, int cellsPerAxis, units::Length cellWidth, int order = 5, Real openingAngle = 0.5,
	physics::BoundaryConditions const& boundaries = {});
}	 // namespace octotigerII::gravity
