/** @file
 * @brief Supported gravity image boundary conditions.
 */
#pragma once
#include "octotigerII/physics/boundary.hpp"

namespace octotigerII::gravity {

/// Analytic transport faces do not specify an exterior gravitational source.
inline void validateBoundaries(physics::BoundaryConditions const& boundaries) {
	boundaries.validate();
	if (boundaries.contains(physics::BoundaryCondition::Analytic))
		throw std::invalid_argument("Gravity supports outflow, inflow, periodic, and reflecting faces; analytic gravity is not implemented");
}

}	 // namespace octotigerII::gravity
