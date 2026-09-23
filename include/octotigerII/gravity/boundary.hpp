/** @file
 * @brief Boundary capability of the current isolated gravity solver.
 */
#pragma once
#include "octotigerII/physics/boundary.hpp"


namespace octotigerII::gravity {


/// Keep gravity capability separate from transport geometry for future extensions.
inline void validateBoundaries(physics::BoundaryConditions const& boundaries) {
	boundaries.validate();
	if (!boundaries.all(physics::BoundaryCondition::Outflow))
		throw std::invalid_argument("Gravity currently supports only outflow boundaries on every face (isolated gravity)");
}


} // namespace octotigerII::gravity
