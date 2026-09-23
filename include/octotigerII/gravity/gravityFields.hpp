/** @file
 * @brief CGS potential and acceleration field ordering.
 * @ingroup numerics
 */
// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#pragma once

#include "octotigerII/mesh.hpp"
#include "octotigerII/units/state.hpp"


namespace octotigerII::gravity {

using State = units::ScalarVectorState<units::VelocitySquared, units::Acceleration>;

using Fields = mesh::PatchData<State>;


}	 // namespace octotigerII::gravity
