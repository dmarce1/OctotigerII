/** @file
 * @brief Time-tagged face-flux values for numerical-kernel consumers.
 * @ingroup mesh
 */
// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#pragma once
#include <vector>
#include "octotigerII/mesh.hpp"


namespace octotigerII {


// Time-tagged face fluxes for consumers of the numerical kernel.
/// Flux values with geometry and a physical time interval.
/// These are fluxes, not time-integrated corrections. A future AMR driver would
/// need to integrate and reconcile them; this type alone does not perform refluxing.
/// @ingroup mesh
template <typename State>
class FieldFluxPacket {
public:

	mesh::BlockLocation location;
	mesh::MeshLayout layout;
	mesh::PhysicalCoordinates lower{};
	units::Length cellWidth = units::Length::from_value(1);
	mesh::TimeInterval interval;
	std::vector<std::vector<State>> fluxes;

	/// Serialize this value with its compile-time quantity types preserved.
	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & location & layout & cellWidth & interval & fluxes;
		for (auto& coordinate : lower)
			archive & coordinate;
	}
};


}	 // namespace octotigerII
