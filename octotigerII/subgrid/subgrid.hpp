/** @file
 * @brief Temporary interior snapshots used for initialization, diagnostics, and output.
 * @ingroup mesh
 */
// Distributed under the Boost Software License, Version 1.0.
#pragma once
#include <vector>
#include "octotigerII/config.hpp"
#include "octotigerII/gravity/gravityFields.hpp"
#include "octotigerII/hydro/hydroSystem.hpp"
#include "octotigerII/radiation/radiationTransport.hpp"
#include "octotigerII/subgrid/topology.hpp"


namespace octotigerII {


/// Owning, interior-only export at one physical time.
/// Used for initialization, diagnostics, and Silo output; it is not persistent
/// subgrid storage or a restart checkpoint format.
/// @ingroup mesh
class Snapshot {
public:

	mesh::BlockLocation location;
	mesh::MeshLayout layout;
	mesh::PhysicalCoordinates lower{};
	units::Length cellWidth = units::Length::from_value(1);
	units::Time time{};
	bool hydroEnabled = false, radiationEnabled = false, gravityEnabled = false;
	hydro::Fields hydro;
	std::vector<mesh::PatchData<units::Density>> species;
	radiation::Fields radiation;
	gravity::Fields gravity;

	mesh::PatchData<units::Density> density;	// Prescribed density for gravity-only problems.
	/// Serialize this value with its compile-time quantity types preserved.
	template <typename Archive>
	void serialize(Archive& a, unsigned) {
		a & location & layout & cellWidth & time;
		for (auto& coordinate : lower)
			a & coordinate;
		a & hydroEnabled & radiationEnabled & gravityEnabled & hydro & radiation & gravity & density & species;
	}
};


// Build a temporary interior-only value for problem initialization and IO tests.
/// Construct a temporary interior-only CGS value for problem initialization.
Snapshot initialSnapshot(Config const& config, mesh::BlockLocation location, bool refinementProbe = false);

}	 // namespace octotigerII
