// Distributed under the Boost Software License, Version 1.0.
#pragma once
#include "octotigerII/config.hpp"
#include "octotigerII/gravity/gravityFields.hpp"
#include "octotigerII/hydro/hydroSystem.hpp"
#include "octotigerII/radiation/radiationTransport.hpp"
#include "octotigerII/subgrid/exchange.hpp"
#include <vector>

namespace octotigerII {
using HydroSnapshot = FieldSnapshot<hydro::ConservedState>;
using RadiationSnapshot = FieldSnapshot<radiation::RadiationSystem::State>;

struct Snapshot {
	mesh::BlockLocation location;
	mesh::MeshLayout layout;
	mesh::PhysicalCoordinates lower{};
	Real cellWidth = 1, time = 0;
	bool hydroEnabled = false, radiationEnabled = false, gravityEnabled = false;
	hydro::Fields hydro;
	radiation::Fields radiation;
	gravity::Fields gravity;
	mesh::PatchData<Real> density; // Prescribed density for gravity-only problems.
	template <class Archive> void serialize(Archive& a, unsigned) {
		a & location & layout & lower & cellWidth & time;
		a & hydroEnabled & radiationEnabled & gravityEnabled & hydro & radiation & gravity &
			density;
	}
};

// The complete per-block application state and operations. It is independent
// of HPX; the runtime adds only component lifetime, placement and action calls.
class Subgrid {
  public:
	Subgrid() = default;
	Subgrid(Config const& config, mesh::BlockLocation location);
	Snapshot snapshot() const;
	Real stableTimestep() const;
	void advance(std::vector<HydroSnapshot> const& hydro,
				 std::vector<RadiationSnapshot> const& radiation, Real dt);
	void setGravity(std::vector<gravity::State> const& fields);
	void kickGravity(Real dt);

  private:
	Config config_;
	Snapshot data_;
};
} // namespace octotigerII
