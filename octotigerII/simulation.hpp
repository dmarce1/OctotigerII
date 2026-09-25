/** @file
 * @brief Global diagnostics and symmetric gravity/transport integration.
 * @ingroup runtime
 */
#pragma once
#include <functional>
#include "octotigerII/config.hpp"
#include "octotigerII/gravity/solver.hpp"
#include "octotigerII/runtime.hpp"


namespace octotigerII {


/// CGS volume integrals and extrema of synchronized interior fields.
/// maximumReducedFlux is the dimensionless magnitude |F|/(cE).
/// @ingroup runtime
class Diagnostics : public ConservedTotals {
public:

	units::Time time{};
	units::Energy kineticEnergy{}, thermalEnergy{}, gasGravityEnergy{};
	BoundaryTransport boundary;
	ConservedTotals norm; ///< Sum of absolute cell contributions, component by component.
	units::Energy gasGravityNorm{};
	units::Density minimumDensity{};
	units::Pressure minimumPressure{};
	units::EnergyDensity minimumRadiationEnergy{};
	Real maximumReducedFlux = 0;
};


/// Integrate synchronized interior fields and check positivity and radiation realizability.
Diagnostics diagnose(std::vector<Snapshot> const& snapshots, Config const& config);


/// Initial/final diagnostics, completed-step count, gravity work, and final snapshots.
/// @ingroup runtime
class RunResult {
public:

	Diagnostics initial, final;
	int steps = 0;
	gravity::Statistics gravityWork;
	std::vector<Snapshot> snapshots;
};


using Observer = std::function<void(std::vector<Snapshot> const&, int, Diagnostics const&)>;

/// Advance the configured problem and call the observer with synchronized snapshots.
/// Gravity uses symmetric source/transport/source composition; see
/// @ref ref_strang1968 "Strang (1968)".
RunResult run(Config const& config, Observer const& observer = {});
}	 // namespace octotigerII
