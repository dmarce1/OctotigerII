/** @file
 * @brief Stage publication and distributed locality work queues.
 * @ingroup runtime
 */
#pragma once
#include <memory>
#include "octotigerII/conservation.hpp"
#include "octotigerII/gravity/solver.hpp"
#include "octotigerII/refinement/criteria.hpp"
#include "octotigerII/subgrid/subgrid.hpp"

namespace octotigerII {

/// Cumulative completed tasks on their owning locality and on a stealing locality.
/// Timestep calculations also count as tasks; these are not network-message counts.
/// @ingroup runtime
class SchedulingStatistics {
public:
	std::uint64_t localTasks = 0;
	std::uint64_t stolenTasks = 0;
	std::vector<std::uint64_t> levelSteps; // Accepted transport steps, indexed by spatial level.

	/// Serialize this value with its compile-time quantity types preserved.
	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & localTasks & stolenTasks & levelSteps;
	}
};

// Owns fields, mesh metadata, and one scheduler per locality. No subgrid
// component owns numerical state. API calls are serialized stage operations.
/// Coordinates field storage and one executor per locality.
/// Public stage operations are serialized by an API mutex. Each update reads an
/// immutable bank and writes the other; publication follows completion of all work
/// and writebacks. Exceptions leave the published bank and time unchanged.
/// See @ref ref_kaiser2014 "Kaiser et al. (2014)" for the HPX execution model.
/// @ingroup runtime
class Runtime {
public:
	explicit Runtime(Config const& config, refinement::Criteria additionalCriteria = {});

	~Runtime();

	Runtime(Runtime const&) = delete;

	Runtime& operator=(Runtime const&) = delete;

	/// Gather interior-only exports of the published bank under the stage lock.
	std::vector<Snapshot> snapshots() const;

	/// Return the next synchronization interval. With eligible mixed-level AMR
	/// this is the coarsest active level's CFL limit; finer levels subcycle.
	/// Runs with uniform external acceleration use the global minimum.
	units::Time stableTimestep() const;

	/// Advance transport only; use advanceGravity for time-refined self-gravity.
	/// Optional AMR time refinement
	/// uses one dyadic step per spatial level and integrated flux correction.
	/// Publish only after all blocks and remote writebacks succeed.
	/// On failure, drain work and preserve the previously published state and time.
	void advance(units::Time dt);

	/// Complete a self-gravitating gas interval, including sources and energy work.
	/// AMR levels subcycle when timestep.refinement is enabled. The hierarchical
	/// mode reconciles provisional physical-state predictors to HOLD shell impulses.
	/// Every refined source ledger closes before publication; failure in the refined
	/// path restores the input. Global steps use the legacy stage sequence.
	gravity::Statistics advanceGravity(units::Time dt);

	/// Cumulative boundary transport from successfully published timesteps.
	BoundaryTransport boundaryTransport() const;

	/// Regrid at synchronization points. Predicted signal travel can exhaust
	/// the buffer before amr.regridEvery. Returns whether leaves changed.
	/// With conserveRegridEnergy, solveGravity() or setGravity() must follow a
	/// changed mesh before a timestep or another regrid. The next field publication
	/// recovers gas energy from conservatively remapped E+rho*phi/2.
	bool regrid(units::Time nextStep, bool force = false);
	std::size_t shadowCellCount() const;

	/// Solve gravity on the partitioned hierarchy and publish its distributed fields.
	/// The input bank and physical time survive a failed solve unchanged.
	gravity::Statistics solveGravity();

	/// Validate a complete gravity field directory, fill the next bank, then publish it.
	void setGravity(std::vector<std::vector<gravity::State>> const& fields);

	/// Apply self gravity plus uniform external acceleration and its kinetic-energy change.
	/// Internal energy is preserved; physical time is unchanged by this substep.
	void kickGravity(units::Time dt);

	/// Retain endpoint gravity and start accounting for the two kinetic-work kicks.
	void beginGravityEnergy();
	/// In mullen mode replace kick work with conservative mass-flux work;
	/// in naive mode retain kick work. Both modes account for boundary potential flux.
	/// Call after transport, the new gravity solve, and the second momentum kick.
	void finishGravityEnergy(units::Time dt);

	/// Return the number of represented elements or blocks.
	std::size_t size() const;

	/// Return the number of completed publications, including gravity assignments and kicks.
	std::uint64_t generation() const;

	/// Return cumulative completed local and stolen tasks from successful phases.
	SchedulingStatistics statistics() const;

	/// Return a human-readable description of the compiled execution backend.
	static char const* backend();

private:
	class Impl;
	std::unique_ptr<Impl> impl_;
};

}	 // namespace octotigerII
