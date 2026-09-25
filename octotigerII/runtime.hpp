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

	/// Serialize this value with its compile-time quantity types preserved.
	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & localTasks & stolenTasks;
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

	/// Reduce the per-block Courant limits after every locality finishes its tasks.
	units::Time stableTimestep() const;

	/// Publish one transport update only after all blocks and remote writebacks succeed.
	/// On failure, drain work and preserve the previously published state and time.
	void advance(units::Time dt);

	/// Cumulative boundary transport from successfully published timesteps.
	BoundaryTransport boundaryTransport() const;

	/// Regrid at synchronization points. Predicted signal travel can exhaust
	/// the buffer before amr.regridEvery. Returns whether leaves changed.
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
	/// Replace those kicks' self-gravity work with conservative mass-flux work.
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
