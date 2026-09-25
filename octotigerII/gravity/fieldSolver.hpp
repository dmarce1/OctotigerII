/** @file
 * @brief Partitioned FMM hierarchy operating directly on distributed fields.
 * @ingroup runtime
 */
#pragma once

#include <memory>
#include "octotigerII/gravity/solver.hpp"
#include "octotigerII/storage/registry.hpp"

namespace octotigerII::gravity {

/// A linear combination of density banks for one hydro block. Setting both
/// weights to zero excludes that block from the source mass distribution.
class DensitySelection {
public:
	unsigned bank = 0, secondBank = 0;
	Real weight = 1, secondWeight = 0;

	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & bank & secondBank & weight & secondWeight;
	}
};

/// A field solve independent of hydro bank publication. All ranges use the
/// constructor's field layout, and block lists use its original block order.
class FieldSolveRequest {
public:
	/// Empty selects the hydro density (or gravity-only density) supplied at construction.
	storage::FieldHandle<units::Density> density;
	/// Empty reads sourceBank on every block; otherwise one selection per block.
	std::vector<DensitySelection> sources;
	unsigned sourceBank = 0;
	/// Empty selects the gravity field supplied at construction.
	storage::ColumnHandle<State> output;
	unsigned outputBank = 0;
	/// Empty publishes all blocks. Unselected output ranges are left untouched.
	std::vector<unsigned char> targets;
	/// Required for signed source fields such as density differences.
	bool allowSignedDensity = false;

	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & density & sources & sourceBank & output & outputBank & targets & allowSignedDensity;
	}
};

/// Persistent hierarchy with one partition per locality and bounded HPX workers.
/// The legacy solve(bank) reads the published bank and fills the other bank.
/// The caller must hold the runtime stage lock and publish only after success.
class FieldSolver {
public:
	FieldSolver(Config const& config, std::vector<Subgrid> const& blocks, FieldDirectory const& fields, std::vector<storage::Locality> const& localities);

	~FieldSolver();

	FieldSolver(FieldSolver const&) = delete;

	FieldSolver& operator=(FieldSolver const&) = delete;

	Statistics solve(unsigned bank);

	/// Evaluate a source selection, writing only the requested gravity output.
	/// The caller holds the stage lock; hydro, radiation and species are untouched.
	/// Pair statistics count directed interactions, since masks need not be symmetric.
	Statistics solve(FieldSolveRequest const& request);

private:
	class Impl;
	std::unique_ptr<Impl> impl_;
};


}	 // namespace octotigerII::gravity
