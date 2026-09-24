/** @file
 * @brief Partitioned FMM hierarchy operating directly on distributed fields.
 * @ingroup runtime
 */
#pragma once

#include <memory>
#include "octotigerII/gravity/solver.hpp"
#include "octotigerII/storage/registry.hpp"

namespace octotigerII::gravity {


/// Persistent hierarchy with one partition per locality and bounded HPX workers.
/// A solve reads the published bank and fills the other bank. The caller must
/// hold the runtime stage lock and publish only after solve() succeeds.
class FieldSolver {
public:
	FieldSolver(Config const& config, std::vector<Subgrid> const& blocks, FieldDirectory const& fields, std::vector<storage::Locality> const& localities);

	~FieldSolver();

	FieldSolver(FieldSolver const&) = delete;

	FieldSolver& operator=(FieldSolver const&) = delete;

	Statistics solve(unsigned bank);

private:
	class Impl;
	std::unique_ptr<Impl> impl_;
};


}	 // namespace octotigerII::gravity
