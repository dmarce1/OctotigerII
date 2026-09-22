#pragma once
#include "octotigerII/subgrid/subgrid.hpp"
#include <memory>

namespace octotigerII {
// Fixed block directory, placement and completion barriers. The numerical
// state stays in Subgrid. Every phase completes before its inputs change.
class Runtime {
  public:
	explicit Runtime(Config const& config);
	~Runtime();
	Runtime(Runtime const&) = delete;
	Runtime& operator=(Runtime const&) = delete;
	std::vector<Snapshot> snapshots() const;
	Real stableTimestep() const;
	void advance(std::vector<Snapshot> const& old, Real dt);
	void setGravity(std::vector<std::vector<gravity::State>> const& fields);
	void kickGravity(Real dt);
	std::size_t size() const;
	static char const* backend();

  private:
	struct Impl;
	std::unique_ptr<Impl> impl_;
};
} // namespace octotigerII
