#pragma once
#include "octotigerII/simulation.hpp"
#include <fstream>

namespace octotigerII {
class Output {
  public:
	explicit Output(Config const& config);
	void operator()(std::vector<Snapshot> const& snapshots, int step,
					Diagnostics const& diagnostics);

  private:
	Config config_;
	std::ofstream diagnostics_, series_;
	int frame_ = 0;
};
} // namespace octotigerII
