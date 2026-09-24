/** @file
 * @brief Silo time-series output from synchronized physical-CGS snapshots.
 * @ingroup runtime
 */
#pragma once
#include <fstream>
#include "octotigerII/simulation.hpp"


namespace octotigerII {


/// Writes zone-centered CGS fields and a VisIt series file.
/// Snapshots must share a physical time. Existing frame files use DB_CLOBBER.
/// @ingroup runtime
class Output {
public:

	explicit Output(Config const& config);

	void operator()(std::vector<Snapshot> const& snapshots, int step, Diagnostics const& diagnostics);

private:

	Config config_;
	std::ofstream series_, conservation_;
	int frame_ = 0;
	Diagnostics initial_;
	bool haveInitial_ = false;
};
}	 // namespace octotigerII
