#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>
#include "octotigerII/mesh.hpp"
#include "octotigerII/output.hpp"
#include "octotigerII/profilingSnapshots.hpp"
#include "octotigerII/verification/analytic.hpp"
#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/hpx_init.hpp>
#endif


namespace {
int application(std::vector<std::string> const& args) {
	try {
		if (std::find(args.begin(), args.end(), "--help") != args.end()) {
			std::cout << octotigerII::helpText();
			return 0;
		}
		auto const config = octotigerII::parseConfig(args);
		std::cout << "OctotigerII 0.1.0 | " << octotigerII::Runtime::backend() << " | " << config.problem << " | " << octotigerII::ndim << "D\n";
		std::cout << std::scientific << std::setprecision(6);
		octotigerII::Output output(config);
		bool progressHeaderPrinted = false;
		int const lastSubgridLevel = config.amr.enabled ? config.amr.maxLevel : config.mesh.level;
		auto const result = octotigerII::run(config, [&](auto const& snapshots, int step, auto const& d) {
			if (step == 0) std::clog << "Initialization: writing initial diagnostics and output\n";
			output(snapshots, step, d);
			if (step == 0) std::clog << "Initialization: initial output complete\n";
			if (step == 0 || step % config.output.every == 0 || d.time >= config.runtime.stopTime) {
				if (!progressHeaderPrinted) {
					std::cout << std::left << std::setw(8) << "step" << std::setw(16) << "time[s]";
					for (int level = 0; level <= lastSubgridLevel; ++level)
						std::cout << std::setw(24) << ("level " + std::to_string(level));
					std::cout << '\n' << std::setw(8) << "" << std::setw(16) << "";
					for (int level = 0; level <= lastSubgridLevel; ++level) {
						std::cout << std::setw(24) << "(total/leafs)";
					}
					std::cout << '\n';
					progressHeaderPrinted = true;
				}
				// Snapshots list active leaves. Add each leaf's ancestors to count all
				// tree nodes once at their own level, including covered subgrids.
				std::vector<std::size_t> leaves(static_cast<std::size_t>(lastSubgridLevel + 1));
				std::vector<std::unordered_set<octotigerII::mesh::BlockLocation,
					octotigerII::mesh::BlockLocationHash>> allSubgrids(static_cast<std::size_t>(lastSubgridLevel + 1));
				for (auto const& block : snapshots) {
					++leaves.at(static_cast<std::size_t>(block.location.level));
					auto node = block.location;
					for (;;) {
						allSubgrids.at(static_cast<std::size_t>(node.level)).insert(node);
						if (node.level == 0) break;
						node = node.parent();
					}
				}
				std::cout << std::right << std::setw(8) << step
						  << std::setw(16) << octotigerII::units::value(d.time);
				for (std::size_t level = 0; level < allSubgrids.size(); ++level)
					std::cout << std::setw(24) << (std::to_string(allSubgrids[level].size()) + "/" + std::to_string(leaves[level]));
				std::cout << '\n' << std::flush;
			}
		});
		auto const comparison = octotigerII::verification::compare(result.snapshots, config);
		comparison.print(std::cout);
		comparison.enforce(config);
		if (config.gravityEnabled()) {
			std::cout << "gravity multipolePairs=" << result.gravityWork.multipolePairs << " directPairs=" << result.gravityWork.directPairs
				<< " workerTasks=" << result.gravityWork.workerTasks << " cellsPerLocality=";
			for (auto count : result.gravityWork.localityCells)
				std::cout << ' ' << count;
			std::cout << '\n';
		}
		std::cout << "Completed " << result.steps << " steps at t=" << octotigerII::units::value(result.final.time) << " s\n";
		return 0;
	} catch (std::exception const& error) {
		std::cerr << "OctotigerII: " << error.what() << '\n';
		return 1;
	}
}
#ifdef OCTOTIGERII_WITH_HPX
std::vector<std::string> applicationArguments;
#endif
}	 // namespace


#ifdef OCTOTIGERII_WITH_HPX
int hpx_main(int, char**) {
	int const status = application(applicationArguments);
	octotigerII::profiling::stopSnapshots();
	hpx::finalize();
	return status;
}
#endif

int main(int argc, char** argv) {
	std::vector<std::string> args(argv + 1, argv + argc);
#ifdef OCTOTIGERII_WITH_HPX
	// HPX consumes only its own arguments; the strict application parser sees
	// exactly the same argument vector in the HPX and serial builds.
	std::vector<char*> runtimeArguments{argv[0]};
	for (int i = 1; i < argc; ++i) {
		if (std::string(argv[i]).starts_with("--hpx:"))
			runtimeArguments.push_back(argv[i]);
		else
			applicationArguments.emplace_back(argv[i]);
	}
	int const runtimeCount = static_cast<int>(runtimeArguments.size());
	runtimeArguments.push_back(nullptr);
	hpx::init_params init;
	init.startup = octotigerII::profiling::startSnapshots;
	return hpx::init(runtimeCount, runtimeArguments.data(), init);
#else
	return application(args);
#endif
}
