#include "octotigerII/output.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
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
		std::cout << "OctotigerII 0.1.0 | " << octotigerII::Runtime::backend() << " | "
				  << config.problem << '\n';
		std::cout << std::scientific << std::setprecision(6);
		octotigerII::Output output(config);
		auto const result =
			octotigerII::run(config, [&](auto const& snapshots, int step, auto const& d) {
				output(snapshots, step, d);
				if (step == 0 || step % config.outputEvery == 0 || d.time >= config.stopTime)
					std::cout << "step=" << step << " time=" << d.time << " mass=" << d.mass
							  << " gasEnergy=" << d.gasEnergy
							  << " radiationEnergy=" << d.radiationEnergy << '\n';
			});
		if (config.gravityEnabled())
			std::cout << "gravity multipolePairs=" << result.gravityWork.multipolePairs
					  << " directPairs=" << result.gravityWork.directPairs << '\n';
		std::cout << "Completed " << result.steps << " steps at t=" << result.final.time << " s\n";
		return 0;
	} catch (std::exception const& error) {
		std::cerr << "OctotigerII: " << error.what() << '\n';
		return 1;
	}
}
#ifdef OCTOTIGERII_WITH_HPX
std::vector<std::string> applicationArguments;
#endif
} // namespace
#ifdef OCTOTIGERII_WITH_HPX
int hpx_main(int, char**) {
	int const status = application(applicationArguments);
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
	return hpx::init(runtimeCount, runtimeArguments.data());
#else
	return application(args);
#endif
}
