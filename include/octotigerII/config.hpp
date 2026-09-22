#pragma once
#include "octotigerII/math/Real.hpp"
#include <string>
#include <vector>

namespace octotigerII {
inline constexpr Real physicalLightSpeed = 2.99792458e10; // cm/s
struct Config {
	std::string problem = "sod";
	int dimensions = 1, cells = 16, level = 1;
	Real lower = 0, upper = 1, stopTime = 0.2, cfl = 0.4;
	Real gamma = 1.4, lightSpeedRatio = 1;
	int multipoleOrder = 5;
	Real openingAngle = 0.5;
	int maxSteps = 100000, outputEvery = 10;
	bool periodic = false, outputEnabled = true;
	std::string outputDirectory = "output", outputFormat = "csv";
	void validate() const;
	bool hydroEnabled() const;
	bool radiationEnabled() const;
	bool gravityEnabled() const;
	template <class Archive> void serialize(Archive& a, unsigned) {
		a & problem & dimensions & cells & level & lower & upper & stopTime & cfl;
		a & gamma & lightSpeedRatio & multipoleOrder & openingAngle & maxSteps;
		a & outputEvery & periodic & outputEnabled & outputDirectory & outputFormat;
	}
};
Config parseConfig(std::vector<std::string> const& arguments);
std::string helpText();
} // namespace octotigerII
