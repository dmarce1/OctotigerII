#include "octotigerII/config.hpp"
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace octotigerII {
namespace {
std::string trim(std::string s) {
	auto first = s.find_first_not_of(" \t\r\n");
	return first == std::string::npos ? ""
									  : s.substr(first, s.find_last_not_of(" \t\r\n") - first + 1);
}
using Settings = std::map<std::string, std::string>;
void setting(Settings& values, std::string line) {
	auto pos = line.find('=');
	if (pos == std::string::npos)
		throw std::invalid_argument("Expected key=value: " + line);
	std::string key = trim(line.substr(0, pos)), value = trim(line.substr(pos + 1));
	if (key.empty() || value.empty())
		throw std::invalid_argument("Empty setting: " + line);
	values[key] = value;
}
Real number(std::string const& value) {
	std::size_t consumed = 0;
	Real x = std::stod(value, &consumed);
	if (consumed != value.size() || !std::isfinite(x))
		throw std::invalid_argument("Invalid number: " + value);
	return x;
}
int integer(std::string const& value) {
	std::size_t consumed = 0;
	int x = std::stoi(value, &consumed);
	if (consumed != value.size())
		throw std::invalid_argument("Invalid integer: " + value);
	return x;
}
bool boolean(std::string const& value) {
	if (value == "on" || value == "true")
		return true;
	if (value == "off" || value == "false")
		return false;
	throw std::invalid_argument("Expected on/off: " + value);
}
} // namespace
bool Config::hydroEnabled() const {
	return problem == "sod" || problem == "kelvin-helmholtz" || problem == "collapse";
}
bool Config::radiationEnabled() const {
	return problem == "streaming" || problem == "radiation-pulse";
}
bool Config::gravityEnabled() const {
	return problem == "gravity-sphere" || problem == "gravity-gaussian" || problem == "collapse";
}
void Config::validate() const {
	if (!hydroEnabled() && !radiationEnabled() && !gravityEnabled())
		throw std::invalid_argument("Unknown problem: " + problem);
	if (dimensions < 1 || dimensions > 3 || cells < 4 || cells > 128 || (cells & (cells - 1)) ||
		level < 0 || level > 6)
		throw std::invalid_argument("mesh: ndim=1..3, cells=power of two in [4,128], level=0..6");
	if (!(upper > lower) || !std::isfinite(upper - lower) || !(stopTime >= 0) ||
		!(cfl > 0 && cfl <= 0.5) || !(gamma > 1) ||
		!(lightSpeedRatio > 0 && lightSpeedRatio <= 1) || maxSteps < 1 || outputEvery < 1)
		throw std::invalid_argument("Invalid domain, timestep, gas, radiation, or output setting");
	if (gravityEnabled() && (dimensions != 3 || periodic))
		throw std::invalid_argument("Gravity currently requires 3D isolated boundaries");
	if (multipoleOrder < 3 || multipoleOrder > 5 ||
		!(openingAngle > 0 && openingAngle < 1 / std::sqrt(3.0)))
		throw std::invalid_argument("Gravity requires order 3..5 and 0<opening_angle<1/sqrt(3)");
	if (problem == "sod" && (dimensions != 1 || periodic))
		throw std::invalid_argument("Sod requires 1D outflow boundaries");
	if (problem == "kelvin-helmholtz" && (dimensions != 2 || !periodic))
		throw std::invalid_argument("Kelvin-Helmholtz requires 2D periodic boundaries");
	if (outputFormat != "csv" && outputFormat != "silo")
		throw std::invalid_argument("output.format must be csv or silo");
#ifndef OCTOTIGERII_WITH_SILO
	if (outputEnabled && outputFormat == "silo")
		throw std::invalid_argument("Silo output requires OCTOTIGERII_WITH_SILO=ON");
#endif
	if (outputDirectory.empty())
		throw std::invalid_argument("Empty output directory");
	if (gravityEnabled() && !hydroEnabled() && stopTime != 0)
		throw std::invalid_argument("Static gravity examples require runtime.stop_time=0");
}

Config parseConfig(std::vector<std::string> const& arguments) {
	Settings values;
	// Config files first, CLI overrides second, independent of argument order.
	for (auto const& arg : arguments)
		if (arg.starts_with("--config=")) {
			std::ifstream in(arg.substr(9));
			if (!in)
				throw std::runtime_error("Cannot open config: " + arg.substr(9));
			std::string line;
			while (std::getline(in, line)) {
				line = trim(line.substr(0, line.find_first_of("#;")));
				if (!line.empty())
					setting(values, line);
			}
		}
	for (auto const& arg : arguments) {
		if (arg.starts_with("--config="))
			continue;
		if (!arg.starts_with("--"))
			throw std::invalid_argument("Expected --key=value: " + arg);
		setting(values, arg.substr(2));
	}
	Config c;
	if (auto it = values.find("problem.name"); it != values.end())
		c.problem = it->second;
	// Small, explicit defaults for each ordinary problem.
	if (c.problem == "kelvin-helmholtz") {
		c.dimensions = 2;
		c.periodic = true;
		c.stopTime = 0.1;
	}
	if (c.radiationEnabled()) {
		c.lower = -3e10;
		c.upper = 3e10;
		c.stopTime = 0.4;
		c.periodic = true;
	}
	if (c.gravityEnabled()) {
		c.dimensions = 3;
		c.cells = 4;
		c.lower = -1e9;
		c.upper = 1e9;
		c.stopTime = 0;
	}
	if (c.problem == "collapse")
		c.stopTime = 1;
	for (auto const& [key, v] : values) {
		if (key == "problem.name")
			c.problem = v;
		else if (key == "mesh.ndim")
			c.dimensions = integer(v);
		else if (key == "mesh.cells")
			c.cells = integer(v);
		else if (key == "mesh.level")
			c.level = integer(v);
		else if (key == "mesh.lower")
			c.lower = number(v);
		else if (key == "mesh.upper")
			c.upper = number(v);
		else if (key == "mesh.periodic")
			c.periodic = boolean(v);
		else if (key == "runtime.stop_time")
			c.stopTime = number(v);
		else if (key == "runtime.max_steps")
			c.maxSteps = integer(v);
		else if (key == "timestep.cfl")
			c.cfl = number(v);
		else if (key == "hydro.gamma")
			c.gamma = number(v);
		else if (key == "radiation.light_speed_ratio")
			c.lightSpeedRatio = number(v);
		else if (key == "gravity.multipole_order")
			c.multipoleOrder = integer(v);
		else if (key == "gravity.opening_angle")
			c.openingAngle = number(v);
		else if (key == "output.enabled")
			c.outputEnabled = boolean(v);
		else if (key == "output.every")
			c.outputEvery = integer(v);
		else if (key == "output.directory")
			c.outputDirectory = v;
		else if (key == "output.format")
			c.outputFormat = v;
		else
			throw std::invalid_argument("Unknown option: " + key);
	}
	c.validate();
	return c;
}
std::string helpText() {
	return "OctotigerII 0.1.0 (cgs)\n"
		   "Usage: octotigerII --config=examples/sod.ini [--key=value ...]\n"
		   "Problems: sod, kelvin-helmholtz, gravity-sphere, gravity-gaussian,\n"
		   "          streaming, radiation-pulse, collapse\n"
		   "Options: problem.name, mesh.ndim, mesh.cells, mesh.level, mesh.lower,\n"
		   "         mesh.upper, mesh.periodic, runtime.stop_time, runtime.max_steps,\n"
		   "         timestep.cfl, hydro.gamma, radiation.light_speed_ratio,\n"
		   "         gravity.multipole_order, gravity.opening_angle, output.enabled,\n"
		   "         output.every, output.directory, output.format (csv|silo)\n"
		   "Booleans: on/off. mesh.cells is cells per block per active axis.\n"
		   "Fixed hierarchy: 2^mesh.level blocks per active axis.\n";
}
} // namespace octotigerII
