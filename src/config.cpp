#include "octotigerII/config.hpp"
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "octotigerII/problems.hpp"

#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/modules/program_options.hpp>
#else
#include <boost/program_options.hpp>
#endif

namespace octotigerII {

namespace {

#ifdef OCTOTIGERII_WITH_HPX
	namespace po = hpx::program_options;
#else
	namespace po = boost::program_options;
#endif

	void addSettings(po::options_description& description) {
		auto options = description.add_options();
		options("randomSeed", po::value<std::int64_t>(), "Global nonnegative random seed (default 5489)");
		options("verification.analytic", po::value<std::string>(), "Analytic comparison: auto (default), on (required), off");
		options("verification.gravityReference", po::value<std::string>(), "Gravity reference: direct (default) or continuum");
		options("verification.directMaxPairs", po::value<std::int64_t>(), "Direct-reference pair budget (default 20000000; at least one complete target)");
		options("verification.directSamples", po::value<std::int64_t>(),
			"Number of direct-reference targets; 0 selects automatically, positive overrides the pair budget");
		options("verification.relativeL1Tolerance", po::value<Real>(), "Maximum relative L1 for each reference field; -1 disables the accuracy gate");
		options("verification.absoluteTolerance", po::value<Real>(), "Maximum absolute Linf for zero-reference fields (CGS), when the gate is enabled");
		options("mesh.cells", po::value<int>(), "Cells per block per active axis");
		options("mesh.level", po::value<int>(), "Uniform block level");
		options("mesh.lower", po::value<Real>(), "Lower domain coordinate (cm)");
		options("mesh.upper", po::value<Real>(), "Upper domain coordinate (cm)");
		options("mesh.periodic", po::value<std::string>(), "Periodic boundaries: on/off");
		options("runtime.stopTime", po::value<Real>(), "Stop time (s)");
		options("runtime.maxSteps", po::value<int>(), "Maximum number of steps");
		options("runtime.workerTasks", po::value<int>(), "Maximum concurrent tasks; 0 uses worker count");
		options("runtime.workStealing", po::value<std::string>(), "Remote work stealing: on/off");
		options("timestep.cfl", po::value<Real>(), "Courant factor");
		options("hydro.gamma", po::value<Real>(), "Ideal-gas adiabatic index");
		options("radiation.lightSpeedRatio", po::value<Real>(), "Radiation transport speed divided by c");
		options("gravity.multipoleOrder", po::value<int>(), "Gravity expansion order (1..10)");
		options("gravity.openingAngle", po::value<Real>(), "Gravity opening angle");
		options("output.enabled", po::value<std::string>(), "Silo output: on/off");
		options("output.every", po::value<int>(), "Output every N steps");
		options("output.directory", po::value<std::string>(), "Silo output directory");
	}

	void addLegacySettings(po::options_description& description) {
		auto options = description.add_options();
		options("runtime.stop_time", po::value<Real>(), "Alias for runtime.stopTime");
		options("runtime.max_steps", po::value<int>(), "Alias for runtime.maxSteps");
		options("runtime.worker_tasks", po::value<int>(), "Alias for runtime.workerTasks");
		options("runtime.work_stealing", po::value<std::string>(), "Alias for runtime.workStealing");
		options("radiation.light_speed_ratio", po::value<Real>(), "Alias for radiation.lightSpeedRatio");
		options("gravity.multipole_order", po::value<int>(), "Alias for gravity.multipoleOrder");
		options("gravity.opening_angle", po::value<Real>(), "Alias for gravity.openingAngle");
	}

	char const* selectedKey(po::variables_map const& values, char const* key, char const* legacy = nullptr) {
		bool const canonicalPresent = values.count(key) != 0;
		bool const legacyPresent = legacy != nullptr && values.count(legacy) != 0;
		if (canonicalPresent && legacyPresent) throw std::invalid_argument(std::string("Conflicting names for ") + key + ": " + legacy);
		if (canonicalPresent) return key;
		if (legacyPresent) return legacy;
		return nullptr;
	}

	template <typename T>
	void readOption(po::variables_map const& values, char const* key, T& target, char const* legacy = nullptr) {
		if (char const* selected = selectedKey(values, key, legacy)) target = values[selected].as<T>();
	}

	void readNumber(po::variables_map const& values, char const* key, Real& target, char const* legacy = nullptr) {
		using std::isfinite;

		if (char const* selected = selectedKey(values, key, legacy)) {
			Real const value = values[selected].as<Real>();
			if (!isfinite(value)) throw std::invalid_argument(std::string("Invalid number for ") + key);
			target = value;
		}
	}

	void readBoolean(po::variables_map const& values, char const* key, bool& target, char const* legacy = nullptr) {
		if (char const* selected = selectedKey(values, key, legacy)) {
			std::string const value = values[selected].as<std::string>();
			if (value == "on" || value == "true")
				target = true;
			else if (value == "off" || value == "false")
				target = false;
			else
				throw std::invalid_argument(std::string("Expected on/off for ") + key + ": " + value);
		}
	}

	void applySettings(Config& config, po::variables_map const& values) {
		if (values.count("problem.name") || values.count("mesh.ndim"))
			throw std::invalid_argument("Problem and dimension are build choices; select OCTOTIGERII_PROBLEM and OCTOTIGERII_NDIM in CMake");

		readOption(values, "randomSeed", config.randomSeed);
		readOption(values, "verification.analytic", config.verification.analytic);
		readOption(values, "verification.gravityReference", config.verification.gravityReference);
		readOption(values, "verification.directMaxPairs", config.verification.directMaxPairs);
		readOption(values, "verification.directSamples", config.verification.directSamples);
		readNumber(values, "verification.relativeL1Tolerance", config.verification.relativeL1Tolerance);
		readNumber(values, "verification.absoluteTolerance", config.verification.absoluteTolerance);
		readOption(values, "mesh.cells", config.mesh.cells);
		readOption(values, "mesh.level", config.mesh.level);
		readBoolean(values, "mesh.periodic", config.mesh.periodic);
		Real lower = units::value(config.mesh.lower), upper = units::value(config.mesh.upper);
		readNumber(values, "mesh.lower", lower);
		readNumber(values, "mesh.upper", upper);
		config.mesh.lower = units::Length::from_value(lower);
		config.mesh.upper = units::Length::from_value(upper);

		Real stopTime = units::value(config.runtime.stopTime);
		readNumber(values, "runtime.stopTime", stopTime, "runtime.stop_time");
		config.runtime.stopTime = units::Time::from_value(stopTime);
		readOption(values, "runtime.maxSteps", config.runtime.maxSteps, "runtime.max_steps");
		readOption(values, "runtime.workerTasks", config.runtime.workerTasks, "runtime.worker_tasks");
		readBoolean(values, "runtime.workStealing", config.runtime.workStealing, "runtime.work_stealing");

		readNumber(values, "timestep.cfl", config.timestep.cfl);
		readNumber(values, "hydro.gamma", config.hydro.gamma);
		readNumber(values, "radiation.lightSpeedRatio", config.radiation.lightSpeedRatio, "radiation.light_speed_ratio");
		readOption(values, "gravity.multipoleOrder", config.gravity.multipoleOrder, "gravity.multipole_order");
		readNumber(values, "gravity.openingAngle", config.gravity.openingAngle, "gravity.opening_angle");
		readBoolean(values, "output.enabled", config.output.enabled);
		readOption(values, "output.every", config.output.every);
		readOption(values, "output.directory", config.output.directory);
	}

}	 // namespace

void Config::validate() const {
	using std::isfinite;
	using std::sqrt;

	if (verification.analytic != "auto" && verification.analytic != "on" && verification.analytic != "off")
		throw std::invalid_argument("verification.analytic must be auto, on, or off");
	if (!isfinite(verification.relativeL1Tolerance) || !isfinite(verification.absoluteTolerance) ||
		(verification.relativeL1Tolerance < 0 && verification.relativeL1Tolerance != -1) || verification.absoluteTolerance < 0 ||
		(verification.analytic == "off" && verification.relativeL1Tolerance >= 0))
		throw std::invalid_argument("Invalid analytic verification tolerance or disabled accuracy gate");
	if (verification.gravityReference != "direct" && verification.gravityReference != "continuum")
		throw std::invalid_argument("verification.gravityReference must be direct or continuum");
	if (verification.directMaxPairs < 1 || verification.directSamples < 0)
		throw std::invalid_argument("Direct pair budget must be positive and directSamples nonnegative");
	if (randomSeed < 0) throw std::invalid_argument("randomSeed must be nonnegative");
	validateProblem(*this);
	if (mesh.cells < 4 || mesh.cells > 128 || (mesh.cells & (mesh.cells - 1)) || mesh.level < 0 || mesh.level > 6)
		throw std::invalid_argument("mesh: cells=power of two in [4,128], level=0..6");
	if (!(mesh.upper > mesh.lower) || !units::finite(mesh.upper - mesh.lower) || !(runtime.stopTime >= units::Time{}) ||
		!(timestep.cfl > 0 && timestep.cfl <= 0.5) || !(hydro.gamma > 1) || !(radiation.lightSpeedRatio > 0 && radiation.lightSpeedRatio <= 1) ||
		runtime.maxSteps < 1 || output.every < 1 || runtime.workerTasks < 0)
		throw std::invalid_argument("Invalid domain, timestep, gas, radiation, or output setting");
	if (gravityEnabled() && mesh.periodic) throw std::invalid_argument("Gravity currently requires 3D isolated boundaries");
	if (gravity.multipoleOrder < 1 || gravity.multipoleOrder > 10 || !(gravity.openingAngle > 0 && gravity.openingAngle < 1 / sqrt(3.0)))
		throw std::invalid_argument("Gravity requires order 1..10 and 0<openingAngle<1/sqrt(3)");
	if (output.directory.empty()) throw std::invalid_argument("Empty output directory");
}

Config parseConfig(std::vector<std::string> const& arguments) {
	po::options_description settings("Simulation options");
	addSettings(settings);
	po::options_description legacy("Legacy aliases");
	addLegacySettings(legacy);
	po::options_description buildChoices("Build-time choices");
	buildChoices.add_options()("problem.name", po::value<std::string>(), "Select OCTOTIGERII_PROBLEM in CMake")(
		"mesh.ndim", po::value<int>(), "Select OCTOTIGERII_NDIM in CMake");
	po::options_description iniOptions;
	iniOptions.add(settings).add(legacy).add(buildChoices);
	po::options_description commandOptions;
	commandOptions.add(iniOptions);
	commandOptions.add_options()("config", po::value<std::vector<std::string>>()->composing(), "INI file (repeatable)");

	po::variables_map commandValues;
	po::store(po::command_line_parser(arguments)
				  .options(commandOptions)
				  .style(po::command_line_style::allow_long | po::command_line_style::long_allow_adjacent)
				  .run(),
		commandValues);
	po::notify(commandValues);

	Config config;
	problemDefaults(config);
	if (commandValues.count("config")) {
		for (std::string const& path : commandValues["config"].as<std::vector<std::string>>()) {
			std::ifstream input(path);
			if (!input) throw std::runtime_error("Cannot open config: " + path);
			po::variables_map iniValues;
			po::store(po::parse_config_file(input, iniOptions), iniValues);
			po::notify(iniValues);
			applySettings(config, iniValues);
		}
	}
	applySettings(config, commandValues);
	config.validate();
	return config;
}

std::string helpText() {
	po::options_description settings("Simulation options");
	addSettings(settings);
	std::ostringstream output;
	output << build::executable << " (CGS, " << ndim << "D)\n"
		   << "Usage: " << build::executable << " [--config=/path/to/bin/problem/inputs] [--key=value ...]\n"
		   << "Problem and dimension are fixed at build time. CLI values override INI files.\n"
		   << "Settings use dotted groups; snake_case names remain supported as aliases.\n"
		   << "Booleans: on/off. mesh.cells is cells per block per active axis.\n"
		   << "Cartesian blocks: 2^mesh.level blocks per axis.\n\n"
		   << settings;
	return output.str();
}

}	 // namespace octotigerII
