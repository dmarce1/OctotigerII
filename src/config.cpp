#include "octotigerII/config.hpp"
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "octotigerII/gravity/boundary.hpp"
#include "octotigerII/problems.hpp"
#include "octotigerII/verification/analytic.hpp"

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
		options("problem.name", po::value<std::string>(), "Problem to run (required in an INI file or on the command line)");
		options("randomSeed", po::value<std::int64_t>(), "Global nonnegative random seed (default 5489)");
		options("verification.analytic", po::value<std::string>(), "Analytic comparison: auto (default), on (required), off");
		options("verification.gravityReference", po::value<std::string>(), "Gravity reference: direct (default) or continuum");
		options("verification.directMaxPairs", po::value<std::int64_t>(), "Direct-reference pair budget (default 20000000; at least one complete target)");
		options("verification.directSamples", po::value<std::int64_t>(),
			"Number of direct-reference targets; 0 selects automatically, positive overrides the pair budget");
		options("verification.relativeL1Tolerance", po::value<Real>(), "Maximum relative L1 for each reference field; -1 disables the accuracy gate");
		options("verification.absoluteTolerance", po::value<Real>(), "Maximum absolute Linf for zero-reference fields (CGS), when the gate is enabled");
		options("mesh.cells", po::value<int>(), "Cells per block per active axis");
		options("mesh.level", po::value<int>(), "Initial block level (also the default AMR minimum)");
		options("amr.enabled", po::value<std::string>(), "Adaptive mesh refinement: on/off");
		options("amr.minLevel", po::value<int>(), "Coarsest block level; -1 uses mesh.level");
		options("amr.maxLevel", po::value<int>(), "Finest allowed block level (default 6)");
		options("amr.regridEvery", po::value<int>(), "Maximum synchronized timesteps between regrids (default 4)");
		options("amr.refineDensity", po::value<Real>(), "Refine cells above this density (g/cm^3); 0 disables density refinement");
		options("amr.maxCellMass", po::value<Real>(), "Maximum mass per active cell (g); 0 disables mass refinement");
		options("amr.shadowTolerance", po::value<Real>(), "Relative fine/shadow difference; 0 disables shadow refinement");
		options("amr.shadowFloor", po::value<Real>(), "Normalization floor as a fraction of each field's maximum magnitude");
		options("amr.coarsenFactor", po::value<Real>(), "Coarsening threshold relative to refinement threshold (default 0.25)");
		options("amr.signalBuffer", po::value<Real>(), "Safety factor for signal travel before the next regrid (at least 1)");
		options("amr.bufferCells", po::value<int>(), "Additional cells around refinement tags (default 1)");
		options("amr.hydro", po::value<std::string>(), "Use hydro fields for shadow refinement: on/off");
		options("amr.radiation", po::value<std::string>(), "Use radiation fields for shadow refinement: on/off");
		options("mesh.lower", po::value<Real>(), "Lower domain coordinate (cm)");
		options("mesh.upper", po::value<Real>(), "Upper domain coordinate (cm)");
		options("mesh.periodic", po::value<std::string>(), "Legacy shorthand: on=all periodic, off=all outflow; per-face settings override");
		for (int axis = 0; axis < ndim; ++axis)
			for (auto side : {"Lower", "Upper"}) {
				auto const key = std::string("mesh.boundary.") + "xyz"[axis] + side;
				options(key.c_str(), po::value<std::string>(), "Face boundary: periodic, reflecting, outflow, inflow, analytic");
			}
		options("runtime.stopTime", po::value<Real>(), "Stop time (s)");
		options("runtime.maxSteps", po::value<int>(), "Maximum number of steps");
		options("runtime.workerTasks", po::value<int>(), "Maximum concurrent tasks; 0 uses worker count");
		options("runtime.workStealing", po::value<std::string>(), "Remote work stealing: on/off");
		options("timestep.cfl", po::value<Real>(), "Courant factor");
		options("timestep.refinement", po::value<std::string>(), "Dyadic level time refinement: on/off (default on; transport and self-gravity without external acceleration)");
		options("massFractions.enabled", po::value<std::string>(), "Material partial densities and massless tracers: on/off (default off)");
		options("massFractions.species", po::value<std::string>(), "Semicolon-separated name:initialFraction:element, A=mass,Z=number, or He=70%,O=30% definitions");
		options("hydro.gamma", po::value<Real>(), "Ideal-gas adiabatic index");
		options("hydro.meanMolecularWeight", po::value<Real>(), "Mean particle mass in atomic mass units, for temperature (default 1)");
		options("hydro.dualEnergy.enabled", po::value<std::string>(), "Dual energy: on/off (default on)");
		options("hydro.dualEnergy.exponent", po::value<Real>(), "Nonzero auxiliary entropy exponent (default 1)");
		options("hydro.dualEnergy.pressureThreshold", po::value<Real>(), "Use total-energy pressure when u/E exceeds this ratio (default 0.001)");
		options("hydro.dualEnergy.syncThreshold", po::value<Real>(), "Reset auxiliary from total energy when u/E exceeds this ratio (default 0.1)");
		for (int axis = 0; axis < ndim; ++axis) {
				auto const key = std::string("hydro.acceleration.") + "xyz"[axis];
				options(key.c_str(), po::value<Real>(), "Uniform external acceleration (cm/s^2)");
			}
		{
			options("rayleighTaylor.densityLower", po::value<Real>(), "Lower-layer density (g/cm^3)");
			options("rayleighTaylor.densityUpper", po::value<Real>(), "Upper-layer density (g/cm^3)");
			options("rayleighTaylor.interfacePressure", po::value<Real>(), "Pressure at the domain midpoint (dyn/cm^2)");
			options("rayleighTaylor.perturbation", po::value<Real>(), "Vertical velocity perturbation amplitude (cm/s)");
		}
		{
			for (int d = 0; d < ndim; ++d) {
				auto const key = std::string("star.center.") + "xyz"[d];
				options(key.c_str(), po::value<Real>(), "Star center coordinate (cm); default is box midpoint");
			}
			options("star.radius", po::value<Real>(), "Lane-Emden surface radius (cm)");
			options("star.centralDensity", po::value<Real>(), "Central density (g/cm^3)");
			options("star.polytropicIndex", po::value<Real>(), "Polytropic index: 0 < n < 5 (default 1.5)");
			options("star.atmosphereFraction", po::value<Real>(), "Ambient density divided by central density (default 1e-8)");
		}
		options("radiation.lightSpeedRatio", po::value<Real>(), "Radiation transport speed divided by c");
		options("gravity.multipoleOrder", po::value<int>(), "Gravity expansion order (1..10)");
		options("gravity.openingAngle", po::value<Real>(), "Gravity opening angle");
		options("gravity.timeIntegration", po::value<std::string>(), "Time-refined self-gravity: hierarchical (default) or conventional");
		options("gravity.energyTreatment", po::value<std::string>(), "Self-gravity energy: mullen (default) or naive kinetic-work kicks");
		options("gravity.conserveRegridEnergy", po::value<std::string>(), "Conserve gas plus gravitational energy during regridding: on (default)/off");
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

	template <typename Quantity>
	void readQuantity(po::variables_map const& values, char const* key, Quantity& target) {
		Real value = units::value(target);
		readNumber(values, key, value);
		target = Quantity::from_value(value);
	}

	void applySettings(Config& config, po::variables_map const& values) {


		readBoolean(values, "massFractions.enabled", config.massFractions.enabled);
		if (values.count("massFractions.species")) config.massFractions.species = composition::parseSpecies(values["massFractions.species"].as<std::string>());
		readOption(values, "randomSeed", config.randomSeed);
		readOption(values, "verification.analytic", config.verification.analytic);
		readOption(values, "verification.gravityReference", config.verification.gravityReference);
		readOption(values, "verification.directMaxPairs", config.verification.directMaxPairs);
		readOption(values, "verification.directSamples", config.verification.directSamples);
		readNumber(values, "verification.relativeL1Tolerance", config.verification.relativeL1Tolerance);
		readNumber(values, "verification.absoluteTolerance", config.verification.absoluteTolerance);
		readOption(values, "mesh.cells", config.mesh.cells);
		readOption(values, "mesh.level", config.mesh.level);
		readBoolean(values, "amr.enabled", config.amr.enabled);
		readBoolean(values, "amr.hydro", config.amr.hydro);
		readBoolean(values, "amr.radiation", config.amr.radiation);
		readOption(values, "amr.minLevel", config.amr.minLevel);
		readOption(values, "amr.maxLevel", config.amr.maxLevel);
		readOption(values, "amr.regridEvery", config.amr.regridEvery);
		readOption(values, "amr.bufferCells", config.amr.bufferCells);
		readQuantity(values, "amr.maxCellMass", config.amr.maxCellMass);
		readQuantity(values, "amr.refineDensity", config.amr.refineDensity);
		readNumber(values, "amr.shadowTolerance", config.amr.shadowTolerance);
		readNumber(values, "amr.shadowFloor", config.amr.shadowFloor);
		readNumber(values, "amr.coarsenFactor", config.amr.coarsenFactor);
		readNumber(values, "amr.signalBuffer", config.amr.signalBuffer);
		if (values.count("mesh.periodic")) {
			bool periodic = false;
			readBoolean(values, "mesh.periodic", periodic);
			config.mesh.boundary = periodic ? physics::BoundaryConditions::periodic() : physics::BoundaryConditions{};
		}
		for (int axis = 0; axis < ndim; ++axis)
			for (bool lower : {true, false}) {
				auto const key = std::string("mesh.boundary.") + "xyz"[axis] + (lower ? "Lower" : "Upper");
				if (values.count(key))
					(lower ? config.mesh.boundary.lower : config.mesh.boundary.upper)[axis] = physics::parseBoundaryCondition(values[key].as<std::string>());
			}
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
		readBoolean(values, "timestep.refinement", config.timestep.refinement);
		readNumber(values, "hydro.gamma", config.hydro.gamma);
		readNumber(values, "hydro.meanMolecularWeight", config.hydro.meanMolecularWeight);
		readBoolean(values, "hydro.dualEnergy.enabled", config.hydro.dualEnergy.enabled);
		readNumber(values, "hydro.dualEnergy.exponent", config.hydro.dualEnergy.exponent);
		readNumber(values, "hydro.dualEnergy.pressureThreshold", config.hydro.dualEnergy.pressureThreshold);
		readNumber(values, "hydro.dualEnergy.syncThreshold", config.hydro.dualEnergy.syncThreshold);
		for (int axis = 0; axis < ndim; ++axis) {
			auto const key = std::string("hydro.acceleration.") + "xyz"[axis];
			readQuantity(values, key.c_str(), config.hydro.acceleration[axis]);
		}
		readQuantity(values, "rayleighTaylor.densityLower", config.rayleighTaylor.densityLower);
		readQuantity(values, "rayleighTaylor.densityUpper", config.rayleighTaylor.densityUpper);
		readQuantity(values, "rayleighTaylor.interfacePressure", config.rayleighTaylor.interfacePressure);
		readQuantity(values, "rayleighTaylor.perturbation", config.rayleighTaylor.perturbation);
		for (int d = 0; d < ndim; ++d) {
			auto const key = std::string("star.center.") + "xyz"[d];
			readQuantity(values, key.c_str(), config.star.center[d]);
		}
		readQuantity(values, "star.radius", config.star.radius);
		readQuantity(values, "star.centralDensity", config.star.centralDensity);
		readNumber(values, "star.polytropicIndex", config.star.polytropicIndex);
		readNumber(values, "star.atmosphereFraction", config.star.atmosphereFraction);
		readNumber(values, "radiation.lightSpeedRatio", config.radiation.lightSpeedRatio, "radiation.light_speed_ratio");
		readOption(values, "gravity.multipoleOrder", config.gravity.multipoleOrder, "gravity.multipole_order");
		readNumber(values, "gravity.openingAngle", config.gravity.openingAngle, "gravity.opening_angle");
		readOption(values, "gravity.timeIntegration", config.gravity.timeIntegration);
		readOption(values, "gravity.energyTreatment", config.gravity.energyTreatment);
		readBoolean(values, "gravity.conserveRegridEnergy", config.gravity.conserveRegridEnergy);
		readBoolean(values, "output.enabled", config.output.enabled);
		readOption(values, "output.every", config.output.every);
		readOption(values, "output.directory", config.output.directory);
	}

}	 // namespace

void Config::validate() const {
	if (gravity.timeIntegration != "hierarchical" && gravity.timeIntegration != "conventional")
		throw std::invalid_argument("gravity.timeIntegration must be hierarchical or conventional");
	if (gravity.energyTreatment != "mullen" && gravity.energyTreatment != "naive")
		throw std::invalid_argument("gravity.energyTreatment must be mullen or naive");
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
	massFractions.validate();
	if (massFractions.enabled && (!build::massFractions || !hydroEnabled()))
		throw std::invalid_argument("massFractions requires a hydro problem and OCTOII_WITH_MASS_FRACTIONS=ON");
	hydro.dualEnergy.validate();
	if (!isfinite(hydro.gamma) || !isfinite(hydro.meanMolecularWeight) || !(hydro.meanMolecularWeight > 0))
		throw std::invalid_argument("Hydro gamma must be finite and meanMolecularWeight finite and positive");
	if (randomSeed < 0) throw std::invalid_argument("randomSeed must be nonnegative");
	for (auto component : hydro.acceleration) {
		if (!units::finite(component)) throw std::invalid_argument("External acceleration must be finite");
	}
	if (!hydroEnabled() && hasExternalAcceleration()) throw std::invalid_argument("External acceleration requires hydro");
	if (ndim != 3 && hasExternalAcceleration()) throw std::invalid_argument("Gravity requires a 3D build");
	mesh.boundary.validate();
	if (gravityEnabled()) gravity::validateBoundaries(mesh.boundary);
	validateProblem(*this);
	if (mesh.cells < 4 || mesh.cells > 128 || (mesh.cells & (mesh.cells - 1)) || mesh.level < 0 || mesh.level > 6)
		throw std::invalid_argument("mesh: cells=power of two in [4,128], level=0..6");
	int const minimumLevel = amr.minLevel < 0 ? mesh.level : amr.minLevel;
	if (amr.minLevel < -1 || minimumLevel > mesh.level || amr.maxLevel < minimumLevel || amr.maxLevel > 16 || (amr.enabled && amr.maxLevel < mesh.level) ||
		amr.regridEvery < 1 || amr.bufferCells < 0 || !(amr.maxCellMass >= units::Mass{}) || !units::finite(amr.maxCellMass) || !(amr.refineDensity >= units::Density{}) || !units::finite(amr.refineDensity) ||
		!isfinite(amr.shadowTolerance) || amr.shadowTolerance < 0 || !isfinite(amr.shadowFloor) || amr.shadowFloor <= 0 || !isfinite(amr.coarsenFactor) ||
		!(amr.coarsenFactor > 0 && amr.coarsenFactor < 1) || !isfinite(amr.signalBuffer) || amr.signalBuffer < 1)
		throw std::invalid_argument("Invalid AMR levels, criteria, buffering, or regrid interval");
	if (amr.enabled && (amr.maxCellMass > units::Mass{} || amr.refineDensity > units::Density{}) && !hydroEnabled() && !gravityEnabled())
		throw std::invalid_argument("Mass/density refinement requires a density field");
	if (!(mesh.upper > mesh.lower) || !units::finite(mesh.upper - mesh.lower) || !(runtime.stopTime >= units::Time{}) ||
		!(timestep.cfl > 0 && timestep.cfl <= 0.5) || !(hydro.gamma > 1) || !(radiation.lightSpeedRatio > 0 && radiation.lightSpeedRatio <= 1) ||
		runtime.maxSteps < 1 || output.every < 1 || runtime.workerTasks < 0)
		throw std::invalid_argument("Invalid domain, timestep, gas, radiation, or output setting");
	if (mesh.boundary.contains(physics::BoundaryCondition::Analytic) && !problemBoundary(*this))
		throw std::invalid_argument(std::string("Analytic boundary is not implemented for problem ") + problem);
	if (gravity.multipoleOrder < 1 || gravity.multipoleOrder > 10 || !(gravity.openingAngle > 0 && gravity.openingAngle < 1 / sqrt(3.0)))
		throw std::invalid_argument("Gravity requires order 1..10 and 0<openingAngle<1/sqrt(3)");
	if (output.directory.empty()) throw std::invalid_argument("Empty output directory");
}

Config parseConfig(std::vector<std::string> const& arguments) {
	po::options_description settings("Simulation options");
	addSettings(settings);
	po::options_description legacy("Legacy aliases");
	addLegacySettings(legacy);
	po::options_description iniOptions;
	iniOptions.add(settings).add(legacy);
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
	std::vector<po::variables_map> files;
	if (commandValues.count("config")) {
		for (auto const& path : commandValues["config"].as<std::vector<std::string>>()) {
			std::ifstream input(path);
			if (!input) throw std::runtime_error("Cannot open config: " + path);
			po::variables_map values;
			po::store(po::parse_config_file(input, iniOptions), values);
			po::notify(values);
			readOption(values, "problem.name", config.problem);
			files.push_back(std::move(values));
		}
	}
	readOption(commandValues, "problem.name", config.problem);
	if (config.problem.empty())
		throw std::invalid_argument("No problem specified. Set problem.name in an INI file or use --problem.name=<name>.");
	problemDefaults(config);
	bool lowerSet = false, upperSet = false, gammaSet = false, densitySet = false;
	std::array<bool, ndim> centerSet{};
	auto apply = [&](po::variables_map const& values) {
		for (int d = 0; d < ndim; ++d) centerSet[d] = centerSet[d] || values.count(std::string("star.center.") + "xyz"[d]);
		lowerSet = lowerSet || values.count("mesh.lower");
		upperSet = upperSet || values.count("mesh.upper");
		gammaSet = gammaSet || values.count("hydro.gamma");
		densitySet = densitySet || values.count("amr.refineDensity");
		applySettings(config, values);
	};
	for (auto const& values : files) apply(values);
	apply(commandValues);
	if (config.problem == "polytrope") {
		if (!lowerSet) config.mesh.lower = -2.0 * config.star.radius;
		if (!upperSet) config.mesh.upper = 2.0 * config.star.radius;
		for (int d = 0; d < ndim; ++d)
			if (!centerSet[d]) config.star.center[d] = (config.mesh.lower + config.mesh.upper) / 2.0;
		if (!gammaSet) config.hydro.gamma = 1 + 1 / config.star.polytropicIndex;
		if (!densitySet) config.amr.refineDensity = 0.01 * config.star.centralDensity;
	}
	config.validate();
	return config;
}

std::string helpText() {
	po::options_description settings("Simulation options");
	addSettings(settings);
	std::ostringstream output;
	output << build::executable << " (CGS, " << ndim << "D)\n"
		   << "Usage: " << build::executable << " [--config=/path/to/bin/problem/inputs] [--problem.name=<name>] [--key=value ...]\n"
		   << "A problem name is required in an INI file or on the command line; dimension is fixed by the executable. CLI values override INI files.\n"
		   << "Settings use dotted groups; snake_case names remain supported as aliases.\n"
		   << "Booleans: on/off. mesh.cells is cells per block per active axis.\n"
		   << "Initial Cartesian mesh: 2^mesh.level blocks per axis.\n\n"
		   << settings << "\nProblems:\n" << problemHelp();
	return output.str();
}

}	 // namespace octotigerII
