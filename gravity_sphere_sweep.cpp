// Standalone OctotigerII gravity-sphere parameter sweep and plotting driver.
//
// Build:
//   g++ -std=c++17 -O2 gravity_sphere_sweep.cpp
//
// Run:
//   ./a.out SPHERE_EXE ORDER_MIN ORDER_MAX ANGLE_MIN ANGLE_MAX ANGLE_STEP
//       [--sweep.repeats=N] [additional OctotigerII options]
//
// Example:
//   ./a.out ./release/octoII-3d
//       2 8 0.10 0.55 0.05 --mesh.cells=8 --mesh.level=2
//       --verification.directSamples=2048 --randomSeed=5489 --hpx:threads=12
//
// Produces gravity-errors.png and a separate gravity-timings.png page.
// Requires an HPX/APEX-enabled Octo-II with profiling enabled; this driver
// itself needs only the C++17 standard library and gnuplot at runtime.
// Reads gravity.solve.wall_ns from APEX's TAU-format user events (no TAU
// installation needed). Setup, reference evaluation and output are excluded;
// communication, waits and first-solve operator-cache construction are included.
// Default: one process/solve per point. --sweep.repeats=N repeats fresh processes,
// NOT warm-cache solves; plot shows median and min/max across these solves.
// APEX aggregate reports cannot discard the first solve or recover a median
// within one process, so this static-sphere driver requires one solve per run.
// Current Octo-II AMR bypasses this sample; missing timing is an error.
// Raw APEX reports and per-repeat timing data are retained for inspection.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>


namespace {


class Errors {
public:
	std::array<double, 3> potential{};
	std::array<double, 3> acceleration{};
};


std::string shellQuote(std::string const& value) {
	std::string result = "'";
	for (char const c : value) {
		if (c == '\'')
			result += "'\\''";
		else
			result += c;
	}
	return result + "'";
}


std::string readFile(std::filesystem::path const& path) {
	std::ifstream input(path);
	if (!input) throw std::runtime_error("Cannot read " + path.string());
	return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}


std::string fieldObject(std::string const& json, std::string const& field) {
	std::string const marker = "\"name\": \"" + field + "\"";
	auto const name = json.find(marker);
	if (name == std::string::npos) throw std::runtime_error("Missing field " + field + " in analytic-errors.json");
	auto const begin = json.rfind('{', name);
	auto const end = json.find('}', name);
	if (begin == std::string::npos || end == std::string::npos)
		throw std::runtime_error("Malformed field " + field + " in analytic-errors.json");
	return json.substr(begin, end - begin + 1);
}


double jsonNumber(std::string const& object, std::string const& key) {
	std::regex const expression("\\\"" + key + "\\\"[[:space:]]*:[[:space:]]*"
		+ "([-+]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][-+]?[0-9]+)?|null)");
	std::smatch match;
	if (!std::regex_search(object, match, expression)) throw std::runtime_error("Missing JSON value " + key);
	if (match[1] == "null") return std::numeric_limits<double>::quiet_NaN();
	return std::stod(match[1]);
}


std::array<double, 3> fieldErrors(std::string const& json, std::string const& field) {
	auto const object = fieldObject(json, field);
	return {jsonNumber(object, "L1"), jsonNumber(object, "L2"), jsonNumber(object, "Linf")};
}


Errors readErrors(std::filesystem::path const& report) {
	auto const json = readFile(report);
	Errors errors;
	errors.potential = fieldErrors(json, "potential");
	for (char const* field : {"accelerationX", "accelerationY", "accelerationZ"}) {
		auto const component = fieldErrors(json, field);
		for (std::size_t norm = 0; norm < errors.acceleration.size(); ++norm)
			errors.acceleration[norm] = std::max(errors.acceleration[norm], component[norm]);
	}
	return errors;
}


// Only read the root locality's completed-solve sample. Worker task totals
// overlap in time and are not a substitute for this elapsed-wall measurement.
double readSolveSeconds(std::filesystem::path const& directory) {
	auto const path = directory / "profile.0.0.0";
	std::ifstream input(path);
	if (!input) throw std::runtime_error("Missing APEX profile: " + path.string()
		+ ". Build HPX with HPX_WITH_APEX=ON and Octo-II with OCTOTIGERII_WITH_PROFILING=ON.");
	bool userEvents = false;
	bool found = false;
	double seconds = 0;
	std::string line;
	while (std::getline(input, line)) {
		if (line.find("# eventname numevents max min mean sumsqr") != std::string::npos) {
			userEvents = true;
			continue;
		}
		if (!userEvents) continue;
		std::istringstream row(line);
		std::string name;
		row >> std::quoted(name);
		if (name != "gravity.solve.wall_ns") continue;
		double count, maximum, minimum, mean, sumSquares;
		if (found || !(row >> count >> maximum >> minimum >> mean >> sumSquares)
			|| !std::isfinite(count) || !std::isfinite(mean) || !std::isfinite(minimum)
			|| !std::isfinite(maximum) || !(minimum > 0) || minimum > maximum
			|| mean < minimum || mean > maximum)
			throw std::runtime_error("Invalid gravity.solve.wall_ns sample in " + path.string());
		if (count != 1)
			throw std::runtime_error("Expected one FMM solve per static-sphere run in " + path.string()
				+ "; multiple aggregated solves cannot supply a median or exclude a warm-up.");
		seconds = mean * 1e-9; // User-event values retain their original nanosecond units.
		found = true;
	}
	if (!found) throw std::runtime_error("No gravity.solve.wall_ns user event in " + path.string()
		+ ". Check profiling is enabled. Current AMR bypasses this measurement; use a uniform grid"
		  " or move the solver's Elapsed guard before its adaptive branch. No total-runtime fallback is used.");
	return seconds;
}


double median(std::vector<double> values) {
	std::sort(values.begin(), values.end());
	auto const middle = values.size() / 2;
	return values.size() % 2 ? values[middle] : (values[middle - 1] + values[middle]) / 2;
}


std::string timeStamp() {
	auto const now = std::chrono::system_clock::now();
	auto const time = std::chrono::system_clock::to_time_t(now);
	std::tm value{};
#ifdef _WIN32
	localtime_s(&value, &time);
#else
	localtime_r(&time, &value);
#endif
	std::ostringstream result;
	result << std::put_time(&value, "%Y%m%d-%H%M%S");
	return result.str();
}


std::string angleTag(double const angle) {
	std::ostringstream result;
	result << std::fixed << std::setprecision(8) << angle;
	auto text = result.str();
	while (!text.empty() && text.back() == '0') text.pop_back();
	if (!text.empty() && text.back() == '.') text.pop_back();
	std::replace(text.begin(), text.end(), '.', '_');
	return text;
}


void writeGnuplot(std::filesystem::path const& root, int const minimumOrder, int const maximumOrder) {
	std::ofstream output(root / "gravity-errors.gnuplot");
	if (!output) throw std::runtime_error("Cannot write gnuplot script");
	output << "set terminal pngcairo size 1800,1000 enhanced font ',14'\n"
		   << "set output 'gravity-errors.png'\n"
		   << "set multiplot layout 2,3 title 'Gravity-sphere FMM error sweep' font ',18'\n"
		   << "set logscale y\n"
		   << "set grid xtics ytics\n"
		   << "set xlabel 'Opening angle'\n"
		   << "set ylabel 'Relative error'\n"
		   << "set key outside right\n"
		   << "set format y '10^{%L}'\n"
		   << "pmin=" << minimumOrder << "\n"
		   << "pmax=" << maximumOrder << "\n"
		   << "file='gravity-errors.dat'\n"
		   << "set title 'Potential L1'\n"
		   << "plot for [p=pmin:pmax] file using ($1==p?$2:1/0):3 with linespoints title sprintf('p=%d',p)\n"
		   << "set title 'Potential L2'\n"
		   << "plot for [p=pmin:pmax] file using ($1==p?$2:1/0):4 with linespoints title sprintf('p=%d',p)\n"
		   << "set title 'Potential Linf'\n"
		   << "plot for [p=pmin:pmax] file using ($1==p?$2:1/0):5 with linespoints title sprintf('p=%d',p)\n"
		   << "set title 'Acceleration L1 (maximum component)'\n"
		   << "plot for [p=pmin:pmax] file using ($1==p?$2:1/0):6 with linespoints title sprintf('p=%d',p)\n"
		   << "set title 'Acceleration L2 (maximum component)'\n"
		   << "plot for [p=pmin:pmax] file using ($1==p?$2:1/0):7 with linespoints title sprintf('p=%d',p)\n"
		   << "set title 'Acceleration Linf (maximum component)'\n"
		   << "plot for [p=pmin:pmax] file using ($1==p?$2:1/0):8 with linespoints title sprintf('p=%d',p)\n"
		   << "unset multiplot\n";
}


void writeTimingGnuplot(std::filesystem::path const& root, int minimumOrder, int maximumOrder, int repeats) {
	std::ofstream output(root / "gravity-timings.gnuplot");
	if (!output) throw std::runtime_error("Cannot write timing gnuplot script");
	output << "set terminal pngcairo size 1400,900 enhanced font ',14'\n"
		   << "set output 'gravity-timings.png'\n"
		   << "set title 'Gravity-sphere FMM timing: median and min/max of " << repeats << " fresh-process solve(s)'\n"
		   << "set xlabel 'Opening angle'\n"
		   << "set ylabel 'FMM solve wall time (s)'\n"
		   << "set logscale y\nset format y '10^{%L}'\n"
		   << "set grid xtics ytics\nset key outside right\n"
		   << "plot for [p=" << minimumOrder << ':' << maximumOrder
		   << "] 'gravity-timings.dat' using ($1==p?$2:1/0):3:4:5 with yerrorlines title sprintf('p=%d',p)\n";
}


void usage(char const* executable) {
	std::cerr << "Usage: " << executable
			  << " SPHERE_EXE ORDER_MIN ORDER_MAX ANGLE_MIN ANGLE_MAX ANGLE_STEP"
				 " [--sweep.repeats=N] [additional OctotigerII options]\n\n"
			  << "Example:\n  " << executable
			  << " ./release/octoII-3d"
				 " 2 8 0.10 0.55 0.05 --mesh.cells=8 --mesh.level=2"
				 " --verification.directSamples=2048 --hpx:threads=12\n\n"
			  << "Writes separate error and FMM timing PNG pages. Repeats default to 1;\n"
				 "repeats are fresh-process solves, not warm-cache solves. Requires APEX.\n";
}


} // namespace


int main(int argc, char** argv) {
	try {
		using std::isfinite;

		if (argc == 2 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
			usage(argv[0]);
			return 0;
		}
		if (argc < 7) {
			usage(argv[0]);
			return 2;
		}

		auto const sphereExecutable = std::filesystem::absolute(argv[1]);
		int const minimumOrder = std::stoi(argv[2]);
		int const maximumOrder = std::stoi(argv[3]);
		double const minimumAngle = std::stod(argv[4]);
		double const maximumAngle = std::stod(argv[5]);
		double const angleStep = std::stod(argv[6]);

		if (!std::filesystem::is_regular_file(sphereExecutable))
			throw std::invalid_argument("Sphere executable does not exist: " + sphereExecutable.string());
		if (minimumOrder < 1 || maximumOrder < minimumOrder)
			throw std::invalid_argument("Require 1 <= ORDER_MIN <= ORDER_MAX");
		if (!isfinite(minimumAngle) || !isfinite(maximumAngle) || !isfinite(angleStep) || minimumAngle <= 0
			|| maximumAngle < minimumAngle || angleStep <= 0)
			throw std::invalid_argument("Require 0 < ANGLE_MIN <= ANGLE_MAX and ANGLE_STEP > 0");

		std::vector<std::string> additionalArguments;
		int repeats = 1;
		for (int i = 7; i < argc; ++i) {
			std::string argument = argv[i];
			if (argument.rfind("--sweep.repeats=", 0) == 0) {
				auto const value = argument.substr(16);
				std::size_t end = 0;
				repeats = std::stoi(value, &end);
				if (end != value.size() || repeats < 1)
					throw std::invalid_argument("--sweep.repeats must be a positive integer");
			} else if (argument.rfind("--sweep.", 0) == 0) {
				throw std::invalid_argument("Unknown sweep option: " + argument);
			} else additionalArguments.push_back(argument);
		}

		auto const root = std::filesystem::absolute("gravity-sphere-sweep-" + timeStamp());
		if (!std::filesystem::create_directory(root))
			throw std::runtime_error("Output directory already exists; retry after one second: " + root.string());
		std::ofstream data(root / "gravity-errors.dat");
		if (!data) throw std::runtime_error("Cannot write gravity-errors.dat");
		data << "# order openingAngle potentialL1 potentialL2 potentialLinf"
				 " accelerationL1Max accelerationL2Max accelerationLinfMax\n";
		data << std::scientific << std::setprecision(17);

		std::ofstream timings(root / "gravity-timings.dat");
		std::ofstream rawTimings(root / "gravity-timing-runs.dat");
		if (!timings || !rawTimings) throw std::runtime_error("Cannot write timing data");
		timings << "# order openingAngle medianSeconds minSeconds maxSeconds repeats\n"
			<< std::scientific << std::setprecision(17);
		rawTimings << "# order openingAngle repeat solveSeconds\n"
			<< std::scientific << std::setprecision(17);

		long long const steps = static_cast<long long>(std::floor((maximumAngle - minimumAngle) / angleStep + 1e-10)) + 1;
		for (int order = minimumOrder; order <= maximumOrder; ++order) {
			for (long long step = 0; step < steps; ++step) {
				double const angle = minimumAngle + step * angleStep;
				if (angle > maximumAngle + 1e-12 * std::max(1.0, std::abs(maximumAngle))) break;

				std::ostringstream runName;
				runName << 'p' << std::setfill('0') << std::setw(2) << order << "-theta-" << angleTag(angle);
				auto const pointDirectory = root / runName.str();
				std::filesystem::create_directories(pointDirectory);
				std::vector<double> solveTimes;
				Errors errors;
				for (int repeat = 0; repeat < repeats; ++repeat) {
					auto const runDirectory = repeat == 0 ? pointDirectory
						: pointDirectory / ("repeat-" + std::to_string(repeat + 1));
					std::filesystem::create_directories(runDirectory);
					auto const apexDirectory = runDirectory / "apex";
					std::filesystem::create_directories(apexDirectory);
					auto const log = runDirectory / "run.log";

					std::ostringstream command;
					command << "APEX_PROFILE_OUTPUT=1 APEX_CSV_OUTPUT=1 APEX_SCREEN_OUTPUT=0 "
						<< "APEX_OUTPUT_FILE_PATH=" << shellQuote(apexDirectory.string()) << ' '
						<< shellQuote(sphereExecutable.string());
					for (auto const& argument : additionalArguments) command << ' ' << shellQuote(argument);
					command << " --problem.name=gravity-sphere --runtime.stopTime=0"
						<< " --verification.analytic=on"
						<< " --verification.gravityReference=direct"
						<< " --verification.relativeL1Tolerance=-1"
						<< " --output.enabled=on"
						<< " --gravity.multipoleOrder=" << order
						<< " --gravity.openingAngle=" << std::setprecision(17) << angle
						<< ' ' << shellQuote("--output.directory=" + runDirectory.string())
						<< " > " << shellQuote(log.string()) << " 2>&1";

					std::cout << "Running p=" << order << " openingAngle=" << angle
						<< " repeat=" << repeat + 1 << '/' << repeats << " ... " << std::flush;
					if (std::system(command.str().c_str()) != 0)
						throw std::runtime_error("OctotigerII failed; see " + log.string());
					auto const runErrors = readErrors(runDirectory / "analytic-errors.json");
					if (repeat == 0) errors = runErrors;
					double const seconds = readSolveSeconds(apexDirectory);
					solveTimes.push_back(seconds);
					rawTimings << order << ' ' << angle << ' ' << repeat + 1 << ' ' << seconds << '\n';
					rawTimings.flush();
					std::cout << "FMM " << std::scientific << std::setprecision(6) << seconds
						<< " s\n" << std::defaultfloat;
				}
				auto const bounds = std::minmax_element(solveTimes.begin(), solveTimes.end());
				timings << order << ' ' << angle << ' ' << median(solveTimes)
					<< ' ' << *bounds.first << ' ' << *bounds.second << ' ' << repeats << '\n';
				timings.flush();

				data << order << ' ' << angle;
				for (double const error : errors.potential) data << ' ' << error;
				for (double const error : errors.acceleration) data << ' ' << error;
				data << '\n';
				data.flush();
			}
		}

		writeGnuplot(root, minimumOrder, maximumOrder);
		writeTimingGnuplot(root, minimumOrder, maximumOrder, repeats);
		std::string const plotCommand = "cd " + shellQuote(root.string()) + " && gnuplot gravity-errors.gnuplot gravity-timings.gnuplot";
		if (std::system(plotCommand.c_str()) != 0)
			throw std::runtime_error("gnuplot failed (install it with: sudo apt install gnuplot); data and script remain in "
				+ root.string());

		std::cout << "\nData: " << root / "gravity-errors.dat"
				  << "\nPlot: " << root / "gravity-errors.png"
				  << "\nTiming data: " << root / "gravity-timings.dat"
				  << "\nTiming plot: " << root / "gravity-timings.png" << '\n';
		return 0;
	} catch (std::exception const& error) {
		std::cerr << "gravity_sphere_sweep: " << error.what() << '\n';
		return 1;
	}
}
