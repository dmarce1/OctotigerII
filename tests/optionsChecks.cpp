#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "octotigerII/config.hpp"

namespace {

void require(bool condition, std::string const& message) {
	if (!condition) throw std::runtime_error(message);
}

void expectFailure(std::vector<std::string> const& arguments) {
	try {
		(void) octotigerII::parseConfig(arguments);
	} catch (std::exception const&) {
		return;
	}
	std::string names;
	for (std::string const& argument : arguments)
		names += " " + argument;
	throw std::runtime_error("Expected an invalid option to be rejected:" + names);
}

}	 // namespace

int main(int argc, char** argv) {
	using octotigerII::parseConfig;

	try {
		require(argc == 2, "Expected a temporary output directory");
		std::filesystem::path const root(argv[1]);
		std::filesystem::create_directories(root);
		auto const first = root / "first.ini";
		auto const second = root / "second.ini";
		{
			std::ofstream file(first);
			file << "[gravity]\nmultipoleOrder=3\nopeningAngle=0.4\n"
				 << "[mesh]\ncells=8\n"
				 << "[runtime]\nworkerTasks=2\nworkStealing=off\n";
			require(bool(file), "Cannot write first INI");
		}
		{
			std::ofstream file(second);
			file << "gravity.multipole_order=4\ngravity.opening_angle=0.45\n";
			require(bool(file), "Cannot write second INI");
		}

		auto const ini = parseConfig({"--config=" + first.string()});
		require(ini.gravity.multipoleOrder == 3 && ini.gravity.openingAngle == 0.4 && ini.mesh.cells == 8 && ini.runtime.workerTasks == 2 &&
				!ini.runtime.workStealing,
			"Hierarchical INI sections");

		auto const fromCli =
			parseConfig({"--gravity.multipoleOrder=5", "--config=" + first.string(), "--gravity.openingAngle=0.5", "--runtime.workStealing=on"});
		require(fromCli.gravity.multipoleOrder == 5 && fromCli.gravity.openingAngle == 0.5 && fromCli.runtime.workStealing,
			"CLI must override INI regardless of argument order");

		auto const legacyCli = parseConfig({"--config=" + first.string(), "--gravity.multipole_order=4"});
		require(legacyCli.gravity.multipoleOrder == 4, "Legacy CLI name must override canonical INI name");

		auto const layered = parseConfig({"--config=" + first.string(), "--config=" + second.string(), "--gravity.multipoleOrder=5"});
		require(layered.gravity.multipoleOrder == 5 && layered.gravity.openingAngle == 0.45, "Later INIs and CLI precedence");

		auto const legacyIni = parseConfig({"--config=" + second.string(), "--gravity.openingAngle=0.4"});
		require(legacyIni.gravity.multipoleOrder == 4 && legacyIni.gravity.openingAngle == 0.4, "Canonical CLI name must override legacy INI name");
		require(parseConfig({"--mesh.cells", "8"}).mesh.cells == 8, "Program_options accepts a value as the next argument");

		expectFailure({"--not.an.option=1"});
		expectFailure({"--problem.name=sod"});
		expectFailure({"--mesh.ndim=3"});
		for (int p = 1; p <= 10; ++p)
			require(parseConfig({"--gravity.multipoleOrder=" + std::to_string(p)}).gravity.multipoleOrder == p, "Orders 1..10 must be accepted");
		expectFailure({"--gravity.multipoleOrder=0"});
		expectFailure({"--gravity.multipoleOrder=11"});
		expectFailure({"--verification.gravityReference=invalid"});
		expectFailure({"--verification.directSamples=-1"});
		expectFailure({"--verification.directMaxPairs=0"});
		expectFailure({"--randomSeed=-1"});
		auto const direct = parseConfig({"--randomSeed=42", "--verification.directSamples=17", "--verification.directMaxPairs=1234"});
		require(
			direct.randomSeed == 42 && direct.verification.directSamples == 17 && direct.verification.directMaxPairs == 1234, "Direct sampling CLI options");
		expectFailure({"--gravity.multipoleOrder=3", "--gravity.multipole_order=3"});
		expectFailure({"--gravity.multipoleO=4"});
		expectFailure({"--config=" + first.string(), "--gravity.openingAngle=nan"});
		expectFailure({"--hydro.gamma=inf"});
		expectFailure({"--runtime.workStealing=yes"});
		require(octotigerII::helpText().find("gravity.multipoleOrder") != std::string::npos, "Generated help must list canonical option names");
		std::cout << "Hierarchical Program_options parsing passed\n";
		return 0;
	} catch (std::exception const& error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
}
