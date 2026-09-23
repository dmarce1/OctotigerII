#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include "octotigerII/gravity/fieldSolver.hpp"
#include "octotigerII/runtime.hpp"
#include "runtimeMain.hpp"
#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/runtime_distributed/find_all_localities.hpp>
#endif

using namespace octotigerII;

namespace {

void require(bool value, char const* message) {
	if (!value) throw std::runtime_error(message);
}

std::size_t globalIndex(Subgrid const& block, mesh::Coordinates c, int cells, int n) {
	for (int d = 0; d < ndim; ++d)
		c[d] += block.location.coordinates[d] * cells;
	return mesh::linearIndex(c, mesh::filledCoordinates(n));
}

std::vector<storage::Locality> localities() {
#ifdef OCTOTIGERII_WITH_HPX
	return hpx::find_all_localities();
#else
	return {0};
#endif
}

class Fixture {
public:
	Fixture(int cells, int level, int order, Real theta, int workers, int pattern)
	  : config(parseConfig({"--mesh.cells=" + std::to_string(cells), "--mesh.level=" + std::to_string(level), "--output.enabled=off"}))
	  , owners(localities())
	  , topology(config, owners.size())
	  , repository(config, topology.storageLayout(), owners) {
		config.gravity.multipoleOrder = order;
		config.gravity.openingAngle = theta;
		config.runtime.workerTasks = workers;
		n = cells * (1 << level);
		h = (config.mesh.upper - config.mesh.lower) / Real(n);
		density.resize(std::size_t(n) * n * n);
		for (std::size_t i = 0; i < density.size(); ++i) {
			Real const value = pattern == 0 ? 1 + Real((i * 37 + 13) % 101) / 73 : (pattern == 1 && i == 0 ? 2 : 0);
			density[i] = units::Density::from_value(value);
		}
		auto const& fields = repository.directory();
		for (auto const& block : topology.blocks()) {
			auto gravity = fields.gravity.output(block.interior, 0);
			for (std::size_t i = 0; i < block.interior.count; ++i)
				gravity.put(i, {});
			fields.gravity.commit(block.interior, 0, gravity);
			if constexpr (build::hydro) {
				auto gas = fields.hydro.output(block.interior, 0);
				block.layout.forEachInterior([&](mesh::Coordinates c, std::size_t i) {
					hydro::ConservedState state{};
					state.density() = density[globalIndex(block, c, cells, n)];
					state.totalEnergy() = units::EnergyDensity::from_value(19 + Real(i));
					gas.put(i, state);
				});
				fields.hydro.commit(block.interior, 0, gas);
			} else {
				auto rho = fields.density.output(block.interior, 0);
				block.layout.forEachInterior([&](mesh::Coordinates c, std::size_t i) { rho.data()[i] = density[globalIndex(block, c, cells, n)]; });
				fields.density.commit(block.interior, 0, rho);
			}
		}
		solver = std::make_unique<gravity::FieldSolver>(config, topology.blocks(), fields, owners);
	}

	std::vector<gravity::State> read(unsigned bank) const {
		std::vector<gravity::State> result(density.size());
		for (auto const& block : topology.blocks()) {
			auto input = repository.directory().gravity.read(block.interior, bank).get();
			block.layout.forEachInterior([&](mesh::Coordinates c, std::size_t i) { result[globalIndex(block, c, config.mesh.cells, n)] = input.at(i); });
		}
		return result;
	}

	void checkCopiedState() const {
		auto const& fields = repository.directory();
		for (auto const& block : topology.blocks()) {
			if constexpr (build::hydro) {
				auto a = fields.hydro.read(block.interior, 0).get();
				auto b = fields.hydro.read(block.interior, 1).get();
				for (std::size_t i = 0; i < block.interior.count; ++i)
					a.at(i).forEach([&](auto f, auto q) { require(q == b.at(i).template get<f>(), "Gravity publication changed gas state"); });
			} else {
				auto a = fields.density.read(block.interior, 0).get();
				auto b = fields.density.read(block.interior, 1).get();
				for (std::size_t i = 0; i < block.interior.count; ++i)
					require(a.data()[i] == b.data()[i], "Gravity publication changed density");
			}
		}
	}

	void firstDensity(units::Density value) {
		auto const range = topology.blocks().front().interior.slice(0, 1);
		auto const& fields = repository.directory();
		if constexpr (build::hydro) {
			auto input = fields.hydro.read(range, 0).get();
			auto state = input.at(0);
			state.density() = value;
			auto output = fields.hydro.output(range, 0);
			output.put(0, state);
			fields.hydro.commit(range, 0, output);
		} else {
			auto output = fields.density.output(range, 0);
			output.data()[0] = value;
			fields.density.commit(range, 0, output);
		}
	}

	Config config;
	std::vector<storage::Locality> owners;
	CartesianTopology topology;
	FieldRepository repository;
	std::unique_ptr<gravity::FieldSolver> solver;
	int n;
	units::Length h;
	std::vector<units::Density> density;
};

void compare(std::vector<gravity::State> const& actual, std::vector<gravity::State> const& expected, Real tolerance) {
	expected.front().forEach([&](auto f, auto q) {
		using std::abs;
		using std::isfinite;
		(void) q;
		Real scale = 0;
		for (auto const& cell : expected)
			scale = std::max(scale, abs(units::value(cell.template get<f>())));
		for (std::size_t i = 0; i < actual.size(); ++i) {
			auto const a = units::value(actual[i].template get<f>());
			auto const b = units::value(expected[i].template get<f>());
			require(isfinite(a) && isfinite(b) && abs(a - b) <= tolerance * scale, "Partitioned FMM field mismatch");
		}
	});
}

std::vector<gravity::State> check(int cells, int level, int order, Real theta, int workers, int pattern) {
	Fixture fixture(cells, level, order, theta, workers, pattern);
	auto const expected = gravity::solve(fixture.density, fixture.n, fixture.h, order, theta);
	auto const work = fixture.solver->solve(0);
	require(work.multipolePairs == expected.statistics.multipolePairs && work.directPairs == expected.statistics.directPairs,
		"Partitioned FMM lost or duplicated interactions");
	require(work.localityCells.size() == fixture.owners.size(), "Missing FMM locality reports");
	require(std::accumulate(work.localityCells.begin(), work.localityCells.end(), std::uint64_t(0)) == fixture.density.size(), "Missing FMM leaf targets");
	if (fixture.density.size() >= fixture.owners.size())
		for (auto count : work.localityCells)
			require(count > 0, "A locality did not evaluate its leaves");
	auto const actual = fixture.read(1);
	compare(actual, expected.fields, 3e-12);
	fixture.checkCopiedState();
	compare(fixture.read(0), std::vector<gravity::State>(actual.size()), 0);
	fixture.solver->solve(1);
	compare(fixture.read(0), actual, 0);
	// Invalid density must be reported across localities, leave the published
	// bank unchanged, drain all tasks, and permit a corrected retry.
	fixture.firstDensity(units::Density::from_value(-1));
	bool failed = false;
	try {
		fixture.solver->solve(0);
	} catch (std::exception const&) {
		failed = true;
	}
	require(failed, "Negative density was accepted");
	compare(fixture.read(0), actual, 0);
	fixture.firstDensity(fixture.density.front());
	fixture.solver->solve(0);
	compare(fixture.read(1), actual, 0);
	std::uint64_t digest = 14695981039346656037ull;
	for (auto const& cell : actual)
		cell.forEach([&](auto, auto q) { digest = (digest ^ std::bit_cast<std::uint64_t>(units::value(q))) * 1099511628211ull; });
	std::cout << "n=" << fixture.n << " p=" << order << " theta=" << theta << " workers=" << workers << " pattern=" << pattern
			  << " multipolePairs=" << work.multipolePairs << " directPairs=" << work.directPairs << " digest=" << digest << " cellsPerLocality=";
	for (auto count : work.localityCells)
		std::cout << ' ' << count;
	std::cout << '\n';
	return actual;
}

int runChecks(int argc, char** argv) {
	try {
		if (argc > 1 && std::string(argv[1]) == "benchmark") {
			int const level = argc > 2 ? std::stoi(argv[2]) : 2;
			int const order = argc > 3 ? std::stoi(argv[3]) : 3;
			Fixture fixture(8, level, order, 0.5, 0, 0);
			unsigned bank = 0;
			for (int iteration = 0; iteration < 4; ++iteration) {
				auto const start = std::chrono::steady_clock::now();
				auto const work = fixture.solver->solve(bank);
				auto const elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
				bank ^= 1;
				std::cout << std::scientific << std::setprecision(6) << "iteration=" << iteration << " seconds=" << elapsed
						  << " cells=" << fixture.density.size() << " workerTasks=" << work.workerTasks << " cellsPerLocality=";
				for (auto count : work.localityCells)
					std::cout << ' ' << count;
				std::cout << '\n';
			}
			return 0;
		}
		auto const single = check(4, 1, 3, 0.5, 1, 0);
		compare(check(4, 1, 3, 0.5, 4, 0), single, 0);
		// Same physical cells, one field block, and potentially several localities.
		compare(check(8, 0, 3, 0.5, 4, 0), single, 0);
		check(4, 1, 1, 0.57, 4, 0);
		check(4, 1, 10, 0.4, 4, 0);
		check(4, 0, 3, 0.1, 4, 1);
		check(4, 0, 3, 0.5, 4, 2);
		check(4, 2, 3, 0.5, 4, 0);
		std::cout << "Partitioned FMM checks passed\n";
		return 0;
	} catch (std::exception const& error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
}

}	 // namespace

int main(int argc, char** argv) {
	return runtimeMain(argc, argv, runChecks);
}
