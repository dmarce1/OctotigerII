#include "octotigerII/gravity/solver.hpp"
#include "octotigerII/subgrid/subgrid.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace octotigerII;
namespace {
void require(bool ok, char const* message) {
	if (!ok)
		throw std::runtime_error(message);
}
void close(Real actual, Real expected, Real relative, char const* message) {
	require(std::isfinite(actual) && std::isfinite(expected), "Nonfinite comparison");
	if (std::abs(actual - expected) > relative * std::max(Real(1), std::abs(expected))) {
		std::cerr << message << ": " << actual << " != " << expected << '\n';
		throw std::runtime_error(message);
	}
}
std::vector<Subgrid> blocks(Config const& c) {
	std::vector<Subgrid> result;
	int const n = 1 << c.level;
	for (int z = 0; z < (c.dimensions > 2 ? n : 1); ++z)
		for (int y = 0; y < (c.dimensions > 1 ? n : 1); ++y)
			for (int x = 0; x < n; ++x)
				result.emplace_back(c, mesh::BlockLocation{c.level, {x, y, z}, c.dimensions});
	return result;
}
Real stable(std::vector<Subgrid> const& blocks) {
	Real dt = std::numeric_limits<Real>::infinity();
	for (auto const& b : blocks)
		dt = std::min(dt, b.stableTimestep());
	return dt;
}
void advance(std::vector<Subgrid>& blocks, Real dt) {
	std::vector<HydroSnapshot> hydro;
	std::vector<RadiationSnapshot> radiation;
	for (auto const& b : blocks) {
		auto const s = b.snapshot();
		if (s.hydroEnabled)
			hydro.push_back({s.location, s.hydro});
		if (s.radiationEnabled)
			radiation.push_back({s.location, s.radiation});
	}
	for (auto& b : blocks)
		b.advance(hydro, radiation, dt);
}
template <class State, class Select>
std::vector<State> flatten(std::vector<Subgrid> const& blocks, Config const& c, Select select) {
	int const n = c.cells * (1 << c.level);
	std::size_t count = 1;
	for (int d = 0; d < c.dimensions; ++d)
		count *= n;
	std::vector<State> result(count);
	for (auto const& b : blocks) {
		auto const s = b.snapshot();
		s.layout.forEachInterior([&](mesh::Coordinates cell, std::size_t index) {
			for (int d = 0; d < c.dimensions; ++d)
				cell[d] += s.location.coordinates[d] * c.cells;
			result[(static_cast<std::size_t>(cell[2]) * n + cell[1]) * n + cell[0]] =
				select(s, index);
		});
	}
	return result;
}
template <class State> State sum(std::vector<State> const& values) {
	State total{};
	for (auto const& value : values)
		total += value;
	return total;
}
void transport() {
	// The same global mesh split into 1 or 2^ndim blocks must give the same
	// solution, including corner halos in the unsplit predictor.
	for (int dimensions = 1; dimensions <= 3; ++dimensions) {
		Config fine =
			parseConfig({"--problem.name=streaming", "--mesh.ndim=" + std::to_string(dimensions),
						 "--mesh.cells=4", "--mesh.level=1", "--runtime.stop_time=0.2",
						 "--output.enabled=off"});
		Config single = fine;
		single.cells = 8;
		single.level = 0;
		auto a = blocks(fine), b = blocks(single);
		auto select = [](Snapshot const& s, std::size_t i) { return s.radiation.values()[i]; };
		using State = radiation::RadiationSystem::State;
		auto initial = flatten<State>(a, fine, select);
		Real time = 0;
		while (time < fine.stopTime) {
			Real dt = std::min({stable(a), stable(b), fine.stopTime - time});
			advance(a, dt);
			advance(b, dt);
			time += dt;
		}
		auto final = flatten<State>(a, fine, select), reference = flatten<State>(b, single, select);
		for (std::size_t i = 0; i < final.size(); ++i) {
			require(radiation::RadiationSystem(physicalLightSpeed).admissible(final[i]),
					"M1 realizability");
			for (int f = 0; f < 4; ++f)
				close(final[i][f], reference[i][f], 3e-12, "Tiled radiation mismatch");
		}
		for (int f = 0; f < 4; ++f)
			close(sum(final)[f], sum(initial)[f], 3e-12, "Periodic radiation conservation");
		std::cout << dimensions
				  << "D radiation: decomposition, conservation, realizability passed\n";
	}
	Config fine =
		parseConfig({"--problem.name=kelvin-helmholtz", "--mesh.cells=8", "--mesh.level=1",
					 "--runtime.stop_time=0.04", "--output.enabled=off"});
	Config single = fine;
	single.cells = 16;
	single.level = 0;
	auto a = blocks(fine), b = blocks(single);
	auto select = [](Snapshot const& s, std::size_t i) { return s.hydro.values()[i]; };
	using State = hydro::ConservedState;
	auto initial = flatten<State>(a, fine, select);
	Real time = 0;
	while (time < fine.stopTime) {
		Real dt = std::min({stable(a), stable(b), fine.stopTime - time});
		advance(a, dt);
		advance(b, dt);
		time += dt;
	}
	auto final = flatten<State>(a, fine, select), reference = flatten<State>(b, single, select);
	for (std::size_t i = 0; i < final.size(); ++i) {
		require(hydro::HydroSystem(fine.gamma).admissible(final[i]), "Hydro positivity");
		for (int f = 0; f < 5; ++f)
			close(final[i][f], reference[i][f], 3e-12, "Tiled hydro mismatch");
	}
	for (int f = 0; f < 5; ++f)
		close(sum(final)[f], sum(initial)[f], 3e-12, "Periodic hydro conservation");
	std::cout << "Hydro: decomposition, conservation, positivity passed\n";
	// A source kick must change momenta and kinetic energy while leaving
	// the gas internal energy exactly unchanged up to rounding.
	Config c = parseConfig({"--problem.name=collapse", "--mesh.level=0", "--output.enabled=off"});
	Subgrid gas(c, {0, {0, 0, 0}, 3});
	std::vector<gravity::State> acceleration(64, gravity::State{0, 2, -3, 4});
	gas.setGravity(acceleration);
	auto before = gas.snapshot();
	gas.kickGravity(0.1);
	auto after = gas.snapshot();
	before.layout.forEachInterior([&](mesh::Coordinates const&, std::size_t i) {
		for (int axis = 0; axis < 3; ++axis)
			close(after.hydro.values()[i].momentum(axis),
				  before.hydro.values()[i].density() * 0.1 * acceleration[0].acceleration(axis),
				  1e-14, "Gravity impulse");
		close(
			hydro::HydroSystem(c.gamma).reconstructionVariables(after.hydro.values()[i]).pressure(),
			hydro::HydroSystem(c.gamma)
				.reconstructionVariables(before.hydro.values()[i])
				.pressure(),
			1e-14, "Kick internal energy");
	});
	std::cout << "Gravity kick internal energy passed\n";
}
void gravityCheck() {
	int constexpr n = 8;
	Real constexpr h = 13;
	std::vector<Real> density(n * n * n);
	for (int z = 0; z < n; ++z)
		for (int y = 0; y < n; ++y)
			for (int x = 0; x < n; ++x)
				density[(z * n + y) * n + x] = 1 + 0.2 * std::sin(0.7 * x + 0.3 * y - 0.9 * z);
	std::vector<gravity::State> reference(density.size());
	for (int z = 0; z < n; ++z)
		for (int y = 0; y < n; ++y)
			for (int x = 0; x < n; ++x) {
				auto& field = reference[(z * n + y) * n + x];
				for (int k = 0; k < n; ++k)
					for (int j = 0; j < n; ++j)
						for (int i = 0; i < n; ++i) {
							if (x == i && y == j && z == k)
								continue;
							std::array<Real, 3> const r{h * (x - i), h * (y - j), h * (z - k)};
							Real const distance = std::hypot(r[0], r[1], r[2]);
							Real const gm = gravity::gravitationalConstant *
											density[(k * n + j) * n + i] * h * h * h;
							field.potential() -= gm / distance;
							for (int d = 0; d < 3; ++d)
								field.acceleration(d) -=
									gm * r[d] / (distance * distance * distance);
						}
			}
	Real previous = std::numeric_limits<Real>::infinity();
	for (int order = 3; order <= 5; ++order) {
		auto const solution = gravity::solve(density, n, h, order, 0.5);
		require(solution.statistics.multipolePairs > 0 && solution.statistics.directPairs > 0,
				"FMM must execute both near and far interactions");
		Real errorPhi = 0, normPhi = 0, errorG = 0, normG = 0;
		for (std::size_t i = 0; i < density.size(); ++i) {
			errorPhi += std::pow(solution.fields[i][0] - reference[i][0], 2);
			normPhi += reference[i][0] * reference[i][0];
			for (int d = 1; d < 4; ++d) {
				errorG += std::pow(solution.fields[i][d] - reference[i][d], 2);
				normG += reference[i][d] * reference[i][d];
			}
		}
		Real const relative = std::sqrt(errorG / normG);
		std::cout << "p=" << order << " relative RMS phi=" << std::sqrt(errorPhi / normPhi)
				  << " g=" << relative << '\n';
		require(relative < previous && relative < 0.02 && std::sqrt(errorPhi / normPhi) < 0.003,
				"FMM direct-reference accuracy/order convergence");
		previous = relative;
		if (order == 5) {
			require(relative < 0.002, "p=5 gravity accuracy");
			auto const scaled = gravity::solve(density, n, 2 * h, order, 0.5);
			for (std::size_t i = 0; i < density.size(); ++i) {
				close(scaled.fields[i][0], 4 * solution.fields[i][0], 1e-13,
					  "Potential length scaling");
				for (int d = 1; d < 4; ++d)
					close(scaled.fields[i][d], 2 * solution.fields[i][d], 1e-13,
						  "Acceleration length scaling");
			}
		}
	}
}
} // namespace
int main(int argc, char** argv) {
	try {
		std::cout << std::scientific << std::setprecision(6);
		if (argc != 2)
			throw std::invalid_argument("Expected transport or gravity");
		if (std::string(argv[1]) == "transport")
			transport();
		else if (std::string(argv[1]) == "gravity")
			gravityCheck();
		else
			throw std::invalid_argument("Unknown check");
		return 0;
	} catch (std::exception const& error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
}
