#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include "octotigerII/gravity/solver.hpp"
#include "octotigerII/runtime.hpp"
#include "runtimeMain.hpp"

using namespace octotigerII;


namespace {
void require(bool ok, char const* message) {
	if (!ok) throw std::runtime_error(message);
}

void close(Real actual, Real expected, Real relative, char const* message) {
	using std::abs;
	using std::isfinite;

	require(isfinite(actual) && isfinite(expected), "Nonfinite comparison");
	if (abs(actual - expected) > relative * std::max(Real(1), abs(expected))) {
		std::cerr << message << ": " << actual << " != " << expected << '\n';
		throw std::runtime_error(message);
	}
}

template <typename U>
void close(boost::units::quantity<U, Real> actual, boost::units::quantity<U, Real> expected, Real relative, char const* message) {
	close(units::value(actual), units::value(expected), relative, message);
}

std::unique_ptr<Runtime> blocks(Config const& c) {
	return std::make_unique<Runtime>(c);
}

units::Time stable(std::unique_ptr<Runtime> const& blocks) {
	return blocks->stableTimestep();
}

void advance(std::unique_ptr<Runtime>& blocks, units::Time dt) {
	blocks->advance(dt);
}

template <typename State, typename Select>
std::vector<State> flatten(std::unique_ptr<Runtime> const& blocks, Config const& c, Select select) {
	int const n = c.mesh.cells * (1 << c.mesh.level);
	std::size_t count = 1;
	for (int d = 0; d < ndim; ++d)
		count *= n;
	std::vector<State> result(count);
	for (auto const& s : blocks->snapshots()) {
		s.layout.forEachInterior([&](mesh::Coordinates cell, std::size_t index) {
			for (int d = 0; d < ndim; ++d)
				cell[d] += s.location.coordinates[d] * c.mesh.cells;
			result[mesh::linearIndex(cell, mesh::filledCoordinates(n))] = select(s, index);
		});
	}
	return result;
}

template <typename State>
State sum(std::vector<State> const& values) {
	State total{};
	for (auto const& value : values)
		total += value;
	return total;
}

void reportScheduling(Runtime const& runtime) {
	auto const statistics = runtime.statistics();
	require(statistics.localTasks + statistics.stolenTasks == 2 * runtime.size() * runtime.generation(), "Transport tasks missing or duplicated");
	std::cout << "  Scheduling: " << statistics.localTasks << " local, " << statistics.stolenTasks << " stolen tasks\n";
}

void streamingProfile() {
	using std::abs;
	using std::exp;
	using std::round;

	// Independent analytic advection: the profile travels at cHat while the
	// stored physical streaming flux remains c E, including at reduced speed.
	for (Real ratio : {Real(1), Real(0.25)}) {
		Config c = parseConfig({"--mesh.cells=32", "--mesh.level=2", "--runtime.stopTime=0.2", "--output.enabled=off"});
		c.radiation.lightSpeedRatio = ratio;
		auto runtime = blocks(c);
		units::Time time{};
		while (time < c.runtime.stopTime) {
			auto const dt = std::min(stable(runtime), c.runtime.stopTime - time);
			advance(runtime, dt);
			time += dt;
		}
		Real error = 0, norm = 0;
		auto const length = c.mesh.upper - c.mesh.lower;
		auto const center = c.mesh.lower + 0.25 * length + ratio * constants::c * time;
		for (auto const& snapshot : runtime->snapshots()) {
			snapshot.layout.forEachInterior([&](mesh::Coordinates const& cell, std::size_t index) {
				auto const position = snapshot.layout.cellCenter(snapshot.lower, snapshot.cellWidth, cell);
				auto distance = position[0] - center;
				distance -= round(distance / length) * length;
				Real const scaled = distance / (0.08 * length);
				Real const expected = 1e-6 + exp(-0.5 * scaled * scaled);
				auto const& state = snapshot.radiation.values()[index];
				error += abs(units::value(state.energy()) - expected);
				norm += expected;
				close(state.radiativeFlux(0), constants::c * state.energy(), 3e-12, "Physical streaming flux must remain c E");
				for (int axis = 1; axis < ndim; ++axis)
					close(state.radiativeFlux(axis), units::EnergyFlux{}, 1e-14, "Spurious transverse radiation flux");
			});
		}
		require(error / norm < 0.003, "Radiation analytic streaming profile");
		std::cout << "Streaming cHat/c=" << ratio << " relative L1=" << error / norm << '\n';
	}
}

template <typename State, typename Select, typename Admissible>
void compareTransport(Select select, Admissible admissible) {
	Config fine = parseConfig({"--mesh.cells=4", "--mesh.level=1", "--output.enabled=off"});
	fine.runtime.workerTasks = 1;
	Config single = fine;
	single.mesh.cells = 8;
	single.mesh.level = 0;
	auto a = blocks(fine), b = blocks(single);
	auto initial = flatten<State>(a, fine, select);
	auto referenceInitial = flatten<State>(b, single, select);
	for (std::size_t i = 0; i < initial.size(); ++i)
		initial[i].forEach([&](auto f, auto q) { close(q, referenceInitial[i].template get<f>(), 3e-12, "Initial decomposition mismatch"); });
	for (int step = 0; step < 3; ++step) {
		auto const dt = std::min(stable(a), stable(b));
		advance(a, dt);
		advance(b, dt);
	}
	auto final = flatten<State>(a, fine, select), reference = flatten<State>(b, single, select);
	for (std::size_t i = 0; i < final.size(); ++i) {
		require(admissible(final[i]), "Transport admissibility");
		final[i].forEach([&](auto f, auto q) { close(q, reference[i].template get<f>(), 3e-12, "Tiled transport mismatch"); });
	}
	if (fine.mesh.periodic) {
		State initialNorm{}, finalNorm{};
		for (auto const& state : initial)
			initialNorm += componentAbs(state);
		for (auto const& state : final)
			finalNorm += componentAbs(state);
		auto const expected = sum(initial);
		// A zero net physical flux can cancel components of order c E. Scale
		// roundoff by their L1 norm, not by the near-zero signed total.
		sum(final).forEach([&](auto f, auto q) {
			auto const scale = std::max(initialNorm.template get<f>(), finalNorm.template get<f>());
			require(units::abs(q - expected.template get<f>()) <= 3e-12 * scale, "Periodic conservation relative to component L1 norm");
		});
	}
	if (std::string(build::problem) == "sod" || std::string(build::problem) == "kelvin-helmholtz") {
		std::size_t const plane = std::string(build::problem) == "sod" ? 8 : 64;
		for (std::size_t i = plane; i < final.size(); ++i)
			final[i].forEach([&](auto f, auto q) { close(q, final[i % plane].template get<f>(), 3e-12, "Extruded problem developed transverse structure"); });
	}
	std::cout << build::problem << ' ' << ndim << "D decomposition and admissibility passed\n";
	reportScheduling(*a);
}

void transport() {
	if (std::string(build::problem) == "streaming" && ndim == 1) streamingProfile();
	if constexpr (build::radiation)
		compareTransport<radiation::RadiationSystem::State>([](Snapshot const& s, std::size_t i) { return s.radiation.values()[i]; },
			[](auto const& state) { return radiation::RadiationSystem(constants::c).admissible(state); });
	if constexpr (build::hydro) {
		auto const c = parseConfig({});
		compareTransport<hydro::ConservedState>([](Snapshot const& s, std::size_t i) { return s.hydro.values()[i]; },
			[&](auto const& state) { return hydro::HydroSystem(c.hydro.gamma).admissible(state); });
	}
	if constexpr (build::gravity && build::hydro) {
		Config c = parseConfig({"--mesh.cells=4", "--mesh.level=0", "--output.enabled=off"});
		Runtime gas(c);
		gravity::State gravity{};
		for (int axis = 0; axis < ndim; ++axis)
			gravity.acceleration(axis) = units::Acceleration::from_value(axis + 2);
		auto const count = gas.snapshots().front().layout.interiorCellCount();
		gas.setGravity({std::vector<gravity::State>(count, gravity)});
		auto before = gas.snapshots().front();
		auto const dt = units::Time::from_value(0.1);
		gas.kickGravity(dt);
		auto after = gas.snapshots().front();
		before.layout.forEachInterior([&](mesh::Coordinates const&, std::size_t i) {
			for (int axis = 0; axis < ndim; ++axis)
				close(after.hydro.values()[i].momentum(axis),
					before.hydro.values()[i].momentum(axis) + before.hydro.values()[i].density() * dt * gravity.acceleration(axis), 1e-14, "Gravity impulse");
			close(hydro::HydroSystem(c.hydro.gamma).reconstructionVariables(after.hydro.values()[i]).pressure(),
				hydro::HydroSystem(c.hydro.gamma).reconstructionVariables(before.hydro.values()[i]).pressure(), 1e-14, "Kick internal energy");
		});
	}
}

#if OCTOTIGERII_GRAVITY
void gravityCheck() {
	using std::sin;

	int constexpr n = 8;
	auto constexpr h = units::Length::from_value(13);
	std::vector<units::Density> density(n * n * n);
	for (int z = 0; z < n; ++z)
		for (int y = 0; y < n; ++y)
			for (int x = 0; x < n; ++x)
				density[(z * n + y) * n + x] = units::Density::from_value(1 + 0.2 * sin(0.7 * x + 0.3 * y - 0.9 * z));
	std::vector<gravity::State> reference(density.size());
	for (int z = 0; z < n; ++z)
		for (int y = 0; y < n; ++y)
			for (int x = 0; x < n; ++x) {
				auto& field = reference[(z * n + y) * n + x];
				for (int k = 0; k < n; ++k)
					for (int j = 0; j < n; ++j)
						for (int i = 0; i < n; ++i) {
							if (x == i && y == j && z == k) continue;
							std::array<units::Length, 3> const r{h * Real(x - i), h * Real(y - j), h * Real(z - k)};
							auto const distance = units::hypot(units::hypot(r[0], r[1]), r[2]);
							auto const gm = constants::G * density[(k * n + j) * n + i] * h * h * h;
							field.potential() -= gm / distance;
							for (int d = 0; d < 3; ++d)
								field.acceleration(d) -= gm * r[d] / (distance * distance * distance);
						}
			}
	Real previous = std::numeric_limits<Real>::infinity();
	for (int order = 1; order <= 10; ++order) {
		auto const solution = gravity::solve(density, n, h, order, 0.5);
		require(solution.statistics.multipolePairs > 0 && solution.statistics.directPairs > 0, "FMM must execute both near and far interactions");
		units::Quantity<4, 0, -4> errorPhi{}, normPhi{};
		units::Quantity<2, 0, -4> errorG{}, normG{};
		for (std::size_t i = 0; i < density.size(); ++i) {
			errorPhi += boost::units::pow<2>(solution.fields[i].potential() - reference[i].potential());
			normPhi += reference[i].potential() * reference[i].potential();
			for (int d = 1; d < 4; ++d) {
				errorG += boost::units::pow<2>(solution.fields[i].acceleration(d - 1) - reference[i].acceleration(d - 1));
				normG += reference[i].acceleration(d - 1) * reference[i].acceleration(d - 1);
			}
		}
		Real const relative = units::sqrt(errorG / normG);
		std::cout << "p=" << order << " relative RMS phi=" << Real(units::sqrt(errorPhi / normPhi)) << " g=" << relative << '\n';
		require(relative < previous && relative < (order >= 3 ? 0.02 : 0.25) && units::sqrt(errorPhi / normPhi) < (order >= 3 ? 0.003 : 0.05),
			"FMM direct-reference accuracy/order convergence");
		previous = relative;
		if (order == 5) {
			require(relative < 0.002, "p=5 gravity accuracy");
			auto const scaled = gravity::solve(density, n, 2.0 * h, order, 0.5);
			for (std::size_t i = 0; i < density.size(); ++i) {
				close(scaled.fields[i].potential(), 4.0 * solution.fields[i].potential(), 1e-13, "Potential length scaling");
				for (int d = 1; d < 4; ++d)
					close(scaled.fields[i].acceleration(d - 1), 2.0 * solution.fields[i].acceleration(d - 1), 1e-13, "Acceleration length scaling");
			}
		}
	}
}

#endif

}	 // namespace


int testMain(int argc, char** argv) {
	try {
		std::cout << std::scientific << std::setprecision(6);
		if (argc != 2) throw std::invalid_argument("Expected transport or gravity");
		if (std::string(argv[1]) == "transport") transport();
#if OCTOTIGERII_GRAVITY
		else if (std::string(argv[1]) == "gravity")
			gravityCheck();
#endif
		else
			throw std::invalid_argument("Unknown check");
		return 0;
	} catch (std::exception const& error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
}

int main(int argc, char** argv) {
	return runtimeMain(argc, argv, testMain);
}
