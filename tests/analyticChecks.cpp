#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include "octotigerII/problems.hpp"
#include "octotigerII/simulation.hpp"
#include "octotigerII/verification/analytic.hpp"
#include "runtimeMain.hpp"

using namespace octotigerII;


namespace {

void require(bool condition, char const* message) {
	if (!condition) throw std::runtime_error(message);
}

void close(Real actual, Real expected, Real tolerance, char const* message) {
	using std::abs;
	using std::isfinite;

	require(isfinite(actual) && abs(actual - expected) <= tolerance * std::max(Real(1), abs(expected)), message);
}

void referenceIdentities() {
	using std::exp;
	using std::sqrt;

	Config c;
	c.mesh.lower = units::Length::from_value(-1);
	c.mesh.upper = units::Length::from_value(1);
	mesh::PhysicalCoordinates x{};
	for (bool gaussian : {false, true}) {
		auto const ref = verification::sphereReference(c, gaussian);
		auto const origin = ref.evaluate(x, {});
		Real const pi = std::numbers::pi_v<Real>, g = units::value(constants::G), rho = 1e4;
		Real const radius = gaussian ? 1 : 0.5, sigma = 0.3;
		Real const centerPhi =
			gaussian ? -4 * pi * g * rho * sigma * sigma * (1 - exp(-radius * radius / (2 * sigma * sigma))) : -2 * pi * g * rho * radius * radius;
		close(units::value(origin.gravity.potential()), centerPhi, 2e-13, "Central analytic potential");
		close(units::value(origin.gravity.acceleration(0)), 0, 0, "Central acceleration must vanish");
		// Independent Simpson integration of shells tests Gaussian mass, potential,
		// and acceleration without using the reference's erf/mass implementation.
		for (Real r : {1e-6, 0.15, 0.5, 1.0, 1.5}) {
			x[0] = units::Length::from_value(r);
			auto integrate = [&](Real lower, Real upper, int power) {
				int constexpr n = 2000;
				Real sum = 0, h = (upper - lower) / n;
				for (int i = 0; i <= n; ++i) {
					Real const s = lower + i * h;
					Real const density = gaussian ? rho * exp(-s * s / (2 * sigma * sigma)) : rho;
					sum += (i == 0 || i == n ? 1 : i % 2 ? 4 : 2) * density * s * (power == 2 ? s : 1);
				}
				return 4 * pi * sum * h / 3;
			};
			Real const mass = integrate(0, std::min(r, radius), 2);
			Real const potential = -g * (mass / r + (r < radius ? integrate(r, radius, 1) : 0));
			auto const state = ref.evaluate(x, {});
			close(units::value(state.gravity.potential()), potential, 2e-11, "Shell-integrated potential");
			close(units::value(state.gravity.acceleration(0)), -g * mass / (r * r), 2e-11, "Shell-integrated acceleration");
		}
		x = {};
	}
	c.mesh.lower = {};
	c.mesh.upper = units::Length::from_value(1);
	c.hydro.gamma = 1.4;
	auto const sod = verification::sodReference(c);
	auto const time = units::Time::from_value(0.2);
	// Independent tabulated star states for standard Sod; both sides of contact.
	for (auto const [position, density] : {std::pair{0.6, 0.4263194281784952}, std::pair{0.75, 0.2655737117053071}}) {
		x[0] = units::Length::from_value(position);
		auto const q = sod.evaluate(x, time).hydro;
		close(units::value(q.density()), density, 3e-12, "Sod star density");
		close(units::value(q.pressure()), 0.3031301780506468, 3e-12, "Sod star pressure");
		close(units::value(q.velocity(0)), 0.9274526200489499, 3e-12, "Sod star velocity");
	}
	// A second gamma checks the shock jump conditions independently.
	c.hydro.gamma = 5.0 / 3;
	auto const monatomic = verification::sodReference(c);
	x[0] = units::Length::from_value(0.75);
	auto const q = monatomic.evaluate(x, time).hydro;
	Real const density = units::value(q.density()), p = units::value(q.pressure()), u = units::value(q.velocity(0));
	Real const shock = density * u / (density - 0.125);
	close(p + density * (shock - u) * (shock - u), 0.1 + 0.125 * shock * shock, 2e-12, "Sod shock momentum jump at gamma=5/3");
	c.mesh.periodic = true;
	c.radiation.lightSpeedRatio = 0.25;
	auto const streaming = verification::streamingReference(c);
	x.fill(units::Length::from_value(0.25));
	auto const travel = (c.mesh.upper - c.mesh.lower) * sqrt(Real(ndim)) / (0.25 * constants::c);
	close(units::value(streaming.evaluate(x, travel).radiation.energy()), 1.000001, 2e-13, "Streaming wraparound at reduced light speed");
}

void normsAndPolicy() {
	auto c = parseConfig({"--mesh.cells=4", "--mesh.level=0", "--runtime.stopTime=0", "--output.enabled=off"});
	c.verification.gravityReference = "continuum";
	auto ref = problemReference(c);
	if (!ref.evaluate) {
		require(verification::compare({initialSnapshot(c, {})}, c).status == "unavailable", "Unsupported problem must report unavailable");
		c.verification.analytic = "on";
		bool rejected = false;
		try {
			verification::reference(c, {});
		} catch (std::invalid_argument const&) {
			rejected = true;
		}
		require(rejected, "Required unavailable reference must be rejected");
		return;
	}
	auto patch = initialSnapshot(c, {});
	patch.layout.forEachInterior([&](mesh::Coordinates const& cell, std::size_t i) {
		auto const q = ref.evaluate(patch.layout.cellCenter(patch.lower, patch.cellWidth, cell), {});
		if constexpr (build::gravity) patch.gravity.values()[i] = q.gravity;
#if OCTOTIGERII_HYDRO
		patch.hydro.values()[i] = hydro::HydroSystem(c.hydro.gamma).conservedState(q.hydro);
#endif
		if constexpr (build::radiation) patch.radiation.values()[i] = q.radiation;
	});
	for (auto const& f : verification::compare({patch}, c).fields)
		close(f.linf, 0, 1e-12, "Exact synthetic field must have zero error");
	patch.layout.forEachInterior([&](mesh::Coordinates const&, std::size_t i) {
		if constexpr (build::gravity) patch.gravity.values()[i].potential() *= 1.1;
		if constexpr (build::hydro) patch.hydro.values()[i].density() *= 1.1;
		if constexpr (build::radiation) patch.radiation.values()[i].energy() *= 1.1;
	});
	auto const comparison = verification::compare({patch}, c);
	auto const& first = comparison.fields.front();
	close(first.l1 / first.referenceL1, 0.1, 2e-12, "Injected 10 percent L1 error");
	close(first.l2 / first.referenceL2, 0.1, 2e-12, "Injected 10 percent L2 error");
	close(first.linf / first.referenceLinf, 0.1, 2e-12, "Injected 10 percent Linf error");
	c.verification.relativeL1Tolerance = 0.01;
	bool rejected = false;
	try {
		comparison.enforce(c);
	} catch (std::runtime_error const&) {
		rejected = true;
	}
	require(rejected, "Accuracy gate must reject injected errors");
	c.verification.relativeL1Tolerance = -1;
	c.verification.analytic = "off";
	require(verification::compare({patch}, c).status == "disabled", "Disabled comparison");
	if (std::string(build::problem) == "sod") {
		c.verification.analytic = "auto";
		require(!verification::reference(c, 2.0 * ref.validUntil).evaluate, "Late Sod reference must be unavailable");
		c.verification.analytic = "on";
		rejected = false;
		try {
			verification::reference(c, 2.0 * ref.validUntil);
		} catch (std::invalid_argument const&) {
			rejected = true;
		}
		require(rejected, "Required late Sod reference must fail");
	}
}

void convergence() {
	std::string const problem = build::problem;
	bool const gravity = problem == "gravity-sphere" || problem == "gravity-gaussian";
	if (!gravity && problem != "sod" && problem != "streaming") return;
	int const coarse = gravity || ndim == 3 ? 16 : 32;
	std::vector<verification::Comparison> results;
	for (int n : {coarse, 2 * coarse}) {
		auto c = parseConfig({"--mesh.cells=" + std::to_string(n), "--mesh.level=0", "--output.enabled=off", "--verification.analytic=on"});
		if (gravity) {
			c.verification.gravityReference = "continuum";
			c.gravity.openingAngle = 0.35;
		} else
			c.runtime.stopTime = units::Time::from_value(problem == "sod" ? 0.15 : 0.2);
		auto const solution = run(c);
		results.push_back(verification::compare(solution.snapshots, c));
		std::cout << "n=" << n << '\n';
		results.back().print(std::cout);
	}
	for (std::size_t f = 0; f < results.front().fields.size(); ++f) {
		auto const& a = results.front().fields[f];
		auto const& b = results.back().fields[f];
		if (b.referenceL1 == 0) {
			require(b.linf < 1e-12, "Zero-reference component drift");
			continue;
		}
		Real const coarseError = a.l1 / a.referenceL1, fineError = b.l1 / b.referenceL1;
		require(fineError < 0.9 * coarseError, "Analytic L1 error must decrease under mesh refinement");
		require(fineError < (gravity ? 0.05 : ndim == 3 ? 0.2 : 0.12), "Analytic reference accuracy");
	}
}

}	 // namespace


int testMain(int argc, char** argv) {
	try {
		std::cout << std::scientific << std::setprecision(6);
		if (argc != 2) throw std::invalid_argument("Expected references or convergence");
		if (std::string(argv[1]) == "references") {
			referenceIdentities();
			normsAndPolicy();
		} else if (std::string(argv[1]) == "convergence")
			convergence();
		else
			throw std::invalid_argument("Unknown analytic check");
		return 0;
	} catch (std::exception const& error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
}

int main(int argc, char** argv) {
	return runtimeMain(argc, argv, testMain);
}
