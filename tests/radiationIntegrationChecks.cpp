#include "testSupport.hpp"
#include "octotigerII/runtime.hpp"
#include "octotigerII/simulation.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <set>

using namespace octotigerII;

namespace {
Config configuration(bool adaptive = false, Real ratio = 1, Real omega = 0) {
	auto c = parseConfig({"--problem.name=radiation-matter", "--mesh.cells=4", "--mesh.level=1", "--output.enabled=off"});
	c.radiation.lightSpeedRatio = ratio;
	c.frame.omega = units::InverseTime::from_value(omega);
	if (omega != 0) c.mesh.boundary = physics::BoundaryConditions{};
	c.amr.enabled = adaptive;
	c.amr.minLevel = 1;
	c.amr.maxLevel = 2;
	c.amr.shadowTolerance = 0;
	c.amr.bufferCells = 0;
	// Fix the test's initial topology; the explicit test driver supplies steps.
	c.runtime.stopTime = {};
	c.validate();
	return c;
}

void budget(Diagnostics const& before, Diagnostics const& after, BoundaryTransport const& boundary, Config const& c, bool momentum = true) {
	Real const weight = 1 / c.radiation.lightSpeedRatio;
	auto const energyFlux = boundary.outward.gasEnergy - boundary.inward.gasEnergy
		+ boundary.outward.potentialEnergy - boundary.inward.potentialEnergy
		+ weight * (boundary.outward.radiationEnergy - boundary.inward.radiationEnergy);
	Real const energyError = Real((after.rslaTotalEnergy + energyFlux - before.rslaTotalEnergy) / before.rslaTotalEnergyNorm);
	EXPECT_NEAR(energyError, 0, 2e-12);
	Real const massError = Real((after.mass + boundary.outward.mass - boundary.inward.mass - before.mass) / before.mass);
	EXPECT_NEAR(massError, 0, 2e-12);
	Real momentumError = 0;
	if (momentum) for (int axis = 0; axis < ndim; ++axis) {
		auto const flux = boundary.outward.momentum[axis] - boundary.inward.momentum[axis]
			+ weight * (boundary.outward.radiationFlux[axis] - boundary.inward.radiationFlux[axis]) / (constants::c * constants::c);
		auto const norm = std::max({before.rslaTotalMomentumNorm[axis], after.rslaTotalMomentumNorm[axis], units::abs(flux)});
		if (norm > units::Momentum{}) {
			Real const error = Real((after.rslaTotalMomentum[axis] + flux - before.rslaTotalMomentum[axis]) / norm);
			momentumError = std::max(momentumError, std::abs(error));
			EXPECT_NEAR(error, 0, 3e-12);
		}
	}
	std::cout << "coupled budget residuals: energy=" << energyError << " mass=" << massError;
	if (momentum) std::cout << " momentum(max)=" << momentumError;
	std::cout << '\n';
}

void transport(bool adaptive, Real ratio, Real omega) {
	auto c = configuration(adaptive, ratio, omega);
	bool refine = adaptive;
	refinement::Criteria criteria{[&](refinement::CellView const& cell) {
		return refine && cell.center[0] < c.mesh.lower + 0.25 * (c.mesh.upper - c.mesh.lower) && cell.level < 2 ? Real(2) : Real(0);
	}};
	Runtime runtime(c, criteria);
	if (adaptive) {
		std::set<int> levels;
		for (auto const& b : runtime.snapshots()) { levels.insert(b.location.level); }
		ASSERT_EQ(levels.size(), 2u);
	}
	auto const before = diagnose(runtime.snapshots(), c);
	for (int i = 0; i < 3; ++i) { runtime.advanceCoupled(0.2 * runtime.stableTimestep()); }
	auto const after = diagnose(runtime.snapshots(), c);
	budget(before, after, runtime.boundaryTransport(), c);
	EXPECT_GT(std::abs(Real((after.radiationEnergy - before.radiationEnergy) / before.radiationEnergy)), 1e-5);
	if (adaptive) {
		refine = false;
		runtime.regrid(units::Time{}, true);
		budget(before, diagnose(runtime.snapshots(), c), runtime.boundaryTransport(), c);
		refine = true;
		runtime.regrid(units::Time{}, true);
		runtime.advanceCoupled(0.1 * runtime.stableTimestep());
		budget(before, diagnose(runtime.snapshots(), c), runtime.boundaryTransport(), c);
	}
}
}

TEST(RadiationIntegration, PhysicalAndReducedSpeedConserveOnPeriodicGrid) {
	for (auto ratio : {Real(1), Real(0.2)}) { transport(false, ratio, 0); }
}

TEST(RadiationIntegration, CoupledTransportConvergesAtSecondOrderInTime) {
	for (bool adaptive : {false, true}) {
		std::cout << "coupled temporal case adaptive=" << adaptive << '\n';
		auto c = configuration(adaptive, 0.2);
		c.mesh.level = adaptive ? 1 : 0;
		c.amr.minLevel = c.mesh.level;
		refinement::Criteria criteria{[&](refinement::CellView const& cell) {
			return adaptive && cell.center[0] < c.mesh.lower + 0.25 * (c.mesh.upper - c.mesh.lower) && cell.level < 2 ? Real(2) : Real(0);
		}};
		Runtime initial(c, criteria);
		auto const interval = 0.4 * initial.stableTimestep();
		auto solve = [&](int steps) {
			Runtime runtime(c, criteria);
			for (int i = 0; i < steps; ++i) { runtime.advanceCoupled(interval / Real(steps)); }
			return runtime.snapshots();
		};
		auto const reference = solve(64);
		std::array<Real, 3> errors{};
		for (int resolution = 0; resolution < 3; ++resolution) {
			auto const result = solve(1 << resolution);
			ASSERT_EQ(result.size(), reference.size());
			Real norm = 0;
			for (std::size_t b = 0; b < result.size(); ++b) {
				for (std::size_t i = 0; i < result[b].radiation.values().size(); ++i) {
					Real const exact = units::value(reference[b].radiation.values()[i].energy());
					errors[resolution] += std::abs(units::value(result[b].radiation.values()[i].energy()) - exact);
					norm += std::abs(exact);
				}
			}
			errors[resolution] /= norm;
		}
		std::cout << "coupled temporal errors adaptive=" << adaptive << ": " << errors[0] << ' ' << errors[1] << ' ' << errors[2] << '\n';
		EXPECT_GT(errors[0] / errors[1], 3.5);
		EXPECT_GT(errors[1] / errors[2], 3.5);
	}
}

TEST(RadiationIntegration, ExternalAccelerationIncludesRadiationMomentumInItsWork) {
	// A wide periodic box makes transport negligible (c*interval/width < 1e-14),
	// while keeping the production initializer and all runtime coupling stages.
	// Each cell can therefore be compared with an independent homogeneous ODE.
	using Number = long double;
	constexpr int momentum = 2, flux = 2 + ndim, work = 2 + 2 * ndim;
	using State = std::array<Number, work + 1>;
	Number const light = units::value(constants::c);
	Number const acceleration = 1e10L;
	for (Real ratio : {Real(1), Real(0.2)}) {
		std::cout << "external-work case ratio=" << ratio << '\n';
		auto c = configuration(false, ratio);
		c.mesh.level = c.amr.minLevel = 0;
		c.mesh.lower = units::Length::from_value(-5e21);
		c.mesh.upper = units::Length::from_value(5e21);
		c.radiation.initialEnergyRatio = 100;
		c.hydro.dualEnergy.pressureThreshold = 0.01;
		c.hydro.acceleration[0] = units::Acceleration::from_value(Real(acceleration));
		c.validate();
		auto const boost = units::Time::from_value(Real(0.005L * light / acceleration));
		auto const interval = units::Time::from_value(Real(0.3L / (ratio * light * 1e-7L * c.radiation.opacity)));
		Runtime initial(c);
		initial.coupleRadiation(units::Time::from_value(1));
		initial.kickGravity(boost);
		auto const initialSnapshots = initial.snapshots();
		ASSERT_EQ(initialSnapshots.size(), 1u);
		auto const& first = initialSnapshots.front();
		hydro::HydroSystem const gas(c.hydro);
		for (auto const& state : first.hydro.values()) {
			auto thermal = state.totalEnergy();
			for (int d = 0; d < ndim; ++d) { thermal -= state.momentum(d) * state.momentum(d) / (2.0 * state.density()); }
			EXPECT_LT(thermal, c.hydro.dualEnergy.pressureThreshold * state.totalEnergy());
		}
		Number const temperatureFactor = (c.hydro.gamma - 1) * c.hydro.meanMolecularWeight
			* units::value(constants::atomicMassUnit / constants::boltzmann);
		auto reference = [&](int steps) {
			Number totalWork = 0;
			for (std::size_t cell = 0; cell < first.hydro.values().size(); ++cell) {
				auto const& g = first.hydro.values()[cell];
				auto const& r = first.radiation.values()[cell];
				Number const rho = units::value(g.density()), chi = rho * c.radiation.opacity;
				State state{};
				state[0] = units::value(g.totalEnergy());
				state[1] = units::value(r.energy());
				Number kinetic = 0;
				for (int d = 0; d < ndim; ++d) {
					state[momentum + d] = units::value(g.momentum(d));
					state[flux + d] = units::value(r.radiativeFlux(d)) / light;
					kinetic += state[momentum + d] * state[momentum + d] / (2 * rho);
				}
				// Preserve the same selected initial heat without borrowing the
				// production source integrator, M1 implementation, or work ledger.
				Number const thermalDefect = units::value(gas.internalEnergy(g)) - (state[0] - kinetic);
				auto rhs = [&](State const& u) {
					Number heat = u[0] + thermalDefect, f2 = 0, betaQ = 0, betaF = 0;
					std::array<Number, ndim> beta{}, f{};
					for (int d = 0; d < ndim; ++d) {
						heat -= u[momentum + d] * u[momentum + d] / (2 * rho);
						beta[d] = u[momentum + d] / (rho * light);
						f[d] = u[flux + d] / u[1];
						f2 += f[d] * f[d];
						betaQ += beta[d] * u[flux + d];
						betaF += beta[d] * f[d];
					}
					Number const temperature = temperatureFactor * heat / rho;
					Number const emission = units::value(constants::radiation) * std::pow(temperature, 4);
					Number const eddington = (3 + 4 * f2) / (5 + 2 * std::sqrt(4 - 3 * f2));
					Number const exchange = chi * light * (u[1] - emission - betaQ);
					State rate{};
					rate[work] = acceleration * u[momentum];
					rate[0] = exchange + rate[work];
					rate[1] = -ratio * exchange;
					for (int d = 0; d < ndim; ++d) {
						Number const pressureBeta = u[1] * ((1 - eddington) * beta[d] / 2
							+ (f2 > 0 ? (3 * eddington - 1) * f[d] * betaF / (2 * f2) : 0));
						Number const force = chi * (u[flux + d] - pressureBeta - emission * beta[d]);
						rate[momentum + d] = force + (d == 0 ? rho * acceleration : 0);
						rate[flux + d] = -ratio * light * force;
					}
					return rate;
				};
				Number const h = units::value(interval) / Number(steps);
				auto add = [](State value, State const& rate, Number dt) {
					for (std::size_t j = 0; j < value.size(); ++j) { value[j] += dt * rate[j]; }
					return value;
				};
				for (int i = 0; i < steps; ++i) {
					auto const k1 = rhs(state), k2 = rhs(add(state, k1, h / 2)),
						k3 = rhs(add(state, k2, h / 2)), k4 = rhs(add(state, k3, h));
					for (int j = 0; j <= work; ++j) { state[j] += h * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6; }
				}
				totalWork += state[work];
			}
			return totalWork;
		};
		Number const exact = reference(1024);
		EXPECT_LT(std::abs((reference(512) - exact) / exact), 1e-13L);
		std::array<Number, 3> errors{};
		for (int resolution = 0; resolution < 3; ++resolution) {
			int const steps = 8 << resolution;
			Runtime runtime(c);
			runtime.coupleRadiation(units::Time::from_value(1));
			runtime.kickGravity(boost);
			for (int i = 0; i < steps; ++i) { runtime.advanceCoupled(interval / Real(steps)); }
			auto const snapshots = runtime.snapshots();
			ASSERT_EQ(snapshots.size(), 1u);
			Number gain = 0;
			for (std::size_t cell = 0; cell < first.hydro.values().size(); ++cell) {
				gain += Number(units::value(snapshots[0].hydro.values()[cell].totalEnergy())) - units::value(first.hydro.values()[cell].totalEnergy());
				gain += (Number(units::value(snapshots[0].radiation.values()[cell].energy())) - units::value(first.radiation.values()[cell].energy())) / ratio;
			}
			errors[resolution] = std::abs((gain - exact) / exact);
		}
		std::cout << "external-work temporal errors ratio=" << ratio << ": " << errors[0] << ' ' << errors[1] << ' ' << errors[2] << '\n';
		EXPECT_GT(errors[0] / errors[1], 3.4L);
		EXPECT_GT(errors[1] / errors[2], 3.4L);
		EXPECT_LT(errors[2], 1e-9L);
	}
}

TEST(RadiationIntegration, MixedLevelSubcyclingAndRegriddingConserve) {
	transport(true, 1, 0);
	transport(true, 0.2, 0);
}

#if OCTOTIGERII_NDIM >= 2
TEST(RadiationIntegration, RotatingGridConservesInertialCombinedBudgets) {
	for (auto omega : {Real(-0.01), Real(0.01)}) { transport(false, 1, omega); }
}
#endif

TEST(RadiationIntegration, ZeroOpacityMatchesUncoupledTransportExactly) {
	auto c = configuration();
	c.radiation.opacity = 0;
	Runtime baseline(c), coupled(c);
	auto const dt = 0.2 * baseline.stableTimestep();
	baseline.advance(dt);
	coupled.advanceCoupled(dt);
	auto const a = baseline.snapshots(), b = coupled.snapshots();
	ASSERT_EQ(a.size(), b.size());
	std::size_t differences = 0;
	for (std::size_t block = 0; block < a.size(); ++block) {
		for (std::size_t i = 0; i < a[block].hydro.values().size(); ++i) {
			a[block].hydro.values()[i].forEach([&](auto field, auto value) { differences += value != b[block].hydro.values()[i].template get<field>(); });
			a[block].radiation.values()[i].forEach([&](auto field, auto value) { differences += value != b[block].radiation.values()[i].template get<field>(); });
		}
	}
	EXPECT_EQ(differences, 0u);
}

TEST(RadiationIntegration, InvalidIntervalLeavesPublishedStateUnchanged) {
	auto c = configuration();
	Runtime runtime(c);
	auto const before = runtime.snapshots();
	auto const generation = runtime.generation();
	EXPECT_THROW(runtime.advanceCoupled(units::Time::from_value(-1)), std::invalid_argument);
	auto const after = runtime.snapshots();
	EXPECT_EQ(generation, runtime.generation());
	ASSERT_EQ(before.size(), after.size());
	for (std::size_t b = 0; b < before.size(); ++b) {
		EXPECT_EQ(before[b].time, after[b].time);
		for (std::size_t i = 0; i < before[b].hydro.values().size(); ++i) {
			before[b].hydro.values()[i].forEach([&](auto field, auto value) { EXPECT_EQ(value, after[b].hydro.values()[i].template get<field>()); });
			before[b].radiation.values()[i].forEach([&](auto field, auto value) { EXPECT_EQ(value, after[b].radiation.values()[i].template get<field>()); });
		}
	}
}

TEST(RadiationIntegration, NumericalFailureRollsBackAndNextStepMatchesCleanRun) {
	auto c = configuration();
	Runtime runtime(c), reference(c);
	auto const step = 0.2 * runtime.stableTimestep();
	auto const before = runtime.snapshots();
	auto const generation = runtime.generation();
	// This interval is positive and finite, so it passes argument validation.
	// Its gross CFL violation forces a numerical-stage failure after the coupled
	// interval has allocated scratch state and begun executing its predictors.
	EXPECT_THROW(runtime.advanceCoupled(1e6 * step), std::runtime_error);
	EXPECT_EQ(runtime.generation(), generation);
	auto expectSame = [](auto const& a, auto const& b) {
		ASSERT_EQ(a.size(), b.size());
		for (std::size_t block = 0; block < a.size(); ++block) {
			EXPECT_EQ(a[block].time, b[block].time);
			for (std::size_t i = 0; i < a[block].hydro.values().size(); ++i) {
				a[block].hydro.values()[i].forEach([&](auto field, auto value) { EXPECT_EQ(value, b[block].hydro.values()[i].template get<field>()); });
				a[block].radiation.values()[i].forEach([&](auto field, auto value) { EXPECT_EQ(value, b[block].radiation.values()[i].template get<field>()); });
			}
		}
	};
	expectSame(before, runtime.snapshots());
	runtime.advanceCoupled(step);
	reference.advanceCoupled(step);
	expectSame(reference.snapshots(), runtime.snapshots());
	budget(diagnose(before, c), diagnose(runtime.snapshots(), c), runtime.boundaryTransport(), c);
}

TEST(RadiationIntegration, AddedEquilibriumRadiationHasNoLocalExchange) {
	auto c = configuration(false, 0.2);
	c.mesh.level = 0;
	c.amr.minLevel = 0;
	c.radiation.initialEnergyRatio = 1;
	Runtime runtime(c);
	auto const before = runtime.snapshots();
	runtime.coupleRadiation(units::Time::from_value(1));
	auto const after = runtime.snapshots();
	ASSERT_EQ(before.size(), after.size());
	for (std::size_t b = 0; b < before.size(); ++b) {
		for (std::size_t i = 0; i < before[b].hydro.values().size(); ++i) {
			test::expectStateNear(before[b].hydro.values()[i], after[b].hydro.values()[i], 3e-13);
			test::expectStateNear(before[b].radiation.values()[i], after[b].radiation.values()[i], 3e-13);
		}
	}
}

#if OCTOTIGERII_GRAVITY
TEST(RadiationIntegration, GasRadiationGravityConserveAcrossAllSchedules) {
	for (auto ratio : {Real(1), Real(0.2)}) { for (auto const* method : {"global", "hierarchical", "conventional", "rotating-hierarchical"}) {
		bool const rotating = std::string(method) == "rotating-hierarchical";
		if (rotating && (ndim < 2 || ratio != 1)) continue;
		std::cout << "coupled gravity case ratio=" << ratio << " method=" << method << '\n';
		auto c = parseConfig({"--problem.name=polytrope", "--mesh.cells=4", "--mesh.level=1", "--amr.enabled=on",
			"--amr.minLevel=1", "--amr.maxLevel=2", "--amr.refineDensity=0", "--amr.shadowTolerance=0", "--amr.bufferCells=0",
			"--radiation.enabled=on", "--radiation.opacity=1e-10", "--radiation.initialEnergyRatio=0.5", "--verification.analytic=off", "--output.enabled=off"});
		c.timestep.refinement = std::string(method) != "global";
		c.gravity.timeIntegration = std::string(method) == "conventional" ? "conventional" : "hierarchical";
		c.radiation.lightSpeedRatio = ratio;
		if (rotating) c.frame.omega = units::InverseTime::from_value(0.01);
		c.validate();
		bool refine = true;
		refinement::Criteria criteria{[&](refinement::CellView const& cell) {
			bool lower = true;
			for (auto x : cell.center) { lower = lower && x < units::Length{}; }
			return refine && lower && cell.level < 2 ? Real(2) : Real(0);
		}};
		Runtime runtime(c, criteria);
		runtime.solveGravity();
		auto const before = diagnose(runtime.snapshots(), c);
		for (int i = 0; i < 2; ++i) { runtime.advanceCoupled(0.1 * runtime.stableTimestep()); }
		budget(before, diagnose(runtime.snapshots(), c), runtime.boundaryTransport(), c, false);
		refine = false;
		runtime.regrid(units::Time{}, true);
		runtime.solveGravity();
		budget(before, diagnose(runtime.snapshots(), c), runtime.boundaryTransport(), c, false);
		refine = true;
		runtime.regrid(units::Time{}, true);
		runtime.solveGravity();
		runtime.advanceCoupled(0.1 * runtime.stableTimestep());
		budget(before, diagnose(runtime.snapshots(), c), runtime.boundaryTransport(), c, false);
	}
	}
}
#endif
