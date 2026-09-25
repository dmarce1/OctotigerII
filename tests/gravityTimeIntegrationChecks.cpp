#include "testSupport.hpp"
#include "octotigerII/runtime.hpp"
#include "octotigerII/simulation.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>

using namespace octotigerII;

namespace {

// GoogleTest's ScopedTrace stack is OS-thread-local. An HPX task can resume
// on another worker after a Runtime call, so a trace spanning such a call
// would pop a different thread's stack. Log case context without scoped TLS.
template <typename Value>
void describeCase(Value const& value) {
	std::cout << "Gravity integration case: " << value << '\n';
}

Config configuration(std::string const& method, bool species = false) {
	// The whole box lies inside the smooth central part of the star. In
	// particular, a surface/atmosphere discontinuity cannot conceal temporal
	// convergence behind a spatial or positivity-limiter error.
	auto c = parseConfig({"--problem.name=polytrope", "--mesh.cells=4", "--mesh.level=1",
		"--mesh.lower=-2e8", "--mesh.upper=2e8", "--amr.enabled=on", "--amr.minLevel=1",
		"--amr.maxLevel=2", "--amr.refineDensity=0", "--amr.maxCellMass=0",
		"--amr.shadowTolerance=0", "--amr.bufferCells=0", "--output.enabled=off",
		"--verification.analytic=off", "--timestep.refinement=on", "--runtime.stopTime=0"});
	c.gravity.timeIntegration = method;
	c.massFractions.enabled = species;
	if (species) c.massFractions.species = composition::parseSpecies("gas:0.7:He;oxygen:0.3:O;dye:0.4:A=0,Z=0");
	return c;
}

refinement::Criterion refinedOctant(bool& enabled) {
	return [&enabled](refinement::CellView const& cell) {
		bool lower = true;
		for (auto x : cell.center) lower = lower && x < units::Length{};
		return enabled && lower && cell.level < 2 ? Real(2) : Real(0);
	};
}

void expectMixed(Runtime const& runtime) {
	std::set<int> levels;
	for (auto const& block : runtime.snapshots()) levels.insert(block.location.level);
	EXPECT_EQ(levels, (std::set<int>{1, 2}));
}

void expectConservation(Runtime const& runtime, Config const& c, Diagnostics const& initial) {
	auto const now = diagnose(runtime.snapshots(), c);
	auto const boundary = runtime.boundaryTransport();
	auto const mass = now.mass + boundary.outward.mass - boundary.inward.mass;
	auto const energy = now.gasGravityEnergy + boundary.outward.gasEnergy - boundary.inward.gasEnergy
		+ boundary.outward.potentialEnergy - boundary.inward.potentialEnergy;
	EXPECT_NEAR(Real((mass - initial.mass) / initial.mass), 0, 5e-13);
	EXPECT_NEAR(Real((energy - initial.gasGravityEnergy) / initial.gasGravityNorm), 0, 8e-13);
}

void exerciseConservation(std::string const& method, bool species) {
	describeCase(method + (species ? " with species" : " without species"));
	auto const c = configuration(method, species);
	bool refine = true;
	Runtime runtime(c, {refinedOctant(refine)});
	runtime.solveGravity();
	expectMixed(runtime);
	auto const initial = diagnose(runtime.snapshots(), c);
	auto const dt = 0.15 * runtime.stableTimestep();
	for (int step = 0; step < 3; ++step) {
		runtime.advanceGravity(dt);
		expectConservation(runtime, c, initial);
		for (auto const& block : runtime.snapshots()) {
			EXPECT_NEAR(Real(block.time / dt), Real(step + 1), 2e-14);
			if (species) for (std::size_t i = 0; i < block.hydro.values().size(); ++i) {
				auto const rho = block.hydro.values()[i].density();
				EXPECT_NEAR(Real(block.species[0].values()[i] / rho), 0.7, 2e-12);
				EXPECT_NEAR(Real(block.species[1].values()[i] / rho), 0.3, 2e-12);
				EXPECT_NEAR(Real(block.species[2].values()[i] / rho), 0.4, 2e-12);
			}
		}
	}
	auto const counts = runtime.statistics().levelSteps;
	ASSERT_GT(counts.size(), 2u);
	EXPECT_EQ(counts[1], 3u);
	EXPECT_GE(counts[2], 2 * counts[1]);
}

// Volume-weighted relative L1 errors. All momentum components share a norm,
// avoiding a zero denominator when a symmetry sets one component to zero.
std::array<Real, 3> error(std::vector<Snapshot> const& actual, std::vector<Snapshot> const& reference) {
	std::array<long double, 3> numerator{}, denominator{};
	if (actual.size() != reference.size()) throw std::logic_error("Temporal comparison changed the leaf mesh");
	for (std::size_t b = 0; b < actual.size(); ++b) {
		auto const& a = actual[b];
		auto const& r = reference[b];
		if (a.location != r.location || a.hydro.values().size() != r.hydro.values().size())
			throw std::logic_error("Temporal comparison changed cell ordering");
		long double const volume = units::value(r.layout.cellMeasure(r.cellWidth));
		for (std::size_t i = 0; i < r.hydro.values().size(); ++i) {
			auto const& av = a.hydro.values()[i];
			auto const& rv = r.hydro.values()[i];
			numerator[0] += volume * std::abs(units::value(av.density() - rv.density()));
			denominator[0] += volume * std::abs(units::value(rv.density()));
			for (int d = 0; d < ndim; ++d) {
				numerator[1] += volume * std::abs(units::value(av.momentum(d) - rv.momentum(d)));
				denominator[1] += volume * std::abs(units::value(rv.momentum(d)));
			}
			numerator[2] += volume * std::abs(units::value(av.totalEnergy() - rv.totalEnergy()));
			denominator[2] += volume * std::abs(units::value(rv.totalEnergy()));
		}
	}
	std::array<Real, 3> result{};
	for (int i = 0; i < 3; ++i) result[i] = Real(numerator[i] / denominator[i]);
	return result;
}

void temporalConvergence(std::string const& method) {
	auto const c = configuration(method);
	bool refine = true;
	Runtime initial(c, {refinedOctant(refine)});
	initial.solveGravity();
	expectMixed(initial);
	auto const stop = 0.6 * initial.stableTimestep();
	auto integrate = [&](int steps) {
		Runtime runtime(c, {refinedOctant(refine)});
		runtime.solveGravity();
		auto const dt = stop / Real(steps);
		for (int i = 0; i < steps; ++i) runtime.advanceGravity(dt);
		return runtime.snapshots();
	};
	auto const reference = integrate(32);
	std::array<std::array<Real, 3>, 3> errors;
	int count = 0;
	for (int steps : {2, 4, 8}) errors[count++] = error(integrate(steps), reference);
	for (int component = 0; component < 3; ++component) {
		std::cout << method << " temporal L1 " << std::array{"density", "momentum", "gas energy"}[component]
			<< ": " << errors[0][component] << ", " << errors[1][component] << ", " << errors[2][component]
			<< "; orders " << std::log2(errors[0][component] / errors[1][component])
			<< ", " << std::log2(errors[1][component] / errors[2][component]) << '\n';
		EXPECT_GT(errors[0][component], 1e-13) << "Test must resolve temporal error above roundoff";
		EXPECT_LT(errors[1][component], 0.35 * errors[0][component]);
		EXPECT_LT(errors[2][component], 0.35 * errors[1][component]);
	}
}

} // namespace

TEST(GravityTimeIntegration, HierarchicalSubcyclesAndConserves) {
	exerciseConservation("hierarchical", false);
	exerciseConservation("hierarchical", true);
}

TEST(GravityTimeIntegration, ConventionalSubcyclesAndConserves) {
	exerciseConservation("conventional", false);
	exerciseConservation("conventional", true);
}

TEST(GravityTimeIntegration, PeriodicAndReflectingShellWorkConserves) {
	for (auto const* method : {"hierarchical", "conventional"})
		for (auto boundary : {physics::BoundaryCondition::Periodic, physics::BoundaryCondition::Reflecting}) {
			describeCase(std::string(method) + (boundary == physics::BoundaryCondition::Periodic ? " periodic" : " reflecting"));
			auto c = configuration(method, true);
			c.mesh.boundary = physics::BoundaryConditions::uniform(boundary);
			bool refine = true;
			Runtime runtime(c, {refinedOctant(refine)});
			runtime.solveGravity();
			expectMixed(runtime);
			auto const before = diagnose(runtime.snapshots(), c);
			runtime.advanceGravity(0.1 * runtime.stableTimestep());
			expectConservation(runtime, c, before);
			auto const steps = runtime.statistics().levelSteps;
			ASSERT_GT(steps.size(), 2u);
			EXPECT_EQ(steps[1], 1u);
			EXPECT_GE(steps[2], 2u);
			// A reflecting wall may exchange momentum with the gas; these tests
			// check the mass and energy identity for each image kernel instead.
		}
}

TEST(GravityTimeIntegration, HierarchicalForcesBalanceMomentumAfterSubcycling) {
	for (Real openingAngle : {Real(0.5), Real(1e-8)}) {
		describeCase(std::string("hierarchical openingAngle=") + ::testing::PrintToString(openingAngle));
		auto c = configuration("hierarchical");
		c.gravity.openingAngle = openingAngle;
		bool refine = true;
		Runtime runtime(c, {refinedOctant(refine)});
		runtime.solveGravity();
		auto const before = diagnose(runtime.snapshots(), c);
		auto const dt = 0.2 * runtime.stableTimestep();
		for (int step = 0; step < 3; ++step) runtime.advanceGravity(dt);
		auto const after = diagnose(runtime.snapshots(), c);
		auto const boundary = runtime.boundaryTransport();
		for (int d = 0; d < ndim; ++d) {
			auto const residual = after.momentum[d] + boundary.outward.momentum[d] - boundary.inward.momentum[d] - before.momentum[d];
			auto const scale = after.norm.momentum[d] + boundary.outward.momentum[d] + boundary.inward.momentum[d];
			ASSERT_GT(scale, units::Momentum{});
			EXPECT_NEAR(Real(residual / scale), 0, 5e-13) << "axis=" << d;
		}
		expectConservation(runtime, c, before);
	}
}

TEST(GravityTimeIntegration, RegridAfterCoupledSubcyclesClosesAllBudgets) {
	for (auto const* method : {"hierarchical", "conventional"}) {
		describeCase(method);
		auto const c = configuration(method, true);
		bool refine = true;
		Runtime runtime(c, {refinedOctant(refine)});
		runtime.solveGravity();
		auto const before = diagnose(runtime.snapshots(), c);
		runtime.advanceGravity(0.1 * runtime.stableTimestep());
		refine = false;
		ASSERT_TRUE(runtime.regrid({}, true));
		runtime.solveGravity();
		expectConservation(runtime, c, before);
		refine = true;
		ASSERT_TRUE(runtime.regrid({}, true));
		runtime.solveGravity();
		expectMixed(runtime);
		runtime.advanceGravity(0.1 * runtime.stableTimestep());
		expectConservation(runtime, c, before);
	}
}

TEST(GravityTimeIntegration, InvalidIntervalPreservesPublishedState) {
	for (auto const* method : {"hierarchical", "conventional"}) {
		auto const c = configuration(method);
		bool refine = true;
		Runtime runtime(c, {refinedOctant(refine)});
		runtime.solveGravity();
		auto const before = runtime.snapshots();
		auto const generation = runtime.generation();
		for (Real invalid : {Real(0), Real(-1), std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()})
			EXPECT_THROW(runtime.advanceGravity(units::Time::from_value(invalid)), std::invalid_argument);
		EXPECT_EQ(runtime.generation(), generation);
		auto const after = runtime.snapshots();
		ASSERT_EQ(before.size(), after.size());
		for (std::size_t b = 0; b < before.size(); ++b) {
			EXPECT_EQ(before[b].time, after[b].time);
			for (std::size_t i = 0; i < before[b].hydro.values().size(); ++i) {
				before[b].hydro.values()[i].forEach([&](auto f, auto value) {
					EXPECT_EQ(value, after[b].hydro.values()[i].template get<f>());
				});
				before[b].gravity.values()[i].forEach([&](auto f, auto value) {
					EXPECT_EQ(value, after[b].gravity.values()[i].template get<f>());
				});
			}
		}
		// A rejected request must not leave an open kick or work interval.
		runtime.advanceGravity(0.05 * runtime.stableTimestep());
	}
}

TEST(GravityTimeIntegration, DisabledTimeRefinementMatchesGlobalComposition) {
	for (auto const* method : {"hierarchical", "conventional"}) {
		describeCase(method);
		auto c = configuration(method);
		c.timestep.refinement = false;
		bool refine = true;
		Runtime coupled(c, {refinedOctant(refine)}), manual(c, {refinedOctant(refine)});
		coupled.solveGravity(); manual.solveGravity();
		auto const dt = 0.1 * coupled.stableTimestep();
		for (int step = 0; step < 2; ++step) {
			coupled.advanceGravity(dt);
			manual.beginGravityEnergy();
			manual.kickGravity(dt / 2.0);
			manual.advance(dt);
			manual.solveGravity();
			manual.kickGravity(dt / 2.0);
			manual.finishGravityEnergy(dt);
		}
		EXPECT_TRUE(coupled.statistics().levelSteps.empty());
		auto const actual = coupled.snapshots(), expected = manual.snapshots();
		ASSERT_EQ(actual.size(), expected.size());
		for (std::size_t b = 0; b < expected.size(); ++b) {
			EXPECT_EQ(actual[b].time, expected[b].time);
			for (std::size_t i = 0; i < expected[b].hydro.values().size(); ++i) {
				expected[b].hydro.values()[i].forEach([&](auto f, auto value) {
					EXPECT_EQ(value, actual[b].hydro.values()[i].template get<f>());
				});
				expected[b].gravity.values()[i].forEach([&](auto f, auto value) {
					EXPECT_EQ(value, actual[b].gravity.values()[i].template get<f>());
				});
			}
		}
	}
}

TEST(GravityTimeIntegration, NaiveSubcyclesAndRegridEnergyOptionRemainsIndependent) {
	for (auto const* method : {"hierarchical", "conventional"}) {
		describeCase(method);
		auto c = configuration(method);
		c.gravity.energyTreatment = "naive";
		c.gravity.conserveRegridEnergy = false;
		bool refine = true;
		Runtime runtime(c, {refinedOctant(refine)});
		runtime.solveGravity();
		auto const initial = diagnose(runtime.snapshots(), c);
		for (int step = 0; step < 2; ++step) runtime.advanceGravity(0.1 * runtime.stableTimestep());
		auto const beforeRegrid = diagnose(runtime.snapshots(), c);
		EXPECT_GT(beforeRegrid.minimumDensity, units::Density{});
		EXPECT_GT(beforeRegrid.minimumPressure, units::Pressure{});
		auto const boundary = runtime.boundaryTransport();
		EXPECT_NEAR(Real((beforeRegrid.mass + boundary.outward.mass - boundary.inward.mass - initial.mass) / initial.mass), 0, 5e-13);
		auto const steps = runtime.statistics().levelSteps;
		ASSERT_GT(steps.size(), 2u);
		EXPECT_EQ(steps[1], 2u);
		EXPECT_GE(steps[2], 2 * steps[1]);
		refine = false;
		ASSERT_TRUE(runtime.regrid({}, true));
		runtime.solveGravity();
		auto const afterRegrid = diagnose(runtime.snapshots(), c);
		EXPECT_NEAR(Real((afterRegrid.gasEnergy - beforeRegrid.gasEnergy) / beforeRegrid.gasGravityNorm), 0, 2e-14);
		EXPECT_GT(std::abs(Real((afterRegrid.gasGravityEnergy - beforeRegrid.gasGravityEnergy) / beforeRegrid.gasGravityNorm)), 1e-10);
	}
}

TEST(GravityTimeIntegration, NestedThreeLevelShellAndFluxRegistersConserve) {
	for (auto const* method : {"hierarchical", "conventional"}) {
		describeCase(method);
		auto c = configuration(method, true);
		c.amr.maxLevel = 3;
		Runtime runtime(c, {[c](refinement::CellView const& cell) {
			bool corner = true;
			for (auto x : cell.center) corner = corner && x < c.mesh.lower + 0.15 * (c.mesh.upper - c.mesh.lower);
			return corner && cell.level < 3 ? Real(2) : Real(0);
		}});
		runtime.solveGravity();
		std::set<int> levels;
		for (auto const& block : runtime.snapshots()) levels.insert(block.location.level);
		ASSERT_EQ(levels, (std::set<int>{1, 2, 3}));
		auto const before = diagnose(runtime.snapshots(), c);
		auto const dt = 0.1 * runtime.stableTimestep();
		runtime.advanceGravity(dt);
		expectConservation(runtime, c, before);
		auto const counts = runtime.statistics().levelSteps;
		ASSERT_GT(counts.size(), 3u);
		EXPECT_EQ(counts[1], 1u); EXPECT_EQ(counts[2], 2u); EXPECT_EQ(counts[3], 4u);
	}
}

TEST(GravityTimeIntegration, FineCflCanRequireMoreThanTwoSubsteps) {
	for (auto const* method : {"hierarchical", "conventional"}) {
		describeCase(method);
		auto c = configuration(method);
		bool refine = true;
		Runtime runtime(c, {refinedOctant(refine)});
		runtime.solveGravity();
		auto const before = diagnose(runtime.snapshots(), c);
		// The refined octant resolves gas nearer the stellar center, where
		// the sound speed is higher. Its limit is smaller than exactly half
		// the coarse limit, requiring the next dyadic fraction.
		runtime.advanceGravity(runtime.stableTimestep());
		auto const counts = runtime.statistics().levelSteps;
		ASSERT_GT(counts.size(), 2u);
		EXPECT_EQ(counts[1], 1u);
		EXPECT_GE(counts[2], 4u);
		expectConservation(runtime, c, before);
	}
}

TEST(GravityTimeIntegration, HierarchicalSmoothTemporalConvergence) { temporalConvergence("hierarchical"); }
TEST(GravityTimeIntegration, ConventionalSmoothTemporalConvergence) { temporalConvergence("conventional"); }
