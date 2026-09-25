#include "testSupport.hpp"
#include "octotigerII/runtime.hpp"
#include "octotigerII/simulation.hpp"
#include "octotigerII/subgrid/view.hpp"
#include "octotigerII/verification/analytic.hpp"
#include <set>
#include <iostream>
using namespace octotigerII;
namespace {
Config config(std::string const& problem = "sod") {
	auto c = parseConfig({"--problem.name=" + problem, "--mesh.cells=8", "--mesh.level=1", "--output.enabled=off"});
	c.amr.enabled = true; c.amr.maxLevel = 2; c.amr.minLevel = 1;
	c.amr.shadowTolerance = 0; c.amr.bufferCells = 0;
	c.mesh.boundary = physics::BoundaryConditions::periodic();
	return c;
}
refinement::Criterion left(Config const& c, bool& enabled) {
	return [c, &enabled](auto const& cell) {
		return enabled && cell.center[0] < (c.mesh.lower + c.mesh.upper) / 2.0 && cell.level < c.amr.maxLevel ? Real(2) : Real(0);
	};
}
void conserved(Diagnostics const& a, Diagnostics const& b) {
	EXPECT_NEAR(units::value(a.mass), units::value(b.mass), 2e-12 * std::max(Real(1), units::value(a.mass)));
	EXPECT_NEAR(units::value(a.gasEnergy), units::value(b.gasEnergy), 2e-12 * std::max(Real(1), units::value(a.gasEnergy)));
	EXPECT_NEAR(units::value(a.radiationEnergy), units::value(b.radiationEnergy), 2e-12 * std::max(Real(1), units::value(a.radiationEnergy)));
	for (int d = 0; d < ndim; ++d)
		EXPECT_NEAR(units::value(a.momentum[d]), units::value(b.momentum[d]), 2e-12 * std::max(Real(1), units::value(a.mass)));
}
}
TEST(TimeRefinement, DefaultOnAndExplicitOff) {
	EXPECT_TRUE(parseConfig({"--problem.name=sod"}).timestep.refinement);
	EXPECT_FALSE(parseConfig({"--problem.name=sod", "--timestep.refinement=off"}).timestep.refinement);
	EXPECT_THROW(parseConfig({"--problem.name=sod", "--timestep.refinement=perhaps"}), std::invalid_argument);
}
TEST(TimeRefinement, HydroUsesFourOrMoreFineStepsAndConserves) {
	auto c = config(); bool enabled = false;
	Runtime runtime(c, {left(c, enabled)});
	enabled = true; ASSERT_TRUE(runtime.regrid({}, true));
	std::set<int> levels;
	for (auto const& b : runtime.snapshots()) levels.insert(b.location.level);
	ASSERT_EQ(levels, (std::set<int>{1, 2}));
	auto initial = diagnose(runtime.snapshots(), c);
	auto const dt = runtime.stableTimestep();
	runtime.advance(dt);
	auto counts = runtime.statistics().levelSteps;
	ASSERT_EQ(counts.size(), 3u);
	EXPECT_EQ(counts[1], 1u);
	EXPECT_GE(counts[2], 4u); // Fine hot gas needs more than the usual two steps.
	for (auto const& b : runtime.snapshots()) EXPECT_EQ(b.time, dt);
	for (int i = 0; i < 3; ++i) runtime.advance(runtime.stableTimestep());
	conserved(initial, diagnose(runtime.snapshots(), c));
}
TEST(TimeRefinement, DisabledPathUsesFineCflAndNoSubcycles) {
	auto c = config(); bool enabled = false;
	Runtime refined(c, {left(c, enabled)});
	c.timestep.refinement = false;
	Runtime global(c, {left(c, enabled)});
	enabled = true; ASSERT_TRUE(refined.regrid({}, true)); ASSERT_TRUE(global.regrid({}, true));
	EXPECT_GT(refined.stableTimestep(), 2.0 * global.stableTimestep());
	auto initial = diagnose(global.snapshots(), c);
	global.advance(global.stableTimestep());
	EXPECT_TRUE(global.statistics().levelSteps.empty());
	conserved(initial, diagnose(global.snapshots(), c));
}
#if OCTOTIGERII_RADIATION
TEST(TimeRefinement, RadiationSubcyclesAndConserves) {
	auto c = config("streaming"); bool enabled = false;
	Runtime runtime(c, {left(c, enabled)});
	enabled = true; ASSERT_TRUE(runtime.regrid({}, true));
	auto const initial = diagnose(runtime.snapshots(), c);
	auto const dt = 0.8 * runtime.stableTimestep();
	runtime.advance(dt);
	auto counts = runtime.statistics().levelSteps;
	ASSERT_EQ(counts.size(), 3u);
	EXPECT_EQ(counts[1], 1u); EXPECT_EQ(counts[2], 2u);
	for (int i = 0; i < 3; ++i) runtime.advance(0.8 * runtime.stableTimestep());
	conserved(initial, diagnose(runtime.snapshots(), c));
}
#endif
TEST(TimeRefinement, RegridAfterSubcyclingPreservesSynchronizedState) {
	auto c = config(); bool enabled = false;
	Runtime runtime(c, {left(c, enabled)});
	enabled = true; ASSERT_TRUE(runtime.regrid({}, true));
	runtime.advance(0.5 * runtime.stableTimestep());
	auto initial = diagnose(runtime.snapshots(), c);
	enabled = false; ASSERT_TRUE(runtime.regrid({}, true));
	conserved(initial, diagnose(runtime.snapshots(), c));
	runtime.advance(0.3 * runtime.stableTimestep());
	conserved(initial, diagnose(runtime.snapshots(), c));
}
TEST(TimeRefinement, HaloUsesEachDonorBankAndCoarseTimeFraction) {
	std::vector<storage::Locality> localities;
#ifdef OCTOTIGERII_WITH_HPX
	localities = {hpx::find_here()};
#else
	localities = {0};
#endif
	storage::PartitionSet store(localities);
	storage::Layout layout({2, 2}, 1);
	storage::Field<units::Density> field(layout, store, 2, "test.temporalHalo");
	for (unsigned bank = 0; bank < 2; ++bank) for (std::size_t block = 0; block < 2; ++block) {
		auto out = field.handle().output(layout.ranges()[block], bank);
		for (int i = 0; i < 2; ++i) out.data()[i] = units::Density::from_value(100 * block + 10 * bank + i);
		field.handle().commit(layout.ranges()[block], bank, out);
	}
	HaloPlan plan; plan.ghostCount = 4;
	plan.reads.push_back({layout.ranges()[0], {{0, 0, 1}, {1, 1, 1}}, 0});
	plan.reads.push_back({layout.ranges()[1], {{0, 2, 1}, {1, 3, 1}}, 1});
	std::vector<units::Density> ghosts;
	readHalo(field.handle(), plan, 0, ghosts, {{0, 0.25}, {1, 0}});
	ASSERT_EQ(ghosts.size(), 4u);
	EXPECT_DOUBLE_EQ(units::value(ghosts[0]), 2.5);
	EXPECT_DOUBLE_EQ(units::value(ghosts[1]), 3.5);
	EXPECT_DOUBLE_EQ(units::value(ghosts[2]), 110);
	EXPECT_DOUBLE_EQ(units::value(ghosts[3]), 111);
}

#if OCTOTIGERII_RADIATION && OCTOTIGERII_NDIM == 1
TEST(TimeRefinement, StreamingRetainsSecondOrderAcrossLevelBoundary) {
	std::vector<Real> errors;
	for (int n : {8, 16, 32}) {
		auto c = config("streaming"); c.mesh.cells = n;
		bool enabled = false;
		Runtime runtime(c, {left(c, enabled)});
		enabled = true; ASSERT_TRUE(runtime.regrid({}, true));
		auto const stop = 0.35 * (c.mesh.upper - c.mesh.lower) / constants::c;
		auto time = units::Time{};
		while (time < stop) {
			auto dt = std::min(runtime.stableTimestep(), stop - time);
			runtime.advance(dt); time += dt;
		}
		Real error = 0;
		for (auto const& block : runtime.snapshots()) block.layout.forEachInterior([&](auto const& cell, std::size_t i) {
			auto const x = block.layout.cellCenter(block.lower, block.cellWidth, cell);
			auto const exact = verification::streamingState(c, x, stop).radiation.energy();
			error += std::abs(units::value(block.radiation.values()[i].energy() - exact)) * units::value(block.cellWidth);
		});
		errors.push_back(error);
	}
	std::cout << "Streaming L1 errors: " << errors[0] << ", " << errors[1] << ", " << errors[2]
		<< "; final observed order " << std::log2(errors[1] / errors[2]) << "\n";
	EXPECT_LT(errors[1], 0.45 * errors[0]);
	EXPECT_LT(errors[2], 0.4 * errors[1]);
}
#endif
TEST(TimeRefinement, NestedThreeLevelRegistersConserve) {
	auto c = config(); c.mesh.cells = 4; c.mesh.level = 2; c.amr.minLevel = 2; c.amr.maxLevel = 4;
	bool enabled = false;
	Runtime runtime(c, {[&](auto const& cell) {
		return enabled && cell.center[0] < c.mesh.lower + 0.08 * (c.mesh.upper - c.mesh.lower) && cell.level < 4 ? Real(2) : Real(0);
	}});
	enabled = true; ASSERT_TRUE(runtime.regrid({}, true));
	std::set<int> levels;
	for (auto const& b : runtime.snapshots()) levels.insert(b.location.level);
	ASSERT_EQ(levels.size(), 3u);
	auto initial = diagnose(runtime.snapshots(), c);
	auto dt = 0.2 * runtime.stableTimestep();
	runtime.advance(dt);
	conserved(initial, diagnose(runtime.snapshots(), c));
	auto const steps = runtime.statistics().levelSteps;
	EXPECT_EQ(steps.at(2), 1u); EXPECT_EQ(steps.at(3), 2u); EXPECT_EQ(steps.at(4), 4u);
	for (auto const& b : runtime.snapshots()) EXPECT_EQ(b.time, dt);
}
TEST(TimeRefinement, OutflowBoundaryLedgerCoversEverySubstep) {
	auto c = config(); c.mesh.cells = 4;
	c.mesh.boundary = physics::BoundaryConditions::uniform(physics::BoundaryCondition::Outflow);
	bool enabled = false;
	Runtime runtime(c, {left(c, enabled)});
	enabled = true; ASSERT_TRUE(runtime.regrid({}, true));
	auto const before = diagnose(runtime.snapshots(), c);
	for (int i = 0; i < 4; ++i) runtime.advance(runtime.stableTimestep());
	auto after = diagnose(runtime.snapshots(), c);
	auto const boundary = runtime.boundaryTransport();
	after.mass += boundary.outward.mass - boundary.inward.mass;
	after.gasEnergy += boundary.outward.gasEnergy - boundary.inward.gasEnergy;
	for (int d = 0; d < ndim; ++d) after.momentum[d] += boundary.outward.momentum[d] - boundary.inward.momentum[d];
	conserved(before, after);
}
