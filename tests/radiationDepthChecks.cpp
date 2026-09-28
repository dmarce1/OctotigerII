#include "testSupport.hpp"
#include "octotigerII/runtime.hpp"
#include "octotigerII/simulation.hpp"
#include <iostream>
#include <limits>
#include <set>

using namespace octotigerII;

namespace {
void runThreeLevels(Config c) {
	c.mesh.cells = 4;
	c.mesh.level = 1;
	c.amr.enabled = true;
	c.amr.minLevel = 1;
	c.amr.maxLevel = 3;
	c.amr.refineDensity = {};
	c.amr.maxCellMass = {};
	c.amr.shadowTolerance = 0;
	c.amr.bufferCells = 0;
	c.runtime.stopTime = {};
	c.timestep.refinement = true;
	c.validate();
	Runtime runtime(c, {[c](refinement::CellView const& cell) {
		bool corner = true;
		for (auto x : cell.center) corner = corner && x < c.mesh.lower + 0.15 * (c.mesh.upper - c.mesh.lower);
		return corner && cell.level < 3 ? Real(2) : Real(0);
	}});
	if (c.gravityEnabled()) runtime.solveGravity();
	std::set<int> levels;
	for (auto const& block : runtime.snapshots()) levels.insert(block.location.level);
	ASSERT_EQ(levels, (std::set<int>{1, 2, 3}));
	auto const before = diagnose(runtime.snapshots(), c);
	for (int step = 0; step < 2; ++step) {
		auto const dt = runtime.stableTimestep();
		units::Time finest = units::Time::from_value(std::numeric_limits<Real>::infinity());
		units::Time coarsest = finest;
		for (auto const& block : runtime.snapshots()) {
			auto const limit = std::min(hydro::Solver(hydro::HydroSystem(c.hydro)).stableTimestep(block.hydro, c.timestep.cfl),
				radiation::Solver(radiation::RadiationSystem(c.radiation.lightSpeedRatio * constants::c)).stableTimestep(block.radiation, c.timestep.cfl));
			finest = std::min(finest, limit);
			if (block.location.level == 1) coarsest = std::min(coarsest, limit);
		}
		EXPECT_LE(dt, 2.0 * finest * (1 + 64 * epsilonR));
		if (step == 0) EXPECT_LT(dt, coarsest);
		std::cout << "Three-level coupled " << c.problem << " step=" << step << " dt=" << dt.value() << '\n';
		ASSERT_NO_THROW(runtime.advanceCoupled(dt));
		auto const after = diagnose(runtime.snapshots(), c);
		auto const boundary = runtime.boundaryTransport();
		auto const energyFlux = boundary.outward.gasEnergy - boundary.inward.gasEnergy
			+ boundary.outward.potentialEnergy - boundary.inward.potentialEnergy
			+ (boundary.outward.radiationEnergy - boundary.inward.radiationEnergy) / c.radiation.lightSpeedRatio;
		EXPECT_NEAR(Real((after.rslaTotalEnergy + energyFlux - before.rslaTotalEnergy) / before.rslaTotalEnergyNorm), 0, 4e-12);
		EXPECT_NEAR(Real((after.mass + boundary.outward.mass - boundary.inward.mass - before.mass) / before.mass), 0, 4e-12);
		for (auto const& block : runtime.snapshots()) {
			for (auto const& gas : block.hydro.values()) EXPECT_TRUE(hydro::HydroSystem(c.hydro).admissible(gas));
			for (auto const& rad : block.radiation.values()) EXPECT_TRUE(radiation::RadiationSystem(c.radiation.lightSpeedRatio * constants::c).admissible(rad));
		}
	}
	auto const steps = runtime.statistics().levelSteps;
	ASSERT_GT(steps.size(), 3u);
	EXPECT_EQ(steps[1], 2u);
	EXPECT_GE(steps[2], 4u);
	EXPECT_GE(steps[3], 8u);
}
}

TEST(RadiationDepth, SmoothMatterUsesReportedCflWithThreeLevels) {
	auto c = parseConfig({"--problem.name=radiation-matter", "--output.enabled=off"});
	// Periodicity connects the nested corner to every level-one block and
	// balance then removes level one. Outflow retains the intended depth gap.
	c.mesh.boundary = finiteVolume::BoundaryConditions{};
	runThreeLevels(c);
}

#if OCTOTIGERII_GRAVITY
TEST(RadiationDepth, StellarSurfaceUsesReportedCflWithThreeLevels) {
	auto c = parseConfig({"--problem.name=polytrope", "--radiation.enabled=on", "--radiation.opacity=1e-10",
		"--radiation.initialEnergyRatio=0.5", "--verification.analytic=off", "--output.enabled=off"});
	// Place stellar structure inside the nested corner refinement, so the
	// deepest predictor is exercised by physical gradients rather than vacuum.
	for (auto& center : c.star.center) center = c.mesh.lower + 0.15 * (c.mesh.upper - c.mesh.lower);
	runThreeLevels(c);
}
#endif
