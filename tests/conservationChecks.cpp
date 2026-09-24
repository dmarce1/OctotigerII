#include "testSupport.hpp"
#include "octotigerII/output.hpp"
#include "octotigerII/runtime.hpp"
#include <fstream>
#include <set>

using namespace octotigerII;

namespace {
void close(Real a, Real b) {
	EXPECT_NEAR(a, b, 2e-12 * std::max({Real(1), std::abs(a), std::abs(b)}));
}
void budget(Diagnostics const& before, Diagnostics const& after, BoundaryTransport const& boundary) {
	auto check = [&](auto a, auto b, auto in, auto out) { close(units::value(a), units::value(b + out - in)); };
	check(before.mass, after.mass, boundary.inward.mass, boundary.outward.mass);
	check(before.gasEnergy, after.gasEnergy, boundary.inward.gasEnergy, boundary.outward.gasEnergy);
	check(before.radiationEnergy, after.radiationEnergy, boundary.inward.radiationEnergy, boundary.outward.radiationEnergy);
	for (int d = 0; d < ndim; ++d) {
		check(before.momentum[d], after.momentum[d], boundary.inward.momentum[d], boundary.outward.momentum[d]);
		check(before.radiationFlux[d], after.radiationFlux[d], boundary.inward.radiationFlux[d], boundary.outward.radiationFlux[d]);
	}
}
Config configuration(std::string problem) {
	return parseConfig({"--problem.name=" + problem, "--mesh.cells=4", "--mesh.level=1", "--output.enabled=off"});
}
void transport(std::string problem, bool adaptive) {
	auto c = configuration(problem);
	if (problem == "streaming") {
		c.mesh.boundary.lower.fill(physics::BoundaryCondition::Analytic);
		c.mesh.boundary.upper.fill(physics::BoundaryCondition::Analytic);
	}
	c.amr.enabled = adaptive;
	c.amr.maxLevel = 2;
	c.amr.shadowTolerance = 0;
	c.amr.refineDensity = {};
	c.amr.bufferCells = 0;
	refinement::Criteria criteria{[&](refinement::CellView const& cell) {
		return cell.center[0] < c.mesh.lower + 0.2 * (c.mesh.upper - c.mesh.lower) && cell.level < 2 ? Real(2) : Real(0);
	}};
	Runtime runtime(c, criteria);
	if (adaptive) {
		std::set<int> levels;
		for (auto const& block : runtime.snapshots()) levels.insert(block.location.level);
		ASSERT_EQ(levels.size(), 2u);
	}
	auto const before = diagnose(runtime.snapshots(), c);
	for (int i = 0; i < 4; ++i) runtime.advance(0.2 * runtime.stableTimestep());
	auto boundary = runtime.boundaryTransport();
	budget(before, diagnose(runtime.snapshots(), c), boundary);
	if (problem == "streaming") {
		EXPECT_GT(boundary.outward.radiationEnergy, units::Energy{});
		EXPECT_GT(boundary.inward.radiationEnergy, units::Energy{});
	} else EXPECT_GT(boundary.outward.momentum[0] + boundary.inward.momentum[0], units::Momentum{});
	// Failed steps may not publish transport, and regridding may not reset it.
	EXPECT_THROW(runtime.advance(units::Time::from_value(-1)), std::invalid_argument);
	if (adaptive) runtime.regrid(units::Time{}, true);
	budget(before, diagnose(runtime.snapshots(), c), runtime.boundaryTransport());
}
}

#if OCTOTIGERII_HYDRO
TEST(Conservation, HydroBoundaryFluxBalancesUniformAndAdaptiveMeshes) {
	transport("sod", false);
	transport("sod", true);
}
TEST(Conservation, ReflectingPressureTractionBalancesMomentum) {
	auto c = configuration("sod");
	c.mesh.boundary.lower.fill(physics::BoundaryCondition::Reflecting);
	c.mesh.boundary.upper.fill(physics::BoundaryCondition::Reflecting);
	Runtime runtime(c);
	auto const before = diagnose(runtime.snapshots(), c);
	runtime.advance(0.2 * runtime.stableTimestep());
	auto const b = runtime.boundaryTransport();
	budget(before, diagnose(runtime.snapshots(), c), b);
	close(units::value(b.inward.mass + b.outward.mass), 0);
	close(units::value(b.inward.gasEnergy + b.outward.gasEnergy), 0);
}
#endif
#if OCTOTIGERII_RADIATION
TEST(Conservation, RadiationBoundaryFluxBalancesUniformAndAdaptiveMeshes) {
	transport("streaming", false);
	transport("streaming", true);
}
#endif
TEST(Conservation, PeriodicBoundaryTransportIsZero) {
	auto c = test::parseConfig({"--mesh.cells=4", "--mesh.level=1", "--output.enabled=off"});
	if (!c.hydroEnabled() && !c.radiationEnabled()) GTEST_SKIP();
	c.mesh.boundary = physics::BoundaryConditions::periodic();
	Runtime runtime(c);
	auto const before = diagnose(runtime.snapshots(), c);
	runtime.advance(0.2 * runtime.stableTimestep());
	auto const b = runtime.boundaryTransport();
	EXPECT_EQ(b.inward.mass + b.outward.mass, units::Mass{});
	EXPECT_EQ(b.inward.gasEnergy + b.outward.gasEnergy, units::Energy{});
	EXPECT_EQ(b.inward.radiationEnergy + b.outward.radiationEnergy, units::Energy{});
	budget(before, diagnose(runtime.snapshots(), c), b);
}
TEST(Conservation, FileIsWrittenEveryStepWithoutSilo) {
	test::TemporaryDirectory directory;
	auto c = test::parseConfig({"--mesh.cells=4", "--mesh.level=0", "--output.enabled=off"});
	c.output.directory = directory.path.string();
	c.output.every = 100;
	c.runtime.stopTime = units::Time::from_value(1e-14);
	if (!c.hydroEnabled() && !c.radiationEnabled()) c.runtime.stopTime = {};
	Output output(c);
	auto result = run(c, [&](auto const& snapshots, int step, auto const& d) { output(snapshots, step, d); });
	std::ifstream file(directory.path / "conservation.csv");
	ASSERT_TRUE(file);
	std::string line;
	std::getline(file, line);
	EXPECT_NE(line.find("step,time_s"), std::string::npos);
	EXPECT_NE(line.find("_norm"), std::string::npos);
	EXPECT_NE(line.find("_drift_scaled"), std::string::npos);
	int rows = 0;
	while (std::getline(file, line)) ++rows;
	EXPECT_EQ(rows, result.steps + 1);
	EXPECT_FALSE(std::filesystem::exists(directory.path / "frames.visit"));
}
#if OCTOTIGERII_HYDRO && OCTOTIGERII_GRAVITY
TEST(Conservation, SelfGravityEnergyHasHalfRhoPhi) {
	auto c = configuration("collapse");
	c.mesh.level = 0;
	Runtime runtime(c);
	runtime.solveGravity();
	auto blocks = runtime.snapshots();
	units::Energy potential{};
	for (auto const& b : blocks) {
		auto volume = b.layout.cellMeasure(b.cellWidth);
		b.layout.forEachInterior([&](auto const&, std::size_t i) {
			potential += 0.5 * volume * b.hydro.values()[i].density() * b.gravity.values()[i].potential();
		});
	}
	auto d = diagnose(blocks, c);
	EXPECT_LT(potential, units::Energy{});
	close(units::value(d.potentialEnergy), units::value(potential));
	close(units::value(d.gasGravityEnergy), units::value(d.gasEnergy + potential));
	close(units::value(d.gasEnergy), units::value(d.kineticEnergy + d.thermalEnergy));
}
#endif

#if OCTOTIGERII_HYDRO
TEST(Conservation, L1NormSurvivesZeroNetMomentum) {
	auto c = configuration("sod");
	c.mesh.level = 0;
	auto block = initialSnapshot(c, {});
	block.layout.forEachInterior([&](auto const& cell, std::size_t i) {
		auto& u = block.hydro.values()[i];
		u.density() = units::Density::from_value(1);
		u.totalEnergy() = units::EnergyDensity::from_value(3);
		for (int d = 0; d < ndim; ++d) u.momentum(d) = {};
		u.momentum(0) = units::MomentumDensity::from_value(cell[0] % 2 ? -1 : 1);
	});
	auto const d = diagnose({block}, c);
	close(units::value(d.momentum[0]), 0);
	EXPECT_GT(d.norm.momentum[0], units::Momentum{});
	close(units::value(d.norm.momentum[0]),
		units::value(block.layout.cellMeasure(block.cellWidth)) * block.layout.interiorCellCount());
}
#endif
