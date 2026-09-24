#include "testSupport.hpp"
#include <gtest/gtest.h>
#include <fstream>
#include "octotigerII/problems.hpp"
#include "octotigerII/problems/laneEmden.hpp"
#include "octotigerII/runtime.hpp"
#include "octotigerII/simulation.hpp"
#include "octotigerII/verification/analytic.hpp"
using namespace octotigerII;

TEST(Polytrope, RadiusCenterAndPressureFollowIniAndCli) {
	test::TemporaryDirectory directory;
	auto const path = directory.path / "star.ini";
	{ std::ofstream out(path); out << "star.radius=2e9\nstar.centralDensity=2e6\nstar.center.x=3e8\nstar.polytropicIndex=1.5\n"; }
	auto const c = test::parseConfig({"--config=" + path.string(), "--star.center.y=-2e8"});
	EXPECT_EQ(c.mesh.upper, units::Length::from_value(4e9));
	EXPECT_EQ(c.mesh.lower, -c.mesh.upper);
	EXPECT_EQ(c.star.center[0], units::Length::from_value(3e8));
	EXPECT_EQ(c.star.center[1], units::Length::from_value(-2e8));
	EXPECT_EQ(c.star.center[2], units::Length{});
	auto const ref = problemReference(c);
	auto const center = ref.evaluate(c.star.center, {});
	EXPECT_EQ(center.hydro.density(), c.star.centralDensity);
	for (int d = 0; d < 3; ++d) EXPECT_EQ(center.gravity.acceleration(d), units::Acceleration{});
	auto x = c.star.center, y = c.star.center;
	x[0] += c.star.radius/2.0;
	y[1] -= c.star.radius/2.0;
	auto const a = ref.evaluate(x, {}), b = ref.evaluate(y, {});
	EXPECT_EQ(a.hydro.density(), b.hydro.density());
	EXPECT_EQ(a.hydro.pressure(), b.hydro.pressure());
	EXPECT_NEAR(units::value(a.gravity.acceleration(0)), -units::value(b.gravity.acceleration(1)), 1e-10);
}

TEST(Polytrope, UnresolvedOffCenterStarTriggersRefinementAndFinalStateIsPhysical) {
	for (auto const* index : {"1.5", "3"}) {
		auto c = test::parseConfig({"--mesh.cells=4", "--mesh.level=0", "--amr.enabled=on", "--amr.maxLevel=4", "--amr.bufferCells=0",
			"--star.radius=2.2e8", std::string("--star.polytropicIndex=") + index,
			"--mesh.lower=-1e9", "--mesh.upper=1e9", "--star.center.x=1.3e8", "--star.center.y=-1.7e8",
			"--runtime.stopTime=0", "--output.enabled=off"});
		auto const coarse = initialSnapshot(c, {0, {}});
		for (auto const& q : coarse.hydro.values())
			EXPECT_EQ(q.density(), c.star.atmosphereFraction*c.star.centralDensity) << "n=" << index;
		Runtime runtime(c);
		auto const reference = problemReference(c);
		units::Density peak{};
		bool finest = false;
		for (auto const& b : runtime.snapshots()) {
			finest = finest || b.location.level == c.amr.maxLevel;
			b.layout.forEachInterior([&](auto const& cell, std::size_t i) {
				auto const expected = reference.evaluate(b.layout.cellCenter(b.lower, b.cellWidth, cell), {});
				EXPECT_EQ(b.hydro.values()[i].density(), expected.hydro.density());
				peak = std::max(peak, b.hydro.values()[i].density());
			});
		}
		EXPECT_TRUE(finest);
		EXPECT_GT(peak, 0.5*c.star.centralDensity);
	}
}

TEST(Polytrope, SmallHydroGravityStepRemainsNearEquilibrium) {
	auto c = test::parseConfig({"--mesh.cells=8", "--mesh.level=1", "--runtime.stopTime=1e-3", "--output.enabled=off", "--verification.analytic=off"});
	auto const result = run(c);
	EXPECT_GT(result.steps, 0);
	EXPECT_NEAR(Real(result.final.mass/result.initial.mass), 1, 1e-7);
	problems::Polytrope const star(c.star.polytropicIndex, c.star.radius, c.star.centralDensity);
	auto const sound = units::sqrt(c.hydro.gamma*star.centralPressure()/c.star.centralDensity);
	units::Velocity maximum{};
	for (auto const& b : result.snapshots)
		for (auto const& q : b.hydro.values())
			if (q.density() > 0.1*c.star.centralDensity)
				for (int d = 0; d < 3; ++d) maximum = std::max(maximum, units::abs(q.momentum(d)/q.density()));
	EXPECT_LT(maximum, 0.01*sound);
}
