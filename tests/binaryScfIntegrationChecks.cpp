#include "testSupport.hpp"
#include "octotigerII/problems/binaryScf.hpp"
#include "octotigerII/output.hpp"
#include "octotigerII/simulation.hpp"
using namespace octotigerII;

TEST(BinaryScfIntegration, RotatingBinaryEntersGravityHydroWithClosedMassAndEnergyBudgets) {
	auto c = test::parseConfig({"--scf.cells=32", "--scf.referenceWidth=2", "--scf.commonPolytropicK=on", "--scf.massRatio=.7", "--scf.donor.fill=1", "--scf.virialTolerance=1",
		"--mesh.cells=8", "--mesh.level=1", "--gravity.multipoleOrder=3", "--output.enabled=off"});
	c.runtime.stopTime = .001 / problems::BinaryScf::get(c)->angularVelocity();
	auto const result = run(c);
	EXPECT_GT(result.steps, 0);
	EXPECT_EQ(result.final.time, c.runtime.stopTime);
	auto const& b = result.final.boundary;
	EXPECT_NEAR(Real((result.final.mass+b.outward.mass-b.inward.mass)/result.initial.mass), 1, 2e-12);
	auto const energy = result.final.gasGravityEnergy + b.outward.gasEnergy + b.outward.potentialEnergy
		- b.inward.gasEnergy - b.inward.potentialEnergy;
	EXPECT_NEAR(Real((energy-result.initial.gasGravityEnergy)/result.initial.gasGravityNorm), 0, 2e-9);
	EXPECT_TRUE(units::finite(result.final.maximumDensity));
	EXPECT_LT(std::abs(Real(result.final.maximumDensity/result.initial.maximumDensity)-1), .01);
}

TEST(BinaryScfIntegration, CroppedDomainCannotProceedToEvolution) {
	auto c = test::parseConfig({"--scf.cells=16", "--scf.primary.fill=.9", "--scf.virialTolerance=1", "--mesh.cells=4",
		"--mesh.level=1", "--mesh.lower=-1e10", "--mesh.upper=1e10", "--output.enabled=off", "--gravity.multipoleOrder=3"});
	EXPECT_THROW(run(c), std::runtime_error);
}

TEST(BinaryScfIntegration, CoarseHydroMeshCannotBypassTheVirialGate) {
	auto c = test::parseConfig({"--scf.cells=32", "--scf.massRatio=.7", "--scf.primary.fill=.9", "--scf.virialTolerance=.2", "--mesh.cells=4",
		"--mesh.level=1", "--output.enabled=off", "--gravity.multipoleOrder=3"});
	EXPECT_THROW(run(c), std::runtime_error);
}

TEST(BinaryScfIntegration, OrbitalOutputCadenceAdvancesAndWritesAtPhysicalTimes) {
	test::TemporaryDirectory directory;
	auto c = test::parseConfig({"--scf.cells=16", "--scf.primary.fill=.9", "--scf.virialTolerance=1",
		"--scf.evolveOrbits=.0002", "--scf.framesPerOrbit=10000", "--mesh.cells=8", "--mesh.level=1",
		"--output.every=100000", "--gravity.multipoleOrder=3"});
	c.output.directory = directory.path.string();
	auto const cadence = problems::BinaryScf::get(c)->orbitalPeriod() / Real(c.scf.framesPerOrbit);
	EXPECT_NEAR(Real(c.runtime.stopTime / cadence), 2, 1e-12);
	Output output(c);
	auto const result = run(c, [&](auto const& snapshots, int step, auto const& d) { output(snapshots, step, d); });
	EXPECT_EQ(result.final.time, c.runtime.stopTime);
	std::ifstream series(directory.path / "frames.visit");
	ASSERT_TRUE(series);
	std::string frame;
	int count = 0;
	while (std::getline(series, frame)) ++count;
	EXPECT_EQ(count, 3);
	EXPECT_TRUE(std::filesystem::exists(directory.path / "frame_000002.silo"));
	EXPECT_FALSE(std::filesystem::exists(directory.path / "frame_000003.silo"));
}
