#include "testSupport.hpp"
#include "octotigerII/problems/binaryScf.hpp"
#include "octotigerII/simulation.hpp"
using namespace octotigerII;

TEST(BinaryScfIntegration, RotatingBinaryEntersGravityHydroWithClosedMassAndEnergyBudgets) {
	auto c = test::parseConfig({"--scf.cells=16", "--scf.primary.fill=.9", "--scf.donor.fill=1", "--scf.virialTolerance=1",
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
