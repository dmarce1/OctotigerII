#include "testSupport.hpp"
#include <gtest/gtest.h>
#include <fstream>
#include "octotigerII/config.hpp"

using namespace octotigerII;

namespace {

class Options : public ::testing::Test {
protected:

	test::TemporaryDirectory directory;
	std::string first, second;

	void SetUp() override {
		first = (directory.path / "first.ini").string();
		second = (directory.path / "second.ini").string();
		std::ofstream a(first), b(second);
		a << "[gravity]\nmultipoleOrder=3\nopeningAngle=0.4\n[mesh]\ncells=8\n[runtime]\nworkerTasks=2\nworkStealing=off\n";
		b << "gravity.multipole_order=4\ngravity.opening_angle=0.45\n";
		ASSERT_TRUE(a && b);
	}
};


TEST_F(Options, HierarchicalIniSections) {
	auto const c = test::parseConfig({"--config=" + first});
	EXPECT_EQ(c.gravity.multipoleOrder, 3);
	EXPECT_DOUBLE_EQ(c.gravity.openingAngle, 0.4);
	EXPECT_EQ(c.mesh.cells, 8);
	EXPECT_EQ(c.runtime.workerTasks, 2);
	EXPECT_FALSE(c.runtime.workStealing);
}

TEST_F(Options, IndependentGravityEnergyControls) {
	auto const defaults = test::parseConfig({});
	EXPECT_EQ(defaults.gravity.energyTreatment, "mullen");
	EXPECT_TRUE(defaults.gravity.conserveRegridEnergy);
	{ std::ofstream out(first); out << "[gravity]\nenergyTreatment=naive\nconserveRegridEnergy=off\n"; }
	auto const legacy = test::parseConfig({"--config=" + first});
	EXPECT_EQ(legacy.gravity.energyTreatment, "naive");
	EXPECT_FALSE(legacy.gravity.conserveRegridEnergy);
	auto const mixed = test::parseConfig({"--config=" + first, "--gravity.conserveRegridEnergy=on"});
	EXPECT_EQ(mixed.gravity.energyTreatment, "naive");
	EXPECT_TRUE(mixed.gravity.conserveRegridEnergy);
	auto const other = test::parseConfig({"--config=" + first, "--gravity.energyTreatment=mullen"});
	EXPECT_EQ(other.gravity.energyTreatment, "mullen");
	EXPECT_FALSE(other.gravity.conserveRegridEnergy);
	EXPECT_THROW(test::parseConfig({"--gravity.energyTreatment=unknown"}), std::invalid_argument);
}


TEST_F(Options, GravityTimeIntegrationControls) {
	auto const defaults = test::parseConfig({});
	EXPECT_EQ(defaults.gravity.timeIntegration, "hierarchical");
	{ std::ofstream out(first); out << "[gravity]\ntimeIntegration=conventional\nenergyTreatment=naive\nconserveRegridEnergy=off\n"; }
	auto const conventional = test::parseConfig({"--config=" + first});
	EXPECT_EQ(conventional.gravity.timeIntegration, "conventional");
	EXPECT_EQ(conventional.gravity.energyTreatment, "naive");
	EXPECT_FALSE(conventional.gravity.conserveRegridEnergy);
	auto const hierarchical = test::parseConfig({"--gravity.timeIntegration=hierarchical", "--config=" + first});
	EXPECT_EQ(hierarchical.gravity.timeIntegration, "hierarchical");
	EXPECT_EQ(hierarchical.gravity.energyTreatment, "naive");
	EXPECT_FALSE(hierarchical.gravity.conserveRegridEnergy);
	auto const global = test::parseConfig({"--gravity.timeIntegration=conventional", "--timestep.refinement=off"});
	EXPECT_EQ(global.gravity.timeIntegration, "conventional");
	EXPECT_FALSE(global.timestep.refinement);
	EXPECT_THROW(test::parseConfig({"--gravity.timeIntegration=unknown"}), std::invalid_argument);
	EXPECT_NE(helpText().find("gravity.timeIntegration"), std::string::npos);
}


TEST_F(Options, CliOverridesIniRegardlessOfArgumentOrder) {
	auto const c = test::parseConfig({"--gravity.multipoleOrder=5", "--config=" + first, "--gravity.openingAngle=0.5", "--runtime.workStealing=on"});
	EXPECT_EQ(c.gravity.multipoleOrder, 5);
	EXPECT_DOUBLE_EQ(c.gravity.openingAngle, 0.5);
	EXPECT_TRUE(c.runtime.workStealing);
}


TEST_F(Options, LegacyAndCanonicalNamesSharePrecedence) {
	EXPECT_EQ(test::parseConfig({"--config=" + first, "--gravity.multipole_order=4"}).gravity.multipoleOrder, 4);
	auto const c = test::parseConfig({"--config=" + second, "--gravity.openingAngle=0.4"});
	EXPECT_EQ(c.gravity.multipoleOrder, 4);
	EXPECT_DOUBLE_EQ(c.gravity.openingAngle, 0.4);
}


TEST_F(Options, LaterIniOverridesEarlierAndCliOverridesBoth) {
	auto const c = test::parseConfig({"--config=" + first, "--config=" + second, "--gravity.multipoleOrder=5"});
	EXPECT_EQ(c.gravity.multipoleOrder, 5);
	EXPECT_DOUBLE_EQ(c.gravity.openingAngle, 0.45);
}


TEST_F(Options, MissingAndMalformedIniAreRejected) {
	EXPECT_THROW(test::parseConfig({"--config=" + first + ".missing"}), std::exception);
	{ std::ofstream file(first); file << "not.an.option=3\n"; }
	EXPECT_THROW(test::parseConfig({"--config=" + first}), std::exception);
}


TEST(OptionValues, SplitArgumentsSamplingAndGeneratedHelp) {
	EXPECT_EQ(test::parseConfig({"--mesh.cells", "8"}).mesh.cells, 8);
	auto const c = test::parseConfig({"--randomSeed=42", "--verification.directSamples=17", "--verification.directMaxPairs=1234"});
	EXPECT_EQ(c.randomSeed, 42u);
	EXPECT_EQ(c.verification.directSamples, 17u);
	EXPECT_EQ(c.verification.directMaxPairs, 1234u);
	EXPECT_NE(helpText().find("gravity.multipoleOrder"), std::string::npos);
}


class MultipoleOption : public ::testing::TestWithParam<int> {};


TEST_P(MultipoleOption, SupportedOrderIsAccepted) {
	EXPECT_EQ(test::parseConfig({"--gravity.multipoleOrder=" + std::to_string(GetParam())}).gravity.multipoleOrder, GetParam());
}


INSTANTIATE_TEST_SUITE_P(Orders, MultipoleOption, ::testing::Range(1, 11));


class InvalidOption : public ::testing::TestWithParam<std::string> {};


TEST_P(InvalidOption, RejectsInvalidValue) {
	EXPECT_THROW(test::parseConfig({GetParam()}), std::exception) << GetParam();
}


INSTANTIATE_TEST_SUITE_P(Validation, InvalidOption, ::testing::Values(
	"--not.an.option=1", "--problem.name=unknown", "--mesh.ndim=3", "--mesh.cells=3", "--mesh.level=-1",
	"--gravity.multipoleOrder=0", "--gravity.multipoleOrder=11", "--gravity.openingAngle=0", "--gravity.openingAngle=0.6",
	"--gravity.openingAngle=nan", "--gravity.multipoleO=4", "--hydro.gamma=inf", "--hydro.gamma=1",
	"--verification.gravityReference=invalid", "--verification.directSamples=-1", "--verification.directMaxPairs=0",
	"--randomSeed=-1", "--runtime.workStealing=yes", "--runtime.stopTime=-1", "--radiation.lightSpeedRatio=0"));


TEST(OptionValues, DuplicateAliasIsRejected) {
	EXPECT_THROW(test::parseConfig({"--gravity.multipoleOrder=3", "--gravity.multipole_order=3"}), std::exception);
}

TEST(OptionValues, DualEnergyDefaultsAndValidation) {
	auto c = test::parseConfig({});
	EXPECT_TRUE(c.hydro.dualEnergy.enabled);
	EXPECT_EQ(c.hydro.dualEnergy.exponent, 1);
	EXPECT_EQ(c.hydro.dualEnergy.pressureThreshold, 0.001);
	EXPECT_EQ(c.hydro.dualEnergy.syncThreshold, 0.1);
	c = test::parseConfig({"--hydro.dualEnergy.enabled=off", "--hydro.dualEnergy.exponent=-0.5", "--hydro.meanMolecularWeight=0.6"});
	EXPECT_FALSE(c.hydro.dualEnergy.enabled);
	EXPECT_EQ(c.hydro.dualEnergy.exponent, -0.5);
	EXPECT_EQ(c.hydro.meanMolecularWeight, 0.6);
	for (auto const* option : {"--hydro.dualEnergy.exponent=0", "--hydro.dualEnergy.exponent=nan", "--hydro.dualEnergy.pressureThreshold=-1",
		"--hydro.dualEnergy.syncThreshold=1", "--hydro.dualEnergy.syncThreshold=0.001", "--hydro.meanMolecularWeight=0"})
		EXPECT_THROW(test::parseConfig({option}), std::exception);
}

TEST_F(Options, DualEnergyIniAndCliPrecedence) {
	{ std::ofstream f(first); f << "[hydro.dualEnergy]\nenabled=off\nexponent=-2\npressureThreshold=0.002\nsyncThreshold=0.2\n"; }
	auto c = test::parseConfig({"--config=" + first, "--hydro.dualEnergy.enabled=on", "--hydro.dualEnergy.exponent=0.5"});
	EXPECT_TRUE(c.hydro.dualEnergy.enabled);
	EXPECT_EQ(c.hydro.dualEnergy.exponent, 0.5);
	EXPECT_EQ(c.hydro.dualEnergy.pressureThreshold, 0.002);
	EXPECT_EQ(c.hydro.dualEnergy.syncThreshold, 0.2);
}

} // namespace
