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

} // namespace
