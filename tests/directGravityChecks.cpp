#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "octotigerII/runtime.hpp"
#include "octotigerII/verification/directGravity.hpp"

using namespace octotigerII;


namespace {


void sampling() {
	auto const a = verification::directTargets(1000, 23, 5489);
	EXPECT_TRUE(a == verification::directTargets(1000, 23, 5489)) << "Repeated seed must reproduce targets";
	EXPECT_TRUE(a != verification::directTargets(1000, 23, 5490)) << "Different seeds must change targets";
	EXPECT_TRUE(a.size() == 23 && std::is_sorted(a.begin(), a.end()) && std::adjacent_find(a.begin(), a.end()) == a.end() && a.back() < 1000) << "Target sample must be sorted, distinct and in bounds";
	EXPECT_TRUE(verification::directTargets(4, 99, 0) == std::vector<std::size_t>({0, 1, 2, 3})) << "Full target census";
	EXPECT_TRUE(verification::directTargets(0, 0, 0).empty()) << "Empty population";
}

void uncertainty() {
	using std::abs;
	using std::sqrt;

	verification::Field field{"test", "1", {1, 1, 1, 2}, {1, 1, 1, 1}};
	auto f = verification::directErrorNorm(field, 100);
	EXPECT_TRUE(f.samplingWarning && f.samplingUncertaintyAvailable) << "Uncertain sparse errors must warn";
	EXPECT_TRUE(abs(f.relativeL1HalfWidth95 - 1.96 * sqrt(0.96 * 0.25 / 4)) < 1e-14) << "L1 ratio uncertainty with finite-population correction";
	f = verification::directErrorNorm(field, 4);
	EXPECT_TRUE(!f.samplingWarning && f.relativeL1HalfWidth95 == 0 && f.relativeL2HalfWidth95 == 0) << "Census has no sampling uncertainty";
	field.numerical = {1.1, 2.2, 3.3, 4.4};
	field.exact = {1, 2, 3, 4};
	f = verification::directErrorNorm(field, 100);
	EXPECT_TRUE(!f.samplingWarning && f.relativeL1HalfWidth95 < 1e-14 && f.relativeL2HalfWidth95 < 1e-14) << "Ratio covariance must cancel proportional error";
	field.numerical = {2};
	field.exact = {1};
	f = verification::directErrorNorm(field, 100);
	EXPECT_TRUE(f.samplingWarning && !f.samplingUncertaintyAvailable) << "One target cannot estimate sample variance";
}

void directReference() {
	using std::sqrt;

	auto c = parseConfig({"--mesh.cells=4", "--mesh.level=0", "--runtime.stopTime=0", "--output.enabled=off", "--verification.analytic=on"});
	auto patch = initialSnapshot(c, {});
	// One point mass in a corner: self potential/force are zero. Every other
	// cell, including vacuum targets, has an independently known field.
	auto const h = patch.cellWidth;
	auto const rho = units::Density::from_value(2);
	patch.layout.forEachInterior([&](mesh::Coordinates const& cell, std::size_t i) {
		bool const source = cell == mesh::Coordinates{};
		if (patch.hydroEnabled)
			patch.hydro.values()[i].density() = source ? rho : units::Density{};
		else
			patch.density.values()[i] = source ? rho : units::Density{};
		auto& field = patch.gravity.values()[i];
		field = {};
		if (!source) {
			Real const r = sqrt(Real(cell[0] * cell[0] + cell[1] * cell[1] + cell[2] * cell[2]));
			field.potential() = -constants::G * rho * h * h / r;
			for (int d = 0; d < ndim; ++d)
				field.acceleration(d) = -constants::G * rho * h * Real(cell[d]) / (r * r * r);
		}
	});
	auto full = verification::compare({patch}, c);
	EXPECT_TRUE(full.cells == 64 && full.totalCells == 64 && full.sourceCells == 1 && full.referenceKind == "direct") << "Default full discrete reference";
	for (auto const& f : full.fields)
		EXPECT_TRUE(f.linf / f.referenceLinf < 1e-14) << "Independent single mass solution and self exclusion";
	c.verification.directSamples = 7;
	c.verification.directMaxPairs = 1;
	auto sampled = verification::compare({patch}, c);
	EXPECT_TRUE(sampled.cells == 7) << "Explicit sample count overrides automatic pair budget";
	c.verification.directSamples = 1000;
	EXPECT_TRUE(verification::compare({patch}, c).cells == 64) << "Explicit count clamps to population";
	c.verification.directSamples = 0;
	EXPECT_TRUE(verification::compare({patch}, c).cells == 1) << "Automatic pair budget must sample";
	c.verification.analytic = "off";
	EXPECT_TRUE(verification::compare({patch}, c).status == "disabled") << "Disabled policy";
	c.verification.analytic = "on";
	c.verification.directSamples = 7;
	patch.layout.forEachInterior([&](mesh::Coordinates const&, std::size_t i) {
		if (patch.hydroEnabled)
			patch.hydro.values()[i].density() = {};
		else
			patch.density.values()[i] = {};
		patch.gravity.values()[i] = {};
	});
	auto const vacuum = verification::compare({patch}, c);
	EXPECT_TRUE(vacuum.sourceCells == 0) << "Vacuum skips all sources";
	for (auto const& f : vacuum.fields)
		EXPECT_TRUE(f.linf == 0 && f.referenceL1 == 0 && f.samplingWarning) << "Zero-reference sample uncertainty is unavailable";
	// Sample selection and accumulation must not depend on snapshot order.
	c.mesh.level = 1;
	Runtime runtime(c);
	auto blocks = runtime.snapshots();
	auto const a = verification::compare(blocks, c);
	std::reverse(blocks.begin(), blocks.end());
	auto const b = verification::compare(blocks, c);
	std::ostringstream ja, jb;
	a.writeJson(ja);
	b.writeJson(jb);
	EXPECT_TRUE(ja.str() == jb.str()) << "Block order must not change direct reference or sampling";
}

}	 // namespace


TEST(DirectReference, ReproducibleDistinctSampling) {
	sampling();
}


TEST(DirectReference, SamplingUncertaintyAndCensusLimits) {
	uncertainty();
}


TEST(DirectReference, SingleMassSelfExclusionSamplingPolicyAndBlockOrder) {
	directReference();
}

TEST(DirectReference, ImageFieldsAndContinuumAvailabilityFollowBoundaries) {
	using B = physics::BoundaryConditions;
	using R = physics::BoundaryCondition;
	for (int mode=0;mode<4;++mode) {
		auto c=parseConfig({"--mesh.cells=4","--mesh.level=0","--runtime.stopTime=0","--output.enabled=off"});
		if(mode==0) c.mesh.boundary.lower[0]=c.mesh.boundary.upper[0]=R::Periodic;
		if(mode==1) c.mesh.boundary.lower[2]=R::Reflecting;
		if(mode==2) c.mesh.boundary=B::periodic();
		if(mode==3) c.mesh.boundary=B::uniform(R::Reflecting);
		c.gravity.openingAngle=0.1;
		c.verification.directSamples=13;
		Runtime runtime(c);
		runtime.solveGravity();
		auto comparison=verification::compare(runtime.snapshots(),c);
		EXPECT_EQ(comparison.status,"available");
		EXPECT_EQ(comparison.cells,13u);
		for(auto const& f:comparison.fields) EXPECT_LE(f.linf,3e-13*f.referenceLinf+1e-20);
		c.verification.gravityReference="continuum";
		EXPECT_FALSE(verification::reference(c,units::Time{}).evaluate);
		c.verification.analytic="on";
		EXPECT_THROW(verification::reference(c,units::Time{}),std::invalid_argument);
	}
}
