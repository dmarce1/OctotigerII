#include "testSupport.hpp"
#include "octotigerII/problems/bipolytrope.hpp"
#include "octotigerII/problems/binaryScf.hpp"
#include "octotigerII/subgrid/subgrid.hpp"
#include "octotigerII/gravity/isolatedPotential.hpp"
#include <numbers>
#include <iostream>
using namespace octotigerII;
namespace {
// Deliberately coarse algorithm tests; the production default rejects their
// large virial errors. Resolution acceptance is tested separately.
Config coarseConfig(std::vector<std::string> args) {
	args.push_back("--scf.virialTolerance=1");
	return test::parseConfig(args);
}
}

TEST(BinaryScfGravity, IsolatedConvolutionMatchesAllPairsIncludingEdgesAndSelfOmission) {
	constexpr int n = 8, count = n*n*n;
	Real const dx = .137;
	std::vector<Real> rho(count);
	for (int i = 0; i < count; ++i) rho[i] = Real((i*37+11)%101)/100;
	gravity::IsolatedPotential solve(n,dx);
	auto const phi = solve(rho);
	for (int i = 0; i < count; ++i) {
		Real reference = 0;
		for (int j = 0; j < count; ++j) if (i != j)
			reference -= rho[j]*dx*dx/std::hypot(Real(i%n-j%n), Real(i/n%n-j/n%n), Real(i/(n*n)-j/(n*n)));
		EXPECT_NEAR(phi[i], reference, 3e-14*std::abs(reference));
	}
	std::fill(rho.begin(), rho.end(), 0); rho.front() = 1;
	auto const single = solve(rho);
	EXPECT_NEAR(single.front(), 0, 1e-16);
	EXPECT_NEAR(single.back(), -dx*dx/(std::sqrt(3.)*(n-1)), 1e-16);
}

TEST(Bipolytrope, SinglePolytropeRecoversAnalyticPressureAndEnthalpy) {
	for (Real n : {.5, 1., 1.5, 3., 5.}) {
		problems::Bipolytrope eos(n, n, .23, 1, 7, 13);
		for (Real rho : {.0001, .01, 1., 2., 7.}) {
			Real const h = 13 * std::pow(rho / 7, 1 / n);
			EXPECT_NEAR(eos.enthalpy(rho), h, 1e-12);
			EXPECT_NEAR(eos.pressure(rho), rho * h / (n+1), 1e-12);
			EXPECT_NEAR(eos.density(h), rho, 1e-12);
		}
	}
}

TEST(Bipolytrope, InterfaceMatchesPressureAndEnthalpyWithIndependentIndices) {
	problems::Bipolytrope eos(3, 1.2, .1, 2, 10, 7);
	EXPECT_DOUBLE_EQ(eos.pressure(1), eos.pressure(.5));
	EXPECT_DOUBLE_EQ(eos.enthalpy(1), eos.enthalpy(.5));
	EXPECT_NEAR(eos.density(eos.interfaceEnthalpy()), .5, 1e-14);
	EXPECT_NEAR(eos.density(eos.interfaceEnthalpy() * (1+1e-10)), 1, 1e-8);
	// Check dh/drho = (1/rho) dP/drho independently by finite differences.
	for (Real rho : {.01, .2, 2., 10.}) {
		Real const eps = 1e-5 * rho;
		Real const dh = (eos.enthalpy(rho+eps)-eos.enthalpy(rho-eps))/(2*eps);
		Real const dp = (eos.pressure(rho+eps)-eos.pressure(rho-eps))/(2*eps);
		EXPECT_NEAR(dh, dp/rho, 1e-8*dh);
		EXPECT_NEAR(eos.density(eos.enthalpy(rho)), rho, 1e-12);
	}
}

TEST(Bipolytrope, MixedInterfaceStateIsAFixedPointOnlyAtTheMatchedEnthalpy) {
	problems::Bipolytrope eos(3,1.5,.1,2,10,7);
	Real const h = eos.interfaceEnthalpy();
	EXPECT_NEAR(eos.updateDensity(h,.75), .75, 1e-14);
	EXPECT_GT(eos.updateDensity(h*1.001,.75), .75);
	EXPECT_LT(eos.updateDensity(h*.999,.75), .75);
	EXPECT_NEAR(eos.updateDensity(eos.enthalpy(.2),.2), .2, 1e-14);
	EXPECT_NEAR(eos.updateDensity(eos.enthalpy(2),2), 2, 1e-14);
	EXPECT_DOUBLE_EQ(eos.updateDensity(-1,.75), 0);
}

TEST(BinaryScf, RejectsInvalidPhysicsAndControls) {
	for (auto option : {"--hydro.gamma=1.4", "--hydro.eos=white-dwarf", "--radiation.enabled=on", "--scf.primary.coreIndex=0",
		"--scf.donor.envelopeIndex=-1", "--scf.donor.densityJump=.9", "--scf.donor.fill=1.01", "--scf.cells=24", "--scf.relaxation=0", "--scf.virialTolerance=-1"})
		EXPECT_THROW(test::parseConfig({option}), std::invalid_argument) << option;
}

TEST(BinaryScf, DetachedEqualStarsConvergeAndFillIsABernoulliCondition) {
	auto c = coarseConfig({"--scf.cells=16", "--scf.primary.fill=.8", "--scf.donor.fill=.8", "--scf.tolerance=1e-4", "--gravity.multipoleOrder=3"});
	auto const model = problems::BinaryScf::get(c);
	auto const& d = model->diagnostics();
	EXPECT_LT(d.densityResidual, c.scf.tolerance);
	EXPECT_LT(d.bernoulliResidual, 2*c.scf.tolerance);
	EXPECT_NEAR(d.mass[1]/d.mass[0], 1, 1e-12);
	EXPECT_NEAR(d.center[0], -d.center[1], 1e-8);
	EXPECT_LT(d.bernoulli[1], d.l1Potential);
	EXPECT_NEAR(Real(model->lengthUnit()*d.separation/c.scf.separation), 1, 1e-14);
}

TEST(BinaryScf, IterationLimitRejectsAnUnfinishedModel) {
	EXPECT_THROW(test::parseConfig({"--scf.cells=16", "--scf.maxIterations=1"}), std::runtime_error);
}

TEST(BinaryScf, SemiDetachedDonorAndConservativeHydroHandoff) {
	auto c = coarseConfig({"--scf.cells=32", "--scf.massRatio=.7", "--scf.primary.fill=.9", "--scf.donor.fill=1", "--mesh.cells=8"});
	auto const model = problems::BinaryScf::get(c);
	auto const& d = model->diagnostics();
	EXPECT_LT(d.densityResidual, c.scf.tolerance);
	EXPECT_NEAR(d.bernoulli[1], d.l1Potential, 1e-14);
	EXPECT_LT(d.bernoulli[0], d.l1Potential);
	EXPECT_EQ(c.frame.omega, model->angularVelocity());
	hydro::HydroSystem gas(c);
	std::array<units::Mass, 5> mass{};
	for (int z = 0; z < (1 << c.mesh.level); ++z) for (int y = 0; y < (1 << c.mesh.level); ++y) for (int x = 0; x < (1 << c.mesh.level); ++x) {
		auto const block = initialSnapshot(c, {c.mesh.level, {x,y,z}});
		auto const dv = block.cellWidth*block.cellWidth*block.cellWidth;
		block.layout.forEachInterior([&](auto const& cell, std::size_t i) {
			auto const position = block.layout.cellCenter(block.lower,block.cellWidth,cell);
			auto const primitive = model->state(model->average(position,block.cellWidth), position);
			auto const& q = block.hydro.values()[i];
			units::EnergyDensity kinetic{};
			for (int a = 0; a < 3; ++a) kinetic += q.momentum(a)*q.momentum(a)/(2.0*q.density());
			EXPECT_NEAR(Real((q.energy()-kinetic)/primitive.pressure()), 1/(c.hydro.gamma-1), 3e-10);
			EXPECT_NEAR(Real(gas.reconstructionVariables(q).pressure()/primitive.pressure()), 1, 3e-10);
			if (c.massFractions.enabled) {
				units::Density sum{};
				for (int s = 0; s < 5; ++s) { sum += block.species[s].values()[i]; mass[s] += block.species[s].values()[i]*dv; }
				EXPECT_EQ(sum, q.density());
			}
		});
	}
	if (c.massFractions.enabled) {
		EXPECT_NEAR(Real((mass[0]+mass[1])/c.scf.primaryMass), 1, 1e-11);
		EXPECT_NEAR(Real((mass[2]+mass[3])/c.scf.primaryMass), c.scf.massRatio, 1e-11);
	}
}

TEST(BinaryScf, IndependentCoreAndEnvelopeIndicesProduceTwoBipolytropes) {
	auto c = coarseConfig({"--scf.cells=32", "--scf.primary.coreIndex=3", "--scf.primary.envelopeIndex=1.5",
		"--scf.donor.coreIndex=2.5", "--scf.donor.envelopeIndex=1.2", "--scf.primary.interfaceFraction=.5", "--scf.donor.interfaceFraction=.5",
		"--scf.primary.densityJump=1.2", "--scf.donor.densityJump=1.1", "--scf.primary.fill=.95", "--scf.donor.fill=1"});
	auto const model = problems::BinaryScf::get(c);
	for (int s = 0; s < 2; ++s) {
		EXPECT_GT(model->diagnostics().coreMass[s], .01*model->diagnostics().mass[s]);
		EXPECT_LT(model->diagnostics().coreMass[s], .9*model->diagnostics().mass[s]);
	}
	EXPECT_DOUBLE_EQ(c.hydro.gamma, Real(5)/3);
}

TEST(BinaryScf, RejectsConvergedButUnresolvedEquilibrium) {
	EXPECT_THROW(test::parseConfig({"--scf.cells=16", "--scf.primary.fill=.8", "--scf.donor.fill=.8"}), std::runtime_error);
}

TEST(BinaryScf, PhysicalRescalingAndExplicitGridFrameAreIndependentOfStructure) {
	auto c = coarseConfig({"--scf.cells=16", "--scf.primary.fill=.8", "--scf.donor.fill=.8", "--scf.tolerance=1e-4", "--frame.omega=0"});
	EXPECT_EQ(c.frame.omega, units::InverseTime{});
	auto a = problems::BinaryScf::get(c);
	c.scf.primaryMass *= 4.; c.scf.separation *= 2.;
	auto b = problems::BinaryScf::get(c);
	EXPECT_NEAR(Real(b->densityUnit()/a->densityUnit()), .5, 1e-14);
	EXPECT_NEAR(Real(b->angularVelocity()/a->angularVelocity()), std::sqrt(.5), 1e-14);
	mesh::PhysicalCoordinates x{};
	x[0] = -.25*c.scf.separation;
	auto const qa = (*a)(x);
	for (auto& v : x) v *= 2.;
	auto const qb = (*b)(x);
	EXPECT_NEAR(Real(qb.density()/qa.density()), .5, 1e-14);
	EXPECT_NEAR(Real(qb.pressure()/qa.pressure()), 1, 1e-14);
	EXPECT_NEAR(Real(qb.velocity(1)/qa.velocity(1)), std::sqrt(2.), 1e-14);
}

// Explicit opt-in: this builds a 128^3 reference and is a resolution study,
// rather than an ordinary unit-test gate. Capture stdout to retain the JSON.
TEST(BinaryScf, DISABLED_ResolutionStudy) {
	Real previous = 1;
	for (int n : {32,64,128}) {
		auto c = coarseConfig({"--scf.cells="+std::to_string(n)});
		auto const model = problems::BinaryScf::get(c);
		model->writeJson(std::cout); std::cout << '\n';
		EXPECT_LT(model->diagnostics().virialResidual, .5*previous);
		previous = model->diagnostics().virialResidual;
	}
}

TEST(BinaryScf, SharedPolytropicConstantProducesLargerRocheFillingDonor) {
	auto c = coarseConfig({"--scf.cells=64", "--scf.referenceWidth=2", "--scf.massRatio=.7", "--scf.commonPolytropicK=on"});
	auto const model = problems::BinaryScf::get(c);
	auto const& d = model->diagnostics();
	EXPECT_NEAR(d.polytropicK[0]/d.polytropicK[1],1,1e-13);
	EXPECT_NEAR(d.mass[1]/d.mass[0],.7,1e-12);
	EXPECT_GT(d.volume[1],d.volume[0]);
	EXPECT_LT(d.fill[0],1);
	EXPECT_NEAR(d.bernoulli[1],d.l1Potential,1e-14);
	EXPECT_LT(d.densityResidual,c.scf.tolerance);
	EXPECT_LT(d.maxStarBernoulliResidual,c.scf.tolerance);
	// Verify P/rho^(5/3) throughout both stars, independently of diagnostics.
	for (auto const& cell : model->cells()) if (cell.density > 1e-8)
		EXPECT_NEAR(cell.pressure/std::pow(cell.density,Real(5)/3)/d.polytropicK[0],1,1e-12);
	model->writeJson(std::cout); std::cout << '\n';
	c.scf.coreIndex[0] = 3;
	EXPECT_THROW(problems::BinaryScf::validate(c),std::invalid_argument);
}

TEST(BinaryScf, DISABLED_SharedKResolutionStudy) {
	Real previous = 1;
	for (int n : {32,64,128}) {
		auto c = coarseConfig({"--scf.cells="+std::to_string(n), "--scf.referenceWidth=2", "--scf.massRatio=.7",
			"--scf.primaryMass=1.193082e33", "--scf.separation=3e9", "--scf.commonPolytropicK=on"});
		auto const model = problems::BinaryScf::get(c);
		auto const& d = model->diagnostics();
		EXPECT_GT(d.volume[1],d.volume[0]);
		EXPECT_NEAR(d.polytropicK[0]/d.polytropicK[1],1,1e-13);
		EXPECT_LT(d.virialResidual,.6*previous);
		previous = d.virialResidual;
		model->writeJson(std::cout); std::cout << '\n';
	}
}
