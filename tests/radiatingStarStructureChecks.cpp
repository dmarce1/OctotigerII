#include <gtest/gtest.h>
#include "octotigerII/problems/radiatingStarStructure.hpp"
#include "octotigerII/problems/laneEmden.hpp"
#include <cmath>
#include <numbers>

using namespace octotigerII;
using namespace octotigerII::problems;

TEST(RadiatingStarStructure, GasRadiationEquationOfStateAndAdiabaticStability) {
	RadiatingStarEos::Parameters parameters;
	RadiatingStarEos eos(parameters);
	for (Real fraction : {1e-16,1e-10,1e-4,0.1,1.0,1e4,1e8}) {
		auto const rho = parameters.centralDensity*fraction;
		auto const q = eos.atDensity(rho);
		EXPECT_NEAR(Real((q.gasPressure+q.radiationEnergy/3.0)/q.pressure),1,3e-14);
		EXPECT_GT(q.temperature,units::Temperature{});
		EXPECT_GT(q.gamma1,Real(4)/3);
		EXPECT_LT(q.nabla,q.nablaAd);
		EXPECT_LT(q.entropyDerivativeLogDensity,0);
		EXPECT_LT(q.radiationPressureSecondDerivative,0);
		EXPECT_NEAR(Real(eos.atIntegralH(q.integralH).density/rho),1,3e-14);
	}
	auto const center = eos.atDensity(parameters.centralDensity);
	EXPECT_NEAR(center.beta,parameters.centralBeta,3e-15);
	EXPECT_NEAR(center.gamma1,1.511111111111111,3e-15);
	EXPECT_NEAR(center.nabla,0.2361111111111111,3e-15);
	EXPECT_NEAR(center.entropyDifference/units::value(constants::boltzmann/constants::atomicMassUnit),0,3e-14);
}

TEST(RadiatingStarStructure, ThermodynamicDerivativesAgreeWithIndependentDifferences) {
	RadiatingStarEos eos({});
	for (Real fraction : {1e-10,1e-4,0.1,1.0,100.0}) {
		auto const rho = eos.parameters().centralDensity*fraction;
		auto const q = eos.atDensity(rho), low = eos.atDensity(rho*std::exp(-1e-4)), high = eos.atDensity(rho*std::exp(1e-4));
		Real const b = Real((high.radiationEnergy-low.radiationEnergy)/(3.0*(high.pressure-low.pressure)));
		Real const db = (high.radiationPressureDerivative-low.radiationPressureDerivative)/units::value(high.integralH-low.integralH);
		Real const ds = (high.entropyDifference-low.entropyDifference)/2e-4;
		EXPECT_NEAR(b/q.radiationPressureDerivative,1,3e-9);
		EXPECT_NEAR(db/q.radiationPressureSecondDerivative,1,3e-8);
		EXPECT_NEAR(ds/q.entropyDerivativeLogDensity,1,3e-8);
	}
}

TEST(RadiatingStarStructure, IsolatedSphericalPoissonLimitConvergesToLaneEmden) {
	LaneEmden spherical(3.5);
	Real previousMassError = 0, previousRadiusError = 0;
	for (int cells : {64,128,256}) {
		RadiatingStarStructure::Parameters parameters;
		parameters.radialCells = cells; parameters.spinFractionOfSphericalBreakup = 0;
		RadiatingStarStructure star(parameters);
		auto const a = star.eos().scaleLength();
		auto const expectedMass = 4*std::numbers::pi*spherical.surfaceMass()*parameters.eos.centralDensity*a*a*a;
		Real const massError = std::abs(Real(star.mass()/expectedMass)-1);
		Real const radiusError = std::abs(Real(star.equatorialRadius()/(a*spherical.surface()))-1);
		EXPECT_LT(massError,0.01);
		EXPECT_LT(radiusError,0.003);
		if (previousMassError > 0) EXPECT_GT(previousMassError/massError,3.8);
		if (previousRadiusError > 0) EXPECT_GT(previousRadiusError/radiusError,3.7);
		previousMassError = massError; previousRadiusError = radiusError;
		EXPECT_NEAR(Real(star.polarRadius()/star.equatorialRadius()),1,1e-12);
		auto const radius = 3.0*star.equatorialRadius();
		auto const exterior = star.sample(radius,units::Length{});
		EXPECT_NEAR(Real(exterior.potential/(-constants::G*star.mass()/radius)),1,3e-14);
		EXPECT_NEAR(Real(exterior.potentialGradient[0]/(constants::G*star.mass()/(radius*radius))),1,3e-14);
		EXPECT_EQ(exterior.thermodynamics.density,units::Density{});
	}
}

TEST(RadiatingStarStructure, RotationIsSelfConsistentAndVirialErrorConverges) {
	Real previousVirial = 0;
	for (int cells : {64,128,256}) {
		RadiatingStarStructure::Parameters parameters;
		parameters.radialCells = cells;
		RadiatingStarStructure star(parameters);
		EXPECT_LT(Real(star.polarRadius()/star.equatorialRadius()),0.99);
		EXPECT_LT(star.diagnostics().densityResidual,parameters.tolerance);
		EXPECT_LT(star.diagnostics().bernoulliResidual,2e-8);
		Real const virial = std::abs(star.diagnostics().virialResidual);
		EXPECT_LT(virial,0.004);
		if (previousVirial > 0) EXPECT_GT(previousVirial/virial,3.8);
		previousVirial = virial;
		auto const r = 0.5*star.sphericalRadius();
		auto const equator = star.sample(r,units::Length{}), pole = star.sample(units::Length{},r);
		EXPECT_GT(equator.thermodynamics.density,1.005*pole.thermodynamics.density);
		for (Real mu : {0.,0.3,0.7,1.}) {
			auto const q = star.sample(r*std::sqrt(1-mu*mu),r*mu);
			auto const rotationPotential = star.angularVelocity()*star.angularVelocity()*r*r*(1-mu*mu)/2.0;
			EXPECT_NEAR(Real((q.thermodynamics.integralH+q.potential-rotationPotential-star.bernoulliConstant())/star.eos().centralIntegralH()),0,3e-14);
		}
	}
}

TEST(RadiatingStarStructure, InteriorHasPositiveHeatingAndThePhysicalDiffusionFlux) {
	RadiatingStarStructure star({});
	for (Real mu : {0.,0.25,0.5,0.75,1.}) { for (Real fraction : {0.,0.2,0.4,0.6,0.8,0.95}) {
		auto const r = fraction*star.equatorialRadius();
		auto const q = star.sample(r*std::sqrt(1-mu*mu),r*mu);
		if (q.thermodynamics.density == units::Density{}) continue;
		auto const d = star.diffusion(q,0.34);
		EXPECT_GT(d.heating,0);
		EXPECT_GE(d.fluxFactor,0);
	}
	}
	for (Real fraction : {0.2,0.4,0.6}) {
		auto const r = fraction*star.equatorialRadius(), h = 1e-5*star.equatorialRadius();
		auto const low = star.sample(r-h,units::Length{}), high = star.sample(r+h,units::Length{}), q = star.sample(r,units::Length{});
		Real const expected = -units::value(constants::c)*(units::value(high.thermodynamics.radiationEnergy-low.thermodynamics.radiationEnergy)/(2*units::value(h)))
			/(3*0.34*units::value(q.thermodynamics.density));
		EXPECT_NEAR(units::value(star.diffusion(q,0.34).flux[0])/expected,1,2e-7);
	}
}
