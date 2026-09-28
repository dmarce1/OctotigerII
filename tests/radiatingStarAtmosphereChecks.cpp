#include <gtest/gtest.h>
#include "octotigerII/problems/radiatingStarAtmosphere.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>

using namespace octotigerII;
using namespace octotigerII::problems;
namespace {
Real value(auto q) { return units::value(q); }
RadiatingStarAtmosphere const& reference() {
	static RadiatingStarAtmosphere const star({});
	return star;
}
}

TEST(RadiatingStarAtmosphere, RegularCoreCrossesM1CriticalPointIntoTransparentEnvelope) {
	auto const& star=reference();auto const center=star.sample(units::Length{});
	EXPECT_EQ(center.radiationFlux,units::EnergyFlux{});
	EXPECT_EQ(center.gravityMagnitude,units::Acceleration{});
	EXPECT_GT(center.heatingCgs,0);
	EXPECT_GT(star.diagnostics().transitionFluxFactor,2*std::sqrt(Real(3))/5);
	EXPECT_LT(star.diagnostics().transitionFluxFactor,1);
	EXPECT_GT(star.diagnostics().centralOpticalDepth,100);
	EXPECT_GT(star.diagnostics().minimumHeatingCgs,0);
	EXPECT_GT(star.diagnostics().envelopeGamma,1);
	EXPECT_LT(star.diagnostics().envelopeGamma,Real(5)/3);
	EXPECT_GT(star.surfaceRadius(),star.transitionRadius());
	auto before=star.sample(star.transitionRadius()*(1-1e-8));
	auto after=star.sample(star.transitionRadius()*(1+1e-8));
	EXPECT_NEAR(Real(before.density/after.density),1,1e-6);
	EXPECT_NEAR(Real(before.gasPressure/after.gasPressure),1,1e-6);
	EXPECT_NEAR(Real(before.radiationEnergy/after.radiationEnergy),1,1e-6);
	EXPECT_NEAR(Real(before.radiationFlux/after.radiationFlux),1,1e-6);
	EXPECT_EQ(after.opacityCgs,0);
	EXPECT_EQ(after.heatingCgs,0);
	EXPECT_EQ(star.sample(star.surfaceRadius()).density,units::Density{});
}

TEST(RadiatingStarAtmosphere, IndependentDifferencesRecoverGravityAndBothMomentumBalances) {
	auto const& star=reference();Real const scale=value(star.eos().scaleLength());
	Real const cut=value(star.transitionRadius())/scale,surface=value(star.surfaceRadius())/scale;
	for(Real const x:{Real(.2),Real(.8),Real(1.6),cut*.85,cut*.98,(cut+surface)/2,surface*1.2,surface*3}) {
		Real const r=x*scale,h=scale*2e-5;
		auto const q=star.sample(units::Length::from_value(r));
		auto const a=star.sample(units::Length::from_value(r-h)),b=star.sample(units::Length::from_value(r+h));
		Real const rho=value(q.density),g=value(q.gravityMagnitude),force=q.opacityCgs*rho*value(q.radiationFlux)/value(constants::c);
		Real const pressureDerivative=(value(b.gasPressure)-value(a.gasPressure))/(2*h);
		Real const hydroScale=std::max({std::abs(pressureDerivative),rho*g,std::abs(force),Real(1e-100)});
		EXPECT_NEAR((pressureDerivative+rho*g-force)/hydroScale,0,3e-6) << x;
		Real const radDerivative=(value(b.radiationPressureRadial)-value(a.radiationPressureRadial))/(2*h);
		Real const curvature=2*value(q.radiationPressureRadial-q.radiationPressureTangential)/r;
		Real const radScale=std::max({std::abs(radDerivative),std::abs(curvature),std::abs(force),Real(1e-100)});
		EXPECT_NEAR((radDerivative+curvature+force)/radScale,0,3e-6) << x;
		Real const potentialDerivative=value(b.potential-a.potential)/(2*h);
		EXPECT_NEAR(potentialDerivative/g,1,3e-7) << x;
		if(rho>0){
			Real const massDerivative=value(b.enclosedMass-a.enclosedMass)/(2*h);
			EXPECT_NEAR(massDerivative/(4*std::numbers::pi_v<Real>*r*r*rho),1,3e-6) << x;
		}
	}
}

TEST(RadiatingStarAtmosphere, PositiveFrozenHeatingSuppliesEscapingLuminosity) {
	auto const& star=reference();constexpr int count=4000;
	Real const cut=value(star.transitionRadius()),dr=cut/count;long double integral=0;
	for(int i=0;i<=count;++i){
		Real const r=i==count?std::nextafter(cut,Real(0)):i*dr;
		auto const q=star.sample(units::Length::from_value(r));
		EXPECT_GE(q.heatingCgs,0);
		EXPECT_LE(q.fluxFactor,1);
		integral+=(i==0 || i==count?1:i%2?4:2)*static_cast<long double>(r)*r*q.heatingCgs;
	}
	Real const generated=4*std::numbers::pi_v<Real>*dr/3*integral;
	EXPECT_NEAR(generated/star.luminosityCgs(),1,2e-6);
	for(Real fraction:{1.01,2.,10.,100.}){
		auto const radius=star.surfaceRadius()*fraction;auto const q=star.sample(radius);
		Real const escaping=4*std::numbers::pi_v<Real>*value(radius)*value(radius)*value(q.radiationFlux);
		EXPECT_NEAR(escaping/star.luminosityCgs(),1,5e-15);
		EXPECT_EQ(q.heatingCgs,0);EXPECT_EQ(q.opacityCgs,0);
	}
}

TEST(RadiatingStarAtmosphere, RejectsUnphysicalOrUnmatchedParameters) {
	RadiatingStarAtmosphere::Parameters p;p.cutoffDensityFraction=0;
	EXPECT_THROW(RadiatingStarAtmosphere{p},std::invalid_argument);
	p.cutoffDensityFraction=.5;
	EXPECT_THROW(RadiatingStarAtmosphere{p},std::runtime_error);
}

TEST(RadiatingStarAtmosphere, CenterSeriesAndRadialToleranceConvergeIndependently) {
	auto const& tight=reference();RadiatingStarAtmosphere::Parameters p;p.tolerance=1e-9;
	RadiatingStarAtmosphere coarse(p);p.centerSeriesRadius*=.5;RadiatingStarAtmosphere smallerCenter(p);
	EXPECT_NEAR(Real(coarse.surfaceRadius()/tight.surfaceRadius()),1,3e-8);
	EXPECT_NEAR(Real(smallerCenter.surfaceRadius()/coarse.surfaceRadius()),1,3e-8);
	EXPECT_NEAR(Real(coarse.mass()/tight.mass()),1,3e-8);
	EXPECT_NEAR(Real(smallerCenter.mass()/coarse.mass()),1,3e-8);
	EXPECT_NEAR(coarse.diagnostics().transitionFluxFactor,tight.diagnostics().transitionFluxFactor,2e-7);
	EXPECT_NEAR(smallerCenter.diagnostics().transitionFluxFactor,coarse.diagnostics().transitionFluxFactor,2e-7);
}
