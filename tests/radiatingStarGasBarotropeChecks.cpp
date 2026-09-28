#include <gtest/gtest.h>
#include "octotigerII/problems/radiatingStarGasBarotrope.hpp"
#include <cmath>
#include <limits>

using namespace octotigerII;
using namespace octotigerII::problems;
namespace {
Real value(auto q) {return units::value(q);}

// Independent integration by parts: integral dp/rho = [p/rho] +
// integral (p/rho) d log rho. This uses pressures, not the table's analytic
// pressure derivatives or its Gauss quadrature.
Real independentIntegral(RadiatingStarGasBarotrope const& gas, units::Density density) {
	int constexpr intervals=4096;
	Real const lower=std::log(Real(gas.cutoffDensity()/gas.parameters().eos.centralDensity));
	Real const upper=std::log(Real(density/gas.parameters().eos.centralDensity));
	Real const step=(upper-lower)/intervals;
	auto thermal=[&](Real x) {
		auto const rho=gas.parameters().eos.centralDensity*std::exp(x);
		return value(gas.eos().atDensity(rho).gasPressure/rho);
	};
	long double integral=0;
	for (int i=0;i<=intervals;++i) { integral+=(i==0 || i==intervals?1:i%2?4:2)*thermal(lower+i*step); }
	return value(gas.cutoffIntegralH())+thermal(upper)-thermal(lower)+Real(integral*step/3);
}
}

TEST(RadiatingStarGasBarotrope, CoreIntegralAgreesWithIndependentPressureQuadrature) {
	for (Real beta:{Real(.2),Real(.8),Real(.99)}) {
		RadiatingStarGasBarotrope::Parameters p;p.eos.centralBeta=beta;
		p.eos.centralDensity=units::Density::from_value(3.7);p.eos.meanMolecularWeight=.71;
		RadiatingStarGasBarotrope const gas(p);
		for (Real ratio:{Real(.031),Real(.37),Real(1),Real(4)}) {
			auto const density=ratio*p.eos.centralDensity;
			EXPECT_NEAR(value(gas.atDensity(density).integralH)/independentIntegral(gas,density),1,3e-12);
		}
	}
	RadiatingStarGasBarotrope const gas({});
	auto const center=gas.atDensity(gas.parameters().eos.centralDensity);
	auto const gasRT=center.gasPressure/center.density;
	// The earlier 6000-point trapezoid prototype gives this approximate value.
	EXPECT_NEAR(Real(center.integralH/gasRT),4.27479908,3e-7);
	EXPECT_GT(std::abs(Real(center.integralH/gasRT)-2.5),1);
	EXPECT_GT(std::abs(Real((center.integralH-gas.eos().centralIntegralH())/gasRT)),1);
}

TEST(RadiatingStarGasBarotrope, TransparentEnvelopeMatchesPressureAndFirstDerivative) {
	RadiatingStarGasBarotrope const gas({});
	auto const cut=gas.atDensity(gas.cutoffDensity());
	auto const core=gas.eos().atDensity(gas.cutoffDensity());
	EXPECT_NEAR(Real(cut.gasPressure/core.gasPressure),1,3e-15);
	EXPECT_NEAR(Real(cut.temperature/core.temperature),1,3e-15);
	EXPECT_NEAR(Real(cut.integralH/(cut.gasPressure/cut.density)),gas.envelopeGamma()/(gas.envelopeGamma()-1),3e-14);
	for (Real ratio:{Real(1e-10),Real(.01),Real(.3)}) {
		auto const q=gas.atDensity(ratio*gas.cutoffDensity());
		EXPECT_NEAR(Real(q.gasPressure/cut.gasPressure)/std::pow(ratio,gas.envelopeGamma()),1,4e-15);
		EXPECT_NEAR(Real(q.temperature/cut.temperature)/std::pow(ratio,gas.envelopeGamma()-1),1,4e-15);
	}
	Real const epsilon=1e-7;
	auto const low=gas.atDensity(gas.cutoffDensity()*std::exp(-epsilon));
	auto const high=gas.atDensity(gas.cutoffDensity()*std::exp(epsilon));
	EXPECT_NEAR(Real(high.densityDerivative/low.densityDerivative),1,2e-7);
	EXPECT_NEAR(Real((high.integralH-low.integralH)*cut.density/(high.gasPressure-low.gasPressure)),1,3e-9);
}

TEST(RadiatingStarGasBarotrope, InverseAndItsDerivativeAreConsistentAndMonotone) {
	RadiatingStarGasBarotrope const gas({});
	units::VelocitySquared previous{};
	for (int i=0;i<=256;++i) {
		Real const ratio=i==256?4.0:std::exp(std::log(1e-12)+(std::log(4.0)-std::log(1e-12))*i/256);
		auto const rho=gas.parameters().eos.centralDensity*ratio;
		auto const q=gas.atDensity(rho),inverse=gas.atIntegralH(q.integralH);
		EXPECT_GT(q.integralH,previous);previous=q.integralH;
		EXPECT_NEAR(Real(inverse.density/rho),1,8e-14);
		EXPECT_NEAR(Real(inverse.densityDerivative/q.densityDerivative),1,8e-14);
		EXPECT_EQ(inverse.integralH,q.integralH);
	}
	for (Real ratio:{Real(1e-8),Real(.002),Real(.016),Real(.07),Real(.2),Real(1),Real(3.9)}) {
		auto const q=gas.atDensity(ratio*gas.parameters().eos.centralDensity);
		auto const delta=q.integralH*1e-4;
		auto const a=gas.atIntegralH(q.integralH-2.0*delta),b=gas.atIntegralH(q.integralH-delta);
		auto const c=gas.atIntegralH(q.integralH+delta),d=gas.atIntegralH(q.integralH+2.0*delta);
		auto const densityDerivative=(a.density-8.0*b.density+8.0*c.density-d.density)/(12.0*delta);
		auto const pressureDerivative=(a.gasPressure-8.0*b.gasPressure+8.0*c.gasPressure-d.gasPressure)/(12.0*delta);
		EXPECT_NEAR(Real(densityDerivative/q.densityDerivative),1,2e-10) << ratio;
		EXPECT_NEAR(Real(pressureDerivative/q.density),1,2e-10) << ratio;
	}
}

TEST(RadiatingStarGasBarotrope, VacuumAndFiniteTableLimitsAreExplicit) {
	RadiatingStarGasBarotrope const gas({});
	EXPECT_EQ(gas.atIntegralH({}).density,units::Density{});
	EXPECT_EQ(gas.atIntegralH(-gas.centralIntegralH()).density,units::Density{});
	EXPECT_EQ(gas.atDensity({}).densityDerivative,(units::Quantity<-5,1,2>{}));
	EXPECT_NEAR(Real(gas.atIntegralH(gas.maximumIntegralH()).density/gas.maximumDensity()),1,3e-15);
	EXPECT_THROW(gas.atIntegralH(1.000001*gas.maximumIntegralH()),std::out_of_range);
	EXPECT_THROW(gas.atDensity(1.000001*gas.maximumDensity()),std::out_of_range);
	EXPECT_THROW(gas.atDensity(-gas.cutoffDensity()),std::invalid_argument);
	EXPECT_THROW(gas.atIntegralH(units::VelocitySquared::from_value(std::numeric_limits<Real>::infinity())),std::invalid_argument);
	EXPECT_THROW(gas.atDensity(units::Density::from_value(std::numeric_limits<Real>::quiet_NaN())),std::invalid_argument);
	for (Real cutoff:{Real(0),Real(-.1),Real(1),std::numeric_limits<Real>::quiet_NaN()}) {
		RadiatingStarGasBarotrope::Parameters p;p.cutoffDensityFraction=cutoff;
		EXPECT_THROW((RadiatingStarGasBarotrope{p}),std::invalid_argument);
	}
}
