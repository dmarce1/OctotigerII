#include <gtest/gtest.h>
#include <limits>
#include "octotigerII/hydro/hydroSystem.hpp"
#include "testSupport.hpp"

using namespace octotigerII;
using namespace octotigerII::hydro;

namespace {

class Hydro : public ::testing::Test {
protected:

	HydroSystem gas{1.4};

	PrimitiveState primitive(Real density = 2, Real pressure = 5) const {
		PrimitiveState p;
		p.density() = units::Density::from_value(density);
		p.pressure() = units::Pressure::from_value(pressure);
		for (int d = 0; d < ndim; ++d) p.velocity(d) = units::Velocity::from_value(d + 1);
		return p;
	}
};


TEST_F(Hydro, PrimitiveRoundTripAndIndependentTotalEnergy) {
	auto p = primitive();
	auto const u = gas.conservedState(p);
	Real speedSquared = 0;
	for (int d = 0; d < ndim; ++d) {
		speedSquared += (d + 1) * (d + 1);
		EXPECT_DOUBLE_EQ(units::value(u.momentum(d)), 2 * (d + 1));
	}
	EXPECT_NEAR(units::value(u.totalEnergy()), 12.5 + speedSquared, 1e-13);
	test::expectStateNear(gas.reconstructionVariables(u), p);
	EXPECT_TRUE(gas.admissible(u));
}


TEST_F(Hydro, EulerFluxAndAcousticSpeedInEveryDirection) {
	auto const p = primitive();
	auto const u = gas.conservedState(p);
	for (int normal = 0; normal < ndim; ++normal) {
		SCOPED_TRACE(normal);
		auto const f = gas.physicalFlux(u, normal);
		Real const vn = normal + 1;
		EXPECT_DOUBLE_EQ(units::value(f.mass()), 2 * vn);
		for (int d = 0; d < ndim; ++d) EXPECT_DOUBLE_EQ(units::value(f.momentum(d)), 2 * (d + 1) * vn + (d == normal ? 5 : 0));
		EXPECT_NEAR(units::value(f.energy()), (units::value(u.totalEnergy()) + 5) * vn, 1e-12);
		EXPECT_NEAR(units::value(gas.maximumSignalSpeed(u, normal)), vn + std::sqrt(3.5), 1e-14);
		test::expectStateNear(gas.riemann(u, u, normal), f);
	}
}


TEST_F(Hydro, HllcResolvesStationaryContactWithoutMassOrEnergyFlux) {
	auto left = primitive(1, 1), right = primitive(0.125, 1);
	for (int d = 0; d < ndim; ++d) left.velocity(d) = right.velocity(d) = {};
	for (int normal = 0; normal < ndim; ++normal) {
		auto const flux = gas.riemann(gas.conservedState(left), gas.conservedState(right), normal);
		EXPECT_NEAR(units::value(flux.mass()), 0, 1e-14);
		EXPECT_NEAR(units::value(flux.energy()), 0, 1e-14);
		for (int d = 0; d < ndim; ++d) EXPECT_NEAR(units::value(flux.momentum(d)), d == normal ? 1 : 0, 1e-14);
	}
}


TEST_F(Hydro, SupersonicRiemannFluxUsesUpwindState) {
	for (int normal = 0; normal < ndim; ++normal) for (Real sign : {-1.0, 1.0}) {
		auto left = primitive(1, 1), right = primitive(2, 0.5);
		left.velocity(normal) = units::Velocity::from_value(sign * 10);
		right.velocity(normal) = units::Velocity::from_value(sign * 11);
		auto const l = gas.conservedState(left), r = gas.conservedState(right);
		test::expectStateNear(gas.riemann(l, r, normal), gas.physicalFlux(sign > 0 ? l : r, normal));
	}
}


TEST_F(Hydro, ReflectionChangesOnlyNormalMomentumAndIsInvolutive) {
	auto const u = gas.conservedState(primitive());
	for (int normal = 0; normal < ndim; ++normal) {
		auto const r = gas.reflected(u, normal);
		EXPECT_EQ(r.density(), u.density());
		EXPECT_EQ(r.totalEnergy(), u.totalEnergy());
		for (int d = 0; d < ndim; ++d) EXPECT_EQ(r.momentum(d), (d == normal ? -1.0 : 1.0) * u.momentum(d));
		test::expectStateNear(gas.reflected(r, normal), u, 0);
	}
}


TEST_F(Hydro, RejectsVacuumNegativePressureAndNonfiniteState) {
	EXPECT_FALSE(gas.admissible({}));
	EXPECT_THROW(gas.reconstructionVariables({}), std::runtime_error);
	auto p = primitive();
	p.pressure() = units::Pressure::from_value(-1);
	EXPECT_THROW(gas.conservedState(p), std::invalid_argument);
	for (Real value : {std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		p = primitive();
		p.velocity(0) = units::Velocity::from_value(value);
		EXPECT_THROW(gas.conservedState(p), std::invalid_argument);
		auto u = gas.conservedState(primitive());
		u.totalEnergy() = units::EnergyDensity::from_value(value);
		EXPECT_FALSE(gas.admissible(u));
	}
}


TEST_F(Hydro, FluxLimiterRepairsDestructiveFluxConservatively) {
	auto const u = gas.conservedState(primitive());
	auto high = gas.physicalFlux(u, 0);
	high.mass() += units::MassFlux::from_value(1e6);
	auto const dtDx = units::TimePerLength::from_value(0.001 / ndim);
	auto const f = gas.limitFlux(u, u, high, 0, dtDx);
	auto const delta = HydroSystem::integratedFlux(f - gas.physicalFlux(u, 0), Real(2 * ndim) * dtDx);
	EXPECT_TRUE(gas.admissible(u - delta));
	EXPECT_TRUE(gas.admissible(u + delta));
	EXPECT_LT(units::value(f.mass()), units::value(high.mass()));
	test::expectStateNear((u - delta) + (u + delta), 2.0 * u);
	test::expectStateNear(gas.limitFlux(u, u, high, 0, {}), high, 0);
}


TEST(HydroValidation, RejectsInvalidEosAndFloors) {
	for (Real gamma : {Real(1), Real(0), std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()})
		EXPECT_THROW(HydroSystem{gamma}, std::invalid_argument);
	for (Real floor : {Real(0), Real(-1), std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		EXPECT_THROW((HydroSystem{1.4, units::Density::from_value(floor)}), std::invalid_argument);
		EXPECT_THROW((HydroSystem{1.4, units::Density::from_value(1e-14), units::Pressure::from_value(floor)}), std::invalid_argument);
	}
}

} // namespace
