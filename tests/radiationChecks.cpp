#include <gtest/gtest.h>
#include <limits>
#include "octotigerII/radiation/radiationTransport.hpp"
#include "testSupport.hpp"

using namespace octotigerII;
using namespace octotigerII::radiation;

namespace {

class LightSpeed : public ::testing::TestWithParam<Real> {};


TEST_P(LightSpeed, IsotropicClosureHasPressureOneThirdAndAcousticWaves) {
	M1::State u;
	u[0] = units::EnergyDensity::from_value(3);
	auto const chat = GetParam() * constants::c;
	for (int normal = 0; normal < ndim; ++normal) {
		auto const waves = M1::physicalFlux(u, normal, chat);
		EXPECT_EQ(waves.flux[0], units::EnergyFlux{});
		for (int d = 0; d < ndim; ++d) EXPECT_NEAR(units::value(waves.flux[d + 1] / chat), d == normal ? 1 : 0, 1e-14);
		EXPECT_NEAR(Real(waves.minus / chat), -1 / std::sqrt(3.0), 1e-14);
		EXPECT_NEAR(Real(waves.plus / chat), 1 / std::sqrt(3.0), 1e-14);
	}
}


TEST_P(LightSpeed, StreamingClosureAndHllUseCorrectPropagationDirection) {
	auto const chat = GetParam() * constants::c;
	RadiationSystem system(chat);
	for (int normal = 0; normal < ndim; ++normal) for (Real sign : {-1.0, 1.0}) {
		RadiationSystem::State u;
		u.energy() = units::EnergyDensity::from_value(2);
		u.radiativeFlux(normal) = sign * constants::c * u.energy();
		auto const waves = M1::physicalFlux(RadiationSystem::toCalculationState(u), normal, chat);
		EXPECT_NEAR(Real(waves.minus / chat), sign, 1e-14);
		EXPECT_NEAR(Real(waves.plus / chat), sign, 1e-14);
		auto const flux = system.physicalFlux(u, normal);
		EXPECT_NEAR(units::value(flux.energy() / chat), sign * 2, 1e-14);
		EXPECT_NEAR(units::value(flux.radiativeFlux(normal) / (chat * constants::c)), 2, 1e-14);
		test::expectStateNear(system.riemann(u, 3.0 * u, normal), system.physicalFlux(sign > 0 ? u : RadiationSystem::State(3.0 * u), normal));
	}
}


TEST_P(LightSpeed, ObliqueClosureIsSymmetricTraceOneAndCausal) {
	auto const chat = GetParam() * constants::c;
	for (Real reduced : {0.0, 0.01, 0.5, 0.99, 1.0}) {
		M1::State u;
		u[0] = units::EnergyDensity::from_value(7);
		Real norm = 0;
		for (int d = 0; d < ndim; ++d) norm += (d + 1) * (d + 1);
		for (int d = 0; d < ndim; ++d) u[d + 1] = u[0] * (reduced * (d + 1) / std::sqrt(norm));
		Real trace = 0;
		for (int a = 0; a < ndim; ++a) {
			auto const f = M1::physicalFlux(u, a, chat);
			trace += Real(f.flux[a + 1] / (chat * u[0]));
			EXPECT_LE(f.minus, f.plus);
			EXPECT_GE(Real(f.minus / chat), -1 - 1e-14);
			EXPECT_LE(Real(f.plus / chat), 1 + 1e-14);
			for (int b = 0; b < ndim; ++b) {
				auto const other = M1::physicalFlux(u, b, chat);
				EXPECT_NEAR(Real(f.flux[b + 1] / (chat * u[0])), Real(other.flux[a + 1] / (chat * u[0])), 1e-14);
			}
			test::expectStateNear(M1::hll(u, u, a, chat), f.flux);
		}
		// The closure is physically 3D even in a reduced-dimensional transport
		// build; each unrepresented diagonal contributes (1-chi)/2.
		Real const chi = (3 + 4 * reduced * reduced) / (5 + 2 * std::sqrt(4 - 3 * reduced * reduced));
		EXPECT_NEAR(trace + (3 - ndim) * (1 - chi) / 2, 1, 2e-14);
	}
}


INSTANTIATE_TEST_SUITE_P(TransportSpeeds, LightSpeed, ::testing::Values(1.0, 0.25));


TEST(Radiation, VacuumFluxAndPhysicalFluxRoundTrip) {
	RadiationSystem system(constants::c);
	auto const flux = system.riemann({}, {}, 0);
	flux.forEach([](auto, auto q) { EXPECT_EQ(q, decltype(q){}); });
	EXPECT_TRUE(system.admissible({}));
	std::array<units::EnergyFlux, ndim> physical{};
	physical[0] = units::EnergyFlux::from_value(1e10);
	auto const u = RadiationSystem::fromPhysical(units::EnergyDensity::from_value(2), physical);
	EXPECT_EQ(RadiationSystem::toPhysicalFlux(u), physical);
	test::expectStateNear(system.conservedState(system.reconstructionVariables(u)), u);
	for (int axis = 0; axis < ndim; ++axis) {
		auto const reflected = system.reflected(u, axis);
		EXPECT_EQ(reflected.energy(), u.energy());
		for (int d = 0; d < ndim; ++d) EXPECT_EQ(reflected.radiativeFlux(d), (d == axis ? -1.0 : 1.0) * u.radiativeFlux(d));
		test::expectStateNear(system.reflected(reflected, axis), u, 0);
	}
}


TEST(Radiation, RealizabilityRejectsNegativeEnergySuperluminalAndNonfiniteFlux) {
	M1::State u;
	u[0] = units::EnergyDensity::from_value(1);
	u[1] = units::EnergyDensity::from_value(1.01);
	EXPECT_FALSE(M1::admissible(u));
	EXPECT_THROW(M1::checkState(u), std::runtime_error);
	u[1] = {};
	u[0] = units::EnergyDensity::from_value(-1);
	EXPECT_FALSE(M1::admissible(u));
	u[0] = units::EnergyDensity::from_value(1);
	for (Real value : {std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		u[1] = units::EnergyDensity::from_value(value);
		EXPECT_FALSE(M1::admissible(u));
	}
}


TEST(Radiation, RoundoffRepairDoesNotHidePhysicalViolations) {
	M1::State u, scale;
	scale[0] = units::EnergyDensity::from_value(1);
	u[0] = units::EnergyDensity::from_value(-epsilonR);
	auto const repaired = M1::roundoffState(u, scale);
	EXPECT_EQ(repaired[0], units::EnergyDensity{});
	u[0] = units::EnergyDensity::from_value(-0.001);
	EXPECT_THROW(M1::roundoffState(u, scale), std::runtime_error);
	u[0] = units::EnergyDensity::from_value(1);
	u[1] = u[0] * (1 + 2 * epsilonR);
	EXPECT_LE(M1::magnitude(M1::canonical(u)), u[0]);
	u[1] = u[0] * 1.01;
	EXPECT_THROW(M1::canonical(u), std::runtime_error);
}


TEST(Radiation, RejectsInvalidSpeedAndAxis) {
	for (Real speed : {Real(0), Real(-1), std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()})
		EXPECT_THROW(RadiationSystem{units::Velocity::from_value(speed)}, std::invalid_argument);
	for (int axis : {-1, ndim}) EXPECT_THROW(M1::physicalFlux({}, axis, constants::c), std::invalid_argument);
}

} // namespace

TEST(Radiation, FluxLimiterPreservesRealizabilityAndPairConservation) {
	RadiationSystem system(0.25 * constants::c);
	RadiationSystem::State u;
	u.energy() = units::EnergyDensity::from_value(1);
	auto const physical = system.physicalFlux(u, 0);
	auto high = physical;
	high.energy() = units::EnergyFlux::from_value(1e20);
	auto const dtDx = Real(0.01 / ndim) / system.reducedLightSpeed();
	auto const limited = system.limitFlux(u, u, high, 0, dtDx);
	auto const delta = RadiationSystem::integratedFlux(limited - physical, Real(2 * ndim) * dtDx);
	// The limiter allows roundoff at the cone boundary, which is canonicalized
	// by the subsequent conservative update before admissibility is enforced.
	auto const left = system.correctRoundoff(u - delta, componentAbs(delta));
	auto const right = system.correctRoundoff(u + delta, componentAbs(delta));
	EXPECT_TRUE(system.admissible(left));
	EXPECT_TRUE(system.admissible(right));
	test::expectStateNear(left + right, RadiationSystem::State(2.0 * u), 2e-12);
	EXPECT_LT(units::value(limited.energy()), units::value(high.energy()));
	test::expectStateNear(system.limitFlux(u, u, high, 0, {}), high, 0);
}
