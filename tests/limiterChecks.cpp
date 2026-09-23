#include <gtest/gtest.h>
#include "octotigerII/physics/finiteVolume.hpp"

using namespace octotigerII;
using namespace octotigerII::physics;

namespace {

class SlopeLimiter : public ::testing::TestWithParam<Limiter> {};


TEST_P(SlopeLimiter, PreservesLinearDataAndFlattensExtrema) {
	for (Real difference : {-3.0, -0.25, 0.0, 0.25, 3.0}) {
		auto const d = units::Density::from_value(difference);
		EXPECT_EQ(limitedSlope(d, d, GetParam()), d);
		EXPECT_EQ(limitedSlope(d, -d, GetParam()), units::Density{});
		EXPECT_EQ(limitedSlope(d, units::Density{}, GetParam()), units::Density{});
	}
}


TEST_P(SlopeLimiter, IsOddSymmetricAndDoesNotExceedTvdBound) {
	for (Real left : {0.01, 0.5, 1.0, 8.0}) for (Real right : {0.01, 0.5, 1.0, 8.0}) {
		auto const a = units::Pressure::from_value(left), b = units::Pressure::from_value(right);
		auto const slope = limitedSlope(a, b, GetParam(), 1.3);
		EXPECT_EQ(limitedSlope(-a, -b, GetParam(), 1.3), -slope);
		EXPECT_EQ(limitedSlope(b, a, GetParam(), 1.3), slope);
		EXPECT_GT(units::value(slope), 0);
		EXPECT_LE(units::value(slope), 2 * std::min(left, right));
	}
}


INSTANTIATE_TEST_SUITE_P(Methods, SlopeLimiter, ::testing::Values(Limiter::Minmod, Limiter::VanLeer, Limiter::MinmodTheta));


TEST(SlopeValues, IndependentUnequalSlopeValues) {
	auto const a = units::Density::from_value(1), b = units::Density::from_value(3);
	EXPECT_DOUBLE_EQ(units::value(limitedSlope(a, b, Limiter::Minmod)), 1);
	EXPECT_DOUBLE_EQ(units::value(limitedSlope(a, b, Limiter::VanLeer)), 1.5);
	EXPECT_DOUBLE_EQ(units::value(limitedSlope(a, b, Limiter::MinmodTheta, 1.3)), 1.3);
	EXPECT_THROW(limitedSlope(a, b, static_cast<Limiter>(-1)), std::logic_error);
}

} // namespace
