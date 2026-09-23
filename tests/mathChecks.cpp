#include <gtest/gtest.h>
#include "octotigerII/math/Vector.hpp"
#include "octotigerII/units/state.hpp"

using namespace octotigerII;

namespace {

TEST(VectorMath, InitializationArithmeticAndDotProduct) {
	Vector<double, 3> zero;
	for (int i = 0; i < 3; ++i) EXPECT_EQ(zero[i], 0);
	Vector<double, 3> a{1, -2, 3}, b{2, 4, -1};
	EXPECT_EQ(a.dot(b), -9);
	EXPECT_EQ(sqr(a), 14);
	EXPECT_EQ(a.min(), -2);
	EXPECT_EQ(a.max(), 3);
	auto const sum = a + b;
	auto const difference = a - b;
	for (int i = 0; i < 3; ++i) {
		EXPECT_EQ(sum[i], a[i] + b[i]);
		EXPECT_EQ(difference[i], a[i] - b[i]);
	}
	EXPECT_NEAR(abs(norm(a)), 1, 2e-15);
	EXPECT_THROW((Vector<double, 2>{1, 2, 3}), std::length_error);
	Vector<double, 3> partial{7};
	EXPECT_EQ(partial[0], 7);
	EXPECT_EQ(partial[1], 0);
	EXPECT_EQ(partial[2], 0);
}


TEST(VectorMath, InPlaceScalarAliasingAndIntegerDivision) {
	Vector<double, 3> a{2, 3, 4};
	a *= a[0];
	EXPECT_EQ(a[0], 4);
	EXPECT_EQ(a[1], 6);
	EXPECT_EQ(a[2], 8);
	a /= a[0];
	EXPECT_EQ(a[0], 1);
	EXPECT_EQ(a[1], 1.5);
	EXPECT_EQ(a[2], 2);
	Vector<int, 3> integers{8, -9, 10};
	integers /= integers[0];
	EXPECT_EQ(integers[0], 1);
	EXPECT_EQ(integers[1], -1);
	EXPECT_EQ(integers[2], 1);
}


TEST(VectorMath, SplitConcatenateAndUnitVectors) {
	Vector<double, 3> a{2, 3, 4};
	auto const [first, rest] = a.split<1>();
	auto const restored = concatenate(first, rest);
	for (int i = 0; i < 3; ++i) EXPECT_EQ(restored[i], a[i]);
	auto const [empty, all] = a.split<0>();
	EXPECT_EQ(empty.size(), 0);
	EXPECT_EQ(all.size(), 3);
	EXPECT_EQ(concatenate(1.0, a)[0], 1);
	EXPECT_EQ(concatenate(a, 5.0)[3], 5);
	for (int axis = 0; axis < 3; ++axis) {
		auto const unit = Vector<double, 3>::unit(axis);
		for (int i = 0; i < 3; ++i) EXPECT_EQ(unit[i], axis == i ? 1 : 0);
	}
}

} // namespace
