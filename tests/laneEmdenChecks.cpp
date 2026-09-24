#include <gtest/gtest.h>
#include <cmath>
#include <numbers>
#include "octotigerII/problems/laneEmden.hpp"
using namespace octotigerII;

TEST(LaneEmden, IndexOneMatchesClosedForm) {
	using std::sin;
	using std::cos;
	problems::LaneEmden const solution(1);
	EXPECT_NEAR(solution.surface(), std::numbers::pi, 2e-10);
	EXPECT_NEAR(solution.surfaceMass(), std::numbers::pi, 2e-9);
	for (Real xi : {0.01, 0.3, 1.0, 2.5, 3.0}) {
		auto const q = solution(xi);
		EXPECT_NEAR(q.theta, sin(xi)/xi, 2e-10);
		EXPECT_NEAR(q.mass, sin(xi)-xi*cos(xi), 2e-9);
	}
	EXPECT_EQ(solution(0).theta, 1);
	EXPECT_EQ(solution(0).mass, 0);
}

TEST(LaneEmden, StandardPolytropeSurfaceAndInvalidIndices) {
	problems::LaneEmden const solution(1.5);
	EXPECT_NEAR(solution.surface(), 3.653753736219, 2e-8);
	EXPECT_NEAR(solution.surfaceMass(), 2.714055120109, 2e-8);
	EXPECT_THROW(problems::LaneEmden(0), std::invalid_argument);
	EXPECT_THROW(problems::LaneEmden(5), std::invalid_argument);
}

TEST(LaneEmden, PhysicalScalingAndHydrostaticBalance) {
	auto const radius = units::Length::from_value(1e9);
	auto const density = units::Density::from_value(1e6);
	problems::Polytrope const star(1.5, radius, density), large(1.5, 2.0*radius, density);
	EXPECT_NEAR(Real(large.mass()/star.mass()), 8, 1e-13);
	EXPECT_NEAR(Real(large.centralPressure()/star.centralPressure()), 4, 1e-13);
	EXPECT_EQ(star(units::Length{}).density, density);
	EXPECT_EQ(star(radius).density, units::Density{});
	for (Real fraction : {0.1, 0.3, 0.7, 0.95}) {
		auto const r = fraction*radius, h = 1e-5*radius;
		auto const q = star(r);
		auto const gradient = (star(r+h).pressure-star(r-h).pressure)/(2.0*h);
		auto const force = -constants::G*q.enclosedMass*q.density/(r*r);
		EXPECT_NEAR(Real(gradient/force), 1, 2e-7);
	}
}
