#include <gtest/gtest.h>
#include <limits>
#include "octotigerII/gravity/solver.hpp"
#include "testSupport.hpp"

using namespace octotigerII;

namespace {

TEST(Gravity, VacuumHasExactlyZeroPotentialAndForce) {
	auto const solution = gravity::solve(std::vector<units::Density>(64), 4, units::Length::from_value(2), 3, 0.5);
	ASSERT_EQ(solution.fields.size(), 64u);
	for (auto const& field : solution.fields) field.forEach([](auto, auto q) { EXPECT_EQ(q, decltype(q){}); });
}


TEST(Gravity, SinglePointMassHasInverseRadiusPotentialAndInverseSquareForce) {
	std::vector<units::Density> density(64);
	density[0] = units::Density::from_value(2);
	auto const h = units::Length::from_value(3);
	auto const solution = gravity::solve(density, 4, h, 5, 0.1);
	ASSERT_EQ(solution.fields.size(), density.size());
	solution.fields[0].forEach([](auto, auto q) { EXPECT_EQ(q, decltype(q){}); });
	for (int z = 0; z < 4; ++z) for (int y = 0; y < 4; ++y) for (int x = 0; x < 4; ++x) {
		if (x == 0 && y == 0 && z == 0) continue;
		Real const radius = std::sqrt(Real(x*x + y*y + z*z));
		auto const& field = solution.fields[(z*4 + y)*4 + x];
		auto const phi = -constants::G * density[0] * h * h / radius;
		EXPECT_NEAR(Real(field.potential() / phi), 1, 2e-13);
		int const coords[] = {x, y, z};
		for (int d = 0; d < 3; ++d) {
			auto const g = -constants::G * density[0] * h * Real(coords[d]) / (radius * radius * radius);
			EXPECT_NEAR(units::value(field.acceleration(d) - g), 0, 2e-13 * units::value(constants::G * density[0] * h));
		}
	}
}


TEST(Gravity, SuperpositionAndDensityScaling) {
	std::vector<units::Density> a(64), b(64), combined(64);
	a[0] = units::Density::from_value(2);
	b[63] = units::Density::from_value(3);
	for (std::size_t i = 0; i < a.size(); ++i) combined[i] = 2.0 * a[i] + b[i];
	auto const h = units::Length::from_value(1);
	auto const fa = gravity::solve(a, 4, h), fb = gravity::solve(b, 4, h), fc = gravity::solve(combined, 4, h);
	for (std::size_t i = 0; i < a.size(); ++i) {
		auto const expected = 2.0 * fa.fields[i] + fb.fields[i];
		expected.forEach([&](auto f, auto q) {
			EXPECT_NEAR(units::value(fc.fields[i].template get<f>() - q), 0, 1e-19);
		});
	}
}


TEST(Gravity, RejectsInvalidGeometryOrderAngleAndDensity) {
	std::vector<units::Density> density(64, units::Density::from_value(1));
	auto const h = units::Length::from_value(1);
	EXPECT_THROW(gravity::solve(density, 3, h), std::exception);
	EXPECT_THROW(gravity::solve({}, 4, h), std::exception);
	for (Real width : {Real(0), Real(-1), std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()})
		EXPECT_THROW(gravity::solve(density, 4, units::Length::from_value(width)), std::exception);
	for (int order : {0, 11}) EXPECT_THROW(gravity::solve(density, 4, h, order), std::exception);
	for (Real theta : {Real(0), Real(-1), Real(0.6), std::numeric_limits<Real>::quiet_NaN()})
		EXPECT_THROW(gravity::solve(density, 4, h, 3, theta), std::exception);
	for (Real rho : {Real(-1), std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		density[0] = units::Density::from_value(rho);
		EXPECT_THROW(gravity::solve(density, 4, h), std::exception);
	}
}

} // namespace
