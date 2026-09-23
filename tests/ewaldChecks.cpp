#include <cmath>
#include <gtest/gtest.h>
#include "octotigerII/gravity/ewald.hpp"
#include "octotigerII/gravity/images.hpp"
#include "octotigerII/gravity/solver.hpp"

using namespace octotigerII;
using namespace octotigerII::gravity;
namespace {
using V = diagonal::Vector;
using B = physics::BoundaryConditions;
using R = physics::BoundaryCondition;
constexpr double pi = 3.14159265358979323846;
std::array<double, 4> pair(V x, V periods) {
	using std::hypot;
	auto k = ewald::derivatives(x, periods, 1);
	double const r = hypot(x[0], x[1], x[2]);
	std::array<double, 4> result{k.values[0], k.at(1, 0, 0), k.at(0, 1, 0), k.at(0, 0, 1)};
	if (r > 0) {
		result[0] -= 1 / r;
		for (int d = 0; d < 3; ++d)
			result[d + 1] += x[d] / (r * r * r);
	}
	return result;
}
TEST(EwaldKernel, SplitAndQuadratureIndependenceIncludingHighDerivatives) {
	using std::abs;
	for (V periods : {V{1, 0, 0}, V{0, 2, 0}, V{1, 2, 0}, V{0, 1, 1}, V{1, 1, 2}}) {
		SCOPED_TRACE(::testing::PrintToString(periods));
		V x{0.17, -0.13, 0.23};
		auto a = ewald::derivatives(x, periods, 12, 2, 64);
		auto b = ewald::derivatives(x, periods, 12, 2.7, 128);
		for (std::size_t i = 0; i < a.values.size(); ++i)
			EXPECT_NEAR(a.values[i], b.values[i], 3e-9 * (1 + abs(a.values[i]))) << i;
	}
}
TEST(EwaldKernel, GradientAndHessianArePotentialDerivatives) {
	for (V periods : {V{1, 0, 0}, V{1, 1, 0}, V{1, 1, 1}}) {
		V x{0.17, -0.13, 0.23};
		auto k = ewald::derivatives(x, periods, 3);
		for (int d = 0; d < 3; ++d) {
			auto a = x, b = x;
			a[d] += 1e-5;
			b[d] -= 1e-5;
			auto plus = ewald::derivatives(a, periods, 2), minus = ewald::derivatives(b, periods, 2);
			for (int j = 0; j < 4; ++j) {
				std::array<int, 3> e{};
				if (j) ++e[j - 1];
				double fd = (plus.at(e[0], e[1], e[2]) - minus.at(e[0], e[1], e[2])) / 2e-5;
				++e[d];
				EXPECT_NEAR(fd, k.at(e[0], e[1], e[2]), 2e-7);
			}
		}
	}
}
TEST(EwaldKernel, SinglePeriodicAxisMatchesLongDirectImageSum) {
	using std::sqrt;
	V x{0.17, 0.31, -0.23};
	for (int axis = 0; axis < 3; ++axis) {
		V periods{};
		periods[axis] = 1;
		auto result = pair(x, periods);
		std::array<long double, 3> sum{};
		for (int image = -20000; image <= 20000; ++image) {
			std::array<long double, 3> r{x[0], x[1], x[2]};
			r[axis] -= image;
			long double distance = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
			for (int d = 0; d < 3; ++d)
				sum[d] += r[d] / (distance * distance * distance);
		}
		for (int d = 0; d < 3; ++d)
			EXPECT_NEAR(result[d + 1], double(sum[d]), 2e-9);
	}
}
TEST(EwaldKernel, SlabMatchesIndependentFourierGreenFunctionAndSheetField) {
	using std::cos;
	using std::exp;
	using std::hypot;
	using std::sin;
	V x{0.21, -0.14, 0.47};
	std::array<double, 4> expected{2 * pi * x[2], 0, 0, 2 * pi};
	for (int a = -20; a <= 20; ++a)
		for (int b = -20; b <= 20; ++b) {
			if (!a && !b) continue;
			double kx = 2 * pi * a, ky = 2 * pi * b, k = hypot(kx, ky), angle = kx * x[0] + ky * x[1];
			double weight = 2 * pi * exp(-k * x[2]);
			expected[0] -= weight / k * cos(angle);
			expected[1] += weight / k * kx * sin(angle);
			expected[2] += weight / k * ky * sin(angle);
			expected[3] += weight * cos(angle);
		}
	auto actual = pair(x, {1, 1, 0});
	for (int d = 0; d < 4; ++d)
		EXPECT_NEAR(actual[d], expected[d], 2e-12);
	EXPECT_NEAR(pair({0.2, 0.1, 4}, {1, 1, 0})[3], 2 * pi, 2e-9);
}
TEST(EwaldKernel, PeriodicSeamsSelfLimitAndCubicMadelungConstant) {
	using std::isfinite;
	for (V period : {V{1, 0, 0}, V{1, 1, 0}, V{1, 1, 1}}) {
		auto self = ewald::derivatives({}, period, 4);
		EXPECT_TRUE(isfinite(self.values[0]));
		for (int d = 0; d < 3; ++d) {
			std::array<int, 3> e{};
			e[d] = 1;
			EXPECT_NEAR(self.at(e[0], e[1], e[2]), 0, 1e-14);
			if (!period[d]) continue;
			V a{0.11, 0.13, 0.17}, b = a;
			a[d] = -0.5;
			b[d] = 0.5;
			auto ka = pair(a, period), kb = pair(b, period);
			for (int j = 0; j < 4; ++j)
				EXPECT_NEAR(ka[j], kb[j], 2e-12);
		}
	}
	auto self = ewald::derivatives({}, {1, 1, 1}, 2);
	EXPECT_NEAR(self.values[0], 2.837297479480619, 3e-13);
	EXPECT_NEAR(self.at(2, 0, 0), -4 * pi / 3, 2e-13);
	EXPECT_NEAR(self.at(0, 2, 0), -4 * pi / 3, 2e-13);
	EXPECT_NEAR(self.at(0, 0, 2), -4 * pi / 3, 2e-13);
}
TEST(EwaldM2L, TraceAndLocalLaplacianAgreeWithIndividualSources) {
	int const p = 6;
	diagonal::Coefficients m(diagonal::coefficientCount(p) + 1, 0), l(m.size(), 0);
	std::array<V, 2> positions{V{0.1, -0.2, 0.27}, V{-0.13, 0.08, -0.24}};
	for (auto s : positions) {
		auto a = imageShiftMultipole({1}, s, 1, p);
		for (std::size_t i = 0; i < m.size(); ++i)
			m[i] += a[i];
	}
	ewald::getOperator(p, {0, 0, 0}, {16, 16, 16})->add(l, m, 1.0 / 16);
	V target{0.13, -0.16, 0.23};
	auto shifted = imageShiftLocal(l, target, 1, p);
	std::array<double, 4> exact{};
	for (auto source : positions) {
		V delta{};
		for (int d = 0; d < 3; ++d)
			delta[d] = (target[d] - source[d]) / 16;
		auto k = ewald::derivatives(delta, {1, 1, 1}, 1);
		for (int j = 0; j < 4; ++j)
			exact[j] += k.values[j] / (j == 0 ? 1 : 16);
	}
	for (int j = 0; j < 4; ++j)
		EXPECT_NEAR(shifted[j], exact[j], 2e-11);
	EXPECT_NEAR(l.back(), -8 * pi / (16 * 16), 2e-15);
}
TEST(GravityImages, GeometryAndReflectionParity) {
	using std::sqrt;
	B bc;
	bc.lower[0] = R::Reflecting;
	bc.upper[1] = R::Reflecting;
	bc.lower[2] = bc.upper[2] = R::Periodic;
	ImageGeometry images(bc);
	ASSERT_EQ(images.images().size(), 4u);
	EXPECT_EQ(images.periods(8), (diagonal::Offset{0, 0, 8}));
	EXPECT_EQ(images.separation({0, 0, 0}, {1, 1, 7}, 8, images.images()[1]), (diagonal::Offset{2, -1, 1}));
	EXPECT_EQ(images.separation({0, 0, 0}, {1, 1, 7}, 8, images.images()[2]), (diagonal::Offset{-1, -14, 1}));
	B walls;
	walls.lower[0] = walls.upper[0] = R::Reflecting;
	EXPECT_EQ(ImageGeometry(walls).periods(8), (diagonal::Offset{16, 0, 0}));
	ImageGeometry periodic(B::periodic());
	EXPECT_NEAR(periodic.nextImageDistance({1, 2, 0}, 8), sqrt(37.0), 1e-14);
	EXPECT_FALSE(periodic.acceptable({4, 0, 0}, 8, 0.5, true));	   // branch cut
	EXPECT_TRUE(periodic.acceptable({0, 0, 0}, 8, 0.5, true));	   // self correction is far
	auto a = imageShiftMultipole({2}, {0.1, -0.2, 0.3}, 1, 5);
	auto b = reflectMultipole(a, 5, 5), expected = imageShiftMultipole({2}, {-0.1, -0.2, -0.3}, 1, 5);
	for (std::size_t i = 0; i < a.size(); ++i)
		EXPECT_NEAR(b[i], expected[i], 1e-15);
}
std::vector<State> directImages(std::vector<units::Density> const& rho, int n, B const& bc) {
	ImageGeometry images(bc);
	std::vector<State> result(rho.size());
	auto xyz = [n](std::size_t i) { return diagonal::Offset{int(i % n), int(i / n % n), int(i / (n * n))}; };
	for (std::size_t i = 0; i < rho.size(); ++i) {
		diagonal::Coefficients l(5, 0);
		for (std::size_t j = 0; j < rho.size(); ++j)
			if (units::value(rho[j]))
				for (auto image : images.images()) {
					auto r = images.separation(xyz(i), xyz(j), n, image);
					if (r != diagonal::Offset{}) diagonal::addDirect(l, units::value(rho[j]), r, 1);
					if (images.periodic()) ewald::addDirect(l, units::value(rho[j]), r, images.periods(n), 1);
				}
		result[i].potential() = units::VelocitySquared::from_value(units::value(constants::G) * l[0]);
		for (int d = 0; d < 3; ++d) {
			std::array<int, 3> e{};
			e[d] = 1;
			result[i].acceleration(d) = units::Acceleration::from_value(-units::value(constants::G) * diagonal::derivative(l, e[0], e[1], e[2]));
		}
	}
	return result;
}
TEST(GravityImages, FmmMatchesDirectWithIndependentImageAcceptance) {
	using std::abs;
	std::vector<B> cases;
	for (int mask = 1; mask < 8; ++mask) {
		B b;
		for (int d = 0; d < 3; ++d)
			if (mask & (1 << d)) b.lower[d] = b.upper[d] = R::Periodic;
		cases.push_back(b);
	}
	B lower;
	lower.lower[0] = R::Reflecting;
	cases.push_back(lower);
	B mixed = lower;
	mixed.upper[1] = R::Reflecting;
	mixed.lower[2] = mixed.upper[2] = R::Periodic;
	cases.push_back(mixed);
	B both;
	both.lower[0] = both.upper[0] = R::Reflecting;
	cases.push_back(both);
	cases.push_back(B::uniform(R::Reflecting));
	int const n = 8;
	std::vector<units::Density> rho(n * n * n);
	for (int i : {0, 17, 201, 456, 511})
		rho[i] = units::Density::from_value(1 + double(i % 7));
	for (std::size_t c = 0; c < cases.size(); ++c) {
		SCOPED_TRACE(c);
		auto exact = directImages(rho, n, cases[c]);
		auto actual = solve(rho, n, units::Length::from_value(1), 8, 0.5, cases[c]);
		double scale = 0, error = 0;
		for (std::size_t i = 0; i < rho.size(); ++i)
			for (int d = 0; d < 3; ++d) {
				scale = std::max(scale, abs(units::value(exact[i].acceleration(d))));
				error = std::max(error, abs(units::value(actual.fields[i].acceleration(d) - exact[i].acceleration(d))));
			}
		EXPECT_LT(error, 2e-5 * scale);
		auto coarse = solve(rho, n, units::Length::from_value(1), 3, 0.5, cases[c]);
		double coarseError = 0;
		for (std::size_t i = 0; i < rho.size(); ++i)
			for (int d = 0; d < 3; ++d)
				coarseError = std::max(coarseError, abs(units::value(coarse.fields[i].acceleration(d) - exact[i].acceleration(d))));
		EXPECT_LT(error, 0.03 * coarseError);
		EXPECT_GT(actual.statistics.multipolePairs, 0u);
		if (ImageGeometry(cases[c]).periodic()) EXPECT_GT(actual.statistics.ewaldPairs, 0u);
	}
}
TEST(GravityImages, ClosedReflectingBoxAndPeriodicUniformLatticeHaveZeroForce) {
	for (auto bc : {B::periodic(), B::uniform(R::Reflecting)}) {
		auto solution = solve(std::vector<units::Density>(64, units::Density::from_value(1)), 4, units::Length::from_value(1), 8, 0.1, bc);
		for (auto const& field : solution.fields)
			for (int d = 0; d < 3; ++d)
				EXPECT_NEAR(units::value(field.acceleration(d)), 0, 1e-19);
	}
}

TEST(EwaldKernel, HighestOrderAndLargeOpenDisplacementsConverge) {
	using std::abs;
	for (auto periods : {V{1, 0, 0}, V{1, 1, 0}, V{1, 1, 1}}) {
		V r{0.49, periods[1] ? 0.41 : 1.8, periods[2] ? 0.39 : -1.9};
		auto a = ewald::derivatives(r, periods, 20, 2, 64);
		auto b = ewald::derivatives(r, periods, 20, 2.3, 128);
		// Scale by the derivative degree; individual symmetry/cancellation zeros
		// can be far smaller than the tensor norm.
		for (int n = 0; n <= 20; ++n) {
			double scale = 1;
			for (int i = n * n; i < (n + 1) * (n + 1); ++i)
				scale = std::max(scale, abs(b.values[i]));
			for (int i = n * n; i < (n + 1) * (n + 1); ++i)
				EXPECT_NEAR(a.values[i], b.values[i], 2e-7 * scale) << n << ":" << i;
		}
	}
}
TEST(GravityImages, ThreeSingleWallsMatchExplicitMirroredMasses) {
	using std::abs;
	int const n = 4, N = 2 * n;
	std::vector<units::Density> original(n * n * n), extended(N * N * N);
	for (int z = 0; z < n; ++z)
		for (int y = 0; y < n; ++y)
			for (int x = 0; x < n; ++x) {
				auto rho = units::Density::from_value(1 + (x + 3 * y + 7 * z) % 11);
				original[(z * n + y) * n + x] = rho;
				for (int mask = 0; mask < 8; ++mask) {
					int a = mask & 1 ? n - 1 - x : n + x, b = mask & 2 ? n - 1 - y : n + y, c = mask & 4 ? n - 1 - z : n + z;
					extended[(c * N + b) * N + a] = rho;
				}
			}
	B bc;
	bc.lower.fill(R::Reflecting);
	auto image = solve(original, n, units::Length::from_value(1), 3, 0.1, bc);
	auto exact = solve(extended, N, units::Length::from_value(1), 3, 0.1);
	for (int z = 0; z < n; ++z)
		for (int y = 0; y < n; ++y)
			for (int x = 0; x < n; ++x) {
				auto const& a = image.fields[(z * n + y) * n + x];
				auto const& b = exact.fields[((z + n) * N + y + n) * N + x + n];
				EXPECT_NEAR(units::value(a.potential() - b.potential()), 0, 2e-18);
				for (int d = 0; d < 3; ++d)
					EXPECT_NEAR(units::value(a.acceleration(d) - b.acceleration(d)), 0, 2e-19);
			}
}
TEST(GravityImages, RootEwaldCorrectionSurvivesDownwardPropagation) {
	using std::abs;
	B const bc = B::uniform(R::Reflecting);
	std::vector<units::Density> rho(64);
	rho[7] = units::Density::from_value(2);
	rho[31] = units::Density::from_value(3);
	auto direct = solve(rho, 4, units::Length::from_value(1), 3, 0.1, bc);
	auto fmm = solve(rho, 4, units::Length::from_value(1), 10, 0.57, bc);
	double scale = 0, error = 0;
	for (std::size_t i = 0; i < rho.size(); ++i)
		for (int d = 0; d < 3; ++d) {
			scale = std::max(scale, abs(units::value(direct.fields[i].acceleration(d))));
			error = std::max(error, abs(units::value(fmm.fields[i].acceleration(d) - direct.fields[i].acceleration(d))));
		}
	EXPECT_LT(error, 1e-4 * scale);
	EXPECT_GT(fmm.statistics.multipolePairs, 0u);
}
}	 // namespace
