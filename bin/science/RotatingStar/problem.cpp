/** @file
 * @brief Original Octo-Tiger oblate SCF rotating-star benchmark.
 */
// Distributed under the Boost Software License, Version 1.0.
#include <array>
#include <cmath>
#include <stdexcept>
#include "octotigerII/problems.hpp"
#include "octotigerII/problems/rotatingStar.hpp"
#include "octotigerII/verification/analytic.hpp"

namespace octotigerII::problems {
namespace {
	constexpr std::array<RotatingStar::Profile, RotatingStar::tableCells * RotatingStar::tableCells> equilibrium{{
#include "equilibrium.inc"
	}};

	// The original initializer extrapolates the first four positive samples at
	// the axis and uses forward four-point interpolation elsewhere. Express its
	// cubic Lagrange coefficients exactly instead of rounded decimal constants.
	std::array<Real, 4> weights(Real x) {
		return {-(x-1)*(x-2)*(x-3)/6, x*(x-2)*(x-3)/2,
			-x*(x-1)*(x-3)/2, x*(x-1)*(x-2)/6};
	}
}

RotatingStar::RotatingStar(units::Length length, units::Density density, Real atmosphere)
  : lengthUnit_(length), densityUnit_(density), atmosphereFraction_(atmosphere) {
	if (!(length > units::Length{}) || !units::finite(length) ||
		!(density > units::Density{}) || !units::finite(density) ||
		!std::isfinite(atmosphere) || !(atmosphere > 0 && atmosphere < 1))
		throw std::invalid_argument("Rotating star requires positive length/density units and an atmosphere fraction between zero and one");
}

RotatingStar::Profile RotatingStar::profile(Real radius, Real z) {
	if (!std::isfinite(radius) || !std::isfinite(z)) throw std::invalid_argument("Nonfinite rotating-star coordinates");
	Real const rcell = std::abs(radius) * tableCells - Real(0.5);
	Real const zcell = std::abs(z) * tableCells - Real(0.5);
	// Check before converting to int, including coordinates far outside the box.
	if (rcell >= tableCells - 3 || zcell >= tableCells - 3) return {};
	int const i = int(rcell), k = int(zcell);
	auto const wr = weights(rcell - i), wz = weights(zcell - k);
	Profile result;
	for (int a = 0; a < 4; ++a) for (int b = 0; b < 4; ++b) {
		auto const& sample = equilibrium[std::size_t((i+a) * tableCells + k+b)];
		Real const w = wr[a] * wz[b];
		result.density += w * sample.density;
		result.internalEnergy += w * sample.internalEnergy;
	}
	result.density = std::max(Real(0), result.density);
	result.internalEnergy = std::max(Real(0), result.internalEnergy);
	return result;
}

units::InverseTime RotatingStar::angularVelocity() const {
	return spin * units::sqrt(constants::G * densityUnit_);
}

units::EnergyDensity RotatingStar::energyUnit() const {
	return constants::G * densityUnit_ * densityUnit_ * lengthUnit_ * lengthUnit_;
}

hydro::PrimitiveState RotatingStar::operator()(mesh::PhysicalCoordinates const& x) const {
	auto const value = profile(Real(units::hypot(x[0], x[1]) / lengthUnit_), Real(x[2] / lengthUnit_));
	hydro::PrimitiveState result;
	result.density() = std::max(value.density, atmosphereFraction_) * densityUnit_;
	result.pressure() = (gamma - 1) * std::max(value.internalEnergy, atmosphereFraction_) * energyUnit();
	// Match the original initializer: momentum is formed before flooring rho,
	// so vacuum atmosphere is initially stationary in the inertial frame.
	Real const massRatio = value.density / std::max(value.density, atmosphereFraction_);
	result.velocity(0) = -massRatio * angularVelocity() * x[1];
	result.velocity(1) = massRatio * angularVelocity() * x[0];
	return result;
}
} // namespace octotigerII::problems

namespace octotigerII::rotatingStar {

ProblemBoundary problemBoundary(Config const&) { return {}; }

verification::Reference problemReference(Config const&) {
	verification::Reference result;
	result.reason = "The tabulated SCF equilibrium has no tabulated gravitational potential; use verification.gravityReference=direct";
	return result;
}

void problemDefaults(Config& c) {
	c.mesh.cells = 4;
	c.mesh.lower = -2.0 * c.star.radius;
	c.mesh.upper = 2.0 * c.star.radius;
	c.hydro.gamma = problems::RotatingStar::gamma;
	c.star.polytropicIndex = 1.5;
	c.star.atmosphereFraction = 1e-10;
	c.amr.refineDensity = 0.01 * c.star.centralDensity;
	c.amr.shadowTolerance = 0;
	c.runtime.stopTime = units::Time::from_value(0.1);
}

void validateProblem(Config const& c) {
	problems::RotatingStar const star(c.star.radius, c.star.centralDensity, c.star.atmosphereFraction);
	for (auto x : c.star.center)
		if (!units::finite(x)) throw std::invalid_argument("Rotating-star center must be finite");
	if (c.star.polytropicIndex != Real(1.5) || std::abs(c.hydro.gamma - problems::RotatingStar::gamma) > 1e-12)
		throw std::invalid_argument("The original rotating-star SCF model requires n=1.5 and hydro.gamma=5/3");
	if (!c.mesh.boundary.all(physics::BoundaryCondition::Outflow))
		throw std::invalid_argument("The isolated rotating-star benchmark requires outflow boundaries");
	if (c.amr.enabled && c.amr.maxLevel >= 0 && c.amr.maxLevel <= 16 && c.mesh.cells > 0 &&
		star.coreLength() < (c.mesh.upper - c.mesh.lower) / Real(c.mesh.cells * (1 << c.amr.maxLevel)))
		throw std::invalid_argument("Rotating-star core is unresolved at amr.maxLevel; increase maxLevel or mesh.cells, or reduce the domain");
	if (c.star.atmosphereFraction * c.star.centralDensity < units::Density::from_value(1e-14) ||
		(problems::RotatingStar::gamma - 1) * c.star.atmosphereFraction * star.energyUnit() < units::Pressure::from_value(1e-14))
		throw std::invalid_argument("Rotating-star atmosphere falls below the hydro density or pressure floor");
}

void initializeProblem(Snapshot& data, Config const& c, bool refinementProbe) {
	problems::RotatingStar star(c.star.radius, c.star.centralDensity, c.star.atmosphereFraction);
	if (refinementProbe) {
		auto const width = initialFeatureWidth(star.coreLength(), data.cellWidth);
		star = problems::RotatingStar(c.star.radius * Real(width / star.coreLength()), c.star.centralDensity, c.star.atmosphereFraction);
	}
	hydro::HydroSystem const gas(c.hydro);
	data.layout.forEachInterior([&](mesh::Coordinates const& cell, std::size_t i) {
		auto relative = data.layout.cellCenter(data.lower, data.cellWidth, cell);
		for (int d = 0; d < ndim; ++d) relative[d] -= c.star.center[d];
		data.hydro.values()[i] = gas.conservedState(star(relative));
	});
}
} // namespace octotigerII::rotatingStar
