#include <algorithm>
#include <cmath>
#include <numbers>
#include "octotigerII/verification/analytic.hpp"


namespace octotigerII::verification {


Reference sodReference(Config const& c) {
	using std::pow;
	using std::sqrt;

	Real const gamma = c.hydro.gamma;
	auto const rhoL = units::Density::from_value(1), rhoR = units::Density::from_value(0.125);
	auto const pL = units::Pressure::from_value(1), pR = units::Pressure::from_value(0.1);
	auto const aL = units::sqrt(gamma * pL / rhoL), aR = units::sqrt(gamma * pR / rhoR);
	// Match the left rarefaction and right shock curves. Bisection brackets the
	// unique star pressure for these Sod states for any supported gamma > 1.
	auto wave = [gamma](units::Pressure p, units::Pressure p0, units::Density rho, units::Velocity a) {
		if (p > p0) return (p - p0) * units::sqrt((2.0 / (gamma + 1) / rho) / (p + (gamma - 1) / (gamma + 1) * p0));
		return 2.0 * a / (gamma - 1) * (pow(Real(p / p0), (gamma - 1) / (2 * gamma)) - 1);
	};
	auto lo = pR, hi = pL;
	for (int iteration = 0; iteration < 100; ++iteration) {
		auto const p = (lo + hi) / 2.0;
		if (wave(p, pL, rhoL, aL) + wave(p, pR, rhoR, aR) > units::Velocity{})
			hi = p;
		else
			lo = p;
	}
	auto const pStar = (lo + hi) / 2.0;
	auto const uStar = (wave(pStar, pR, rhoR, aR) - wave(pStar, pL, rhoL, aL)) / 2.0;
	Real const ratioL = pStar / pL, ratioR = pStar / pR;
	auto const rhoStarL = rhoL * pow(ratioL, 1 / gamma);
	Real const beta = (gamma - 1) / (gamma + 1);
	auto const rhoStarR = rhoR * (ratioR + beta) / (beta * ratioR + 1);
	auto const aStarL = aL * pow(ratioL, (gamma - 1) / (2 * gamma));
	auto const shock = aR * sqrt((gamma + 1) / (2 * gamma) * ratioR + (gamma - 1) / (2 * gamma));
	auto const middle = (c.mesh.lower + c.mesh.upper) / 2.0;
	Reference result;
	result.name = "Sod exact Riemann solution";
	result.reason = "Sod reference is valid only before the first wave reaches a domain boundary";
	result.validUntil = (c.mesh.upper - c.mesh.lower) / (2.0 * std::max(aL, shock));
	result.evaluate = [=](mesh::PhysicalCoordinates const& position, units::Time time) {
		using std::pow;

		ExactState state;
		auto& q = state.hydro;
		if (time == units::Time{}) {
			q.density() = position[0] < middle ? rhoL : rhoR;
			q.pressure() = position[0] < middle ? pL : pR;
			return state;
		}
		auto const xi = (position[0] - middle) / time;
		if (xi <= -aL) {
			q.density() = rhoL;
			q.pressure() = pL;
		} else if (xi < uStar - aStarL) {
			auto const velocity = 2.0 / (gamma + 1) * (aL + xi);
			auto const sound = 2.0 / (gamma + 1) * (aL - (gamma - 1) / 2.0 * xi);
			q.density() = rhoL * pow(Real(sound / aL), 2 / (gamma - 1));
			q.pressure() = pL * pow(Real(sound / aL), 2 * gamma / (gamma - 1));
			q.velocity(0) = velocity;
		} else if (xi < uStar) {
			q.density() = rhoStarL;
			q.pressure() = pStar;
			q.velocity(0) = uStar;
		} else if (xi < shock) {
			q.density() = rhoStarR;
			q.pressure() = pStar;
			q.velocity(0) = uStar;
		} else {
			q.density() = rhoR;
			q.pressure() = pR;
		}
		return state;
	};
	return result;
}


Reference sphereReference(Config const& c, bool gaussian) {
	using std::erf;
	using std::exp;
	using std::expm1;
	using std::sqrt;

	auto const length = c.mesh.upper - c.mesh.lower;
	auto const center = (c.mesh.lower + c.mesh.upper) / 2.0;
	auto const radius = (gaussian ? 0.5 : 0.25) * length;
	auto const sigma = 0.15 * length;
	auto const density = units::Density::from_value(1e4);
	Real constexpr pi = std::numbers::pi_v<Real>;
	// Shell theorem: phi(r)=-G[M(<r)/r + 4 pi integral_r^R rho(s) s ds].
	// A series avoids cancellation in the enclosed Gaussian mass near r=0.
	auto mass = [=](units::Length r) -> units::Mass {
		using std::erf;
		using std::exp;
		using std::sqrt;

		if (!gaussian) return (4 * pi / 3) * density * r * r * r;
		Real const q = r / sigma, q2 = q * q;
		Real const integral =
			q < 0.01 ? q * q2 * (1.0 / 3 - q2 / 10 + q2 * q2 / 56 - q2 * q2 * q2 / 432) : sqrt(pi / 2) * erf(q / sqrt(2.0)) - q * exp(-q2 / 2);
		return 4 * pi * density * sigma * sigma * sigma * integral;
	};
	Reference result;
	result.name = gaussian ? "Spherically truncated Gaussian gravity" : "Uniform sphere gravity";
	result.evaluate = [=](mesh::PhysicalCoordinates const& position, units::Time) {
		using std::exp;
		using std::expm1;

		ExactState state;
		units::Length r{};
		for (auto const& x : position)
			r = units::hypot(r, x - center);
		auto const enclosed = mass(std::min(r, radius));
		if (r >= radius)
			state.gravity.potential() = -constants::G * enclosed / r;
		else {
			units::Quantity<-1, 1, 0> shells{};
			if (gaussian) {
				Real const q = r / sigma, cutoff = radius / sigma;
				shells = 4 * pi * density * sigma * sigma * exp(-q * q / 2) * (-expm1(-(cutoff * cutoff - q * q) / 2));
			} else
				shells = 2 * pi * density * (radius * radius - r * r);
			state.gravity.potential() = -constants::G * (shells + (r == units::Length{} ? units::Quantity<-1, 1, 0>{} : enclosed / r));
		}
		if (r != units::Length{})
			for (int axis = 0; axis < ndim; ++axis)
				state.gravity.acceleration(axis) = -constants::G * enclosed / (r * r) * ((position[axis] - center) / r);
		return state;
	};
	return result;
}


Reference streamingReference(Config const& c) {
	using std::sqrt;

	Reference result;
	result.name = "Periodic streaming translation";
	if (!c.mesh.periodic) {
		result.reason = "Streaming reference requires periodic boundaries";
		return result;
	}
	auto const length = c.mesh.upper - c.mesh.lower;
	auto const speed = c.radiation.lightSpeedRatio * constants::c / sqrt(Real(ndim));
	result.evaluate = [=](mesh::PhysicalCoordinates const& position, units::Time time) {
		using std::exp;
		using std::round;
		using std::sqrt;

		Real r2 = 0;
		for (auto const& x : position) {
			auto distance = x - (c.mesh.lower + 0.25 * length + speed * time);
			distance -= round(distance / length) * length;
			Real const q = distance / (0.08 * length);
			r2 += q * q;
		}
		ExactState state;
		state.radiation.energy() = units::EnergyDensity::from_value(1e-6 + exp(-0.5 * r2));
		for (int axis = 0; axis < ndim; ++axis)
			state.radiation.radiativeFlux(axis) = constants::c * state.radiation.energy() / sqrt(Real(ndim));
		return state;
	};
	return result;
}

}	 // namespace octotigerII::verification
