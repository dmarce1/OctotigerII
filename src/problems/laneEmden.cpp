#include "octotigerII/problems/laneEmden.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <mutex>
#include <shared_mutex>
#include <numbers>
#include <stdexcept>
#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/synchronization/shared_mutex.hpp>
#endif

namespace octotigerII::problems {
namespace {
	using Value = LaneEmden::Value;
	Value derivative(Real xi, Value y, Real n) {
		using std::pow;
		return {-y.mass / (xi * xi), xi * xi * pow(std::max(Real(0), y.theta), n)};
	}
	Value add(Value a, Value b, Real h) { return {a.theta + h * b.theta, a.mass + h * b.mass}; }
	Value step(Real xi, Value y, Real h, Real n) {
		auto const a = derivative(xi, y, n);
		auto const b = derivative(xi + h / 2, add(y, a, h / 2), n);
		auto const c = derivative(xi + h / 2, add(y, b, h / 2), n);
		auto const d = derivative(xi + h, add(y, c, h), n);
		return {y.theta + h / 6 * (a.theta + 2 * b.theta + 2 * c.theta + d.theta),
			y.mass + h / 6 * (a.mass + 2 * b.mass + 2 * c.mass + d.mass)};
	}
	LaneEmden const& solution(Real n) {
		if (!std::isfinite(n) || !(n > 0 && n < 5)) throw std::invalid_argument("Lane-Emden requires 0 < n < 5");
#ifdef OCTOTIGERII_WITH_HPX
		static hpx::shared_mutex mutex;
#else
		static std::shared_mutex mutex;
#endif
		static std::map<Real, std::unique_ptr<LaneEmden>> cache;
		{
			std::shared_lock lock(mutex);
			if (auto it = cache.find(n); it != cache.end()) return *it->second;
		}
		auto value = std::make_unique<LaneEmden>(n);
		std::unique_lock lock(mutex);
		auto const [it, inserted] = cache.try_emplace(n, std::move(value));
		return *it->second;
	}
}

LaneEmden::LaneEmden(Real n) : index_(n) {
	using std::isfinite;
	if (!isfinite(n) || !(n > 0 && n < 5)) throw std::invalid_argument("Lane-Emden requires 0 < n < 5");
	points_.push_back({0, {1, 0}});
	// Regular central series avoids the removable singularity at xi=0.
	Real xi = 1e-5;
	Value y{1 - xi * xi / 6 + n * xi * xi * xi * xi / 120,
		xi * xi * xi / 3 - n * xi * xi * xi * xi * xi / 30};
	points_.push_back({xi, y});
	for (int count = 0; count < 100000; ++count) {
		// Resolve the regular center before taking the normal table steps.
		Real const h = std::min(5e-4 * std::max(Real(1), xi), 0.05 * xi);
		auto next = step(xi, y, h, n);
		if (next.theta <= 0) {
			Real lower = 0, upper = h;
			for (int i = 0; i < 60; ++i) {
				Real const middle = (lower + upper) / 2;
				if (step(xi, y, middle, n).theta > 0) lower = middle;
				else upper = middle;
			}
			Real const offset = (lower + upper) / 2;
			next = step(xi, y, offset, n);
			next.theta = 0;
			points_.push_back({xi + offset, next});
			return;
		}
		xi += h;
		y = next;
		points_.push_back({xi, y});
	}
	throw std::runtime_error("Lane-Emden surface was not reached");
}

LaneEmden::Value LaneEmden::operator()(Real xi) const {
	using std::isfinite;
	if (!isfinite(xi) || xi < 0) throw std::invalid_argument("Invalid Lane-Emden radius");
	if (xi == 0) return {1, 0};
	if (xi >= surface()) return {0, surfaceMass()};
	auto const upper = std::upper_bound(points_.begin(), points_.end(), xi, [](Real x, Point const& p) { return x < p.xi; });
	auto const& a = *(upper - 1);
	auto const& b = *upper;
	Real const h = b.xi - a.xi, t = (xi - a.xi) / h;
	auto const da = a.xi == 0 ? Value{0, 0} : derivative(a.xi, a.value, index_);
	auto const db = derivative(b.xi, b.value, index_);
	auto interpolate = [&](Real u, Real v, Real du, Real dv) {
		return (2*t*t*t - 3*t*t + 1)*u + (t*t*t - 2*t*t + t)*h*du + (-2*t*t*t + 3*t*t)*v + (t*t*t - t*t)*h*dv;
	};
	return {std::max(Real(0), interpolate(a.value.theta, b.value.theta, da.theta, db.theta)),
		std::clamp(interpolate(a.value.mass, b.value.mass, da.mass, db.mass), a.value.mass, b.value.mass)};
}

Polytrope::Polytrope(Real n, units::Length radius, units::Density density)
  : solution_(&solution(n)), radius_(radius), scale_(radius / solution_->surface()), centralDensity_(density) {
	if (!(radius > units::Length{}) || !units::finite(radius) || !(density > units::Density{}) || !units::finite(density))
		throw std::invalid_argument("Polytrope radius and central density must be positive and finite");
}

units::Mass Polytrope::mass() const {
	return 4 * std::numbers::pi_v<Real> * scale_ * scale_ * scale_ * centralDensity_ * solution_->surfaceMass();
}

units::Pressure Polytrope::centralPressure() const {
	return 4 * std::numbers::pi_v<Real> * constants::G * scale_ * scale_ * centralDensity_ * centralDensity_ / (solution_->index() + 1);
}

Polytrope::State Polytrope::operator()(units::Length r) const {
	using std::pow;
	if (!(r >= units::Length{}) || !units::finite(r)) throw std::invalid_argument("Invalid polytrope radius");
	if (r >= radius_) return {{}, {}, mass(), -constants::G * mass() / r};
	auto const value = (*solution_)(Real(r / scale_));
	return {centralDensity_ * pow(value.theta, solution_->index()), centralPressure() * pow(value.theta, solution_->index() + 1),
		4 * std::numbers::pi_v<Real> * scale_ * scale_ * scale_ * centralDensity_ * value.mass,
		-constants::G * mass() / radius_ - (solution_->index() + 1) * centralPressure() / centralDensity_ * value.theta};
}

} // namespace octotigerII::problems
