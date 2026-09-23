#include "octotigerII/gravity/images.hpp"
#include <algorithm>
#include <limits>

namespace octotigerII::gravity {
ImageGeometry::ImageGeometry(physics::BoundaryConditions const& bc) {
	validateBoundaries(bc);
	using Rule = physics::BoundaryCondition;
	for (int axis = 0; axis < 3; ++axis) {
		if (bc.periodic(axis)) periods_[axis] = 1;
		bool const low = bc.lower[axis] == Rule::Reflecting;
		bool const high = bc.upper[axis] == Rule::Reflecting;
		if (low && high) periods_[axis] = 2;
		if (low || high) {
			auto const size = images_.size();
			for (std::size_t i = 0; i < size; ++i) {
				auto image = images_[i];
				image.mask |= 1u << axis;
				image.translation[axis] = low ? 0 : 2;
				images_.push_back(image);
			}
		}
	}
}

diagonal::Offset ImageGeometry::periods(int count) const {
	auto result = periods_;
	for (auto& x : result)
		x *= count;
	return result;
}

diagonal::Offset ImageGeometry::separation(diagonal::Offset a, diagonal::Offset b, int count, SourceImage const& image) const {
	for (int d = 0; d < 3; ++d) {
		if (image.mask & (1u << d)) b[d] = -b[d] - 1 + image.translation[d] * count;
		a[d] -= b[d];
		int const period = periods_[d] * count;
		// Keep the sign at a half-period tie: antisymmetry of the pair split.
		if (period) {
			a[d] %= period;
			if (2 * a[d] > period) a[d] -= period;
			if (2 * a[d] < -period) a[d] += period;
		}
	}
	return a;
}

double ImageGeometry::nextImageDistance(diagonal::Offset r, int count) const {
	using std::abs;
	using std::sqrt;
	double const radius2 = double(r[0]) * r[0] + double(r[1]) * r[1] + double(r[2]) * r[2];
	double next2 = std::numeric_limits<double>::infinity();
	for (int d = 0; d < 3; ++d)
		if (periods_[d]) {
			double const p = double(periods_[d]) * count;
			next2 = std::min(next2, radius2 + p * p - 2 * p * abs(r[d]));
		}
	return sqrt(next2);
}

bool ImageGeometry::acceptable(diagonal::Offset r, int count, double theta, bool correction) const {
	using std::abs;
	using std::hypot;
	for (int d = 0; d < 3; ++d)
		if (periods_[d] && 2 * (abs(r[d]) + 1) > periods_[d] * count) return false;
	double const distance = correction ? nextImageDistance(r, count) : hypot(double(r[0]), double(r[1]), double(r[2]));
	return distance * theta > 1;
}

diagonal::Coefficients imageShiftMultipole(diagonal::Coefficients const& a, diagonal::Vector const& s, double ratio, int p) {
	int const count = diagonal::coefficientCount(p);
	bool const extended = a.size() == std::size_t(count + 1);
	diagonal::Coefficients base(a.begin(), a.begin() + (extended ? count : a.size()));
	auto b = diagonal::shiftMultipole(base, s, ratio, p);
	double const mz = base.size() > 1 ? base[3] : 0;
	b.push_back((extended ? a.back() * ratio * ratio : 0) + ratio * s[2] * mz + 0.5 * s[2] * s[2] * a[0]);
	return b;
}

diagonal::Coefficients imageShiftLocal(diagonal::Coefficients const& a, diagonal::Vector const& s, double ratio, int p) {
	int const count = diagonal::coefficientCount(p);
	diagonal::Coefficients base(a.begin(), a.begin() + count);
	auto b = diagonal::shiftLocal(base, s, ratio, p);
	double const laplace = a.size() > std::size_t(count) ? a.back() : 0;
	b[0] += 0.5 * laplace * s[2] * s[2];
	b[3] += ratio * laplace * s[2];
	b.push_back(laplace * ratio * ratio);
	return b;
}

diagonal::Coefficients reflectMultipole(diagonal::Coefficients const& a, unsigned mask, int p) {
	auto result = a;
	if (!mask || a.size() == 1) return result;
	auto const& ix = diagonal::indices(p);
	for (std::size_t i = 0; i < ix.size(); ++i) {
		int const parity = ((mask & 1) ? ix[i].x : 0) + ((mask & 2) ? ix[i].y : 0) + ((mask & 4) ? ix[i].z : 0);
		if (parity & 1) result[i] = -result[i];
	}
	return result;
}
}	 // namespace octotigerII::gravity
