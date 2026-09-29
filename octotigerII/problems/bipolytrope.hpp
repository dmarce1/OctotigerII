#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "octotigerII/units/cgs.hpp"

namespace octotigerII::problems {

/// Structural barotrope, independent of the adiabatic index used by hydro.
/// Units may be dimensionless or consistent physical units. Pressure and
/// enthalpy are continuous at a possible downward density jump into the envelope.
class Bipolytrope {
public:
	Bipolytrope(Real coreIndex, Real envelopeIndex, Real interfaceFraction, Real densityJump,
		Real centralDensity, Real centralEnthalpy)
	  : nc_(coreIndex), ne_(envelopeIndex), rc_(interfaceFraction * centralDensity),
		re_(rc_ / densityJump) {
		for (auto v : {coreIndex, envelopeIndex, interfaceFraction, densityJump, centralDensity, centralEnthalpy})
			if (!std::isfinite(v) || v <= 0) throw std::invalid_argument("Bipolytrope parameters must be finite and positive");
		if (interfaceFraction > 1 || densityJump < 1)
			throw std::invalid_argument("Bipolytrope requires interfaceFraction <= 1 and densityJump >= 1");
		if (interfaceFraction == 1 && densityJump > 1)
			throw std::invalid_argument("A nonzero density jump requires a core with interfaceFraction < 1");
		pi_ = centralEnthalpy / ((ne_ + 1) / re_ + (nc_ + 1) / rc_ * (std::pow(1 / interfaceFraction, 1 / nc_) - 1));
		hi_ = (ne_ + 1) * pi_ / re_;
		if (!std::isfinite(pi_) || pi_ <= 0) throw std::invalid_argument("Unrepresentable bipolytropic pressure scale");
	}
	Real density(Real h) const {
		if (h <= 0) return 0;
		return h <= hi_ ? re_ * std::pow(h / hi_, ne_)
			: rc_ * std::pow(1 + (h - hi_) * rc_ / ((nc_ + 1) * pi_), nc_);
	}
	/// Invert the continuous, monotone enthalpy graph, including its interface
	/// plateau. A proximal step avoids jumping back and forth across the density
	/// discontinuity. Its fixed points satisfy H(rho)=h, even for mixed cells.
	Real updateDensity(Real h, Real previous) const {
		Real const target = density(h);
		if (h <= 0 || previous == 0 || rc_ == re_ || (previous <= re_ && target <= re_) || (previous >= rc_ && target >= rc_)) return target;
		// A five-percent interface-enthalpy scale preconditions the graph; it
		// changes the iteration, not the fixed-point barotrope or acceptance test.
		Real const lambda = .05*hi_/(rc_-re_), rhs = h+lambda*previous;
		if (rhs >= hi_+lambda*re_ && rhs <= hi_+lambda*rc_)
			return (rhs-hi_)/lambda;
		Real lo = std::min(previous,target), hi = std::max(previous,target);
		for (int j = 0; j < 48; ++j) {
			Real const mid = (lo+hi)/2;
			if (enthalpy(mid)+lambda*mid < rhs) lo = mid; else hi = mid;
		}
		return (lo+hi)/2;
	}
	Real pressure(Real rho) const {
		if (rho <= 0) return 0;
		if (rho < re_) return pi_ * std::pow(rho / re_, 1 + 1 / ne_);
		if (rho <= rc_) return pi_; // mixed interface cells during iteration/remapping
		return pi_ * std::pow(rho / rc_, 1 + 1 / nc_);
	}
	Real enthalpy(Real rho) const {
		if (rho <= 0) return 0;
		if (rho < re_) return hi_ * std::pow(rho / re_, 1 / ne_);
		if (rho <= rc_) return hi_;
		return hi_ + (nc_ + 1) * pi_ / rc_ * (std::pow(rho / rc_, 1 / nc_) - 1);
	}
	Real coreFraction(Real rho) const {
		if (rho <= re_) return 0;
		if (rho >= rc_) return 1;
		// Volume mixture of the two limiting states at the same interface pressure.
		return rc_ * (rho - re_) / ((rc_ - re_) * rho);
	}
	Real interfacePressure() const { return pi_; }
	Real interfaceEnthalpy() const { return hi_; }
private:
	Real nc_, ne_, rc_, re_, pi_, hi_;
};
} // namespace octotigerII::problems
