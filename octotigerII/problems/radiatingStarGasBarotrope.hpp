#pragma once
#include "octotigerII/problems/radiatingStarStructure.hpp"
#include <vector>

namespace octotigerII::problems {

/// Gas pressure on the prescribed gas+radiation structural EOS, joined to a
/// transparent gas polytrope. integralH is integral(dp_g/rho), with its zero at
/// vacuum; it is not thermodynamic enthalpy or an assumption of constant entropy.
/// Construction builds an immutable monotone table; const queries need no locks.
class RadiatingStarGasBarotrope {
public:
	struct Parameters {
		RadiatingStarEos::Parameters eos;
		Real cutoffDensityFraction = 0.015;
	};
	struct State {
		units::Density density{};
		units::Temperature temperature{};
		units::Pressure gasPressure{};
		units::VelocitySquared integralH{};
		/// d rho / d integralH, in g s^2 / cm^5.
		units::Quantity<-5,1,2> densityDerivative{};
	};
	explicit RadiatingStarGasBarotrope(Parameters parameters);
	/// Nonpositive finite h gives vacuum. Values above maximumIntegralH throw.
	State atIntegralH(units::VelocitySquared h) const;
	/// Supports 0 <= rho <= 4 rho_c; rejects overshoots rather than extrapolating.
	State atDensity(units::Density density) const;
	Parameters const& parameters() const { return parameters_; }
	RadiatingStarEos const& eos() const { return eos_; }
	units::Density cutoffDensity() const;
	units::Density maximumDensity() const;
	units::VelocitySquared cutoffIntegralH() const;
	units::VelocitySquared centralIntegralH() const;
	units::VelocitySquared maximumIntegralH() const;
	Real envelopeGamma() const { return envelopeGamma_; }
private:
	struct Node { Real x{}, h{}, derivative{}; };
	struct Interpolation { Real h{}, derivative{}; };
	Real coreDerivative(Real logDensityRatio) const;
	Interpolation interpolate(std::size_t cell, Real logDensityRatio) const;
	State coreState(units::Density density, Real h, Real derivative) const;
	Parameters parameters_;
	RadiatingStarEos eos_;
	Real hScale_{}, gasConstant_{}, cutoffH_{}, cutoffPressure_{}, envelopeGamma_{};
	std::vector<Node> table_;
};

} // namespace octotigerII::problems
