/** @file
 * @brief Zero-temperature electron degeneracy for carbon/oxygen white dwarfs.
 */
#pragma once
#include "octotigerII/units/constants.hpp"
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace octotigerII::hydro {

/// Cold, ideal Fermi electrons; ions supply mass but no degeneracy pressure.
/// Thermal-ion pressure and energy are added separately by HydroSystem.
class WhiteDwarfEos {
public:
	explicit WhiteDwarfEos(Real meanMassPerElectron = 2) : muElectron_(meanMassPerElectron) {
		if (!(muElectron_ > 0) || !std::isfinite(muElectron_))
			throw std::invalid_argument("White-dwarf electron molecular weight must be finite and positive");
	}

	Real meanMassPerElectron() const { return muElectron_; }

	units::Pressure pressure(units::Density density) const {
		Real const x = fermiMomentum(density);
		if (x == 0) return {};
		// The closed expression loses significant digits at small x.
		Real integral;
		if (x < Real(0.5)) {
			Real const x2 = x*x;
			Real coefficient = 1, power = x*x*x*x*x;
			integral = power/5;
			for (int k=1; k<=32; ++k) {
				coefficient *= -Real(2*k-1)/Real(2*k);
				power *= x2;
				Real const term = coefficient*power/Real(5+2*k);
				integral += term;
				if (std::abs(term) < 1e-16*std::abs(integral)) break;
			}
		} else {
			integral = (x*(2*x*x-3)*std::sqrt(1+x*x)+3*std::asinh(x))/8;
		}
		return units::Pressure::from_value(pressureScale()*integral);
	}

	units::VelocitySquared enthalpy(units::Density density) const {
		Real const x = fermiMomentum(density), root = std::hypot(Real(1),x);
		return units::VelocitySquared::from_value(enthalpyScale()*x*x/(root+1));
	}

	units::EnergyDensity internalEnergy(units::Density density) const {
		return density*enthalpy(density)-pressure(density);
	}

	units::VelocitySquared pressureDerivative(units::Density density) const {
		Real const x = fermiMomentum(density);
		return units::VelocitySquared::from_value(enthalpyScale()*x*x/(3*std::hypot(Real(1),x)));
	}

	units::Density densityFromEnthalpy(units::VelocitySquared value) const {
		if (!(value >= units::VelocitySquared{}) || !units::finite(value))
			throw std::invalid_argument("White-dwarf enthalpy must be finite and nonnegative");
		Real const q=units::value(value)/enthalpyScale();
		Real const x2=q*(q+2);
		return units::Density::from_value(densityScale()*std::pow(x2,Real(1.5)));
	}

private:
	static constexpr Real electronMass = 9.1093837139e-28; // g, CODATA 2022
	Real muElectron_;

	Real densityScale() const {
		Real const c=units::value(constants::c);
		Real const hbar=units::value(constants::planck)/(2*std::numbers::pi_v<Real>);
		return muElectron_*units::value(constants::atomicMassUnit)*
			std::pow(electronMass*c/hbar,3)/(3*std::numbers::pi_v<Real>*std::numbers::pi_v<Real>);
	}
	Real enthalpyScale() const {
		Real const c=units::value(constants::c);
		return electronMass*c*c/(muElectron_*units::value(constants::atomicMassUnit));
	}
	Real pressureScale() const {
		Real const c=units::value(constants::c);
		Real const hbar=units::value(constants::planck)/(2*std::numbers::pi_v<Real>);
		return std::pow(electronMass,4)*std::pow(c,5)/
			(3*std::numbers::pi_v<Real>*std::numbers::pi_v<Real>*std::pow(hbar,3));
	}
	Real fermiMomentum(units::Density density) const {
		if (!(density >= units::Density{}) || !units::finite(density))
			throw std::invalid_argument("White-dwarf density must be finite and nonnegative");
		return std::cbrt(units::value(density)/densityScale());
	}
};
} // namespace octotigerII::hydro
