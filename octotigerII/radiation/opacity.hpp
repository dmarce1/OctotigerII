/** @file
 * @brief Gray thermal and flux opacities for constant and fully ionized gas.
 */
#pragma once
#include "octotigerII/hydro/hydroSystem.hpp"
#include <cmath>
#include <stdexcept>

namespace octotigerII::radiation {
using Opacity = units::Quantity<2, -1, 0>;

struct GrayOpacities {
	Opacity planckAbsorption{};
	Opacity fluxExtinction{};
	Opacity rosselandAbsorption{};
	Opacity scattering{};
};

struct OpacityLaw {
	bool ionizedGas = false;
	Opacity constantAbsorption{}, constantScattering{};
	Real hydrogenFraction = 0.7, metalFraction = 0.02;

	GrayOpacities evaluate(Real density, Real temperature, Real electronFraction = -1) const {
		if (!(density > 0) || !(temperature > 0) || !std::isfinite(density) || !std::isfinite(temperature))
			throw std::invalid_argument("Opacity requires positive finite density and temperature");
		if (!ionizedGas) return {constantAbsorption, constantAbsorption + constantScattering,
			constantAbsorption, constantScattering};
		// Fully ionized, nondegenerate free-free approximation. The factor 37
		// distinguishes the thermal Planck mean from the Rosseland absorption
		// mean; it is not an appropriate fit to bound-free or line opacity.
		Real const rosseland = 3.68e22 * (1 + hydrogenFraction) * (1 - metalFraction)
			* density * std::pow(temperature, -3.5);
		Real const ye = electronFraction >= 0 ? electronFraction : (1 + hydrogenFraction) / 2;
		Real const scattering = 0.4 * ye;
		if (!std::isfinite(37 * rosseland) || !std::isfinite(rosseland + scattering)
			|| scattering < 0 || rosseland < 0)
			throw std::runtime_error("Analytic gray opacity is not representable at this state");
		return {Opacity::from_value(37 * rosseland), Opacity::from_value(rosseland + scattering),
			Opacity::from_value(rosseland), Opacity::from_value(scattering)};
	}

	GrayOpacities evaluate(hydro::ConservedState const& gas, hydro::HydroSystem const& system) const {
		if (!ionizedGas) return {constantAbsorption, constantAbsorption + constantScattering,
			constantAbsorption, constantScattering};
		Real const ye = gas.electrons() > units::Density{} ? Real(gas.electrons() / gas.density()) : -1;
		return evaluate(units::value(gas.density()), units::value(system.temperature(gas)), ye);
	}
};

inline OpacityLaw opacityLaw(Config const& config, Opacity prescribedAbsorption) {
	return {config.radiation.opacityModel == "ionized-gas", prescribedAbsorption,
		Opacity::from_value(config.radiation.scatteringOpacity),
		config.radiation.hydrogenFraction, config.radiation.metalFraction};
}

inline bool radiationCouplingEnabled(Config const& config) {
	return config.radiation.opacityModel == "ionized-gas" || config.radiation.opacity > 0 || config.radiation.scatteringOpacity > 0;
}
} // namespace octotigerII::radiation
