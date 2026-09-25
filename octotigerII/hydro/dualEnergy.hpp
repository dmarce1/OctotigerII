// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#pragma once
#include <cmath>
#include <stdexcept>
#include "octotigerII/math/Real.hpp"

namespace octotigerII::hydro {

/// Controls A = rho [(u/u0)/(rho/rho0)^gamma]^exponent, with fixed CGS
/// references rho0=1 g/cm^3 and u0=1 erg/cm^3. A therefore has density units.
class DualEnergyOptions {
public:

	bool enabled = true;
	Real exponent = 1;
	Real pressureThreshold = 0.001;
	Real syncThreshold = 0.1;

	void validate() const {
		if (!std::isfinite(exponent) || exponent == 0)
			throw std::invalid_argument("hydro.dualEnergy.exponent must be finite and nonzero");
		if (!std::isfinite(pressureThreshold) || !std::isfinite(syncThreshold) ||
			!(pressureThreshold >= 0 && pressureThreshold < syncThreshold && syncThreshold < 1))
			throw std::invalid_argument("Dual energy requires 0 <= pressureThreshold < syncThreshold < 1");
	}

	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & enabled & exponent & pressureThreshold & syncThreshold;
	}
};

} // namespace octotigerII::hydro
