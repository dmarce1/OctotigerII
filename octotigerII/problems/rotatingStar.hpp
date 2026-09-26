/** @file
 * @brief Original Octo-Tiger oblate rotating-star SCF equilibrium, scaled to CGS.
 */
// Distributed under the Boost Software License, Version 1.0.
#pragma once
#include "octotigerII/hydro/hydroSystem.hpp"

namespace octotigerII::problems {

/// A fixed n=3/2 SCF model with independent physical length and density units.
/// The original table has G=1, polar/equatorial radii approximately 0.605/0.905,
/// and uniform spin. No files are needed at runtime.
class RotatingStar {
public:

	struct Profile {
		Real density{}, internalEnergy{};
	};

	static constexpr Real spin = 0.5155532816213834;
	static constexpr Real gamma = Real(5) / 3;
	static constexpr int tableCells = 100;

	RotatingStar(units::Length lengthUnit, units::Density densityUnit, Real atmosphereFraction = 1e-10);

	/// Dimensionless original SCF profile; returns zero outside table coverage.
	static Profile profile(Real cylindricalRadius, Real z);

	/// Inertial gas state at a position relative to the stellar center.
	hydro::PrimitiveState operator()(mesh::PhysicalCoordinates const& relativePosition) const;
	units::InverseTime angularVelocity() const;
	units::EnergyDensity energyUnit() const;
	units::Length coreLength() const { return 0.2 * lengthUnit_; }

private:
	units::Length lengthUnit_;
	units::Density densityUnit_;
	Real atmosphereFraction_;
};

} // namespace octotigerII::problems
