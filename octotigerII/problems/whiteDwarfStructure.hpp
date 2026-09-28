/** @file
 * @brief Uniformly rotating carbon/oxygen white-dwarf reference.
 */
#pragma once
#include "octotigerII/hydro/whiteDwarfEos.hpp"
#include <array>
#include <vector>

namespace octotigerII::problems {

/// Cold-electron/thermal-ion barotrope with an axisymmetric Poisson SCF solve.
/// The ion pressure is a small fixed fraction of degeneracy pressure in the
/// reference only; subsequent evolution uses a separate ideal-ion energy.
class WhiteDwarfStructure {
public:
	struct Parameters {
		units::Density centralDensity = units::Density::from_value(1e9);
		Real meanMassPerElectron = 2, meanMassPerIon = 14;
		Real thermalPressureFraction = 1e-4;
		Real spinFractionOfSphericalBreakup = 0.2;
		int radialCells = 256, angularPoints = 32, maxMultipole = 12, maxIterations = 400;
		Real tolerance = 1e-9, relaxation = 0.5;
	};
	struct State {
		units::Density density{};
		units::Pressure pressure{}, ionPressure{};
		units::EnergyDensity radiationEnergy{};
		units::Temperature temperature{};
		units::VelocitySquared potential{};
		std::array<units::Acceleration,2> potentialGradient{}, effectivePotentialGradient{};
	};
	struct Diagnostics {
		int iterations{};
		Real densityResidual{}, bernoulliResidual{}, virialResidual{};
	};
	explicit WhiteDwarfStructure(Parameters);
	State sample(units::Length cylindricalRadius, units::Length z) const;
	std::array<units::EnergyFlux,2> diffusionFlux(State const&, Real opacityCgs) const;
	Parameters const& parameters() const { return parameters_; }
	Diagnostics const& diagnostics() const { return diagnostics_; }
	units::Length scaleLength() const { return scaleLength_; }
	units::Length equatorialRadius() const { return equatorialRadius_*scaleLength_; }
	units::Length polarRadius() const { return polarRadius_*scaleLength_; }
	units::Length sphericalRadius() const { return sphericalRadius_*scaleLength_; }
	units::Mass mass() const;
	units::InverseTime angularVelocity() const;
	units::Velocity centralSoundSpeed() const;
private:
	struct Potential { Real value{}, radialDerivative{}, muDerivative{}; };
	Real densityRatio(Real normalizedEnthalpy) const;
	void poisson(std::vector<Real> const& density);
	Potential potential(Real radius, Real mu) const;
	Real surface(Real mu) const;
	Parameters parameters_;
	hydro::WhiteDwarfEos cold_;
	Diagnostics diagnostics_;
	units::VelocitySquared centralEnthalpy_{};
	units::Length scaleLength_{};
	Real step_{}, outerRadius_{}, omegaSquared_{}, centralPotential_{}, massIntegral_{};
	Real equatorialRadius_{}, polarRadius_{}, sphericalRadius_{}, sphericalMass_{};
	std::vector<Real> mu_, angularWeights_, legendre_, coefficients_, coefficientDerivatives_, exteriorMoments_;
};
} // namespace octotigerII::problems
