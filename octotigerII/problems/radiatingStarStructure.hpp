#pragma once
#include "octotigerII/units/constants.hpp"
#include <array>
#include <vector>

namespace octotigerII::problems {

/// Gas plus LTE radiation on a prescribed, nonisentropic structural polytrope.
/// The evolution gas equation of state is still an ideal gas with gamma=5/3.
class RadiatingStarEos {
public:
	struct Parameters {
		Real index = 3.5;
		units::Density centralDensity = units::Density::from_value(1);
		Real centralBeta = 0.8, meanMolecularWeight = 0.6;
	};
	struct State {
		units::Density density{};
		units::Pressure pressure{}, gasPressure{};
		units::EnergyDensity radiationEnergy{};
		units::Temperature temperature{};
		/// Integral dP/rho, not the thermodynamic specific enthalpy.
		units::VelocitySquared integralH{};
		Real beta{}, radiationPressureDerivative{}, gamma1{}, nabla{}, nablaAd{};
		/// db/dH in s^2/cm^2 and specific entropy differences in erg/(g K).
		Real radiationPressureSecondDerivative{}, entropyDifference{}, entropyDerivativeLogDensity{};
	};
	explicit RadiatingStarEos(Parameters parameters);
	State atDensity(units::Density density) const;
	State atIntegralH(units::VelocitySquared h) const;
	Parameters const& parameters() const { return parameters_; }
	units::Pressure centralPressure() const;
	units::Temperature centralTemperature() const;
	units::VelocitySquared centralIntegralH() const;
	units::Length scaleLength() const;
private:
	Parameters parameters_;
	Real gasConstant_{}, centralTemperature_{}, centralPressure_{}, centralH_{}, scaleLength_{};
};

/// Isolated axisymmetric, uniformly rotating diffusion-interior SCF reference.
/// Construction performs the Poisson iteration. Subsequent queries are const,
/// thread safe, and O(maxMultipole), with no mutable caches or synchronization.
/// Its formal surface is not a physical radiation atmosphere.
class RadiatingStarStructure {
public:
	struct Parameters {
		RadiatingStarEos::Parameters eos;
		/// Omega / sqrt(G M0/R0^3), with M0,R0 from the nonrotating polytrope.
		Real spinFractionOfSphericalBreakup = 0.2;
		int radialCells = 256, angularPoints = 32, maxMultipole = 12, maxIterations = 400;
		Real tolerance = 1e-10, relaxation = 0.5;
	};
	struct State {
		RadiatingStarEos::State thermodynamics;
		units::VelocitySquared potential{};
		/// Cylindrical components (R,z) of grad Phi and grad Psi, respectively.
		std::array<units::Acceleration, 2> potentialGradient{}, effectivePotentialGradient{};
	};
	struct Diffusion {
		std::array<units::EnergyFlux, 2> flux{};
		/// Required external photon heating in erg/(cm^3 s), at leading order v/c.
		Real heating{};
		Real fluxFactor{};
	};
	struct Diagnostics {
		int iterations{};
		Real densityResidual{}, bernoulliResidual{};
		/// Scalar virial residual (2 Trot + W + 3 integral P)/|W|.
		Real virialResidual{};
	};
	explicit RadiatingStarStructure(Parameters parameters);
	State sample(units::Length cylindricalRadius, units::Length z) const;
	/// Constant-opacity diffusion formula only; does not impose a flux limiter
	/// or extend the interior solution through its transport atmosphere.
	Diffusion diffusion(State const& state, Real opacityCgs) const;
	RadiatingStarEos const& eos() const { return eos_; }
	Parameters const& parameters() const { return parameters_; }
	Diagnostics const& diagnostics() const { return diagnostics_; }
	units::Mass mass() const;
	units::Length equatorialRadius() const;
	units::Length polarRadius() const;
	units::Length sphericalRadius() const;
	units::InverseTime angularVelocity() const;
	units::VelocitySquared bernoulliConstant() const;
private:
	struct Potential { Real value{}, radialDerivative{}, muDerivative{}; };
	Potential potential(Real dimensionlessRadius, Real mu) const;
	void poisson(std::vector<Real> const& density);
	Real surface(Real mu) const;
	Parameters parameters_;
	RadiatingStarEos eos_;
	Diagnostics diagnostics_;
	Real step_{}, outerRadius_{}, omegaSquared_{}, centralPotential_{}, massIntegral_{};
	Real equatorialRadius_{}, polarRadius_{}, sphericalRadius_{};
	std::vector<Real> mu_, angularWeights_, legendre_, coefficients_, coefficientDerivatives_, exteriorMoments_;
};

} // namespace octotigerII::problems
