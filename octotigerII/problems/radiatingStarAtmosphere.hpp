#pragma once
#include "octotigerII/problems/radiatingStarStructure.hpp"
#include <array>
#include <vector>

namespace octotigerII::problems {

/// Spherical, static full-M1 reference with prescribed positive photon heating.
/// A density-dependent opacity tapers to zero before a transparent, gas-only
/// polytropic envelope. This is an intentionally prescribed-opacity benchmark.
/// Heating is calculated once from the reference, never from an evolved state.
class RadiatingStarAtmosphere {
public:
	struct Parameters {
		RadiatingStarEos::Parameters eos;
		/// kappa_0 rho_c alpha; kappa = kappa_0 max(0,1-rho_cut/rho)^2.
		Real opacityLengthScale = 240;
		Real cutoffDensityFraction = 0.016;
		Real tolerance = 1e-10, maximumRadialStep = 0.002;
		/// Dimensionless radius at which the regular central series starts the ODE.
		Real centerSeriesRadius = 0.001;
		int maximumSteps = 200000;
	};
	struct State {
		units::Density density{};
		units::Temperature temperature{};
		units::Pressure gasPressure{}, radiationPressureRadial{}, radiationPressureTangential{};
		units::EnergyDensity radiationEnergy{};
		units::EnergyFlux radiationFlux{};
		units::Mass enclosedMass{};
		units::VelocitySquared potential{};
		units::Acceleration gravityMagnitude{};
		/// CGS derivatives with respect to physical radius, and external heating.
		Real densityDerivative{}, gasPressureDerivative{}, radiationEnergyDerivative{}, radiationFluxDerivative{};
		Real opacityCgs{}, heatingCgs{}, fluxFactor{};
	};
	struct Diagnostics {
		int radialSteps{};
		Real transitionFluxFactor{}, centralOpticalDepth{}, minimumHeatingCgs{};
		Real envelopeGamma{}, coreMassFraction{};
	};
	explicit RadiatingStarAtmosphere(Parameters parameters);
	State sample(units::Length radius) const;
	RadiatingStarEos const& eos() const { return eos_; }
	Parameters const& parameters() const { return parameters_; }
	Diagnostics const& diagnostics() const { return diagnostics_; }
	units::Length transitionRadius() const;
	units::Length surfaceRadius() const;
	units::Mass mass() const;
	Real luminosityCgs() const;
	/// The frozen opacity profile, independent of any subsequently evolved gas.
	Real opacityCgs(units::Length radius) const { return sample(radius).opacityCgs; }
private:
	struct Node {
		Real radius{}, potential{};
		std::array<Real,3> value{}, derivative{};
	};
	struct Thermodynamics {
		Real temperature{}, energy{}, energyDerivative{}, pressure{}, pressureDerivative{};
		Real energySecondDerivative{}, pressureSecondDerivative{};
	};
	Thermodynamics thermodynamics(Real density) const;
	Real opacity(Real density) const;
	std::array<Real,3> coreRhs(Real radius, std::array<Real,3> const& state) const;
	std::array<Real,3> envelopeRhs(Real radius, std::array<Real,3> const& state) const;
	Real vacuumFluxFactor(Real radius) const;
	std::array<Real,3> interpolate(std::vector<Node> const& nodes, Real radius) const;
	std::array<Real,3> interpolateDerivative(std::vector<Node> const& nodes, Real radius) const;
	Real potential(std::vector<Node> const& nodes, Real radius) const;
	Parameters parameters_;
	RadiatingStarEos eos_;
	Diagnostics diagnostics_;
	Real transitionRadius_{}, surfaceRadius_{}, mass_{}, luminosity_{}, transitionFluxFactor_{};
	Real envelopeGamma_{}, envelopeK_{}, centerDensitySecondDerivative_{}, centerFluxDerivative_{};
	std::vector<Node> core_, envelope_;
};

} // namespace octotigerII::problems
