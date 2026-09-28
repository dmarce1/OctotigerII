// Deferred research prototype: not part of the build and not application-tested.
// See docs/research/radiating-star-2026-09/README.md before reviving this reader.
#pragma once
#include "octotigerII/problems/radiatingStarAtmosphere.hpp"
#include "octotigerII/problems/radiatingStarGasBarotrope.hpp"
#include <array>
#include <string_view>
#include <vector>

namespace octotigerII::problems {

/// Immutable reader for an offline, axisymmetric reference about the spherical
/// full-M1 solution. Separate core/exterior interpolants preserve the opacity
/// interface. This evaluator supplies LTE, zeroth-order meridional fields and
/// the supplied azimuthal radiation flux; it does not by itself impose the
/// moving-gas thermal balance Q-v.G=0 or certify the table's steady PDE residual.
class RadiatingStarReference {
public:
	struct Parameters {
		units::Density centralDensity = units::Density::from_value(1);
		Real meanMolecularWeight = 0.6;
		Real atmosphereTolerance = 1e-9;
	};
	struct Metadata {
		Real index{}, centralBeta{}, opacity{}, cutoffDensityFraction{};
		Real spin{}, extent{}, lightSpeed{};
	};
	struct State {
		units::Density density{};
		units::Temperature temperature{};
		units::Pressure gasPressure{};
		units::EnergyDensity radiationEnergy{};
		units::VelocitySquared integralH{}, potential{};
		units::EnergyFlux radialFlux{}, thetaFlux{}, cylindricalFlux{}, verticalFlux{}, azimuthalFlux{};
		units::InverseTime angularVelocity{};
		std::array<units::Acceleration,2> gravity{}, potentialGradient{}, integralHGradient{};
		std::array<units::Quantity<-4,1,0>,2> densityGradient{};
		std::array<units::Quantity<-1,0,0,1>,2> temperatureGradient{};
		std::array<units::Quantity<-2,1,-2>,2> gasPressureGradient{}, radiationEnergyGradient{};
		/// Spherical derivatives use physical r and signed mu=z/r. q is defined
		/// by F_theta=-c E sqrt(1-mu^2) q and is odd across the equator.
		Real mu{}, radialFluxFactor{}, polarFluxFactor{};
		units::Quantity<-1,0,0> radialFactorRadialDerivative{}, polarFactorRadialDerivative{};
		Real radialFactorMuDerivative{}, polarFactorMuDerivative{};
		/// Regular azimuthal coefficient W=F_phi/(c p_unit sqrt(1-mu^2)).
		Real azimuthalFluxCoefficient{}, azimuthalCoefficientMuDerivative{};
		units::Quantity<-1,0,0> azimuthalCoefficientRadialDerivative{};
		units::Quantity<-2,1,-2> energyRadialDerivative{};
		units::EnergyDensity energyMuDerivative{};
		/// Frozen opacity and div(F_p), in cm^2/g and erg/(cm^3 s).
		/// The transparent exterior has the explicitly source-free prescription.
		Real opacityCgs{}, photonHeatingCgs{};
	};

	/// The text format is versioned, whitespace separated, and locale independent:
	/// OCTOII_RADIATING_STAR 1; seven metadata numbers; coreRows exteriorRows
	/// nMu 6; mu axis; core r axis; exterior r axis; then, for each grid, nodes
	/// in radial-major order and six fields (value,dr,dmu,drdmu) per node.
	/// Fields are delta h_g, delta ln E, delta f_r, delta q, delta Phi, W.
	/// W=F_phi/(c p_unit sqrt(1-mu^2)); lengths use 4*pi*G=rho_c=R_gas*T_c=1.
	RadiatingStarReference(Parameters parameters, std::string_view text);
	State sample(units::Length cylindricalRadius, units::Length z) const;
	Metadata const& metadata() const { return metadata_; }
	Parameters const& parameters() const { return parameters_; }
	RadiatingStarAtmosphere const& baseline() const { return baseline_; }
	RadiatingStarGasBarotrope const& gasBarotrope() const { return barotrope_; }
	RadiatingStarEos const& eos() const { return baseline_.eos(); }
	units::Length lengthScale() const;
	units::Pressure pressureScale() const;
	units::VelocitySquared velocitySquaredScale() const;
	units::Length maximumRadius() const;
	units::Length transitionRadius() const;

private:
	static constexpr std::size_t fieldCount = 6;
	struct Jet { Real value{}, radial{}, angular{}, mixed{}, radialSecond{}; };
	using Node = std::array<Jet,fieldCount>;
	struct Grid {
		std::vector<Real> radius;
		std::vector<Node> nodes;
	};
	struct Parsed {
		Metadata metadata;
		std::vector<Real> mu;
		Grid core, exterior;
	};
	static Parsed parse(std::string_view text);
	RadiatingStarReference(Parameters parameters, Parsed parsed);
	std::array<Jet,fieldCount> interpolate(Grid const& grid, Real radius, Real mu) const;
	Parameters parameters_;
	Metadata metadata_;
	RadiatingStarAtmosphere baseline_;
	RadiatingStarGasBarotrope barotrope_;
	std::vector<Real> mu_;
	Grid core_, exterior_;
	Real lengthScale_{}, pressureScale_{}, hScale_{}, surfacePotential_{};
};

} // namespace octotigerII::problems
