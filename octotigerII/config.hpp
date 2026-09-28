/** @file
 * @brief Validated CGS application options and command-line parsing.
 * @ingroup runtime
 */
#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <string>
#include <vector>
#include "octotigerII/buildConfig.hpp"
#include "octotigerII/hydro/dualEnergy.hpp"
#include "octotigerII/composition/species.hpp"
#include "octotigerII/physics/boundary.hpp"
#include "octotigerII/units/constants.hpp"

namespace octotigerII {

/// Problem and runtime options. Lengths and times are physical CGS quantities.
/// Dimension and available physics are build choices; problem identity is a runtime choice.
/// @ingroup runtime
class Config {
public:
	composition::Options massFractions;
	std::string problem;
	std::int64_t randomSeed = 5489;

	class MeshOptions {
	public:
		int cells = 16, level = 1;
		units::Length lower{}, upper = units::Length::from_value(1);
		physics::BoundaryConditions boundary;

		template <typename Archive>
		void serialize(Archive& archive, unsigned) {
			archive & cells & level & lower & upper & boundary;
		}
	} mesh;

	class FrameOptions {
	public:
		/// Constant angular velocity about the inertial z axis, in radians/s.
		units::InverseTime omega{};
		template <typename Archive>
		void serialize(Archive& archive, unsigned) { archive & omega; }
	} frame;

	class AmrOptions {
	public:
		bool enabled = false, hydro = true, radiation = true;
		int minLevel = -1, maxLevel = 6, regridEvery = 4, bufferCells = 1;
		units::Mass maxCellMass{};
		units::Density refineDensity{};
		Real shadowTolerance = 0.05, shadowFloor = 1e-8, coarsenFactor = 0.25, signalBuffer = 1;

		template <typename Archive>
		void serialize(Archive& archive, unsigned) {
			archive & enabled & hydro & radiation & minLevel & maxLevel & regridEvery & bufferCells;
			archive & maxCellMass & refineDensity & shadowTolerance & shadowFloor & coarsenFactor & signalBuffer;
		}
	} amr;

	class RuntimeOptions {
	public:
		units::Time stopTime = units::Time::from_value(0.2);
		int maxSteps = 100000, workerTasks = 0;
		bool workStealing = true;

		template <typename Archive>
		void serialize(Archive& archive, unsigned) {
			archive & stopTime & maxSteps & workerTasks & workStealing;
		}
	} runtime;

	class TimestepOptions {
	public:
		Real cfl = 0.4;
		bool refinement = true;

		template <typename Archive>
		void serialize(Archive& archive, unsigned) {
			archive & cfl & refinement;
		}
	} timestep;

	class HydroOptions {
	public:
		Real gamma = 1.4;
		Real meanMolecularWeight = 1;
		std::string eos = "ideal";
		Real meanMassPerElectron = 2;
		std::string helmholtzTable;
		Real temperatureFloor = 1000;
		hydro::DualEnergyOptions dualEnergy;
		std::array<units::Acceleration, ndim> acceleration{};

		template <typename Archive>
		void serialize(Archive& archive, unsigned) {
			archive & gamma & meanMolecularWeight & eos & meanMassPerElectron & dualEnergy & helmholtzTable & temperatureFloor;
			for (auto& component : acceleration)
				archive & component;
		}
	} hydro;

	class RayleighTaylorOptions {
	public:
		units::Density densityLower = units::Density::from_value(1), densityUpper = units::Density::from_value(2);
		units::Pressure interfacePressure = units::Pressure::from_value(2.5);
		units::Velocity perturbation = units::Velocity::from_value(0.01);

		template <typename Archive>
		void serialize(Archive& archive, unsigned) {
			archive & densityLower & densityUpper & interfacePressure & perturbation;
		}
	} rayleighTaylor;

	class StarOptions {
	public:
		units::Length radius = units::Length::from_value(1e9);
		std::array<units::Length, ndim> center{};
		units::Density centralDensity = units::Density::from_value(1e6);
		Real polytropicIndex = 1.5, atmosphereFraction = 1e-8;
		template <typename Archive>
		void serialize(Archive& archive, unsigned) {
			archive & radius & centralDensity & polytropicIndex & atmosphereFraction;
			for (auto& x : center) archive & x;
		}
	} star;

	/// Controls of the constructed gas+radiation stellar reference. The gas EOS
	/// remains gamma=5/3; star.polytropicIndex describes only the structure.
	class RadiatingStarOptions {
	public:
		Real centralGasFraction = 0.8;
		/// Omega / sqrt(G M_spherical / R_spherical^3).
		Real rotationFraction = 0.2;
		/// kappa*rho_c*alpha when radiation.opacity is not explicitly supplied.
		Real opticalDepthScale = 240;
		/// Reference density / rho_c below which the benchmark envelope is transparent.
		Real opacityCutoffFraction = 0.016;
		int radialCells = 256, angularPoints = 32, multipoles = 12;
		Real structureTolerance = 1e-9;

		template <typename Archive>
		void serialize(Archive& archive, unsigned) {
			archive & centralGasFraction & rotationFraction & opticalDepthScale & opacityCutoffFraction
				& radialCells & angularPoints & multipoles & structureTolerance;
		}
	} radiatingStar;

	class WhiteDwarfOptions {
	public:
		Real thermalPressureFraction = 1e-4;
		Real rotationFraction = 0.2;
		int radialCells = 256, angularPoints = 32, multipoles = 12;
		Real structureTolerance = 1e-9;
		template <typename Archive>
		void serialize(Archive& archive, unsigned) {
			archive & thermalPressureFraction & rotationFraction & radialCells & angularPoints
				& multipoles & structureTolerance;
		}
	} whiteDwarf;

	class RadiationOptions {
	public:
		/// Add radiation to a hydro problem; radiation-only problems enable it themselves.
		bool enabled = false;
		/// Zero radiation energy transport at a fixed, otherwise free box boundary.
		/// This is an insulating numerical constraint, not a moving mirror.
		bool closedBoundary = false;
		Real lightSpeedRatio = 1;
		/// Equal gray absorption/emission opacity in cm^2/g; zero disables exchange.
		Real opacity = 0;
		/// Fixed physical length in cm for nonfatal RSLA diagnostics; zero uses box width.
		Real diagnosticLength = 0;
		/// Initial comoving radiation energy relative to a*T^4 for added radiation.
		Real initialEnergyRatio = 1;

		template <typename Archive>
		void serialize(Archive& archive, unsigned) {
			archive & enabled & lightSpeedRatio & opacity & diagnosticLength & initialEnergyRatio & closedBoundary;
		}
	} radiation;

	class GravityOptions {
	public:
		int multipoleOrder = 5;
		Real openingAngle = 0.5;
		std::string timeIntegration = "hierarchical";
		std::string energyTreatment = "mullen";
		bool conserveRegridEnergy = true;

		template <typename Archive>
		void serialize(Archive& archive, unsigned) {
			archive & multipoleOrder & openingAngle & timeIntegration & energyTreatment & conserveRegridEnergy;
		}
	} gravity;

	class OutputOptions {
	public:
		int every = 10;
		bool enabled = true;
		std::string directory = "output";

		template <typename Archive>
		void serialize(Archive& archive, unsigned) {
			archive & every & enabled & directory;
		}
	} output;

	class VerificationOptions {
	public:
		std::string analytic = "auto";
		std::string gravityReference = "direct";
		std::int64_t directMaxPairs = 20000000, directSamples = 0;
		Real relativeL1Tolerance = -1, absoluteTolerance = 1e-12;

		template <typename Archive>
		void serialize(Archive& archive, unsigned) {
			archive & analytic & gravityReference & directMaxPairs & directSamples & relativeL1Tolerance & absoluteTolerance;
		}
	} verification;

	/// Reject invalid runtime capacities and physical configuration values.
	void validate() const;

	bool hasExternalAcceleration() const {
		return std::any_of(hydro.acceleration.begin(), hydro.acceleration.end(), [](auto component) { return component != units::Acceleration{}; });
	}

	bool hydroEnabled() const;

	bool radiationEnabled() const;

	bool gravityEnabled() const;

	/// Serialize this value with its compile-time quantity types preserved.
	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & massFractions & problem & randomSeed & mesh & frame & amr & runtime & timestep & hydro & rayleighTaylor & star & radiatingStar & radiation & gravity & output & verification;
	}
};

/// Parse command-line/INI options into typed CGS values and validate the result.
Config parseConfig(std::vector<std::string> const& arguments);

/// Return the command-line option summary.
std::string helpText();
}	 // namespace octotigerII
