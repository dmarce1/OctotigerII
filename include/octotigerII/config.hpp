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
#include "octotigerII/physics/boundary.hpp"
#include "octotigerII/units/constants.hpp"

namespace octotigerII {

/// Problem and runtime options. Lengths and times are physical CGS quantities.
/// Problem identity, dimension, and enabled physics are immutable build choices.
/// @ingroup runtime
class Config {
public:
	static constexpr auto problem = build::problem;
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

		template <typename Archive>
		void serialize(Archive& archive, unsigned) {
			archive & cfl;
		}
	} timestep;


	class HydroOptions {
	public:
		Real gamma = 1.4;
		std::array<units::Acceleration, ndim> acceleration{};

		template <typename Archive>
		void serialize(Archive& archive, unsigned) {
			archive & gamma;
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

	class RadiationOptions {
	public:
		Real lightSpeedRatio = 1;

		template <typename Archive>
		void serialize(Archive& archive, unsigned) {
			archive & lightSpeedRatio;
		}
	} radiation;


	class GravityOptions {
	public:
		int multipoleOrder = 5;
		Real openingAngle = 0.5;

		template <typename Archive>
		void serialize(Archive& archive, unsigned) {
			archive & multipoleOrder & openingAngle;
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

	static constexpr bool hydroEnabled() {
		return build::hydro;
	}

	static constexpr bool radiationEnabled() {
		return build::radiation;
	}

	static constexpr bool gravityEnabled() {
		return build::gravity;
	}

	/// Serialize this value with its compile-time quantity types preserved.
	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & randomSeed & mesh & runtime & timestep & hydro & rayleighTaylor & radiation & gravity & output & verification;
	}
};

/// Parse command-line/INI options into typed CGS values and validate the result.
Config parseConfig(std::vector<std::string> const& arguments);

/// Return the command-line option summary.
std::string helpText();
}	 // namespace octotigerII
