/** @file
 * @brief Rotating full-M1 reference with prescribed opacity and photon heating.
 */
#include "octotigerII/problems.hpp"
#include "octotigerII/problems/radiatingStarReference.hpp"
#include "octotigerII/verification/analytic.hpp"
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <mutex>
#include <shared_mutex>
#include <stdexcept>
#include <string_view>
#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/synchronization/shared_mutex.hpp>
#endif

namespace octotigerII::radiatingStar {
namespace {
using Model = problems::RadiatingStarReference;
constexpr std::string_view referenceData =
#include "reference.inc"
;

std::shared_ptr<Model const> model(Config const& config) {
	Model::Parameters parameters;
	parameters.centralDensity = config.star.centralDensity;
	parameters.meanMolecularWeight = config.hydro.meanMolecularWeight;
	parameters.atmosphereTolerance = config.radiatingStar.structureTolerance;
	using Key = std::array<Real, 3>;
	Key const key{units::value(parameters.centralDensity), parameters.meanMolecularWeight,
		parameters.atmosphereTolerance};
#ifdef OCTOTIGERII_WITH_HPX
	static hpx::shared_mutex mutex;
#else
	static std::shared_mutex mutex;
#endif
	static std::map<Key, std::shared_ptr<Model const>> cache;
	{
		std::shared_lock lock(mutex);
		if (auto found = cache.find(key); found != cache.end()) return found->second;
	}
	auto value = std::make_shared<Model const>(parameters, referenceData);
	std::unique_lock lock(mutex);
	if (cache.size() >= 8) cache.clear();
	return cache.try_emplace(key, std::move(value)).first->second;
}

mesh::PhysicalCoordinates relativePosition(Config const& config, mesh::PhysicalCoordinates position) {
	for (int axis = 0; axis < ndim; ++axis) { position[axis] -= config.star.center[axis]; }
	return position;
}

verification::ExactState state(Config const& config, Model const& star,
	mesh::PhysicalCoordinates const& position, Real probeScale = 1) {
	auto const x = relativePosition(config, position);
	auto const cylindricalRadius = units::hypot(x[0], x[1]);
	auto const profile = star.sample(cylindricalRadius / probeScale, x[2] / probeScale);
	auto const floor = star.gasBarotrope().atDensity(config.star.atmosphereFraction * config.star.centralDensity);
	verification::ExactState result;
	result.hydro.density() = std::max(profile.density, floor.density);
	result.hydro.pressure() = std::max(profile.gasPressure, floor.gasPressure);
	Real const materialFraction = Real(profile.density / result.hydro.density());
	result.hydro.velocity(0) = -materialFraction * profile.angularVelocity * x[1];
	result.hydro.velocity(1) = materialFraction * profile.angularVelocity * x[0];
	result.radiation.energy() = profile.radiationEnergy;
	result.radiation.radiativeFlux(2) = profile.verticalFlux;
	result.gravity.potential() = profile.potential;
	result.gravity.acceleration(2) = profile.gravity[1];
	if (cylindricalRadius > units::Length{}) {
		Real const nx = Real(x[0] / cylindricalRadius), ny = Real(x[1] / cylindricalRadius);
		result.radiation.radiativeFlux(0) = nx * profile.cylindricalFlux - ny * profile.azimuthalFlux;
		result.radiation.radiativeFlux(1) = ny * profile.cylindricalFlux + nx * profile.azimuthalFlux;
		result.gravity.acceleration(0) = nx * profile.gravity[0];
		result.gravity.acceleration(1) = ny * profile.gravity[0];
	}
	return result;
}
} // namespace

ProblemBoundary problemBoundary(Config const& config) {
	auto star = model(config);
	return [config, star](auto const& position, units::Time) { return state(config, *star, position); };
}

ProblemRadiationMaterial problemRadiationMaterial(Config const& config) {
	auto star = model(config);
	return [config, star](auto const& position, units::Time) {
		auto const x = relativePosition(config, position);
		auto const cylindricalRadius = units::hypot(x[0], x[1]);
		// The transparent source-free exterior also covers boundary ghost cells
		// outside the finite reference table.
		if (units::hypot(cylindricalRadius, x[2]) >= star->transitionRadius()) return RadiationMaterial{};
		auto const value = star->sample(cylindricalRadius, x[2]);
		return RadiationMaterial{units::Quantity<2, -1, 0>::from_value(value.opacityCgs),
			units::Quantity<-1, 1, -3>::from_value(value.photonHeatingCgs)};
	};
}

verification::Reference problemReference(Config const& config) {
	verification::Reference result;
	result.name = "Rotating full-M1 balanced reference (fixed photon heater; transparent envelope)";
	result.evaluate = radiatingStar::problemBoundary(config);
	return result;
}

void problemDefaults(Config& config) {
	config.star.centralDensity = units::Density::from_value(1);
	config.star.polytropicIndex = 3.5;
	config.star.atmosphereFraction = 1e-12;
	config.radiatingStar.centralGasFraction = 0.8;
	config.radiatingStar.rotationFraction = 0.1;
	config.radiatingStar.opticalDepthScale = 100 * std::sqrt(5.625);
	config.radiatingStar.opacityCutoffFraction = 0.015;
	config.hydro.gamma = Real(5) / 3;
	config.hydro.meanMolecularWeight = 0.6;
	config.mesh.cells = 8;
	config.mesh.level = 0;
	config.mesh.boundary = physics::BoundaryConditions::uniform(physics::BoundaryCondition::Outflow);
	config.amr.enabled = true;
	config.amr.minLevel = 0;
	config.amr.maxLevel = 5;
	config.amr.refineDensity = 0.001 * config.star.centralDensity;
	config.amr.shadowTolerance = 0;
	config.radiation.lightSpeedRatio = 1;
	// The manifest enables radiation; the initializer already supplies M1 moments.
	config.radiation.enabled = false;
}

void validateProblem(Config const& config) {
	if (std::abs(config.hydro.gamma - Real(5) / 3) > 1e-12 || config.radiation.enabled)
		throw std::invalid_argument("The rotating radiating star requires gas gamma=5/3 and its explicit M1 initializer");
	if (!(config.star.atmosphereFraction > 0 && config.star.atmosphereFraction < 1e-3)
		|| !std::isfinite(config.star.atmosphereFraction))
		throw std::invalid_argument("Invalid radiating-star numerical atmosphere fraction");
	if (!config.mesh.boundary.all(physics::BoundaryCondition::Outflow))
		throw std::invalid_argument("The rotating radiating star requires outflow transport boundaries and isolated gravity");
	for (auto center : config.star.center) {
		if (!units::finite(center)) throw std::invalid_argument("Nonfinite stellar center");
	}
	if (config.frame.omega != units::InverseTime{}
		&& (config.star.center[0] != units::Length{} || config.star.center[1] != units::Length{}))
		throw std::invalid_argument("The prescribed rotating-star material must be centered on the grid rotation axis");
	auto star = model(config);
	auto const& metadata = star->metadata();
	auto const close = [](Real a, Real b) { return std::abs(a - b) <= 1e-12 * std::max(Real(1), std::abs(b)); };
	Real const opacityScale = config.radiation.opacity * units::value(config.star.centralDensity)
		* units::value(star->lengthScale());
	if (!close(config.star.polytropicIndex, metadata.index)
		|| !close(config.radiatingStar.centralGasFraction, metadata.centralBeta)
		|| !close(config.radiatingStar.rotationFraction, metadata.spin)
		|| !close(config.radiatingStar.opacityCutoffFraction, metadata.cutoffDensityFraction)
		|| !close(opacityScale, metadata.opacity))
		throw std::invalid_argument("The selected structural parameters do not match the generated rotating reference; regenerate the table for a different model");
	units::Length farthest{};
	for (int axis = 0; axis < ndim; ++axis) {
		auto const distance = std::max(units::abs(config.mesh.lower - config.star.center[axis]),
			units::abs(config.mesh.upper - config.star.center[axis]));
		farthest = units::hypot(farthest, distance);
	}
	if (farthest > star->maximumRadius())
		throw std::invalid_argument("The rotating reference table must cover every corner of the simulation box");
	for (int axis = 0; axis < ndim; ++axis) { for (auto face : {config.mesh.lower, config.mesh.upper}) {
		auto point = config.star.center;
		point[axis] = face;
		auto const x = relativePosition(config, point);
		if (star->sample(units::hypot(x[0], x[1]), x[2]).density > units::Density{})
			throw std::invalid_argument("The complete rotating star must fit inside the domain");
	}
	}
	if (config.star.atmosphereFraction * config.star.centralDensity < units::Density::from_value(1e-14))
		throw std::invalid_argument("The stellar numerical atmosphere falls below the hydro density floor");
}

void initializeProblem(Snapshot& data, Config const& config, bool refinementProbe) {
	auto star = model(config);
	Real const scale = refinementProbe ? Real(initialFeatureWidth(star->eos().scaleLength(), data.cellWidth)
		/ star->eos().scaleLength()) : 1;
	hydro::HydroSystem const gas(config.hydro);
	data.layout.forEachInterior([&](auto const& cell, std::size_t index) {
		auto const position = data.layout.cellCenter(data.lower, data.cellWidth, cell);
		auto const value = state(config, *star, position, scale);
		data.hydro.values()[index] = gas.conservedState(value.hydro);
		data.radiation.values()[index] = value.radiation;
	});
}
} // namespace octotigerII::radiatingStar
