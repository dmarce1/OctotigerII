/** @file
 * @brief Opaque, source-free, uniformly rotating gas+radiation SCF star.
 */
#include "octotigerII/problems.hpp"
#include "octotigerII/problems/radiatingStarStructure.hpp"
#include "octotigerII/radiation/opacity.hpp"
#include "octotigerII/verification/analytic.hpp"
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <mutex>
#include <shared_mutex>
#include <stdexcept>
#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/synchronization/shared_mutex.hpp>
#endif

namespace octotigerII::radiatingStar {
namespace {
using Model = problems::RadiatingStarStructure;

std::shared_ptr<Model const> model(Config const& config) {
	Model::Parameters parameters;
	parameters.eos = {config.star.polytropicIndex, config.star.centralDensity,
		config.radiatingStar.centralGasFraction, config.hydro.meanMolecularWeight};
	parameters.spinFractionOfSphericalBreakup = config.radiatingStar.rotationFraction;
	parameters.radialCells = config.radiatingStar.radialCells;
	parameters.angularPoints = config.radiatingStar.angularPoints;
	parameters.maxMultipole = config.radiatingStar.multipoles;
	parameters.tolerance = config.radiatingStar.structureTolerance;
	using Key = std::array<Real, 9>;
	Key const key{parameters.eos.index, units::value(parameters.eos.centralDensity),
		parameters.eos.centralBeta, parameters.eos.meanMolecularWeight,
		parameters.spinFractionOfSphericalBreakup, Real(parameters.radialCells),
		Real(parameters.angularPoints), Real(parameters.maxMultipole), parameters.tolerance};
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
	auto value = std::make_shared<Model const>(parameters);
	std::unique_lock lock(mutex);
	if (cache.size() >= 8) cache.clear();
	return cache.try_emplace(key, std::move(value)).first->second;
}

verification::ExactState state(Config const& config, Model const& star,
	mesh::PhysicalCoordinates position, Real probeScale = 1) {
	for (int axis = 0; axis < ndim; ++axis) { position[axis] -= config.star.center[axis]; }
	auto const R = units::hypot(position[0], position[1]);
	auto const profile = star.sample(R / probeScale, position[2] / probeScale);
	auto const floor = star.eos().atDensity(config.star.atmosphereFraction * config.star.centralDensity);
	auto const& matter = profile.thermodynamics.density > floor.density ? profile.thermodynamics : floor;
	verification::ExactState result;
	result.hydro.density() = matter.density;
	result.hydro.pressure() = matter.gasPressure;
	// As in rotatingStar, momentum belongs to the stellar material before the
	// density floor is applied. The exterior atmosphere is inertially at rest.
	Real const materialFraction = Real(profile.thermodynamics.density / matter.density);
	result.hydro.velocity(0) = -materialFraction * star.angularVelocity() * position[1];
	result.hydro.velocity(1) = materialFraction * star.angularVelocity() * position[0];

	std::array<units::EnergyFlux, ndim> comovingFlux{};
	Real const nx = R > units::Length{} ? Real(position[0] / R) : 0;
	Real const ny = R > units::Length{} ? Real(position[1] / R) : 0;
	if (profile.thermodynamics.density > floor.density) {
		auto const fluxOpacity = radiation::opacityLaw(config, radiation::Opacity::from_value(config.radiation.opacity))
			.evaluate(units::value(profile.thermodynamics.density), units::value(profile.thermodynamics.temperature)).fluxExtinction;
		auto const diffusion = star.diffusion(profile, units::value(fluxOpacity));
		comovingFlux[0] = nx * diffusion.flux[0];
		comovingFlux[1] = ny * diffusion.flux[0];
		comovingFlux[2] = diffusion.flux[1];
	}
	// Diffusion ceases to be valid in the final thin surface layer. Limit only
	// the comoving flux there; never add the formal div(F) as a photon heater.
	auto const E0 = matter.radiationEnergy;
	units::EnergyFlux norm{};
	for (auto flux : comovingFlux) { norm = units::hypot(norm, flux); }
	auto const maximumFlux = (1 - 32 * epsilonR) * constants::c * E0;
	if (norm > maximumFlux) {
		for (auto& flux : comovingFlux) { flux *= Real(maximumFlux / norm); }
		norm = maximumFlux;
	}
	Real const f = Real(norm / (constants::c * E0));
	Real const f2 = std::min(Real(1), f * f);
	auto const transversePressure = E0 * ((1 - f2) / (1 + std::sqrt(4 - 3 * f2)));
	Real beta2 = 0;
	for (int axis = 0; axis < ndim; ++axis) {
		Real const beta = Real(result.hydro.velocity(axis) / constants::c);
		beta2 += beta * beta;
	}
	if (!(beta2 < 1)) throw std::invalid_argument("Radiating-star rotation must remain subluminal");
	Real const boost2 = 1 / (1 - beta2), boost = std::sqrt(boost2);
	// The azimuthal velocity is perpendicular to the meridional M1 flux.
	// This exact boost keeps the retained mixed-frame local thermal source zero.
	result.radiation.energy() = boost2 * (E0 + beta2 * transversePressure);
	for (int axis = 0; axis < ndim; ++axis) {
		result.radiation.radiativeFlux(axis) = boost * comovingFlux[axis]
			+ boost2 * (E0 + transversePressure) * result.hydro.velocity(axis);
	}
	result.gravity.potential() = profile.potential;
	result.gravity.acceleration(0) = -nx * profile.potentialGradient[0];
	result.gravity.acceleration(1) = -ny * profile.potentialGradient[0];
	result.gravity.acceleration(2) = -profile.potentialGradient[1];
	return result;
}
} // namespace

ProblemBoundary problemBoundary(Config const& config) {
	auto star = model(config);
	return [config, star](auto const& position, units::Time) { return state(config, *star, position); };
}

verification::Reference problemReference(Config const& config) {
	verification::Reference result;
	result.name = "Initial opaque rotating SCF gas+radiation star (zero photon source)";
	result.evaluate = radiatingStar::problemBoundary(config);
	return result;
}

void problemDefaults(Config& config) {
	config.star.centralDensity = units::Density::from_value(1);
	config.star.polytropicIndex = 3.5;
	config.star.atmosphereFraction = 1e-12;
	config.radiatingStar.centralGasFraction = 0.8;
	config.radiatingStar.rotationFraction = 0.2;
	config.hydro.gamma = Real(5) / 3;
	config.hydro.meanMolecularWeight = 0.6;
	config.frame.omega = {};
	config.mesh.cells = 8;
	config.mesh.level = 0;
	config.mesh.boundary = physics::BoundaryConditions::uniform(physics::BoundaryCondition::Outflow);
	config.amr.enabled = true;
	config.amr.minLevel = 0;
	config.amr.maxLevel = 4;
	config.amr.refineDensity = 0.001 * config.star.centralDensity;
	config.amr.shadowTolerance = 0;
	config.radiation.lightSpeedRatio = 1;
	config.radiation.opacity = 0.34;
	config.radiation.closedBoundary = true;
	// Radiation is enabled by the manifest; this flag adds radiation to a
	// hydro-only initializer and must not replace the explicit moments above.
	config.radiation.enabled = false;
}

void validateProblem(Config const& config) {
	if (std::abs(config.hydro.gamma - Real(5) / 3) > 1e-12 || config.radiation.enabled)
		throw std::invalid_argument("The opaque rotating star requires gas gamma=5/3 and its explicit radiation initializer");
	if (config.radiation.lightSpeedRatio != 1 || !config.radiation.closedBoundary)
		throw std::invalid_argument("The opaque rotating star requires physical light speed and radiation.closedBoundary=on");
	if (config.radiation.opacityModel == "constant" && !(config.radiation.opacity + config.radiation.scatteringOpacity > 0))
		throw std::invalid_argument("The opaque rotating star requires positive flux opacity");
	if (config.frame.omega != units::InverseTime{}
		|| !config.mesh.boundary.all(physics::BoundaryCondition::Outflow))
		throw std::invalid_argument("The insulating rotating-star benchmark requires a fixed grid with outflow gas boundaries and isolated gravity");
	if (!(config.star.atmosphereFraction > 0 && config.star.atmosphereFraction < 1e-3)
		|| !std::isfinite(config.star.atmosphereFraction))
		throw std::invalid_argument("Invalid radiating-star numerical atmosphere fraction");
	for (auto center : config.star.center) {
		if (!units::finite(center)) throw std::invalid_argument("Nonfinite stellar center");
	}
	auto const star = model(config);
	for (int axis = 0; axis < ndim; ++axis) {
		auto const radius = axis == 2 ? star->polarRadius() : star->equatorialRadius();
		if (config.star.center[axis] - radius <= config.mesh.lower
			|| config.star.center[axis] + radius >= config.mesh.upper)
			throw std::invalid_argument("The complete rotating SCF star must fit inside the simulation box");
	}
	auto const floor = star->eos().atDensity(config.star.atmosphereFraction * config.star.centralDensity);
	if (floor.density < units::Density::from_value(1e-14) || floor.gasPressure < units::Pressure::from_value(1e-14))
		throw std::invalid_argument("The stellar numerical atmosphere falls below the hydro density or pressure floor");
}

void initializeProblem(Snapshot& data, Config const& config, bool refinementProbe) {
	auto star = model(config);
	Real const scale = refinementProbe ? Real(initialFeatureWidth(star->eos().scaleLength(), data.cellWidth)
		/ star->eos().scaleLength()) : 1;
	hydro::HydroSystem const gas(config);
	data.layout.forEachInterior([&](auto const& cell, std::size_t index) {
		auto const position = data.layout.cellCenter(data.lower, data.cellWidth, cell);
		auto const value = state(config, *star, position, scale);
		data.hydro.values()[index] = gas.conservedState(value.hydro);
		data.radiation.values()[index] = value.radiation;
	});
}
} // namespace octotigerII::radiatingStar
