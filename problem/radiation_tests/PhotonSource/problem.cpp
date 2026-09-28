/** @file
 * @brief Homogeneous gray matter heated by a fixed isotropic photon source.
 */
#include "octotigerII/problems.hpp"
#include "octotigerII/subgrid/subgrid.hpp"
#include "octotigerII/verification/analytic.hpp"

namespace octotigerII::photon_source {

ProblemBoundary problemBoundary(Config const&) { return {}; }
verification::Reference problemReference(Config const&) { return {}; }
void validateProblem(Config const&) {}

void problemDefaults(Config& config) {
	config.mesh.cells = 4;
	config.mesh.level = 0;
	config.mesh.lower = units::Length::from_value(-5e11);
	config.mesh.upper = units::Length::from_value(5e11);
	config.mesh.boundary = finiteVolume::BoundaryConditions::periodic();
	config.hydro.gamma = Real(5) / 3;
	config.radiation.enabled = true;
	config.radiation.opacity = 1;
	config.radiation.initialEnergyRatio = 1;
	config.amr.shadowTolerance = 0;
	config.runtime.stopTime = units::Time::from_value(0.1);
}

ProblemRadiationMaterial problemRadiationMaterial(Config const& config) {
	auto const opacity = units::Quantity<2, -1, 0>::from_value(config.radiation.opacity);
	return [opacity](mesh::PhysicalCoordinates const&, units::Time) {
		// Physical photon emissivity: exactly 1 erg/(cm^3 s), independent of
		// the evolving gas or radiation. Runtime applies the RSLA factor chat/c.
		return RadiationMaterial{opacity, units::Quantity<-1, 1, -3>::from_value(1)};
	};
}

void initializeProblem(Snapshot& data, Config const& config, bool) {
	hydro::HydroSystem const system(config);
	hydro::PrimitiveState primitive;
	primitive.density() = units::Density::from_value(1e-10);
	auto const temperature = units::Temperature::from_value(1e4);
	primitive.pressure() = primitive.density() * constants::boltzmann * temperature
		/ (config.hydro.meanMolecularWeight * constants::atomicMassUnit);
	auto const gas = system.conservedState(primitive);
	radiation::RadiationSystem::State rad;
	rad.energy() = config.radiation.initialEnergyRatio * constants::radiation * boost::units::pow<4>(temperature);
	for (auto& state : data.hydro.values()) state = gas;
	for (auto& state : data.radiation.values()) state = rad;
}

} // namespace octotigerII::photon_source
