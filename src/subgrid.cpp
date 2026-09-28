#include "octotigerII/subgrid/subgrid.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include "octotigerII/problems.hpp"
#include "octotigerII/radiation/diffusionFlux.hpp"

namespace octotigerII {
Snapshot initialSnapshot(Config const& c, mesh::BlockLocation location, bool refinementProbe) {
	using std::ldexp;

	Snapshot data_;
	c.validate();
	if (location.level < 0 || location.level > 16) throw std::invalid_argument("Invalid subgrid level");
	data_.location = location;
	data_.layout = mesh::MeshLayout(c.mesh.cells);
	auto const blockWidth = (c.mesh.upper - c.mesh.lower) * ldexp(Real(1), -location.level);
	data_.cellWidth = blockWidth / Real(c.mesh.cells);
	for (int axis = 0; axis < ndim; ++axis) {
		if (location.coordinates[axis] < 0 || location.coordinates[axis] >= (1 << location.level))
			throw std::invalid_argument("Subgrid coordinate out of bounds");
		data_.lower[axis] = c.mesh.lower + Real(location.coordinates[axis]) * blockWidth;
	}
	data_.hydroEnabled = c.hydroEnabled();
	data_.radiationEnabled = c.radiationEnabled();
	data_.gravityEnabled = c.gravityEnabled();
	if (build::hydro && c.hydroEnabled()) data_.hydro = hydro::Fields(data_.layout, data_.cellWidth, data_.lower);
	if (build::radiation && c.radiationEnabled()) data_.radiation = radiation::Fields(data_.layout, data_.cellWidth, data_.lower);
	if (build::gravity && c.gravityEnabled()) {
		data_.gravity = gravity::Fields(data_.layout, data_.cellWidth, data_.lower);
		if (!data_.hydroEnabled) data_.density = mesh::PatchData<units::Density>(data_.layout, data_.cellWidth, data_.lower);
	}
	initializeProblem(data_, c, refinementProbe);
	if constexpr (build::hydro && build::radiation) if (c.radiation.enabled) {
		hydro::HydroSystem const gas(c.hydro);
		for (std::size_t i = 0; i < data_.hydro.values().size(); ++i) {
			auto const& material = data_.hydro.values()[i];
			auto const temperature = gas.temperature(material);
			auto const thermalEnergy = c.radiation.initialEnergyRatio * constants::radiation * boost::units::pow<4>(temperature);
			radiation::MaterialVelocity velocity{};
			for (int d = 0; d < ndim; ++d) velocity[d] = material.momentum(d) / material.density();
			// Match the retained mixed-frame source equations. At ratio one both
			// local exchange rates vanish; this does not imply spatial equilibrium.
			auto const unitEnergy = units::EnergyDensity::from_value(1);
			auto const unit = radiation::materialEquilibriumMoments(unitEnergy, velocity);
			Real workFraction = 0;
			for (int d = 0; d < ndim; ++d)
				workFraction += Real(velocity[d] * unit.radiativeFlux(d) / (constants::c * constants::c * unitEnergy));
			data_.radiation.values()[i] = radiation::materialEquilibriumMoments(thermalEnergy / (1 - workFraction), velocity);
		}
	}
	if (c.massFractions.enabled) {
		if (data_.species.empty()) for (auto const& species : c.massFractions.species) {
			data_.species.emplace_back(data_.layout, data_.cellWidth, data_.lower);
			for (std::size_t i = 0; i < data_.hydro.values().size(); ++i)
				data_.species.back().values()[i] = species.initialFraction * data_.hydro.values()[i].density();
		}
		if (data_.species.size() != c.massFractions.species.size()) throw std::invalid_argument("Problem species field count mismatch");
		for (auto const& species : data_.species) {
			if (species.values().size() != data_.hydro.values().size()) throw std::invalid_argument("Problem species field size mismatch");
			for (auto value : species.values())
				if (!units::finite(value) || value < units::Density{}) throw std::invalid_argument("Problem species density must be finite and nonnegative");
		}
		for (std::size_t i = 0; i < data_.hydro.values().size(); ++i) {
			units::Density rho{};
			for (std::size_t s = 0; s < data_.species.size(); ++s) if (!c.massFractions.species[s].tracer()) rho += data_.species[s].values()[i];
			if (!(rho > units::Density{})) throw std::invalid_argument("Problem material density must be positive");
			data_.hydro.values()[i].density() = rho;
		}
	}
	return data_;
}

}	 // namespace octotigerII
