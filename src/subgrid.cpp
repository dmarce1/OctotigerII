#include "octotigerII/subgrid/subgrid.hpp"
#include "octotigerII/problems.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace octotigerII {
Subgrid::Subgrid(Config const& c, mesh::BlockLocation location) : config_(c) {
	c.validate();
	if (location.dimensionCount != c.dimensions || location.level != c.level)
		throw std::invalid_argument("Subgrid location/config mismatch");
	data_.location = location;
	data_.layout = mesh::MeshLayout(c.dimensions, c.cells, 2);
	Real const blockWidth = std::ldexp(c.upper - c.lower, -c.level);
	data_.cellWidth = blockWidth / c.cells;
	for (int axis = 0; axis < c.dimensions; ++axis) {
		if (location.coordinates[axis] < 0 || location.coordinates[axis] >= (1 << c.level))
			throw std::invalid_argument("Subgrid coordinate out of bounds");
		data_.lower[axis] = c.lower + location.coordinates[axis] * blockWidth;
	}
	data_.hydroEnabled = c.hydroEnabled();
	data_.radiationEnabled = c.radiationEnabled();
	data_.gravityEnabled = c.gravityEnabled();
	if (data_.hydroEnabled)
		data_.hydro = hydro::Fields(data_.layout, data_.cellWidth, data_.lower);
	if (data_.radiationEnabled)
		data_.radiation = radiation::Fields(data_.layout, data_.cellWidth, data_.lower);
	if (data_.gravityEnabled) {
		data_.gravity = gravity::Fields(data_.layout, data_.cellWidth, data_.lower);
		if (!data_.hydroEnabled)
			data_.density = mesh::PatchData<Real>(data_.layout, data_.cellWidth, data_.lower);
	}
	initializeProblem(data_, c);
}
Snapshot Subgrid::snapshot() const {
	return data_;
}

Real Subgrid::stableTimestep() const {
	Real result = std::numeric_limits<Real>::infinity();
	if (data_.hydroEnabled) {
		result = transportExchange::stableStep(data_.hydro, hydro::HydroSystem(config_.gamma),
											   config_.cfl);
		if (data_.gravityEnabled) {
			Real acceleration = 0;
			data_.layout.forEachInterior([&](mesh::Coordinates const&, std::size_t i) {
				Real norm = 0;
				for (int axis = 0; axis < 3; ++axis)
					norm += std::abs(data_.gravity.values()[i].acceleration(axis));
				acceleration = std::max(acceleration, norm);
			});
			if (acceleration > 0) {
				// Bound the extra half-kick velocity in the transport CFL.
				Real const b = config_.cfl / result;
				Real const a = 0.5 * acceleration / data_.cellWidth;
				result = 2 * config_.cfl / (b + std::sqrt(b * b + 4 * a * config_.cfl));
				result = std::min(result, 0.2 * std::sqrt(data_.cellWidth / acceleration));
			}
		}
	}
	if (data_.radiationEnabled)
		result = std::min(
			result, transportExchange::stableStep(
						data_.radiation,
						radiation::RadiationSystem(config_.lightSpeedRatio * physicalLightSpeed),
						config_.cfl));
	return result;
}
void Subgrid::advance(std::vector<HydroSnapshot> const& hydroSources,
					  std::vector<RadiationSnapshot> const& radiationSources, Real dt) {
	if (!(dt > 0) || !std::isfinite(dt))
		throw std::invalid_argument("Invalid step size");
	ExchangeDomain const domain{config_.dimensions, config_.lower, config_.upper, config_.periodic};
	if (data_.hydroEnabled)
		transportExchange::advancePatch(data_.location, data_.hydro, hydroSources, domain,
										hydro::HydroSystem(config_.gamma), dt);
	if (data_.radiationEnabled)
		transportExchange::advancePatch(
			data_.location, data_.radiation, radiationSources, domain,
			radiation::RadiationSystem(config_.lightSpeedRatio * physicalLightSpeed), dt);
	data_.time += dt;
}
void Subgrid::setGravity(std::vector<gravity::State> const& fields) {
	if (!data_.gravityEnabled || fields.size() != data_.layout.interiorCellCount())
		throw std::invalid_argument("Gravity assignment size/field mismatch");
	std::size_t source = 0;
	data_.layout.forEachInterior([&](mesh::Coordinates const&, std::size_t i) {
		data_.gravity.values()[i] = fields[source++];
	});
	data_.gravity.timeState().time = data_.time;
}
void Subgrid::kickGravity(Real dt) {
	if (!data_.hydroEnabled || !data_.gravityEnabled || !(dt > 0) || !std::isfinite(dt) ||
		!transportExchange::sameTime(data_.gravity.timeState().time, data_.time))
		throw std::logic_error("Gravity kick requires synchronized gas and gravity fields");
	data_.layout.forEachInterior([&](mesh::Coordinates const&, std::size_t i) {
		auto& u = data_.hydro.values()[i];
		Real work = 0;
		for (int axis = 0; axis < 3; ++axis) {
			Real const old = u.momentum(axis),
					   impulse = dt * u.density() * data_.gravity.values()[i].acceleration(axis);
			u.momentum(axis) += impulse;
			work += impulse * (old + 0.5 * impulse) / u.density();
		}
		// Preserve internal energy during the kick. No gravitational field
		// energy variable, potential-energy flux, or global energy correction.
		u.totalEnergy() += work;
		if (!hydro::HydroSystem(config_.gamma).admissible(u))
			throw std::runtime_error("Invalid gravity kick state");
	});
}
} // namespace octotigerII
