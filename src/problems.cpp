// Distributed under the Boost Software License, Version 1.0.
#include "octotigerII/problems.hpp"
#include <cmath>
#include "octotigerII/subgrid/subgrid.hpp"


namespace octotigerII {
/** @brief Fill physical CGS initial states on one temporary interior patch.
 * The shock-tube ratios follow @ref ref_sod1978 "Sod (1978)".
 * Other fixtures are project regression configurations.
 */
void initializeProblem(Snapshot& data, Config const& c) {
	using std::exp;
	using std::round;
	using std::sin;
	using std::sqrt;
	using std::tanh;

	hydro::HydroSystem gas(c.hydro.gamma);
	data.layout.forEachInterior([&](mesh::Coordinates const& cell, std::size_t i) {
		auto const point = data.layout.cellCenter(data.lower, data.cellWidth, cell);
		auto const length = c.mesh.upper - c.mesh.lower;
		Real radius2 = 0;
		for (int d = 0; d < c.dimensions; ++d) {
			Real const r = (point[d] - (c.mesh.lower + c.mesh.upper) / 2.0) / length;
			radius2 += r * r;
		}
		Real const gaussian = exp(-radius2 / (2 * 0.15 * 0.15));
		if (data.hydroEnabled) {
			hydro::PrimitiveState primitive{};
			primitive.density() = units::Density::from_value(1);
			primitive.pressure() = units::Pressure::from_value(1);
			if (c.problem == "sod") {
				bool const left = point[0] < (c.mesh.lower + c.mesh.upper) / 2.0;
				primitive.density() = units::Density::from_value(left ? 1 : 0.125);
				primitive.pressure() = units::Pressure::from_value(left ? 1 : 0.1);
			} else if (c.problem == "kelvin-helmholtz") {
				Real const x = (point[0] - c.mesh.lower) / length, y = (point[1] - c.mesh.lower) / length;
				Real const layer = 0.5 * (tanh((y - 0.25) / 0.025) - tanh((y - 0.75) / 0.025));
				primitive.density() = units::Density::from_value(1 + layer);
				primitive.pressure() = units::Pressure::from_value(2.5);
				primitive.velocity(0) = units::Velocity::from_value(layer - 0.5);
				primitive.velocity(1) = units::Velocity::from_value(0.01 * sin(4 * piR * x));
			} else if (c.problem == "collapse") {
				primitive.density() = units::Density::from_value(1e4 * (0.01 + gaussian));
				primitive.pressure() = units::Pressure::from_value(1e12);
			}
			data.hydro.values()[i] = gas.conservedState(primitive);
		}
		if (data.radiationEnabled) {
			Real distance2 = 0;
			for (int axis = 0; axis < c.dimensions; ++axis) {
				auto distance = point[axis] - (c.mesh.lower + 0.25 * length);
				if (c.mesh.boundary.periodic(axis)) distance -= round(distance / length) * length;
				distance2 += distance * distance / (0.08 * length * 0.08 * length);
			}
			auto& state = data.radiation.values()[i];
			state.energy() = units::EnergyDensity::from_value(1e-6 + exp(-0.5 * distance2));
			if (c.problem == "streaming")
				for (int axis = 0; axis < c.dimensions; ++axis)
					state.radiativeFlux(axis) = constants::c * state.energy() / sqrt(Real(c.dimensions));
		}
		if (data.gravityEnabled && !data.hydroEnabled)
			data.density.values()[i] = units::Density::from_value(1e4 * (c.problem == "gravity-sphere" ? Real(radius2 < 0.25 * 0.25) : gaussian));
	});
}
}	 // namespace octotigerII
