// Distributed under the Boost Software License, Version 1.0.
#include "octotigerII/problems.hpp"
#include "octotigerII/subgrid/subgrid.hpp"
#include <cmath>

namespace octotigerII {
void initializeProblem(Snapshot& data, Config const& c) {
	hydro::HydroSystem gas(c.gamma);
	data.layout.forEachInterior([&](mesh::Coordinates const& cell, std::size_t i) {
		auto const point = data.layout.cellCenter(data.lower, data.cellWidth, cell);
		Real const length = c.upper - c.lower;
		Real radius2 = 0;
		for (int d = 0; d < c.dimensions; ++d) {
			Real const r = (point[d] - (c.lower + c.upper) / 2) / length;
			radius2 += r * r;
		}
		Real const gaussian = std::exp(-radius2 / (2 * 0.15 * 0.15));
		if (data.hydroEnabled) {
			hydro::PrimitiveState primitive{};
			primitive.density() = 1;
			primitive.pressure() = 1;
			if (c.problem == "sod") {
				bool const left = point[0] < (c.lower + c.upper) / 2;
				primitive.density() = left ? 1 : 0.125;
				primitive.pressure() = left ? 1 : 0.1;
			} else if (c.problem == "kelvin-helmholtz") {
				Real const x = (point[0] - c.lower) / length, y = (point[1] - c.lower) / length;
				Real const layer =
					0.5 * (std::tanh((y - 0.25) / 0.025) - std::tanh((y - 0.75) / 0.025));
				primitive.density() = 1 + layer;
				primitive.pressure() = 2.5;
				primitive.velocity(0) = layer - 0.5;
				primitive.velocity(1) = 0.01 * std::sin(4 * piR * x);
			} else if (c.problem == "collapse") {
				primitive.density() = 1e4 * (0.01 + gaussian);
				primitive.pressure() = 1e12;
			}
			data.hydro.values()[i] = gas.conservedState(primitive);
		}
		if (data.radiationEnabled) {
			Real distance2 = 0;
			for (int axis = 0; axis < c.dimensions; ++axis) {
				Real distance = point[axis] - (c.lower + 0.25 * length);
				if (c.periodic)
					distance -= std::round(distance / length) * length;
				distance2 += distance * distance / (0.08 * length * 0.08 * length);
			}
			auto& state = data.radiation.values()[i];
			state[0] = 1e-6 + std::exp(-0.5 * distance2);
			if (c.problem == "streaming")
				for (int axis = 0; axis < c.dimensions; ++axis)
					state[axis + 1] = state[0] / std::sqrt(Real(c.dimensions));
		}
		if (data.gravityEnabled && !data.hydroEnabled)
			data.density.values()[i] =
				1e4 * (c.problem == "gravity-sphere" ? Real(radius2 < 0.25 * 0.25) : gaussian);
	});
}
} // namespace octotigerII
