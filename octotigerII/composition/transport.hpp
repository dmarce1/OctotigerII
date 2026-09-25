/** @file
 * @brief Conservative material/tracer transport driven by hydro's mass flux.
 */
#pragma once
#include "octotigerII/composition/species.hpp"
#include "octotigerII/mesh.hpp"
#include <algorithm>
#include <array>
#include <stdexcept>

namespace octotigerII::composition {
using Fluxes = std::array<std::vector<units::MassFlux>, ndim>;

/// Upwind concentrations use the sign of the *limited numerical mass flux*.
/// All material fluxes partition that one flux; tracers never participate in
/// normalization. Donor-cell concentration transport is first order at contacts.
/// Input values are partial densities (or rho*C for a massless tracer).
template <typename Read, typename MassFlux>
std::vector<Fluxes> fluxes(Options const& options, mesh::MeshLayout const& layout, Read const& read, MassFlux const& massFlux) {
	std::size_t const count = options.species.size();
	std::vector<Fluxes> result(count);
	for (int d = 0; d < ndim; ++d) {
		for (auto& s : result) s[d].resize(layout.faceCount(d));
		mesh::forEachCoordinate(layout.faceExtents(d), [&](auto const& face) {
			auto const f = massFlux(d, face);
			auto donor = face;
			if (f >= units::MassFlux{}) --donor[d];
			std::vector<units::Density> values(count);
			units::Density rho{};
			std::size_t closure = count;
			for (std::size_t s = 0; s < count; ++s) {
				values[s] = read(s, donor);
				if (!units::finite(values[s]) || values[s] < units::Density{}) throw std::runtime_error("Invalid species donor density");
				if (!options.species[s].tracer()) {
					rho += values[s];
					if (closure == count || values[s] > values[closure]) closure = s;
				}
			}
			if (!(rho > units::Density{})) throw std::runtime_error("Species transport requires positive material density");
			auto remaining = f;
			auto const index = layout.faceIndex(d, face);
			for (std::size_t s = 0; s < count; ++s) if (s != closure) {
				auto const sf = Real(values[s] / rho) * f;
				result[s][d][index] = sf;
				if (!options.species[s].tracer()) remaining -= sf;
			}
			result[closure][d][index] = remaining;
		});
	}
	return result;
}

inline units::Density update(units::Density value, Fluxes const& flux, mesh::MeshLayout const& layout,
	mesh::Coordinates const& cell, units::TimePerLength dtOverDx) {
	units::MassFlux divergence{};
	for (int d = 0; d < ndim; ++d) {
		auto upper = cell; ++upper[d];
		divergence += flux[d][layout.faceIndex(d, cell)] - flux[d][layout.faceIndex(d, upper)];
	}
	auto const next = value + dtOverDx * divergence;
	if (!units::finite(next) || next < units::Density{})
		throw std::runtime_error("Species update is negative/nonfinite; reduce timestep (no species clipping is applied)");
	return next;
}
} // namespace octotigerII::composition
