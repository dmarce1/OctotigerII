/** @file
 * @brief Conservative gravitational work from accepted finite-volume mass fluxes.
 * @ingroup gravity
 */
// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#pragma once

#include "octotigerII/conservation.hpp"
#include "octotigerII/gravity/gravityFields.hpp"
#include "octotigerII/storage/columns.hpp"
#include "octotigerII/subgrid/topology.hpp"

#include <exception>
#include <optional>
#include <stdexcept>
#include <vector>

namespace octotigerII::gravity {

/// Cell energy-density increments and the corresponding physical-boundary ledger.
struct FluxWork {
	std::vector<units::EnergyDensity> work;
	BoundaryTransport boundary;
};

/// Compute work using the mean of two endpoint fields without changing storage.
/// massFlux holds interval-averaged fluxes, so dt supplies their time integral.
/// Canonical work uses the plan's fine subface flux at coarse/fine interfaces.
/// Provisional work (canonical=false) uses this block's own flux on those faces.
/// The endpoint fields may be masked interaction components. Call on both sides
/// of each interface, including cells outside the component's active set.
inline FluxWork fluxWork(Subgrid const& block, std::vector<GravityWorkFace> const& plan,
	storage::ColumnHandle<State> const& oldGravity, unsigned oldBank,
	storage::ColumnHandle<State> const& newGravity, unsigned newBank,
	storage::FieldHandle<units::MassFlux> const& massFlux, unsigned fluxBank, units::Time dt, bool canonical = true) {
	if (!(dt >= units::Time{}) || !units::finite(dt)) throw std::invalid_argument("Gravity flux-work interval must be finite and nonnegative");
	if (!(block.cellWidth > units::Length{}) || !units::finite(block.cellWidth))
		throw std::invalid_argument("Gravity flux work requires a positive finite cell width");
	if (block.layout.ghostWidth() != 0 || block.layout.cellCount() != block.interior.count)
		throw std::invalid_argument("Gravity flux work requires a compact interior block layout");
	FluxWork result{std::vector<units::EnergyDensity>(block.interior.count), {}};
	if (dt == units::Time{}) return result;

	using PotentialRead = storage::Future<storage::Buffer<units::VelocitySquared>>;
	using MassRead = storage::Future<storage::Buffer<units::MassFlux>>;
	struct PendingFace {
		std::optional<PotentialRead> oldPotential, newPotential;
		std::optional<MassRead> mass;
	};
	std::optional<storage::PendingColumns<State>> oldPending, newPending;
	std::optional<MassRead> massPending;
	std::vector<PendingFace> pending(plan.size());
	std::exception_ptr error;
	// Launch independent transfers before awaiting them. If setup or a transfer
	// fails, drain every already-issued read before propagating its exception.
	try {
		oldPending = oldGravity.read(block.interior, oldBank);
		newPending = newGravity.read(block.interior, newBank);
		massPending = massFlux.read(block.massFlux, fluxBank);
		for (std::size_t i = 0; i < plan.size(); ++i) {
			auto const& face = plan[i];
			if (canonical) pending[i].mass = massFlux.read(face.massFlux, fluxBank);
			if (!face.physical) {
				pending[i].oldPotential = std::get<0>(oldGravity.fields).read(face.neighbor, oldBank);
				pending[i].newPotential = std::get<0>(newGravity.fields).read(face.neighbor, newBank);
			}
		}
	} catch (...) { error = std::current_exception(); }
	auto consume = [&](auto& read, auto&& accept) {
		if (!read) return;
		try { accept(read->get()); }
		catch (...) { if (!error) error = std::current_exception(); }
	};
	storage::Columns<State> old, now;
	storage::Buffer<units::MassFlux> mass;
	consume(oldPending, [&](auto value) { old = std::move(value); });
	consume(newPending, [&](auto value) { now = std::move(value); });
	consume(massPending, [&](auto value) { mass = std::move(value); });
	std::vector<units::MassFlux> faceMass(plan.size());
	std::vector<units::VelocitySquared> neighborPotential(plan.size());
	for (std::size_t i = 0; i < pending.size(); ++i) {
		consume(pending[i].mass, [&](auto value) { faceMass[i] = value.data()[0]; });
		consume(pending[i].oldPotential, [&](auto value) { neighborPotential[i] += 0.5 * value.data()[0]; });
		consume(pending[i].newPotential, [&](auto value) { neighborPotential[i] += 0.5 * value.data()[0]; });
	}
	if (error) std::rethrow_exception(error);

	std::vector<units::VelocitySquared> potential(block.interior.count);
	for (std::size_t i = 0; i < potential.size(); ++i) potential[i] = 0.5 * (old.at(i).potential() + now.at(i).potential());
	for (int axis = 0; axis < ndim; ++axis) {
		auto extents = block.layout.faceExtents(axis);
		extents[axis] -= 2;
		mesh::forEachCoordinate(extents, [&](auto face) {
			++face[axis];
			auto left = face;
			--left[axis];
			auto const l = block.layout.index(left), r = block.layout.index(face);
			auto const amount = (dt / block.cellWidth) * mass.data()[allFaceIndex(block.layout, axis, face)];
			auto const energy = -0.5 * amount * (potential[r] - potential[l]);
			result.work[l] += energy;
			result.work[r] += energy;
		});
	}

	auto const volume = block.layout.cellMeasure(block.cellWidth);
	for (std::size_t i = 0; i < plan.size(); ++i) {
		auto const& face = plan[i];
		if (!canonical) {
			auto index = face.cell;
			mesh::Coordinates coordinate{};
			for (int axis = 0; axis < ndim; ++axis) {
				coordinate[axis] = int(index % std::size_t(block.layout.interiorExtent(axis)));
				index /= std::size_t(block.layout.interiorExtent(axis));
			}
			if (face.sign > 0) ++coordinate[face.axis];
			faceMass[i] = mass.data()[allFaceIndex(block.layout, face.axis, coordinate)];
		}
		auto const outward = Real(face.sign) * face.areaFraction * (dt / block.cellWidth) * faceMass[i];
		units::VelocitySquared difference{};
		if (face.physical) {
			auto const acceleration = 0.5 * (old.at(face.cell).acceleration(face.axis) + now.at(face.cell).acceleration(face.axis));
			difference = -Real(face.sign) * 0.5 * block.cellWidth * acceleration;
			auto const transported = volume * outward * (potential[face.cell] + difference);
			if (transported < units::Energy{}) result.boundary.inward.potentialEnergy -= transported;
			else result.boundary.outward.potentialEnergy += transported;
		} else {
			difference = face.potentialFraction * (neighborPotential[i] - potential[face.cell]);
		}
		result.work[face.cell] -= outward * difference;
	}
	return result;
}

/// Compute work from a single supplied potential/acceleration field.
inline FluxWork fluxWork(Subgrid const& block, std::vector<GravityWorkFace> const& plan,
	storage::ColumnHandle<State> const& gravity, unsigned gravityBank,
	storage::FieldHandle<units::MassFlux> const& massFlux, unsigned fluxBank, units::Time dt, bool canonical = true) {
	return fluxWork(block, plan, gravity, gravityBank, gravity, gravityBank, massFlux, fluxBank, dt, canonical);
}

} // namespace octotigerII::gravity
