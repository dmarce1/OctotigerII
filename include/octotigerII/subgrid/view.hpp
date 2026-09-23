/** @file
 * @brief Temporary patch views joining interiors to compact halo storage.
 * @ingroup mesh
 */
// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#pragma once

#include "octotigerII/profiling.hpp"

#include "octotigerII/storage/columns.hpp"
#include "octotigerII/subgrid/topology.hpp"


namespace octotigerII {


// An execution-time patch: interior columns are direct local storage views
// (or received buffers for stolen work). Only ghost values occupy workspace.
/// Execution-time patch combining direct interior columns with temporary ghosts.
/// The block, plan, and ghost vector must outlive this view. On stolen work the
/// interior columns are received buffers rather than local allocation views.
/// @ingroup mesh
template <typename State>
class PatchView {
public:

	PatchView(Subgrid const& block, storage::Columns<State> interior, HaloPlan const& plan, std::vector<State> const& ghosts)
	  : block_(block)
	  , interior_(std::move(interior))
	  , plan_(plan)
	  , ghosts_(ghosts)
	  , padded_(block.layout.cellsPerActiveDimension(), 2) {}

	/// Return the indexing layout represented by this object.
	mesh::MeshLayout const& layout() const {
		return padded_;
	}

	/// Return the physical width of one active cell in centimeters.
	units::Length cellWidth() const {
		return block_.cellWidth;
	}

	/// Read a cell using interior coordinates without a ghost offset.
	State atInterior(mesh::Coordinates const& cell) const {
		return interior_.at(block_.layout.index(cell));
	}

	/// Read a cell using the coordinates of the padded storage layout.
	State atStorage(mesh::Coordinates const& cell) const {
		if (padded_.isInterior(cell)) return atInterior(padded_.interiorCoordinates(cell));
		return ghosts_.at(plan_.ghostIndices.at(padded_.index(cell)));
	}

private:

	Subgrid const& block_;
	storage::Columns<State> interior_;
	HaloPlan const& plan_;
	std::vector<State> const& ghosts_;
	mesh::MeshLayout padded_;
};


/// Launch independent column reads, scatter them to compact ghost storage,
/// and drain all transfers before returning or propagating an exception.
template <typename State>
void readHalo(storage::ColumnHandle<State> const& fields, HaloPlan const& plan, unsigned bank, std::vector<State>& ghosts) {
	profiling::Elapsed profile("transport.halo.wall_ns");
	std::vector<storage::PendingColumns<State>> pending;
	pending.reserve(plan.reads.size());
	// Launch all independent reads before awaiting any result.
	for (auto const& read : plan.reads)
		pending.push_back(fields.read(read.range, bank));
	ghosts.resize(plan.ghostCount);
	std::exception_ptr error;
	for (std::size_t i = 0; i < pending.size(); ++i) {
		try {
			auto values = pending[i].get();
			for (auto const& copy : plan.reads[i].copies)
				ghosts[copy.destination] = values.at(copy.source);
		} catch (...) {
			if (!error) error = std::current_exception();
		}
	}
	if (error) std::rethrow_exception(error);
}

/// Complete physical boundary values after all donor reads have finished.
/// Analytic corners use the original position and take precedence over other faces.
template <typename System>
void applyHaloBoundaries(HaloPlan const& plan, std::vector<typename System::State>& ghosts, System const& system,
	units::Time time, physics::AnalyticBoundary<typename System::State> const& analytic = {}) {
	for (std::size_t i = 0; i < plan.reflectionMasks.size(); ++i)
		ghosts.at(i) = physics::transformBoundary(ghosts.at(i), plan.reflectionMasks[i], plan.outflowLowerMasks.at(i), plan.outflowUpperMasks.at(i), system);
	for (auto const& ghost : plan.analyticGhosts)
		ghosts.at(ghost.destination) = physics::evaluateBoundary(analytic, ghost.position, time, system);
}


}	 // namespace octotigerII
