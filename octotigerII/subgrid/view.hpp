/** @file
 * @brief Temporary patch views joining interiors to compact halo storage.
 * @ingroup mesh
 */
// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#pragma once

#include <optional>
#include "octotigerII/profiling.hpp"

#include "octotigerII/amr/interpolation.hpp"
#include "octotigerII/storage/columns.hpp"
#include "octotigerII/subgrid/topology.hpp"

namespace octotigerII {

struct HaloTime {
	unsigned bank = 0;
	Real fraction = 0;
	unsigned nextBank = ~0u; // Default selects bank^1; coupled gravity uses a separate predictor.
	template <typename Archive>
	void serialize(Archive& archive, unsigned) { archive & bank & fraction & nextBank; }
};

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

/// Read the donor's own bank and, for an advancing coarse donor, interpolate
/// its two endpoints before spatial prolongation. Drain every launched read.
template <typename Field, typename State, typename Access>
void readTimedHalo(Field const& field, HaloPlan const& plan, unsigned bank, std::vector<State>& ghosts,
	std::vector<HaloTime> const& times, Access access) {
	profiling::Elapsed profile("transport.halo.wall_ns");
	using Pending = decltype(field.read(storage::Range{}, bank));
	std::vector<Pending> old, next;
	std::vector<HaloTime> selected;
	for (auto const& read : plan.reads) {
		auto const t = times.empty() ? HaloTime{bank, 0} : times.at(read.level);
		selected.push_back(t);
		old.push_back(field.read(read.range, t.bank));
		// A zero fraction needs no new endpoint, which may still be unwritten.
		if (t.fraction != 0) next.push_back(field.read(read.range, t.nextBank == ~0u ? (t.bank ^ 1) : t.nextBank));
	}
	ghosts.assign(plan.valueCount ? plan.valueCount : plan.ghostCount, State{});
	std::exception_ptr error;
	std::size_t j = 0;
	for (std::size_t i = 0; i < old.size(); ++i) {
		std::optional<decltype(old[i].get())> values, future;
		try { values = old[i].get(); } catch (...) { if (!error) error = std::current_exception(); }
		if (selected[i].fraction != 0) {
			try { future = next[j].get(); } catch (...) { if (!error) error = std::current_exception(); }
			++j;
		}
		if (!values || (selected[i].fraction != 0 && !future)) continue;
		for (auto const& copy : plan.reads[i].copies) {
			auto value = access(*values, copy.source);
			if (future) value = (1 - selected[i].fraction) * value + selected[i].fraction * access(*future, copy.source);
			ghosts[copy.destination] += copy.weight * value;
		}
	}
	if (error) std::rethrow_exception(error);
}

template <typename State>
void readHalo(storage::ColumnHandle<State> const& fields, HaloPlan const& plan, unsigned bank, std::vector<State>& ghosts,
	std::vector<HaloTime> const& times = {}) {
	readTimedHalo(fields, plan, bank, ghosts, times, [](auto const& v, std::size_t i) { return v.at(i); });
}

template <typename T>
void readHalo(storage::FieldHandle<T> const& field, HaloPlan const& plan, unsigned bank, std::vector<T>& ghosts,
	std::vector<HaloTime> const& times = {}) {
	readTimedHalo(field, plan, bank, ghosts, times, [](auto const& v, std::size_t i) { return v.data()[i]; });
}

/// Complete physical boundary values after all donor reads have finished.
/// Analytic corners use the original position and take precedence over other faces.
template <typename System>
void applyHaloBoundaries(HaloPlan const& plan, std::vector<typename System::State>& ghosts, System const& system, units::Time time,
	physics::AnalyticBoundary<typename System::State> const& analytic = {}) {
	for (auto const& ghost : plan.analyticGhosts)
		ghosts.at(ghost.destination) = physics::evaluateBoundary(analytic, ghost.position, time, system);
	auto transform = [&](std::size_t i) {
		ghosts.at(i) = physics::transformBoundary(ghosts.at(i), plan.reflectionMasks[i], plan.outflowLowerMasks.at(i), plan.outflowUpperMasks.at(i), system);
	};
	for (std::size_t i = plan.ghostCount; i < plan.reflectionMasks.size(); ++i)
		transform(i);
	for (auto const& interpolation : plan.prolongations) {
		std::array<typename System::State, ndim> slopes;
		auto const& center = ghosts.at(interpolation.center);
		for (int d = 0; d < ndim; ++d)
			slopes[d] = amr::slope(ghosts.at(interpolation.left[d]), center, ghosts.at(interpolation.right[d]));
		ghosts.at(interpolation.destination) = amr::interpolate(center, slopes, interpolation.offset, system);
	}
	for (std::size_t i = 0; i < plan.ghostCount; ++i)
		transform(i);
}

}	 // namespace octotigerII
