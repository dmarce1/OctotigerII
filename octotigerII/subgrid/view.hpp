/** @file
 * @brief Temporary patch views joining interiors to compact halo storage.
 * @ingroup mesh
 */
// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#pragma once

#include <atomic>
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

/// Cumulative donor requests issued by readHalo on this locality. Values count
/// scalar field elements, including repeated requests and AMR time endpoints.
/// Remote payload bytes count numerical buffers only, excluding parcel metadata,
/// requests, reflux, gravity solves, migration, and all non-halo field reads.
/// Failed asynchronous reads remain requests and are included in these totals.
struct HaloReadStatistics {
	std::uint64_t haloCalls = 0;
	std::uint64_t fieldReads = 0;
	std::uint64_t values = 0;
	std::uint64_t remoteFieldReads = 0;
	std::uint64_t remoteValues = 0;
	std::uint64_t remotePayloadBytes = 0;

	HaloReadStatistics& operator+=(HaloReadStatistics const& other) {
		haloCalls += other.haloCalls;
		fieldReads += other.fieldReads;
		values += other.values;
		remoteFieldReads += other.remoteFieldReads;
		remoteValues += other.remoteValues;
		remotePayloadBytes += other.remotePayloadBytes;
		return *this;
	}

	friend HaloReadStatistics operator-(HaloReadStatistics const& after, HaloReadStatistics const& before) {
		return {after.haloCalls - before.haloCalls, after.fieldReads - before.fieldReads, after.values - before.values,
			after.remoteFieldReads - before.remoteFieldReads, after.remoteValues - before.remoteValues,
			after.remotePayloadBytes - before.remotePayloadBytes};
	}

	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & haloCalls & fieldReads & values & remoteFieldReads & remoteValues & remotePayloadBytes;
	}
};

namespace detail {
struct HaloReadCounters {
	std::atomic<std::uint64_t> haloCalls{0}, fieldReads{0}, values{0};
	std::atomic<std::uint64_t> remoteFieldReads{0}, remoteValues{0}, remotePayloadBytes{0};
};
inline HaloReadCounters haloReadCounters;

/// Expand derived columns into the scalar buffers that their reads request.
template <typename T>
HaloReadStatistics haloReadPayload(storage::FieldHandle<T> const& field, storage::Range const& range) {
	HaloReadStatistics result;
	if (!field.sumSources.empty()) {
		for (auto const& source : field.sumSources) { result += haloReadPayload(source, range); }
	} else {
		result.fieldReads = 1;
		result.values = range.count;
		if (!field.local(range)) {
			result.remoteFieldReads = 1;
			result.remoteValues = range.count;
			result.remotePayloadBytes = range.count * sizeof(T);
		}
	}
	return result;
}

template <typename State>
HaloReadStatistics haloReadPayload(storage::ColumnHandle<State> const& fields, storage::Range const& range) {
	HaloReadStatistics result;
	std::apply([&](auto const&... field) { ((result += haloReadPayload(field, range)), ...); }, fields.fields);
	return result;
}

// Aggregate per halo call, so scalar reads do not each contend on atomics.
class HaloReadMeter {
public:
	HaloReadStatistics counts{1};
	~HaloReadMeter() {
		haloReadCounters.haloCalls.fetch_add(counts.haloCalls, std::memory_order_relaxed);
		haloReadCounters.fieldReads.fetch_add(counts.fieldReads, std::memory_order_relaxed);
		haloReadCounters.values.fetch_add(counts.values, std::memory_order_relaxed);
		haloReadCounters.remoteFieldReads.fetch_add(counts.remoteFieldReads, std::memory_order_relaxed);
		haloReadCounters.remoteValues.fetch_add(counts.remoteValues, std::memory_order_relaxed);
		haloReadCounters.remotePayloadBytes.fetch_add(counts.remotePayloadBytes, std::memory_order_relaxed);
		profiling::sample("transport.halo.calls", counts.haloCalls);
		profiling::sample("transport.halo.field_reads", counts.fieldReads);
		profiling::sample("transport.halo.values", counts.values);
		profiling::sample("transport.halo.remote_field_reads", counts.remoteFieldReads);
		profiling::sample("transport.halo.remote_values", counts.remoteValues);
		profiling::sample("transport.halo.remote_payload_bytes", counts.remotePayloadBytes);
	}
};
} // namespace detail

/// Take differences around completed operations; readings during active halo
/// calls need not form a coherent snapshot. Counters are never reset, so one
/// measurement cannot erase another. Distributed totals require all localities.
inline HaloReadStatistics haloReadStatistics() {
	auto const& c = detail::haloReadCounters;
	return {c.haloCalls.load(std::memory_order_relaxed), c.fieldReads.load(std::memory_order_relaxed),
		c.values.load(std::memory_order_relaxed), c.remoteFieldReads.load(std::memory_order_relaxed),
		c.remoteValues.load(std::memory_order_relaxed), c.remotePayloadBytes.load(std::memory_order_relaxed)};
}

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

	/// Logical lower corner; the rotating frame maps it to inertial coordinates.
	mesh::PhysicalCoordinates const& lower() const { return block_.lower; }

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
	detail::HaloReadMeter meter;
	using Pending = decltype(field.read(storage::Range{}, bank));
	std::vector<Pending> old, next;
	std::vector<HaloTime> selected;
	for (auto const& read : plan.reads) {
		auto const t = times.empty() ? HaloTime{bank, 0} : times.at(read.level);
		selected.push_back(t);
		old.push_back(field.read(read.range, t.bank));
		meter.counts += detail::haloReadPayload(field, read.range);
		// A zero fraction needs no new endpoint, which may still be unwritten.
		if (t.fraction != 0) {
			next.push_back(field.read(read.range, t.nextBank == ~0u ? (t.bank ^ 1) : t.nextBank));
			meter.counts += detail::haloReadPayload(field, read.range);
		}
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
	physics::AnalyticBoundary<typename System::State> const& analytic = {}, physics::RotatingFrame const& frame = physics::RotatingFrame{}) {
	for (auto const& ghost : plan.analyticGhosts)
		ghosts.at(ghost.destination) = physics::evaluateBoundary(analytic, ghost.position, time, system);
	auto transform = [&](std::size_t i) {
		ghosts.at(i) = frame.active()
			? physics::transformBoundary(ghosts.at(i), plan.reflectionMasks[i], plan.outflowLowerMasks.at(i), plan.outflowUpperMasks.at(i),
				system, frame, plan.boundaryPositions.at(i), time)
			: physics::transformBoundary(ghosts.at(i), plan.reflectionMasks[i], plan.outflowLowerMasks.at(i), plan.outflowUpperMasks.at(i), system);
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
