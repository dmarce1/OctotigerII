/** @file
 * @brief Cartesian indexing, physical geometry, and patch-local time metadata.
 * @ingroup mesh
 */
// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#pragma once

#include "octotigerII/buildConfig.hpp"
#include "octotigerII/units/cgs.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace octotigerII::mesh {

using Coordinates = std::array<int, ndim>;
using PhysicalCoordinates = std::array<units::Length, ndim>;

/// Construct equal extents in exactly ndim directions.
inline Coordinates filledCoordinates(int value) {
	Coordinates result{};
	result.fill(value);
	return result;
}

/// Flatten coordinates with x contiguous, without inactive axes.
inline std::size_t linearIndex(Coordinates const& cell, Coordinates const& extents) {
	std::size_t index = 0, stride = 1;
	for (int axis = 0; axis < ndim; ++axis) {
		if (cell[axis] < 0 || cell[axis] >= extents[axis]) throw std::out_of_range("Coordinate outside extent");
		index += static_cast<std::size_t>(cell[axis]) * stride;
		stride *= extents[axis];
	}
	return index;
}

/// Visit the half-open coordinate box in x-contiguous order.
/// No singleton inactive dimensions are represented.
template <typename Function>
void forEachCoordinate(Coordinates const& begin, Coordinates const& end, Function&& function) {
	for (int axis = 0; axis < ndim; ++axis)
		if (end[axis] <= begin[axis]) return;
	Coordinates cell = begin;
	while (true) {
		function(cell);
		int axis = 0;
		for (; axis < ndim; ++axis) {
			if (++cell[axis] < end[axis]) break;
			cell[axis] = begin[axis];
		}
		if (axis == ndim) return;
	}
}

template <typename Function>
void forEachCoordinate(Coordinates const& extents, Function&& function) {
	forEachCoordinate(Coordinates{}, extents, std::forward<Function>(function));
}

/// Indexing for an ndim-dimensional Cartesian patch with x contiguous.
/// CGS integrated totals use unit transverse measure in reduced dimensions;
/// that convention does not add coordinate axes, cells, or vector components.
/// @ingroup mesh
class MeshLayout {
public:
	MeshLayout()
	  : MeshLayout(2) {}

	explicit MeshLayout(int cellsPerActiveDimension, int ghostWidth = 0);

	int cellsPerActiveDimension() const;

	int ghostWidth() const;

	int interiorExtent(int axis) const;

	int extent(int axis) const;

	Coordinates interiorExtents() const;

	Coordinates extents() const;

	std::size_t interiorCellCount() const;

	std::size_t cellCount() const;

	int childCount() const;

	bool isActiveChild(int childSlot) const;

	std::vector<int> activeChildSlots() const;

	std::size_t index(Coordinates const& storageCoordinates) const;

	Coordinates storageCoordinates(Coordinates const& interiorCoordinates) const;

	Coordinates interiorCoordinates(Coordinates const& storageCoordinates) const;

	bool isInterior(Coordinates const& storageCoordinates) const;

	Coordinates faceExtents(int normal) const;

	std::size_t faceCount(int normal) const;

	std::size_t faceIndex(int normal, Coordinates const& faceCoordinates) const;

	units::Volume cellMeasure(units::Length cellWidth) const;

	units::Quantity<2, 0, 0> faceMeasure(units::Length cellWidth) const;

	PhysicalCoordinates cellCenter(PhysicalCoordinates const& lower, units::Length cellWidth, Coordinates const& interiorCoordinates) const;

	template <typename Function>
	void forEachInterior(Function&& function) const {
		forEachCoordinate(interiorExtents_, [&](Coordinates const& cell) { function(cell, index(storageCoordinates(cell))); });
	}

	/// Serialize this value with its compile-time quantity types preserved.
	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & cellsPerActiveDimension_;
		archive & ghostWidth_;
		rebuild();
	}

private:
	int cellsPerActiveDimension_ = 2;
	int ghostWidth_ = 0;
	Coordinates interiorExtents_{};
	Coordinates extents_{};
	Coordinates strides_{};

	void rebuild();

	void validateAxis(int axis) const;

	static std::size_t product(Coordinates const& extents);
};

/// Topological identity (level and ndim integer coordinates).
/// This value carries no ownership of numerical field arrays.
/// @ingroup mesh
class BlockLocation {
public:
	int level = 0;
	Coordinates coordinates{};

	bool isRoot() const;

	BlockLocation parent() const;

	int childSlot() const;

	BlockLocation child(int slot) const;

	/// Serialize this value with its compile-time quantity types preserved.
	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & level;
		// Inline metadata is deserialized before its enclosing value may move.
		for (auto& coordinate : coordinates)
			archive & coordinate;
	}

	friend bool operator==(BlockLocation const&, BlockLocation const&) = default;
};

/// Hash of the block location for topology directories.
/// @ingroup mesh
class BlockLocationHash {
public:
	std::size_t operator()(BlockLocation const& location) const noexcept;
};

/// Morton order on the common finest lattice; ancestors precede descendants.
/// Sixteen block levels fit in 48 bits in three dimensions.
inline std::uint64_t mortonKey(BlockLocation const& location) {
	if (location.level < 0 || location.level > 16) throw std::invalid_argument("Morton block level outside [0,16]");
	std::uint64_t key = 0;
	for (int d = 0; d < ndim; ++d) {
		auto const coordinate = std::uint64_t(location.coordinates[d]) << (16 - location.level);
		for (int bit = 0; bit < 16; ++bit)
			key |= ((coordinate >> bit) & 1) << (ndim * bit + d);
	}
	return key;
}

// Time is patch-local even while the first implementation advances every
// level synchronously. This deliberately leaves room for temporal refinement:
// boundary and reflux operations can be tagged by physical time intervals,
// rather than assuming one global cycle number.
/// Physical time and integration metadata for a patch.
/// The current runtime advances all blocks synchronously; temporalLevel/substep
/// are metadata for future temporal refinement, not an implemented subcycling driver.
/// @ingroup mesh
class TimeState {
public:
	units::Time time{};
	units::Time stepSize{};
	std::uint64_t step = 0;
	int temporalLevel = 0;
	int substep = 0;

	units::Time nextTime() const;

	bool synchronizedWith(TimeState const& other, Real tolerance = 64 * epsilonR) const;

	void completeStep(units::Time completedStepSize);

	/// Serialize this value with its compile-time quantity types preserved.
	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & time;
		archive & stepSize;
		archive & step;
		archive & temporalLevel;
		archive & substep;
	}
};

/// Closed physical time interval used to label fluxes and boundary data.
/// @ingroup mesh
class TimeInterval {
public:
	units::Time begin{};
	units::Time end{};

	units::Time duration() const;

	bool contains(units::Time sampleTime, Real tolerance = 64 * epsilonR) const;

	/// Serialize this value with its compile-time quantity types preserved.
	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & begin;
		archive & end;
	}
};

/// Owning patch value used for temporary snapshots and standalone kernel tests.
/// Persistent distributed interiors live in storage::StoragePartition instead.
/// Values use MeshLayout indexing; this type may include ghosts when used in tests.
/// @ingroup mesh
template <typename State>
class PatchData {
public:
	PatchData() = default;

	PatchData(MeshLayout layout, units::Length cellWidth, PhysicalCoordinates lower = {})
	  : layout_(std::move(layout))
	  , cellWidth_(cellWidth)
	  , lower_(lower)
	  , values_(layout_.cellCount()) {
		if (!(cellWidth_ > units::Length{}) || !units::finite(cellWidth_)) {
			throw std::invalid_argument("Patch cell width must be positive and finite");
		}
		for (int axis = 0; axis < ndim; ++axis) {
			if (!units::finite(lower_[axis])) {
				throw std::invalid_argument("Patch bounds must be finite");
			}
		}
	}

	/// Return the indexing layout represented by this object.
	MeshLayout const& layout() const {
		return layout_;
	}

	/// Return the physical width of one active cell in centimeters.
	units::Length cellWidth() const {
		return cellWidth_;
	}

	/// Return the physical lower corner in centimeters.
	PhysicalCoordinates const& lower() const {
		return lower_;
	}

	TimeState const& timeState() const {
		return timeState_;
	}

	TimeState& timeState() {
		return timeState_;
	}

	std::vector<State> const& values() const {
		return values_;
	}

	std::vector<State>& values() {
		return values_;
	}

	/// Read a cell using the coordinates of the padded storage layout.
	State const& atStorage(Coordinates const& coordinates) const {
		return values_.at(layout_.index(coordinates));
	}

	/// Read a cell using the coordinates of the padded storage layout.
	State& atStorage(Coordinates const& coordinates) {
		return values_.at(layout_.index(coordinates));
	}

	/// Read a cell using interior coordinates without a ghost offset.
	State const& atInterior(Coordinates const& coordinates) const {
		return atStorage(layout_.storageCoordinates(coordinates));
	}

	/// Read a cell using interior coordinates without a ghost offset.
	State& atInterior(Coordinates const& coordinates) {
		return atStorage(layout_.storageCoordinates(coordinates));
	}

	/// Serialize this value with its compile-time quantity types preserved.
	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & layout_;
		archive & cellWidth_;
		for (auto& coordinate : lower_)
			archive & coordinate;
		archive & timeState_;
		archive & values_;
	}

private:
	MeshLayout layout_;
	units::Length cellWidth_ = units::Length::from_value(1);
	PhysicalCoordinates lower_{};
	TimeState timeState_{};
	std::vector<State> values_;
};

}	 // namespace octotigerII::mesh
