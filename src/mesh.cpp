// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.

#include "octotigerII/mesh.hpp"

#include <algorithm>
#include <functional>
#include <limits>


namespace octotigerII::mesh {

MeshLayout::MeshLayout(int cellsPerActiveDimension, int ghostWidth)
  : cellsPerActiveDimension_(cellsPerActiveDimension)
  , ghostWidth_(ghostWidth) {
	rebuild();
}

int MeshLayout::cellsPerActiveDimension() const {
	return cellsPerActiveDimension_;
}

int MeshLayout::ghostWidth() const {
	return ghostWidth_;
}

int MeshLayout::interiorExtent(int axis) const {
	validateAxis(axis);
	return interiorExtents_[axis];
}

int MeshLayout::extent(int axis) const {
	validateAxis(axis);
	return extents_[axis];
}

Coordinates MeshLayout::interiorExtents() const {
	return interiorExtents_;
}

Coordinates MeshLayout::extents() const {
	return extents_;
}

std::size_t MeshLayout::interiorCellCount() const {
	return product(interiorExtents_);
}

std::size_t MeshLayout::cellCount() const {
	return product(extents_);
}

int MeshLayout::childCount() const {
	return 1 << ndim;
}

bool MeshLayout::isActiveChild(int childSlot) const {
	return childSlot >= 0 && childSlot < childCount();
}

std::vector<int> MeshLayout::activeChildSlots() const {
	std::vector<int> slots(childCount());
	for (int slot = 0; slot < childCount(); ++slot) {
		slots[slot] = slot;
	}
	return slots;
}

std::size_t MeshLayout::index(Coordinates const& storageCoordinates) const {
	std::size_t result = 0;
	for (int axis = 0; axis < ndim; ++axis) {
		if (storageCoordinates[axis] < 0 || storageCoordinates[axis] >= extents_[axis]) {
			throw std::out_of_range("Mesh storage coordinate is outside the allocated extent");
		}
		result += static_cast<std::size_t>(storageCoordinates[axis]) * strides_[axis];
	}
	return result;
}

Coordinates MeshLayout::storageCoordinates(Coordinates const& interiorCoordinates) const {
	Coordinates result{};
	for (int axis = 0; axis < ndim; ++axis) {
		if (interiorCoordinates[axis] < 0 || interiorCoordinates[axis] >= interiorExtents_[axis]) {
			throw std::out_of_range("Mesh interior coordinate is outside the active extent");
		}
		result[axis] = interiorCoordinates[axis] + ghostWidth_;
	}
	return result;
}

Coordinates MeshLayout::interiorCoordinates(Coordinates const& storageCoordinates) const {
	Coordinates result{};
	if (!isInterior(storageCoordinates)) {
		throw std::out_of_range("Mesh storage coordinate is not an interior cell");
	}
	for (int axis = 0; axis < ndim; ++axis) {
		result[axis] = storageCoordinates[axis] - ghostWidth_;
	}
	return result;
}

bool MeshLayout::isInterior(Coordinates const& storageCoordinates) const {
	for (int axis = 0; axis < ndim; ++axis) {
		if (storageCoordinates[axis] < 0 || storageCoordinates[axis] >= extents_[axis]) {
			return false;
		}
		if (storageCoordinates[axis] < ghostWidth_ || storageCoordinates[axis] >= ghostWidth_ + cellsPerActiveDimension_) {
			return false;
		}
	}
	return true;
}

Coordinates MeshLayout::faceExtents(int normal) const {
	validateAxis(normal);
	Coordinates result = interiorExtents_;
	++result[normal];
	return result;
}

std::size_t MeshLayout::faceCount(int normal) const {
	return product(faceExtents(normal));
}

std::size_t MeshLayout::faceIndex(int normal, Coordinates const& faceCoordinates) const {
	Coordinates const dimensions = faceExtents(normal);
	std::size_t stride = 1;
	std::size_t result = 0;
	for (int axis = 0; axis < ndim; ++axis) {
		if (faceCoordinates[axis] < 0 || faceCoordinates[axis] >= dimensions[axis]) {
			throw std::out_of_range("Face coordinate is outside the face array");
		}
		result += static_cast<std::size_t>(faceCoordinates[axis]) * stride;
		stride *= static_cast<std::size_t>(dimensions[axis]);
	}
	return result;
}

units::Volume MeshLayout::cellMeasure(units::Length width) const {
	if (!(width > units::Length{}) || !units::finite(width)) throw std::invalid_argument("Cell width must be positive and finite");
	// Lower-dimensional totals use one centimeter along each inactive axis.
	auto const cm = units::Length::from_value(1);
	return boost::units::pow<ndim>(width) * boost::units::pow<3 - ndim>(cm);
}

units::Quantity<2, 0, 0> MeshLayout::faceMeasure(units::Length width) const {
	return cellMeasure(width) / width;
}

PhysicalCoordinates MeshLayout::cellCenter(PhysicalCoordinates const& lower, units::Length cellWidth, Coordinates const& interiorCoordinates) const {
	PhysicalCoordinates result{};
	for (int axis = 0; axis < ndim; ++axis) {
		if (interiorCoordinates[axis] < 0 || interiorCoordinates[axis] >= interiorExtents_[axis]) {
			throw std::out_of_range("Cell-center coordinate is outside the interior");
		}
		result[axis] = lower[axis] + (interiorCoordinates[axis] + Real(0.5)) * cellWidth;
	}
	return result;
}

void MeshLayout::rebuild() {
	if (cellsPerActiveDimension_ < 2 || cellsPerActiveDimension_ % 2 != 0) {
		throw std::invalid_argument("Active mesh extent must be a positive even number");
	}
	if (ghostWidth_ < 0) {
		throw std::invalid_argument("Mesh ghost width cannot be negative");
	}
	for (int axis = 0; axis < ndim; ++axis) {
		interiorExtents_[axis] = cellsPerActiveDimension_;
		extents_[axis] = cellsPerActiveDimension_ + 2 * ghostWidth_;
	}
	strides_[0] = 1;
	for (int axis = 1; axis < ndim; ++axis) {
		strides_[axis] = strides_[axis - 1] * extents_[axis - 1];
	}
}

void MeshLayout::validateAxis(int axis) const {
	if (axis < 0 || axis >= ndim) {
		throw std::out_of_range("Mesh axis must lie in [0, ndim)");
	}
}

std::size_t MeshLayout::product(Coordinates const& extents) {
	std::size_t result = 1;
	for (int extent : extents) {
		result *= static_cast<std::size_t>(extent);
	}
	return result;
}

bool BlockLocation::isRoot() const {
	return level == 0;
}

BlockLocation BlockLocation::parent() const {
	if (isRoot()) {
		throw std::logic_error("The root block has no parent");
	}
	if (level < 0) {
		throw std::invalid_argument("Invalid block location");
	}
	BlockLocation result = *this;
	--result.level;
	for (int axis = 0; axis < ndim; ++axis) {
		result.coordinates[axis] = coordinates[axis] / 2;
	}
	return result;
}

int BlockLocation::childSlot() const {
	if (isRoot()) {
		throw std::logic_error("The root block is not a child");
	}
	if (level < 0) {
		throw std::invalid_argument("Invalid block location");
	}
	int result = 0;
	for (int axis = 0; axis < ndim; ++axis) {
		result |= (coordinates[axis] & 1) << axis;
	}
	return result;
}

BlockLocation BlockLocation::child(int slot) const {
	if (level < 0 || level >= std::numeric_limits<unsigned>::digits - 2 || slot < 0 || slot >= (1 << ndim)) {
		throw std::invalid_argument("Invalid child slot for block dimensionality");
	}
	BlockLocation result = *this;
	++result.level;
	for (int axis = 0; axis < ndim; ++axis) {
		result.coordinates[axis] = 2 * coordinates[axis] + ((slot >> axis) & 1);
	}
	return result;
}

std::size_t BlockLocationHash::operator()(BlockLocation const& location) const noexcept {
	std::size_t result = std::hash<int>{}(location.level);
	for (int coordinate : location.coordinates) {
		result ^= std::hash<int>{}(coordinate) + 0x9e3779b9 + (result << 6) + (result >> 2);
	}
	return result;
}

units::Time TimeState::nextTime() const {
	return time + stepSize;
}

bool TimeState::synchronizedWith(TimeState const& other, Real tolerance) const {
	auto const scale = std::max({units::Time::from_value(1), units::abs(time), units::abs(other.time)});
	return units::abs(time - other.time) <= tolerance * scale;
}

void TimeState::completeStep(units::Time completedStepSize) {
	if (!(completedStepSize > units::Time{}) || !units::finite(completedStepSize)) {
		throw std::invalid_argument("Completed timestep must be positive and finite");
	}
	stepSize = completedStepSize;
	time += completedStepSize;
	++step;
	++substep;
}

units::Time TimeInterval::duration() const {
	if (end < begin) {
		throw std::logic_error("Time interval ends before it begins");
	}
	return end - begin;
}

bool TimeInterval::contains(units::Time sampleTime, Real tolerance) const {
	auto const scale = std::max({units::Time::from_value(1), units::abs(begin), units::abs(end), units::abs(sampleTime)});
	return sampleTime >= begin - tolerance * scale && sampleTime <= end + tolerance * scale;
}

}	 // namespace octotigerII::mesh
