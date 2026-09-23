#include <gtest/gtest.h>
#include <limits>
#include <numeric>
#include <set>
#include "octotigerII/mesh.hpp"
#include "octotigerII/storage/layout.hpp"

using namespace octotigerII;
using namespace octotigerII::mesh;

namespace {

class MeshGhostWidth : public ::testing::TestWithParam<int> {};


TEST_P(MeshGhostWidth, InteriorStorageRoundTripAndTraversal) {
	MeshLayout layout(4, GetParam());
	std::set<std::size_t> visited;
	layout.forEachInterior([&](Coordinates const& cell, std::size_t index) {
		auto const storage = layout.storageCoordinates(cell);
		EXPECT_EQ(layout.interiorCoordinates(storage), cell);
		EXPECT_TRUE(layout.isInterior(storage));
		EXPECT_EQ(layout.index(storage), index);
		EXPECT_TRUE(visited.insert(index).second);
	});
	EXPECT_EQ(visited.size(), layout.interiorCellCount());
	std::size_t expected = 0;
	forEachCoordinate(layout.extents(), [&](Coordinates const& cell) { EXPECT_EQ(layout.index(cell), expected++); });
	EXPECT_EQ(expected, layout.cellCount());
}


INSTANTIATE_TEST_SUITE_P(Halos, MeshGhostWidth, ::testing::Values(0, 1, 2));


TEST(Mesh, FaceIndicesCoverEachFaceExactlyOnce) {
	MeshLayout layout(4, 2);
	for (int axis = 0; axis < ndim; ++axis) {
		SCOPED_TRACE(axis);
		std::size_t index = 0;
		forEachCoordinate(layout.faceExtents(axis), [&](Coordinates const& face) { EXPECT_EQ(layout.faceIndex(axis, face), index++); });
		EXPECT_EQ(index, layout.faceCount(axis));
		EXPECT_EQ(index, 5 * layout.interiorCellCount() / 4);
		EXPECT_THROW(layout.faceIndex(axis, layout.faceExtents(axis)), std::out_of_range);
	}
}


TEST(Mesh, RejectsInvalidExtentsCoordinatesAndAxes) {
	for (int n : {-2, 0, 1, 3}) EXPECT_THROW(MeshLayout{n}, std::invalid_argument);
	EXPECT_THROW(MeshLayout(4, -1), std::invalid_argument);
	MeshLayout layout(4, 2);
	for (int axis : {-1, ndim}) {
		EXPECT_THROW(layout.extent(axis), std::out_of_range);
		EXPECT_THROW(layout.faceExtents(axis), std::out_of_range);
	}
	EXPECT_THROW(layout.index(filledCoordinates(-1)), std::out_of_range);
	EXPECT_THROW(layout.index(filledCoordinates(8)), std::out_of_range);
	EXPECT_THROW(layout.storageCoordinates(filledCoordinates(4)), std::out_of_range);
	EXPECT_THROW(layout.interiorCoordinates({}), std::out_of_range);
	EXPECT_THROW(linearIndex(filledCoordinates(-1), filledCoordinates(4)), std::out_of_range);
	int visited = 0;
	forEachCoordinate(Coordinates{}, [&](auto const&) { ++visited; });
	EXPECT_EQ(visited, 0);
}


TEST(Mesh, CgsGeometryAndPatchBounds) {
	MeshLayout layout(4, 2);
	PhysicalCoordinates lower;
	lower.fill(units::Length::from_value(-2));
	auto const width = units::Length::from_value(0.5);
	auto const center = layout.cellCenter(lower, width, filledCoordinates(3));
	for (auto x : center) EXPECT_DOUBLE_EQ(units::value(x), -0.25);
	EXPECT_DOUBLE_EQ(units::value(layout.cellMeasure(width)), std::pow(0.5, ndim));
	EXPECT_DOUBLE_EQ(units::value(layout.faceMeasure(width)), std::pow(0.5, ndim - 1));
	for (Real invalid : {Real(0), Real(-1), std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		EXPECT_THROW(layout.cellMeasure(units::Length::from_value(invalid)), std::invalid_argument);
		EXPECT_THROW((PatchData<units::Density>(layout, units::Length::from_value(invalid))), std::invalid_argument);
	}
	PatchData<units::Density> patch(layout, width, lower);
	patch.atInterior({}) = units::Density::from_value(7);
	EXPECT_EQ(patch.atStorage(filledCoordinates(2)), units::Density::from_value(7));
	EXPECT_THROW(patch.atInterior(filledCoordinates(4)), std::out_of_range);
	lower[0] = units::Length::from_value(std::numeric_limits<Real>::quiet_NaN());
	EXPECT_THROW((PatchData<units::Density>(layout, width, lower)), std::invalid_argument);
}


TEST(BlockTopology, ChildBitsParentAndInvalidSlots) {
	BlockLocation root;
	EXPECT_TRUE(root.isRoot());
	EXPECT_THROW(root.parent(), std::logic_error);
	EXPECT_THROW(root.childSlot(), std::logic_error);
	for (int slot = 0; slot < (1 << ndim); ++slot) {
		auto const child = root.child(slot);
		EXPECT_EQ(child.parent(), root);
		EXPECT_EQ(child.childSlot(), slot);
		for (int d = 0; d < ndim; ++d) EXPECT_EQ(child.coordinates[d], (slot >> d) & 1);
		for (int second = 0; second < (1 << ndim); ++second) EXPECT_EQ(child.child(second).parent(), child);
	}
	EXPECT_THROW(root.child(-1), std::invalid_argument);
	EXPECT_THROW(root.child(1 << ndim), std::invalid_argument);
}


TEST(TimeMetadata, StepPublicationSynchronizationAndIntervals) {
	TimeState a, b;
	a.completeStep(units::Time::from_value(0.125));
	a.completeStep(units::Time::from_value(0.25));
	EXPECT_DOUBLE_EQ(units::value(a.time), 0.375);
	EXPECT_DOUBLE_EQ(units::value(a.nextTime()), 0.625);
	EXPECT_EQ(a.step, 2u);
	EXPECT_EQ(a.substep, 2);
	EXPECT_FALSE(a.synchronizedWith(b));
	b.time = a.time;
	EXPECT_TRUE(a.synchronizedWith(b));
	for (Real dt : {Real(0), Real(-1), std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		EXPECT_THROW(a.completeStep(units::Time::from_value(dt)), std::invalid_argument);
		EXPECT_EQ(a.time, b.time);
		EXPECT_EQ(a.step, 2u);
	}
	TimeInterval interval{units::Time::from_value(1), units::Time::from_value(2)};
	EXPECT_EQ(interval.duration(), units::Time::from_value(1));
	EXPECT_TRUE(interval.contains(interval.begin));
	EXPECT_TRUE(interval.contains(interval.end));
	EXPECT_FALSE(interval.contains(units::Time::from_value(0.9)));
	std::swap(interval.begin, interval.end);
	EXPECT_THROW(interval.duration(), std::logic_error);
}


TEST(StorageLayout, RecordsAreContiguousNonoverlappingAndCapacityIsExact) {
	storage::Layout layout({3, 7, 1, 9, 2}, 3);
	std::vector<std::size_t> sizes(3);
	for (auto range : layout.ranges()) {
		ASSERT_LT(range.partition, sizes.size());
		EXPECT_EQ(range.offset, sizes[range.partition]);
		sizes[range.partition] += range.count;
		EXPECT_NO_THROW(layout.validate(range));
	}
	EXPECT_EQ(std::accumulate(sizes.begin(), sizes.end(), std::size_t(0)), 22u);
	for (std::size_t p = 0; p < sizes.size(); ++p) EXPECT_EQ(layout.capacity(p), sizes[p]);
	auto const addresses = layout.addressSpace();
	EXPECT_TRUE(addresses.ranges().empty());
	EXPECT_EQ(addresses.capacity(), layout.capacity());
	EXPECT_THROW(layout.validate({3, 0, 1}), std::out_of_range);
	EXPECT_THROW(layout.validate({0, layout.capacity(0), 1}), std::out_of_range);
}


TEST(StorageLayout, EmptyPartitionsOverflowAndSlices) {
	EXPECT_THROW((storage::Layout({1}, 0)), std::invalid_argument);
	EXPECT_THROW((storage::Layout({0}, 1)), std::invalid_argument);
	EXPECT_THROW((storage::Layout({std::numeric_limits<std::size_t>::max(), 1}, 1)), std::overflow_error);
	storage::Layout empty({}, 3);
	EXPECT_EQ(empty.capacity(), 0u);
	EXPECT_EQ(empty.partitionCount(), 3u);
	storage::Range range{2, 7, 4};
	auto const slice = range.slice(1, 2);
	EXPECT_EQ(slice.partition, 2u);
	EXPECT_EQ(slice.offset, 8u);
	EXPECT_EQ(slice.count, 2u);
	EXPECT_NO_THROW(range.slice(4, 0));
	EXPECT_THROW(range.slice(3, 2), std::out_of_range);
	EXPECT_THROW(range.slice(std::numeric_limits<std::size_t>::max(), 1), std::out_of_range);
}

} // namespace
