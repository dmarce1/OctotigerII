#include "octotigerII/subgrid/topology.hpp"
#include "octotigerII/profiling.hpp"
#include <algorithm>
#include <limits>
#include <map>
#include <tuple>


namespace octotigerII {

CartesianTopology::CartesianTopology(Config const& config, std::size_t partitions)
  : config_(config) {
	config.validate();
	int const n = 1 << config.mesh.level;
	mesh::MeshLayout const cells(config.mesh.cells);
	mesh::forEachCoordinate(mesh::filledCoordinates(n), [&](mesh::Coordinates const& coordinate) {
		Subgrid block;
		block.id = blocks_.size();
		block.location = {config.mesh.level, coordinate};
		block.layout = cells;
		block.cellWidth = (config.mesh.upper - config.mesh.lower) / Real(n * config.mesh.cells);
		for (int axis = 0; axis < ndim; ++axis)
			block.lower[axis] = config.mesh.lower + Real(coordinate[axis] * config.mesh.cells) * block.cellWidth;
		blocks_.push_back(block);
	});
	layout_ = storage::Layout(std::vector<std::size_t>(blocks_.size(), cells.interiorCellCount()), partitions);
	for (std::size_t i = 0; i < blocks_.size(); ++i)
		blocks_[i].interior = layout_.ranges()[i];
}

HaloPlan makeHaloPlan(Config const& config_, std::vector<Subgrid> const& blocks_, std::size_t index) {
	profiling::Region profile("mesh.halo_plan");
	auto const& block = blocks_.at(index);
	int const n = 1 << config_.mesh.level;
	int const globalCells = n * config_.mesh.cells;
	mesh::MeshLayout const padded(config_.mesh.cells, 2);
	HaloPlan plan;
	plan.ghostIndices.resize(padded.cellCount(), std::numeric_limits<std::size_t>::max());
	// Group identical source cells and then join adjacent addresses. This is
	// a geometry plan; transport message coalescing remains HPX's responsibility.
	using Address = std::pair<std::size_t, std::size_t>;
	std::map<Address, std::vector<std::size_t>> donors;
	auto const extents = padded.extents();
	mesh::forEachCoordinate(extents, [&](mesh::Coordinates cell) {
		if (padded.isInterior(cell)) return;
		auto const ghost = plan.ghostCount++;
		plan.ghostIndices[padded.index(cell)] = ghost;
		mesh::Coordinates sourceBlock{};
		for (int axis = 0; axis < ndim; ++axis) {
			int coordinate = block.location.coordinates[axis] * config_.mesh.cells + cell[axis] - 2;
			if (config_.mesh.periodic)
				coordinate = (coordinate % globalCells + globalCells) % globalCells;
			else
				coordinate = std::clamp(coordinate, 0, globalCells - 1);
			sourceBlock[axis] = coordinate / config_.mesh.cells;
			cell[axis] = coordinate % config_.mesh.cells;
		}
		auto const sourceId = mesh::linearIndex(sourceBlock, mesh::filledCoordinates(n));
		auto const& source = blocks_.at(sourceId);
		auto const offset = source.interior.offset + source.layout.index(cell);
		donors[{source.interior.partition, offset}].push_back(ghost);
	});
	for (auto const& [address, destinations] : donors) {
		auto const [partition, offset] = address;
		if (plan.reads.empty() || plan.reads.back().range.partition != partition || plan.reads.back().range.offset + plan.reads.back().range.count != offset)
			plan.reads.push_back({{partition, offset, 0}, {}});
		auto& read = plan.reads.back();
		for (auto const destination : destinations)
			read.copies.push_back({read.range.count, destination});
		++read.range.count;
	}
	return plan;
}

HaloPlan CartesianTopology::halo(std::size_t index) const {
	return makeHaloPlan(config_, blocks_, index);
}

}	 // namespace octotigerII
