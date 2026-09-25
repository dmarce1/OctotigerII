#include "octotigerII/subgrid/topology.hpp"
#include <algorithm>
#include <bit>
#include <limits>
#include <map>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include "octotigerII/profiling.hpp"

namespace octotigerII {
namespace {
	using Location = mesh::BlockLocation;
	std::vector<Location> uniformLeaves(Config const& c) {
		std::vector<Location> result;
		mesh::forEachCoordinate(mesh::filledCoordinates(1 << c.mesh.level), [&](auto const& coordinate) { result.push_back({c.mesh.level, coordinate}); });
		return result;
	}
	class CellAddress {
	public:
		std::size_t block;
		mesh::Coordinates cell;
		int level;
	};
	class Directory {
	public:
		Directory(Config const& config, std::vector<Subgrid> const& blocks)
		  : blocks_(blocks)
		  , bits_(std::countr_zero(unsigned(config.mesh.cells))) {
			for (std::size_t i = 0; i < blocks.size(); ++i)
				leaves_.emplace(blocks[i].location, i);
		}
		std::optional<CellAddress> locate(Location const& cell) const {
			if (cell.level < bits_) return {};
			Location block{cell.level - bits_, {}};
			for (int d = 0; d < ndim; ++d)
				block.coordinates[d] = cell.coordinates[d] >> bits_;
			for (;;) {
				if (auto found = leaves_.find(block); found != leaves_.end()) {
					mesh::Coordinates local{};
					for (int d = 0; d < ndim; ++d)
						local[d] = (cell.coordinates[d] >> (cell.level - block.level - bits_)) - (block.coordinates[d] << bits_);
					return CellAddress{found->second, local, block.level + bits_};
				}
				if (block.isRoot()) return {};
				block = block.parent();
			}
		}
		template <typename F>
		void average(Location const& cell, Real weight, F const& f) const {
			if (auto source = locate(cell)) {
				f(*source, weight);
				return;
			}
			if (cell.level >= 24) throw std::logic_error("Incomplete AMR leaf directory");
			for (int slot = 0; slot < (1 << ndim); ++slot)
				average(cell.child(slot), weight / Real(1 << ndim), f);
		}
		int bits() const {
			return bits_;
		}

	private:
		std::vector<Subgrid> const& blocks_;
		int bits_;
		std::unordered_map<Location, std::size_t, mesh::BlockLocationHash> leaves_;
	};
}	 // namespace

std::size_t allFaceCount(int cells) {
	return ndim * mesh::MeshLayout(cells).faceCount(0);
}
std::size_t allFaceIndex(mesh::MeshLayout const& layout, int axis, mesh::Coordinates const& face) {
	return axis * layout.faceCount(0) + layout.faceIndex(axis, face);
}
std::size_t boundaryFluxCount(int cells) {
	std::size_t face = 1;
	for (int d = 1; d < ndim; ++d)
		face *= cells;
	return 2 * ndim * face;
}
std::size_t boundaryFluxIndex(int cells, int axis, bool upper, mesh::Coordinates const& cell) {
	std::size_t index = 0, stride = 1;
	for (int d = 0; d < ndim; ++d)
		if (d != axis) {
			index += stride * cell[d];
			stride *= cells;
		}
	return index + (2 * axis + int(upper)) * stride;
}

CartesianTopology::CartesianTopology(Config const& config, std::size_t partitions)
  : CartesianTopology(config, partitions, uniformLeaves(config)) {}

CartesianTopology::CartesianTopology(Config const& config, std::size_t partitions, std::vector<Location> const& leaves)
  : config_(config) {
	config.validate();
	if (leaves.empty()) throw std::invalid_argument("An AMR mesh cannot be empty");
	std::unordered_set<Location, mesh::BlockLocationHash> nodes, leafSet;
	Real volume = 0;
	for (auto location : leaves) {
		if (location.level < 0 || location.level > 16 || !leafSet.insert(location).second) throw std::invalid_argument("Invalid or duplicate AMR leaf");
		for (auto x : location.coordinates)
			if (x < 0 || x >= (1 << location.level)) throw std::invalid_argument("AMR block outside domain");
		Real v = 1;
		for (int d = 0; d < ndim; ++d)
			v /= Real(1 << location.level);
		volume += v;
		while (!location.isRoot()) {
			location = location.parent();
			nodes.insert(location);
		}
	}
	for (auto const& leaf : leaves)
		if (nodes.contains(leaf)) throw std::invalid_argument("Overlapping AMR leaves");
	if (volume < 1 - 64 * epsilonR || volume > 1 + 64 * epsilonR) throw std::invalid_argument("AMR leaves do not cover the domain");
	mesh::MeshLayout const cells(config.mesh.cells);
	auto ordered = leaves;
	if (config.amr.enabled) std::sort(ordered.begin(), ordered.end(), [](auto const& a, auto const& b) { return mesh::mortonKey(a) < mesh::mortonKey(b); });
	for (auto const& location : ordered) {
		Subgrid block;
		block.id = blocks_.size();
		block.location = location;
		block.layout = cells;
		block.cellWidth = (config.mesh.upper - config.mesh.lower) / Real((1 << location.level) * config.mesh.cells);
		for (int d = 0; d < ndim; ++d)
			block.lower[d] = config.mesh.lower + Real(location.coordinates[d] * config.mesh.cells) * block.cellWidth;
		blocks_.push_back(block);
	}
	layout_ = storage::Layout(std::vector<std::size_t>(blocks_.size(), cells.interiorCellCount()), partitions);
	storage::Layout const fluxLayout(std::vector<std::size_t>(blocks_.size(), boundaryFluxCount(config.mesh.cells)), partitions);
	storage::Layout const massFluxLayout(std::vector<std::size_t>(blocks_.size(), allFaceCount(config.mesh.cells)), partitions);
	for (std::size_t i = 0; i < blocks_.size(); ++i) {
		blocks_[i].interior = layout_.ranges()[i];
		blocks_[i].boundaryFlux = fluxLayout.ranges()[i];
		blocks_[i].massFlux = massFluxLayout.ranges()[i];
	}
}

HaloPlan makeHaloPlan(Config const& config, std::vector<Subgrid> const& blocks, std::size_t index) {
	profiling::Region profile("mesh.halo_plan");
	auto const& block = blocks.at(index);
	Directory const directory(config, blocks);
	int const level = block.location.level + directory.bits();
	mesh::MeshLayout const padded(config.mesh.cells, 2);
	HaloPlan plan;
	plan.ghostIndices.resize(padded.cellCount(), std::numeric_limits<std::size_t>::max());
	mesh::forEachCoordinate(padded.extents(), [&](auto const& cell) {
		if (!padded.isInterior(cell)) plan.ghostIndices[padded.index(cell)] = plan.ghostCount++;
	});
	plan.valueCount = plan.ghostCount;
	using Address = std::pair<std::size_t, std::size_t>;
	std::map<Address, std::vector<HaloCopy>> donors;
	std::map<Address, int> donorLevels;
	auto describe = [&](Location const& cell, std::size_t destination, bool interpolate) {
		auto const mapped = config.mesh.boundary.map(cell.coordinates, 1 << cell.level);
		plan.reflectionMasks.resize(std::max(plan.reflectionMasks.size(), destination + 1));
		plan.outflowLowerMasks.resize(plan.reflectionMasks.size());
		plan.outflowUpperMasks.resize(plan.reflectionMasks.size());
		plan.reflectionMasks[destination] = mapped.analytic ? 0 : mapped.reflectionMask;
		plan.outflowLowerMasks[destination] = mapped.analytic ? 0 : mapped.outflowLowerMask;
		plan.outflowUpperMasks[destination] = mapped.analytic ? 0 : mapped.outflowUpperMask;
		if (mapped.analytic) {
			mesh::PhysicalCoordinates position{};
			auto const h = (config.mesh.upper - config.mesh.lower) / Real(1 << cell.level);
			for (int d = 0; d < ndim; ++d)
				position[d] = config.mesh.lower + (cell.coordinates[d] + 0.5) * h;
			plan.analyticGhosts.push_back({destination, position});
			return std::optional<Location>{};
		}
		Location const source{cell.level, mapped.source};
		if (interpolate)
			if (auto donor = directory.locate(source); donor && donor->level < source.level) return std::optional<Location>{source};
		directory.average(source, 1, [&](CellAddress const& donor, Real weight) {
			auto const& b = blocks[donor.block];
			Address const address{b.interior.partition, b.interior.offset + b.layout.index(donor.cell)};
			donors[address].push_back({0, destination, weight});
			donorLevels[address] = b.location.level;
		});
		return std::optional<Location>{};
	};
	mesh::forEachCoordinate(padded.extents(), [&](auto const& cell) {
		if (padded.isInterior(cell)) return;
		auto const destination = plan.ghostIndices[padded.index(cell)];
		Location global{level, {}};
		for (int d = 0; d < ndim; ++d)
			global.coordinates[d] = block.location.coordinates[d] * config.mesh.cells + cell[d] - 2;
		auto const fine = describe(global, destination, true);
		if (!fine) return;
		auto const donor = directory.locate(*fine).value();
		Location coarse = *fine;
		while (coarse.level > donor.level)
			coarse = coarse.parent();
		HaloPlan::Prolongation p;
		p.destination = destination;
		p.center = plan.valueCount++;
		describe(coarse, p.center, false);
		Real const ratio = Real(1 << (fine->level - coarse.level));
		for (int d = 0; d < ndim; ++d) {
			auto left = coarse, right = coarse;
			--left.coordinates[d];
			++right.coordinates[d];
			p.left[d] = plan.valueCount++;
			p.right[d] = plan.valueCount++;
			describe(left, p.left[d], false);
			describe(right, p.right[d], false);
			p.offset[d] = (fine->coordinates[d] + 0.5) / ratio - (coarse.coordinates[d] + 0.5);
		}
		plan.prolongations.push_back(p);
	});
	for (auto const& [address, copies] : donors) {
		auto const [partition, offset] = address;
		int const donorLevel = donorLevels.at(address);
		if (plan.reads.empty() || plan.reads.back().level != donorLevel || plan.reads.back().range.partition != partition || plan.reads.back().range.offset + plan.reads.back().range.count != offset)
			plan.reads.push_back({{partition, offset, 0}, {}, donorLevel});
		auto& read = plan.reads.back();
		for (auto copy : copies) {
			copy.source = read.range.count;
			read.copies.push_back(copy);
		}
		++read.range.count;
	}
	return plan;
}

std::vector<FluxCorrection> makeRefluxPlan(Config const& c, std::vector<Subgrid> const& blocks, std::size_t id) {
	Directory const directory(c, blocks);
	auto const& block = blocks.at(id);
	int const level = block.location.level + directory.bits(), n = c.mesh.cells;
	std::vector<FluxCorrection> result;
	for (int axis = 0; axis < ndim; ++axis)
		for (bool upper : {false, true}) {
			auto extents = mesh::filledCoordinates(n);
			extents[axis] = 1;
			mesh::forEachCoordinate(extents, [&](auto cell) {
				cell[axis] = upper ? n - 1 : 0;
				Location neighbor{level, {}};
				for (int d = 0; d < ndim; ++d)
					neighbor.coordinates[d] = block.location.coordinates[d] * n + cell[d];
				neighbor.coordinates[axis] += upper ? 1 : -1;
				if ((neighbor.coordinates[axis] < 0 || neighbor.coordinates[axis] >= (1 << level)) && !c.mesh.boundary.periodic(axis)) return;
				neighbor.coordinates = c.mesh.boundary.map(neighbor.coordinates, 1 << level).source;
				if (directory.locate(neighbor)) return;
				FluxCorrection correction;
				correction.cell = block.layout.index(cell);
				correction.coarseFlux = boundaryFluxIndex(n, axis, upper, cell);
				correction.sign = upper ? -1 : 1;
				for (int slot = 0; slot < (1 << ndim); ++slot) {
					if (bool(slot & (1 << axis)) != !upper) continue;
					auto const fine = neighbor.child(slot);
					auto const donor = directory.locate(fine);
					if (!donor || donor->level != fine.level) throw std::logic_error("Refluxing requires a 2:1 balanced mesh");
					auto const& source = blocks.at(donor->block);
					correction.fineFluxes.push_back(source.boundaryFlux.slice(boundaryFluxIndex(n, axis, !upper, donor->cell), 1));
				}
				result.push_back(std::move(correction));
			});
		}
	return result;
}

std::vector<GravityWorkFace> makeGravityWorkPlan(Config const& c, std::vector<Subgrid> const& blocks, std::size_t id) {
	Directory const directory(c, blocks);
	auto const& block = blocks.at(id);
	int const level = block.location.level + directory.bits(), n = c.mesh.cells;
	std::vector<GravityWorkFace> result;
	for (int axis = 0; axis < ndim; ++axis) for (bool upper : {false, true}) {
		auto extents = mesh::filledCoordinates(n); extents[axis] = 1;
		mesh::forEachCoordinate(extents, [&](auto cell) {
			cell[axis] = upper ? n - 1 : 0;
			auto face = cell; if (upper) ++face[axis];
			GravityWorkFace entry;
			entry.cell = block.layout.index(cell); entry.axis = axis; entry.sign = upper ? 1 : -1;
			entry.massFlux = block.massFlux.slice(allFaceIndex(block.layout, axis, face), 1);
			Location neighbor{level, {}};
			for (int d = 0; d < ndim; ++d) neighbor.coordinates[d] = block.location.coordinates[d] * n + cell[d];
			neighbor.coordinates[axis] += entry.sign;
			if ((neighbor.coordinates[axis] < 0 || neighbor.coordinates[axis] >= (1 << level)) && !c.mesh.boundary.periodic(axis)) {
				entry.physical = true; result.push_back(entry); return;
			}
			neighbor.coordinates = c.mesh.boundary.map(neighbor.coordinates, 1 << level).source;
			if (auto donor = directory.locate(neighbor)) {
				auto const& source = blocks.at(donor->block);
				entry.neighbor = source.interior.slice(source.layout.index(donor->cell), 1);
				entry.potentialFraction = Real(block.cellWidth / (block.cellWidth + source.cellWidth));
				result.push_back(entry);
			} else {
				for (int slot = 0; slot < (1 << ndim); ++slot) {
					if (bool(slot & (1 << axis)) != !upper) continue;
					auto const fine = neighbor.child(slot);
					auto const donor = directory.locate(fine);
					if (!donor || donor->level != fine.level) throw std::logic_error("Gravity work requires a 2:1 balanced mesh");
					auto const& source = blocks.at(donor->block);
					auto f = donor->cell; if (!upper) ++f[axis];
					auto part = entry;
					part.neighbor = source.interior.slice(source.layout.index(donor->cell), 1);
					part.massFlux = source.massFlux.slice(allFaceIndex(source.layout, axis, f), 1);
					part.areaFraction = Real(1) / (1 << (ndim - 1));
					part.potentialFraction = Real(block.cellWidth / (block.cellWidth + source.cellWidth));
					result.push_back(part);
				}
			}
		});
	}
	return result;
}

HaloPlan CartesianTopology::halo(std::size_t index) const {
	return makeHaloPlan(config_, blocks_, index);
}
}	 // namespace octotigerII
