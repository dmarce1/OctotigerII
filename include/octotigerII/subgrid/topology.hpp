/** @file
 * @brief Lightweight block geometry and reusable halo plans.
 * @ingroup mesh
 */
// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#pragma once

#include "octotigerII/config.hpp"
#include "octotigerII/mesh.hpp"
#include "octotigerII/storage/layout.hpp"


namespace octotigerII {


// A persistent block is geometry and an address. It owns no field arrays.
/// Persistent block descriptor containing geometry and an interior address only.
/// @ingroup mesh
class Subgrid {
public:

	std::uint64_t id = 0;
	mesh::BlockLocation location;
	mesh::MeshLayout layout;
	mesh::PhysicalCoordinates lower{};
	units::Length cellWidth{};
	storage::Range interior;

	/// Serialize this value with its compile-time quantity types preserved.
	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & id & location & layout & cellWidth & interior;
		for (auto& coordinate : lower)
			archive & coordinate;
	}
};


/// One element mapping from a fetched contiguous range to a compact ghost buffer.
/// @ingroup mesh
class HaloCopy {
public:

	std::size_t source = 0;
	std::size_t destination = 0;

	/// Serialize this value with its compile-time quantity types preserved.
	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & source & destination;
	}
};


/// Contiguous source range plus its destination mappings in the halo.
/// @ingroup mesh
class HaloRead {
public:

	storage::Range range;
	std::vector<HaloCopy> copies;

	/// Serialize this value with its compile-time quantity types preserved.
	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & range & copies;
	}
};


/// A prescribed ghost has no donor read and is evaluated at this physical position.
class AnalyticGhost {
public:

	std::size_t destination = 0;
	mesh::PhysicalCoordinates position{};

	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & destination;
		for (auto& coordinate : position)
			archive & coordinate;
	}
};


// Read-only geometry, reused for every stage. Each read is a contiguous run
// of required source cells. Halos do not include the block's interior.
/// Reusable geometry plan for a block's halo.
/// ghostIndices maps padded indices to compact halo slots; interior indices use
/// a sentinel. The plan contains no stage values and may be reused between updates.
/// @ingroup mesh
class HaloPlan {
public:

	std::vector<HaloRead> reads;
	// Padded index -> compact ghost index; interiors contain the sentinel.
	std::vector<std::size_t> ghostIndices;
	std::size_t ghostCount = 0;
	std::vector<unsigned> reflectionMasks;
	std::vector<AnalyticGhost> analyticGhosts;

	/// Serialize this value with its compile-time quantity types preserved.
	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & reads & ghostIndices & ghostCount & reflectionMasks & analyticGhosts;
	}
};


/// Map the fixed-level Cartesian halo to contiguous source runs.
/// Physical boundaries use per-face wrapping, clamping, reflection, or analytic data.
HaloPlan makeHaloPlan(Config const& config, std::vector<Subgrid> const& blocks, std::size_t index);


// Current mesh adapter: fixed, uniform Cartesian blocks. Only this adapter
// interprets mesh.level. The field store accepts arbitrary range lengths.
/// Current fixed-level Cartesian mesh adapter.
/// Supplies descriptors and halo dependencies to an otherwise topology-independent
/// store. Mixed-level AMR and refinement operations are not implemented here.
/// @ingroup mesh
class CartesianTopology {
public:

	CartesianTopology(Config const& config, std::size_t partitions);

	/// Return the current ordered block descriptors.
	std::vector<Subgrid> const& blocks() const {
		return blocks_;
	}

	/// Return the record placement corresponding to the block directory.
	storage::Layout const& storageLayout() const {
		return layout_;
	}

	/// Build the required source ranges and destination mappings for one block.
	HaloPlan halo(std::size_t block) const;

private:

	Config config_;
	storage::Layout layout_;
	std::vector<Subgrid> blocks_;
};


}	 // namespace octotigerII
