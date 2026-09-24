/** @file
 * @brief Indivisible record ranges and independently sized storage partitions.
 * @ingroup storage
 */
// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>


namespace octotigerII::storage {


// Logical addresses contain no mesh coordinates, refinement level, or topology.
/// Half-open [offset, offset+count) range within one storage partition.
/// Offsets count typed elements, not bytes. The owning FieldHandle supplies identity.
/// @ingroup storage
class Range {
public:

	std::size_t partition = 0;
	std::size_t offset = 0;
	std::size_t count = 0;

	/// Return a bounds-checked subrange relative to this range.
	Range slice(std::size_t begin, std::size_t size) const {
		if (begin > count || size > count - begin) throw std::out_of_range("Field slice exceeds its allocation");
		return {partition, offset + begin, size};
	}

	/// Serialize this value with its compile-time quantity types preserved.
	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & partition & offset & count;
	}
};


// Independently sized partitions contain arbitrary, indivisible record ranges.
// Storage bank count is a field policy, independent of this address layout.
/// Immutable record placement independent of any mesh or particle representation.
/// Records are indivisible, assigned in contiguous groups, and need not have equal
/// lengths. Partition capacities are their actual element totals, without equal-size padding.
/// @ingroup storage
class Layout {
public:

	Layout() = default;

	Layout(std::vector<std::size_t> const& counts, std::size_t partitions) {
		if (partitions == 0) throw std::invalid_argument("Storage needs at least one partition");
		std::vector<std::size_t> sizes(partitions, 0);
		for (std::size_t i = 0; i < counts.size(); ++i) {
			if (!counts[i]) throw std::invalid_argument("Empty record allocation");
			// Contiguous chunks preserve the order chosen by the topology adapter.
			auto const owner = std::min(partitions - 1, i * partitions / counts.size());
			if (counts[i] > std::numeric_limits<std::size_t>::max() - sizes[owner]) throw std::overflow_error("Field allocation overflow");
			ranges_.push_back({owner, sizes[owner], counts[i]});
			sizes[owner] += counts[i];
		}
		capacities_ = std::move(sizes);
	}

	std::vector<Range> const& ranges() const {
		return ranges_;
	}

	/// Return the number of storage partitions.
	std::size_t partitionCount() const {
		return capacities_.size();
	}

	/// Return a partition capacity in elements; the no-argument overload returns the maximum.
	std::size_t capacity() const {
		return capacities_.empty() ? 0 : *std::max_element(capacities_.begin(), capacities_.end());
	}

	/// Return a partition capacity in elements; the no-argument overload returns the maximum.
	std::size_t capacity(std::size_t partition) const {
		return capacities_.at(partition);
	}

	/// Copy partition capacities without copying the record directory.
	Layout addressSpace() const {
		Layout result;
		result.capacities_ = capacities_;
		return result;
	}

	/// Reject invalid dimensions, capacities, or physical configuration values.
	void validate(Range const& range) const {
		auto const size = capacity(range.partition);
		if (range.offset > size || range.count > size - range.offset) throw std::out_of_range("Invalid field range");
	}

	/// Serialize this value with its compile-time quantity types preserved.
	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & capacities_ & ranges_;
	}

private:

	std::vector<std::size_t> capacities_;
	std::vector<Range> ranges_;
};


}	 // namespace octotigerII::storage
