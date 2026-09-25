/** @file
 * @brief Named field identities and shared partition handles, independent of topology.
 * @ingroup storage
 */
// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#pragma once

#include <exception>
#include <limits>
#include "octotigerII/storage/layout.hpp"
#include "octotigerII/storage/partition.hpp"
#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/runtime_distributed/find_here.hpp>
#include <hpx/synchronization/mutex.hpp>
#endif


namespace octotigerII::storage {


/// Shared component directory. Several storage partitions may occupy one locality.
/// Fields allocated through copies of this object share a non-reusing identity counter.
/// @ingroup storage
class PartitionSet {
public:

#ifdef OCTOTIGERII_WITH_HPX
	using Partition = hpx::id_type;
#else
	using Partition = std::shared_ptr<StoragePartition>;
#endif

	explicit PartitionSet(std::vector<Locality> owners)
	  : state_(std::make_shared<State>()) {
		if (owners.empty()) throw std::invalid_argument("Storage needs at least one partition");
		state_->owners = std::move(owners);
		for (auto const& owner : state_->owners) {
#ifdef OCTOTIGERII_WITH_HPX
			state_->partitions.push_back(hpx::new_<StoragePartition>(owner).get());
#else
			(void) owner;
			state_->partitions.push_back(std::make_shared<StoragePartition>());
#endif
		}
	}

	/// Return the locality assigned to each partition; repetitions are permitted.
	std::vector<Locality> const& owners() const {
		return state_->owners;
	}

	/// Return component handles shared by the fields in this partition set.
	std::vector<Partition> const& partitions() const {
		return state_->partitions;
	}

	/// Reserve a fresh field identity. Retired IDs are never reused.
	FieldId newField() const {
		std::lock_guard lock(state_->mutex);
		if (state_->next == std::numeric_limits<FieldId>::max()) throw std::overflow_error("Storage field identity exhausted");
		return state_->next++;
	}

private:

	/// Shared partition directory and monotonic field-ID allocator.
	class State {
	public:

		std::vector<Locality> owners;
		std::vector<Partition> partitions;
#ifdef OCTOTIGERII_WITH_HPX
		hpx::mutex mutex;
#else
		std::mutex mutex;
#endif
		FieldId next = 1;
	};


	std::shared_ptr<State> state_;
};


/// Copyable typed address directory; copying it does not copy numerical data.
/// Local reads retain a bank allocation, remote reads return received buffers.
/// The caller enforces immutable inputs and non-overlapping writers. Low-level
/// buffers are writable; this interface does not implement stage scheduling.
/// @ingroup storage
template <typename T>
class FieldHandle {
public:

	Layout layout;
	std::vector<Locality> owners;
	std::vector<PartitionSet::Partition> partitions;
	FieldId id = 0;
	unsigned banks = 0;
	/// Read-only sum, evaluated into temporary buffers on demand. No backing allocation.
	std::vector<FieldHandle<T>> sumSources;

	/// Report whether a range belongs to the current locality.
	bool local(Range const& range) const {
		if (!sumSources.empty()) return sumSources.front().local(range);
#ifdef OCTOTIGERII_WITH_HPX
		return owners.at(range.partition) == hpx::find_here();
#else
		(void) range;
		return true;
#endif
	}

	/// Return a direct retained local view or asynchronously fetch a contiguous remote range.
	Future<Buffer<T>> read(Range range, unsigned bank) const {
		if (!sumSources.empty()) {
			std::vector<Future<Buffer<T>>> pending;
			for (auto const& source : sumSources) pending.push_back(source.read(range, bank));
			auto sum = [range](auto reads) {
				Buffer<T> result(range.count);
				std::fill_n(result.data(), range.count, T{});
				std::exception_ptr error;
				for (auto& read : reads) {
					try {
						auto values = read.get();
						for (std::size_t i = 0; i < range.count; ++i) result.data()[i] += values.data()[i];
					} catch (...) { if (!error) error = std::current_exception(); }
				}
				if (error) std::rethrow_exception(error);
				return result;
			};
#ifdef OCTOTIGERII_WITH_HPX
			return hpx::async(std::move(sum), std::move(pending));
#else
			return Future<Buffer<T>>(sum(std::move(pending)));
#endif
		}
		validate(range, bank);
#ifdef OCTOTIGERII_WITH_HPX
		if (local(range)) return hpx::make_ready_future(readField<T>(partitions.at(range.partition), id, bank, range.offset, range.count));
		return hpx::async<ReadAction<T>>(owners.at(range.partition), partitions.at(range.partition), id, bank, range.offset, range.count);
#else
		return Future<Buffer<T>>(partitions.at(range.partition)->template read<T>(id, bank, range.offset, range.count));
#endif
	}

	/// Acquire a direct local output view or allocate a temporary remote result buffer.
	Buffer<T> output(Range range, unsigned bank) const {
		if (!sumSources.empty()) { sumSources.front().validate(range, bank); return Buffer<T>(range.count); }
		validate(range, bank);
		return local(range) ? read(range, bank).get() : Buffer<T>(range.count);
	}

	/// Finish writing the selected output range before returning to the stage scheduler.
	void commit(Range range, unsigned bank, Buffer<T> const& output) const {
		if (!sumSources.empty()) {
			sumSources.front().validate(range, bank);
			if (output.size() != range.count) throw std::invalid_argument("Derived field size mismatch");
			return; // Writes of a temporary reconstructed component are deliberately discarded.
		}
		validate(range, bank);
		if (output.size() != range.count) throw std::invalid_argument("Field commit size mismatch");
#ifdef OCTOTIGERII_WITH_HPX
		if (!local(range))
			hpx::async<WriteAction<T>>(owners.at(range.partition), partitions.at(range.partition), id, bank, range.offset, output).get();
		else
#endif
		{
			// Also accept caller-owned buffers, while leaving direct local
			// output views untouched. Both paths validate the field identity.
			auto destination = read(range, bank).get();
			if (destination.data() != output.data() && range.count) std::copy_n(output.data(), range.count, destination.data());
		}
	}

	/// Serialize this value with its compile-time quantity types preserved.
	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & layout & owners & partitions & id & banks & sumSources;
	}

private:

	/// Reject an out-of-bounds range, invalid identity, or nonexistent bank.
	void validate(Range const& range, unsigned bank) const {
		layout.validate(range);
		if (!id || bank >= banks) throw std::out_of_range("Invalid field identity or stage bank");
	}
};


/// Named scalar column distributed over a shared PartitionSet.
/// Capacity is independent per partition and bank count is configurable.
/// A replacement allocation receives a fresh field ID. Old fields survive until
/// retire() or destruction of the last partition handle; this is not an in-place
/// resize API. Application stage publication remains a Runtime responsibility.
/// @ingroup storage
template <typename T>
class Field {
public:

	Field(Layout const& layout, std::vector<Locality> const& owners, unsigned banks = 2, std::string name = {})
	  : Field(layout, PartitionSet(owners), banks, std::move(name)) {}

	Field(Layout const& layout, PartitionSet const& store, unsigned banks = 2, std::string name = {}) {
		if (store.owners().size() != layout.partitionCount() || banks == 0) throw std::invalid_argument("Invalid field placement or bank count");
		handle_.layout = layout.addressSpace();
		handle_.owners = store.owners();
		handle_.partitions = store.partitions();
		handle_.id = store.newField();
		handle_.banks = banks;
		try {
			for (std::size_t i = 0; i < handle_.partitions.size(); ++i) {
#ifdef OCTOTIGERII_WITH_HPX
				hpx::async<CreateAction<T>>(handle_.owners[i], handle_.partitions[i], handle_.id, layout.capacity(i), banks, name).get();
#else
				handle_.partitions[i]->template create<T>(handle_.id, layout.capacity(i), banks, name);
#endif
			}
		} catch (...) {
			auto error = std::current_exception();
			try {
				retire();
			} catch (...) {
			}
			std::rethrow_exception(error);
		}
	}

	/// Return lightweight field identities and placement without copying values.
	FieldHandle<T> const& handle() const {
		return handle_;
	}

	/// Remove the field from subsequent lookup after writers/transfers drain.
	/// Already acquired views retain their allocation; stale handle lookups fail.
	void retire() {
		std::exception_ptr error;
		for (std::size_t i = 0; i < handle_.partitions.size(); ++i) {
			try {
#ifdef OCTOTIGERII_WITH_HPX
				hpx::async<RetireAction>(handle_.owners[i], handle_.partitions[i], handle_.id).get();
#else
				handle_.partitions[i]->retire(handle_.id);
#endif
			} catch (...) {
				if (!error) error = std::current_exception();
			}
		}
		if (error) std::rethrow_exception(error);
	}

private:

	FieldHandle<T> handle_;
};


}	 // namespace octotigerII::storage
