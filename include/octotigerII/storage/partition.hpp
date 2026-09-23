/** @file
 * @brief Generic HPX storage component and typed range transfer actions.
 * @ingroup storage
 */
// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#pragma once

#include <algorithm>
#include <cstdint>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include "octotigerII/storage/buffer.hpp"

#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/include/actions.hpp>
#include <hpx/include/components.hpp>
#include <hpx/parcelset/coalescing_message_handler_registration.hpp>
#endif


namespace octotigerII::storage {

using FieldId = std::uint64_t;

/// HPX component owning multiple typed columns, independent of geometry.
/// Each bank has a stable allocation. A short metadata lock protects lookup;
/// numerical reads/writes obey the caller's immutable-input/single-writer policy.
/// Retained buffers keep their bank alive through transmission and retirement.
/// See @ref ref_biddiscombe2017 "Biddiscombe et al. (2017)" for HPX buffer transfers.
/// @ingroup storage
class StoragePartition
#ifdef OCTOTIGERII_WITH_HPX
  : public hpx::components::component_base<StoragePartition>
#endif
{

public:

	StoragePartition() = default;

	/// Allocate a new typed column before publishing its metadata.
	/// @param id Fresh identity, unique within the partition set.
	/// @param capacity Number of elements in each bank of this partition.
	/// @param banks Number of independently allocated versions, at least one.
	/// @param name Descriptive field name; it does not determine the C++ type.
	template <typename T>
	void create(FieldId id, std::size_t capacity, unsigned banks, std::string name) {
		if (banks == 0) throw std::invalid_argument("A field needs at least one bank");
		auto column = std::make_shared<Column<T>>(capacity, banks, std::move(name));
		std::lock_guard lock(mutex_);
		if (!columns_.emplace(id, std::move(column)).second) throw std::invalid_argument("Field identity already exists in partition");
	}

	/// Retire a field after its tasks drain. Existing views retain their bank.
	void retire(FieldId id) {
		std::lock_guard lock(mutex_);
		columns_.erase(id);
	}

	/// Retain a typed bank range after checking identity, type, bank, and bounds.
	/// The HPX reference buffer keeps the allocation alive during serialization.
	/// Callers must keep its contents immutable until readers and sends finish.
	template <typename T>
	Buffer<T> read(FieldId id, unsigned bank, std::size_t begin, std::size_t count) const {
		auto data = allocation<T>(id, bank, begin, count);
#ifdef OCTOTIGERII_WITH_HPX
		auto pointer = data->data();
		if (begin) pointer += begin;
		return Buffer<T>(pointer, count, Buffer<T>::reference, [data = std::move(data)](T*) {});
#else
		return Buffer<T>(std::move(data), begin, count);
#endif
	}

	/// Copy a received range into an exclusively owned output allocation.
	/// This final placement copy remains even with zero-copy receive serialization.
	template <typename T>
	void write(FieldId id, unsigned bank, std::size_t begin, Buffer<T> const& values) {
		auto data = allocation<T>(id, bank, begin, values.size());
		if (values.size()) std::copy_n(values.data(), values.size(), data->data() + begin);
	}

private:

	/// Type-erased metadata; concrete columns retain their actual quantity type.
	class ColumnBase {
	public:

		explicit ColumnBase(std::string name)
		  : name(std::move(name)) {}

		virtual ~ColumnBase() = default;

		std::string name;
	};


	template <typename T>
	class Column : public ColumnBase {
	public:

		Column(std::size_t capacity, unsigned count, std::string name)
		  : ColumnBase(std::move(name)) {
			for (unsigned bank = 0; bank < count; ++bank)
				banks.push_back(std::make_shared<std::vector<T>>(capacity));
		}

		std::vector<std::shared_ptr<std::vector<T>>> banks;
	};


	/// Resolve under a short metadata lock and retain the validated bank.
	template <typename T>
	std::shared_ptr<std::vector<T>> allocation(FieldId id, unsigned bank, std::size_t begin, std::size_t count) const {
		std::shared_ptr<Column<T>> column;
		{
			std::lock_guard lock(mutex_);
			column = std::dynamic_pointer_cast<Column<T>>(columns_.at(id));
		}
		if (!column) throw std::invalid_argument("Field quantity type mismatch");
		auto data = column->banks.at(bank);
		if (begin > data->size() || count > data->size() - begin) throw std::out_of_range("Field range exceeds partition capacity");
		return data;
	}

	mutable std::mutex mutex_;
	std::map<FieldId, std::shared_ptr<ColumnBase>> columns_;
};


#ifdef OCTOTIGERII_WITH_HPX
template <typename T>
void createField(hpx::id_type id, FieldId field, std::size_t capacity, unsigned banks, std::string name) {
	hpx::get_ptr<StoragePartition>(id).get()->create<T>(field, capacity, banks, std::move(name));
}


template <typename T>
Buffer<T> readField(hpx::id_type id, FieldId field, unsigned bank, std::size_t begin, std::size_t count) {
	return hpx::get_ptr<StoragePartition>(id).get()->read<T>(field, bank, begin, count);
}


template <typename T>
void writeField(hpx::id_type id, FieldId field, unsigned bank, std::size_t begin, Buffer<T> values) {
	hpx::get_ptr<StoragePartition>(id).get()->write<T>(field, bank, begin, values);
}


void retireField(hpx::id_type id, FieldId field);
using RetireAction = hpx::actions::make_action<decltype(&retireField), &retireField>::type;

template <typename T>
using CreateAction = typename hpx::actions::make_action<decltype(&createField<T>), &createField<T>>::type;

template <typename T>
using ReadAction = typename hpx::actions::make_action<decltype(&readField<T>), &readField<T>>::type;

template <typename T>
using WriteAction = typename hpx::actions::make_action<decltype(&writeField<T>), &writeField<T>>::type;
#endif

}	 // namespace octotigerII::storage


#ifdef OCTOTIGERII_WITH_HPX
HPX_REGISTER_ACTION_DECLARATION(octotigerII::storage::RetireAction, octoRetireField)
#define OCTO_STORAGE_TYPE(Type, Name)                                                                                                              \
	HPX_REGISTER_ACTION_DECLARATION(octotigerII::storage::CreateAction<Type>, Name##Create)                                                        \
	HPX_REGISTER_ACTION_DECLARATION(octotigerII::storage::ReadAction<Type>, Name##Read)                                                            \
	HPX_REGISTER_ACTION_DECLARATION(octotigerII::storage::WriteAction<Type>, Name##Write)                                                          \
	HPX_ACTION_USES_MESSAGE_COALESCING_NOTHROW_DECLARATION(octotigerII::storage::ReadAction<Type>, #Name "Read", std::size_t(-1), std::size_t(-1)) \
	HPX_ACTION_USES_MESSAGE_COALESCING_NOTHROW_DECLARATION(octotigerII::storage::WriteAction<Type>, #Name "Write", std::size_t(-1), std::size_t(-1))
#include "octotigerII/storage/types.def"
#undef OCTO_STORAGE_TYPE
#endif
