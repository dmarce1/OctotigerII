/** @file
 * @brief Reference-backed typed buffers and serial/HPX future abstractions.
 * @ingroup storage
 */
// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#pragma once

#include <memory>
#include <type_traits>
#include <utility>
#include <vector>
#include "octotigerII/units/cgs.hpp"

#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/include/async.hpp>
#include <hpx/include/serialization.hpp>
#include <hpx/serialization/serialize_buffer.hpp>

// Explicit opt-in for the concrete CGS scalar representations. These Boost
// quantities have user-defined copies, so HPX cannot infer the optimization.
#define OCTO_STORAGE_TYPE(Type, Name)                                               \
	static_assert(sizeof(Type) == sizeof(Real) && std::is_standard_layout_v<Type>); \
	HPX_IS_BITWISE_SERIALIZABLE(Type)
#include "octotigerII/storage/types.def"
#undef OCTO_STORAGE_TYPE
#endif


namespace octotigerII::storage {

#ifdef OCTOTIGERII_WITH_HPX
using Locality = hpx::id_type;

template <typename T>
using Future = hpx::future<T>;

template <typename T>
using Buffer = hpx::serialization::serialize_buffer<T>;
#else
using Locality = std::size_t;


/// Already-completed value used by the serial backend of the field API.
/// @ingroup storage
template <typename T>
class Future {
public:

	explicit Future(T value)
	  : value_(std::move(value)) {}

	/// Move the already-completed serial result to the caller.
	T get() {
		return std::move(value_);
	}

private:

	T value_;
};


/// Shared ownership of a typed contiguous range in the serial backend.
/// The HPX backend uses serialize_buffer with a retained allocation instead.
/// @ingroup storage
template <typename T>
class Buffer {
public:

	Buffer() = default;

	explicit Buffer(std::size_t size)
	  : Buffer(std::make_shared<std::vector<T>>(size), 0, size) {}

	Buffer(std::shared_ptr<std::vector<T>> owner, std::size_t begin, std::size_t size)
	  : owner_(std::move(owner))
	  , data_(owner_->data())
	  , size_(size) {
		if (begin) data_ += begin;
	}

	/// Access the contiguous range. The owning buffer must remain alive.
	T* data() {
		return data_;
	}

	/// Access the contiguous range. The owning buffer must remain alive.
	T const* data() const {
		return data_;
	}

	/// Return the number of represented elements or blocks.
	std::size_t size() const {
		return size_;
	}

private:

	std::shared_ptr<std::vector<T>> owner_;
	T* data_ = nullptr;
	std::size_t size_ = 0;
};
#endif

}	 // namespace octotigerII::storage
