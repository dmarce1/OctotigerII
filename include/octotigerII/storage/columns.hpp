/** @file
 * @brief Structure-of-arrays storage for heterogeneous physical states.
 * @ingroup storage
 */
// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#pragma once

#include <exception>
#include "octotigerII/storage/field.hpp"
#include "octotigerII/units/state.hpp"

namespace octotigerII::storage {

template <typename State>
class Columns;
template <typename State>
class ColumnHandle;
template <typename State>
class ColumnFields;
template <typename State>
class PendingColumns;

template <typename... Q>
class Columns<units::State<Q...>> {
public:
    using State = units::State<Q...>;
    std::tuple<Buffer<Q>...> buffers;

    State at(std::size_t i) const {
        return std::apply([i](auto const&... buffer) { return State(buffer.data()[i]...); }, buffers);
    }

    void put(std::size_t i, State const& state) {
        state.forEach([&](auto field, auto const& quantity) { std::get<field>(buffers).data()[i] = quantity; });
    }
};

template <typename... Q>
class PendingColumns<units::State<Q...>> {
public:
    std::tuple<Future<Buffer<Q>>...> pending;

    /// Return the component selected by its compile-time tuple index.
    Columns<units::State<Q...>> get() {
        Columns<units::State<Q...>> result;
        std::exception_ptr error;
        units::State<Q...>{}.forEach([&](auto field, auto const&) {
            try {
                std::get<field>(result.buffers) = std::get<field>(pending).get();
            } catch (...) {
                if (!error) error = std::current_exception();
            }
        });
        if (error) std::rethrow_exception(error);
        return result;
    }
};

template <typename... Q>
class ColumnHandle<units::State<Q...>> {
public:
    using State = units::State<Q...>;
    std::tuple<FieldHandle<Q>...> fields;

    PendingColumns<State> read(Range range, unsigned bank) const {
        return {std::apply([&](auto const&... field) { return std::tuple{field.read(range, bank)...}; }, fields)};
    }

    /// Acquire a direct local output view or allocate a temporary remote result buffer.
    Columns<State> output(Range range, unsigned bank) const {
        return {std::apply([&](auto const&... field) { return std::tuple{field.output(range, bank)...}; }, fields)};
    }

    /// Finish writing the selected output range before returning to the stage scheduler.
    void commit(Range range, unsigned bank, Columns<State> const& columns) const {
        State{}.forEach([&](auto i, auto const&) { std::get<i>(fields).commit(range, bank, std::get<i>(columns.buffers)); });
    }

    /// Serialize this value with its compile-time quantity types preserved.
    template <typename Archive>
    void serialize(Archive& archive, unsigned) {
        std::apply([&](auto&... field) { ((archive & field), ...); }, fields);
    }
};

template <typename... Q>
class ColumnFields<units::State<Q...>> {
public:
    ColumnFields(Layout const& layout, PartitionSet const& store, std::string const& name)
      : fields_(makeFields(layout, store, name, std::index_sequence_for<Q...>{})) {}

    /// Return lightweight field identities and placement without copying values.
    ColumnHandle<units::State<Q...>> handle() const {
        return {std::apply([](auto const&... field) { return std::tuple{field->handle()...}; }, fields_)};
    }

private:
    template <std::size_t... I>
    static auto makeFields(Layout const& layout, PartitionSet const& store, std::string const& name, std::index_sequence<I...>) {
        return std::tuple{std::make_unique<Field<Q>>(layout, store, 2, name + "." + std::to_string(I))...};
    }

    std::tuple<std::unique_ptr<Field<Q>>...> fields_;
};

}    // namespace octotigerII::storage
