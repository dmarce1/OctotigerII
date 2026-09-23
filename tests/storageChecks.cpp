#include <iostream>
#include <limits>
#include <stdexcept>
#include "octotigerII/runtime.hpp"
#include "octotigerII/storage/field.hpp"
#include "runtimeMain.hpp"
#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/runtime_distributed/find_all_localities.hpp>
#endif

using namespace octotigerII;


namespace {

void require(bool ok, char const* message) {
	if (!ok) throw std::runtime_error(message);
}

void fields() {
	std::vector<storage::Locality> owners;
#ifdef OCTOTIGERII_WITH_HPX
	owners = hpx::find_all_localities();
#else
	owners = {0};
#endif
	owners.push_back(owners.front());	 // Multiple partitions may share a locality.
	storage::PartitionSet store(owners);
	// Deliberately not mesh-shaped: these can equally be particle attributes.
	storage::Layout layout({3, 8192, 7, 19, 5}, owners.size());
	storage::Field<units::Density> field(layout, store, 3, "particles.density");
	auto const& handle = field.handle();
	for (auto range : layout.ranges()) {
		auto output = handle.output(range, 0);
		for (std::size_t i = 0; i < range.count; ++i)
			output.data()[i] = units::Density::from_value(100 * range.partition + range.offset + i);
		handle.commit(range, 0, output);
		auto input = handle.read(range, 0).get();
		require(!handle.local(range) || input.data() == output.data(), "Local checkout copied field data");
		for (std::size_t i = 0; i < range.count; ++i)
			require(input.data()[i] == units::Density::from_value(100 * range.partition + range.offset + i), "Field read/write roundtrip");
		auto next = handle.output(range, 1);
		for (std::size_t i = 0; i < range.count; ++i)
			next.data()[i] = 2.0 * input.data()[i];
		handle.commit(range, 1, next);
		auto current = handle.read(range, 0).get();
		auto future = handle.read(range, 1).get();
		for (std::size_t i = 0; i < range.count; ++i)
			require(current.data()[i] == input.data()[i] && future.data()[i] == 2.0 * input.data()[i], "Stage banks alias");
		auto slice = handle.read(range.slice(1, range.count - 1), 0).get();
		require(slice.data()[0] == input.data()[1], "Subrange read");
	}
	// Two fields with the same quantity type must remain separate allocations.
	storage::Field<units::Density> other(layout, store, 1, "other.density");
	require(other.handle().id != handle.id && other.handle().partitions == handle.partitions, "Fields did not share a partition with distinct identities");
	for (auto range : layout.ranges()) {
		auto independent = other.handle().read(range, 0).get();
		for (std::size_t i = 0; i < range.count; ++i)
			require(independent.data()[i] == units::Density{}, "Distinct fields alias");
		auto third = handle.output(range, 2);
		for (std::size_t i = 0; i < range.count; ++i)
			third.data()[i] = units::Density::from_value(42);
		handle.commit(range, 2, third);
		require(handle.read(range, 2).get().data()[0] == units::Density::from_value(42), "Configurable bank count failed");
	}
	// A second indexing layout uses the same partitions without mesh coupling.
	storage::Layout fluxLayout({37, 8201, 11}, owners.size());
	storage::Field<units::EnergyFlux> flux(fluxLayout, store, 1, "radiation.flux.x");
	for (auto range : fluxLayout.ranges()) {
		auto output = flux.handle().output(range, 0);
		for (std::size_t i = 0; i < range.count; ++i)
			output.data()[i] = units::EnergyFlux::from_value(2.99792458e10 * (i + 1));
		flux.handle().commit(range, 0, output);
		auto input = flux.handle().read(range, 0).get();
		for (std::size_t i = 0; i < range.count; ++i)
			require(input.data()[i] == units::EnergyFlux::from_value(2.99792458e10 * (i + 1)), "Physical flux transfer changed units or values");
	}
	bool rejected = false;
	try {
		handle.read({0, layout.capacity(), 1}, 0).get();
	} catch (std::out_of_range const&) {
		rejected = true;
	}
	require(rejected, "Invalid field address accepted");
	rejected = false;
	try {
		other.handle().read(layout.ranges()[0], 1).get();
	} catch (std::out_of_range const&) {
		rejected = true;
	}
	require(rejected, "Invalid bank accepted");
	// A forged handle must not reinterpret a different Boost quantity type.
	storage::FieldHandle<units::EnergyDensity> wrong;
	wrong.layout = handle.layout;
	wrong.owners = handle.owners;
	wrong.partitions = handle.partitions;
	wrong.id = handle.id;
	wrong.banks = handle.banks;
	rejected = false;
	try {
		wrong.read(layout.ranges().back(), 0).get();
	} catch (std::exception const&) {
		rejected = true;
	}
	require(rejected, "Mismatched quantity type accepted");
#ifdef OCTOTIGERII_WITH_HPX
	// Exercise the exact typed buffer representation used by field reads.
	// Large eligible quantities must become pointer chunks, not packed values.
	auto checkChunks = [&](auto const& buffer) {
		std::vector<char> bytes;
		std::vector<hpx::serialization::serialization_chunk> chunks;
		{
			hpx::serialization::output_archive archive(bytes, 0, &chunks);
			archive & buffer;
		}
		bool pointerChunk = false;
		for (auto const& chunk : chunks)
			pointerChunk |= chunk.type_ == hpx::serialization::chunk_type::chunk_type_pointer;
		require(pointerChunk, "Typed field buffer lost zero-copy serialization");
	};
	checkChunks(handle.read(layout.ranges()[1], 0).get());
	checkChunks(flux.handle().read(fluxLayout.ranges()[1], 0).get());
	std::cout << "Density and physical radiation-flux zero-copy pointer chunks verified\n";
#endif
	auto retained = flux.handle().read(fluxLayout.ranges()[0], 0).get();
	flux.retire();
	require(retained.data()[0] == units::EnergyFlux::from_value(2.99792458e10), "Retirement invalidated a retained buffer");
	rejected = false;
	try {
		flux.handle().read(fluxLayout.ranges().back(), 0).get();
	} catch (std::exception const&) {
		rejected = true;
	}
	require(rejected, "Retired field handle accepted");
	storage::Field<units::EnergyFlux> replacement(fluxLayout, store, 1);
	require(replacement.handle().id != flux.handle().id, "Retired field identity reused");
	std::cout << "Shared partitions, independent layouts/fields, physical flux, views, banks, type validation and retirement passed\n";
}

void same(std::vector<Snapshot> const& a, std::vector<Snapshot> const& b) {
	require(a.size() == b.size(), "Snapshot count changed");
	for (std::size_t block = 0; block < a.size(); ++block) {
		require(a[block].time == b[block].time, "Failed stage changed time");
		require(a[block].layout.ghostWidth() == 0, "Persistent snapshot contains ghosts");
		for (std::size_t i = 0; i < a[block].radiation.values().size(); ++i)
			a[block].radiation.values()[i].forEach(
				[&](auto field, auto q) { require(q == b[block].radiation.values()[i].template get<field>(), "Published radiation input changed"); });
		for (std::size_t i = 0; i < a[block].hydro.values().size(); ++i)
			a[block].hydro.values()[i].forEach(
				[&](auto field, auto q) { require(q == b[block].hydro.values()[i].template get<field>(), "Published input changed"); });
	}
}

void stages() {
	auto config = parseConfig({"--mesh.cells=4", "--mesh.level=1", "--output.enabled=off"});
	Runtime runtime(config), reference(config);
	auto initial = runtime.snapshots();
	bool rejected = false;
	try {
		if (std::string(build::problem) == "sod")
			runtime.advance(units::Time::from_value(10));
		else
			runtime.advance(units::Time::from_value(-1));
	} catch (std::exception const&) {
		rejected = true;
	}
	require(rejected && runtime.generation() == 0, "Invalid numerical stage was published");
	same(initial, runtime.snapshots());
	auto dt = runtime.stableTimestep();
	runtime.advance(dt);
	reference.advance(dt);
	same(runtime.snapshots(), reference.snapshots());
	require(runtime.generation() == 1, "Completed stage did not publish exactly once");
	auto const stats = runtime.statistics();
	require(stats.localTasks + stats.stolenTasks == 2 * runtime.size(), "Successful tasks missing or duplicated");
	std::cout << "Failed stage rollback, retry, publication and task accounting passed (" << stats.localTasks << " local, " << stats.stolenTasks
			  << " stolen tasks)\n";
}

}	 // namespace


int testMain(int, char**) {
	try {
		fields();
		if constexpr (build::hydro || build::radiation) stages();
		return 0;
	} catch (std::exception const& error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
}

int main(int argc, char** argv) {
	return runtimeMain(argc, argv, testMain);
}
