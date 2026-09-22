#include "octotigerII/output.hpp"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <silo.h>
#include <stdexcept>

using namespace octotigerII;
namespace {
void require(bool ok, char const* message) {
	if (!ok)
		throw std::runtime_error(message);
}
void check(std::string const& problem, int dimensions, std::filesystem::path const& root) {
	auto c = parseConfig(
		{"--problem.name=" + problem, "--mesh.ndim=" + std::to_string(dimensions), "--mesh.level=0",
		 "--mesh.cells=4", "--output.format=silo", "--runtime.stop_time=0",
		 "--output.directory=" + (root / (problem + std::to_string(dimensions))).string()});
	Subgrid block(c, {0, {0, 0, 0}, dimensions});
	auto patch = block.snapshot();
	if (patch.gravityEnabled) {
		block.setGravity(std::vector<gravity::State>(patch.layout.interiorCellCount(),
													 gravity::State{-12, 1, 2, 3}));
		patch = block.snapshot();
	}
	std::vector<Snapshot> patches{patch};
	// Write twice to the same name: exercise the DB_CLOBBER requirement.
	for (int pass = 0; pass < 2; ++pass) {
		Output output(c);
		output(patches, 0, diagnose(patches, c));
	}
	std::string const filename =
		(std::filesystem::path(c.outputDirectory) / "frame_000000.silo").string();
	std::unique_ptr<DBfile, decltype(&DBClose)> file(DBOpen(filename.c_str(), DB_HDF5, DB_READ),
													 &DBClose);
	require(bool(file), "Silo open");
	std::unique_ptr<DBmultimesh, decltype(&DBFreeMultimesh)> multi(
		DBGetMultimesh(file.get(), "mesh"), &DBFreeMultimesh);
	require(multi && multi->nblocks == 1, "Silo multimesh");
	std::unique_ptr<DBquadmesh, decltype(&DBFreeQuadmesh)> mesh(
		DBGetQuadmesh(file.get(), "block0/mesh"), &DBFreeQuadmesh);
	require(mesh && mesh->ndims == dimensions, "Silo dimensional geometry");
	for (int axis = 0; axis < dimensions; ++axis) {
		require(mesh->dims[axis] == 5, "Silo coordinate count");
		auto const* values = static_cast<double const*>(mesh->coords[axis]);
		require(values[0] == patch.lower[axis] &&
					values[4] == patch.lower[axis] + 4 * patch.cellWidth,
				"Silo coordinates");
	}
	auto field = [&](char const* name, auto expected) {
		std::unique_ptr<DBquadvar, decltype(&DBFreeQuadvar)> var(
			DBGetQuadvar(file.get(), ("block0/" + std::string(name)).c_str()), &DBFreeQuadvar);
		require(var && var->datatype == DB_DOUBLE && var->centering == DB_ZONECENT &&
					var->nels == int(patch.layout.interiorCellCount()),
				"Silo field metadata");
		require(var->dtime == 0 && var->cycle == 0, "Silo time/cycle");
		auto const* values = static_cast<double const*>(var->vals[0]);
		std::size_t j = 0;
		patch.layout.forEachInterior([&](mesh::Coordinates const&, std::size_t i) {
			Real const value = expected(i);
			require(std::isfinite(values[j]) &&
						std::abs(values[j] - value) <= 1e-14 * std::max(Real(1), std::abs(value)),
					"Silo field roundtrip");
			++j;
		});
	};
	if (patch.hydroEnabled) {
		field("density", [&](std::size_t i) { return patch.hydro.values()[i].density(); });
		field("gasEnergy", [&](std::size_t i) { return patch.hydro.values()[i].totalEnergy(); });
	}
	if (patch.radiationEnabled) {
		field("radiationEnergy", [&](std::size_t i) { return patch.radiation.values()[i][0]; });
		field("radiationFluxX",
			  [&](std::size_t i) { return physicalLightSpeed * patch.radiation.values()[i][1]; });
	}
	if (patch.gravityEnabled) {
		field("potential", [](std::size_t) { return -12.0; });
		field("accelerationZ", [](std::size_t) { return 3.0; });
	}
	std::cout << problem << ' ' << dimensions << "D Silo roundtrip passed\n";
}
} // namespace
int main(int argc, char** argv) {
	try {
		if (argc != 2)
			throw std::invalid_argument("Expected output directory");
		std::filesystem::path const root(argv[1]);
		check("sod", 1, root);
		check("kelvin-helmholtz", 2, root);
		check("collapse", 3, root);
		for (int d = 1; d <= 3; ++d)
			check("streaming", d, root);
		check("gravity-gaussian", 3, root);
		return 0;
	} catch (std::exception const& error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
}
