#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <silo.h>
#include <stdexcept>
#include "octotigerII/output.hpp"
#include "octotigerII/verification/analytic.hpp"

using namespace octotigerII;


namespace {
void require(bool ok, char const* message) {
	if (!ok) throw std::runtime_error(message);
}

void check(std::string const& problem, int dimensions, std::filesystem::path const& root) {
	using std::abs;
	using std::isfinite;

	auto c = parseConfig({"--mesh.level=0", "--mesh.cells=4", "--runtime.stopTime=0", "--output.directory=" + root.string()});
	c.verification.gravityReference = "continuum";
	auto patch = initialSnapshot(c, {0, {}});
	if (patch.gravityEnabled) {
		gravity::State state{};
		state.potential() = units::VelocitySquared::from_value(-12);
		for (int axis = 0; axis < ndim; ++axis)
			state.acceleration(axis) = units::Acceleration::from_value(axis + 1);
		patch.gravity.values().assign(patch.layout.interiorCellCount(), state);
	}
	std::vector<Snapshot> patches{patch};
	// Write twice to the same name: exercise the DB_CLOBBER requirement.
	for (int pass = 0; pass < 2; ++pass) {
		Output output(c);
		output(patches, 0, diagnose(patches, c));
	}
	std::string const filename = (std::filesystem::path(c.output.directory) / "frame_000000.silo").string();
	std::unique_ptr<DBfile, decltype(&DBClose)> file(DBOpen(filename.c_str(), DB_HDF5, DB_READ), &DBClose);
	require(bool(file), "Silo open");
	std::unique_ptr<DBmultimesh, decltype(&DBFreeMultimesh)> multi(DBGetMultimesh(file.get(), "mesh"), &DBFreeMultimesh);
	require(multi && multi->nblocks == 1, "Silo multimesh");
	std::unique_ptr<DBquadmesh, decltype(&DBFreeQuadmesh)> mesh(DBGetQuadmesh(file.get(), "block0/mesh"), &DBFreeQuadmesh);
	require(mesh && mesh->ndims == dimensions, "Silo dimensional geometry");
	for (int axis = 0; axis < dimensions; ++axis) {
		require(mesh->dims[axis] == 5, "Silo coordinate count");
		require(mesh->units[axis] && std::string(mesh->units[axis]) == "cm", "Silo coordinate units");
		auto const* values = static_cast<double const*>(mesh->coords[axis]);
		require(values[0] == units::value(patch.lower[axis]) && values[4] == units::value(patch.lower[axis] + 4.0 * patch.cellWidth), "Silo coordinates");
	}
	auto field = [&](char const* name, char const* units, auto expected) {
		std::unique_ptr<DBquadvar, decltype(&DBFreeQuadvar)> var(DBGetQuadvar(file.get(), ("block0/" + std::string(name)).c_str()), &DBFreeQuadvar);
		require(
			var && var->datatype == DB_DOUBLE && var->centering == DB_ZONECENT && var->nels == int(patch.layout.interiorCellCount()), "Silo field metadata");
		require(var->units && std::string(var->units) == units, "Silo CGS field units");
		require(var->dtime == 0 && var->cycle == 0, "Silo time/cycle");
		auto const* values = static_cast<double const*>(var->vals[0]);
		std::size_t j = 0;
		patch.layout.forEachInterior([&](mesh::Coordinates const&, std::size_t i) {
			Real const value = units::value(expected(i));
			require(isfinite(values[j]) && abs(values[j] - value) <= 1e-14 * std::max(Real(1), abs(value)), "Silo field roundtrip");
			++j;
		});
	};
	if (patch.hydroEnabled) {
		field("density", "g/cm^3", [&](std::size_t i) { return patch.hydro.values()[i].density(); });
		field("gasEnergy", "erg/cm^3", [&](std::size_t i) { return patch.hydro.values()[i].totalEnergy(); });
	}
	if (patch.radiationEnabled) {
		field("radiationEnergy", "erg/cm^3", [&](std::size_t i) { return patch.radiation.values()[i].energy(); });
		field("radiationFluxX", "erg/(cm^2 s)", [&](std::size_t i) { return patch.radiation.values()[i].radiativeFlux(0); });
	}
	if (patch.gravityEnabled) {
		field("potential", "cm^2/s^2", [](std::size_t) { return units::VelocitySquared::from_value(-12); });
		field("accelerationZ", "cm/s^2", [](std::size_t) { return units::Acceleration::from_value(3); });
	}
	for (auto const& f : verification::sample(patch, c, verification::reference(c, patch.time))) {
		std::unique_ptr<DBquadvar, decltype(&DBFreeQuadvar)> exact(DBGetQuadvar(file.get(), ("block0/" + f.name + "Exact").c_str()), &DBFreeQuadvar);
		std::unique_ptr<DBquadvar, decltype(&DBFreeQuadvar)> error(DBGetQuadvar(file.get(), ("block0/" + f.name + "Error").c_str()), &DBFreeQuadvar);
		require(exact && error && exact->nels == int(f.exact.size()) && error->nels == exact->nels, "Analytic Silo arrays");
		require(exact->units && f.units == exact->units && error->units && f.units == error->units, "Analytic Silo CGS units");
		for (std::size_t i = 0; i < f.exact.size(); ++i) {
			require(static_cast<double const*>(exact->vals[0])[i] == f.exact[i], "Analytic Silo value");
			require(static_cast<double const*>(error->vals[0])[i] == f.numerical[i] - f.exact[i], "Signed Silo error");
		}
		std::unique_ptr<DBmultivar, decltype(&DBFreeMultivar)> multiExact(DBGetMultivar(file.get(), (f.name + "Exact").c_str()), &DBFreeMultivar);
		require(multiExact && multiExact->nvars == 1, "Analytic Silo multivar");
	}
	auto* toc = DBGetToc(file.get());
	require(toc != nullptr, "Silo table of contents");
	for (int axis = ndim; axis < 3; ++axis) {
		for (int i = 0; i < toc->nmultivar; ++i) {
			std::string const name = toc->multivar_names[i];
			require(
				name != std::string("momentum") + "XYZ"[axis] && name != std::string("radiationFlux") + "XYZ"[axis], "Silo contains inactive vector component");
		}
	}
	std::cout << problem << ' ' << dimensions << "D Silo roundtrip passed\n";
}

}	 // namespace


int main(int argc, char** argv) {
	try {
		if (argc != 2) throw std::invalid_argument("Expected output directory");
		std::filesystem::path const root(argv[1]);
		check(build::problem, ndim, root);
		return 0;
	} catch (std::exception const& error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
}
