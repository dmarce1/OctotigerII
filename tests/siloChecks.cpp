#include "testSupport.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <iostream>
#include <memory>
#include <set>
#include <silo.h>
#include <stdexcept>
#include "octotigerII/output.hpp"
#include "octotigerII/verification/analytic.hpp"

using namespace octotigerII;

namespace {

void checkVectorExpression(DBfile* file, std::string const& stem, std::string const& suffix, int blocks) {
	std::unique_ptr<DBdefvars, decltype(&DBFreeDefvars)> expressions(DBGetDefvars(file, "expressions"), &DBFreeDefvars);
	ASSERT_TRUE(expressions);
	auto const name = stem + suffix;
	int index = -1;
	for (int i = 0; i < expressions->ndefs; ++i) {
		if (name == expressions->names[i]) {
			ASSERT_EQ(index, -1) << "Duplicate expression: " << name;
			index = i;
		}
	}
	ASSERT_GE(index, 0) << "Missing vector expression: " << name;
	EXPECT_EQ(expressions->types[index], DB_VARTYPE_VECTOR);
	EXPECT_TRUE(!expressions->guihides || expressions->guihides[index] == 0);
	std::array<std::string, 3> components;
	components.fill("zonal_constant(<mesh>,0)");
	for (int axis = 0; axis < ndim; ++axis) {
		components[axis] = stem + "XYZ"[axis] + suffix;
		std::unique_ptr<DBmultivar, decltype(&DBFreeMultivar)> multi(DBGetMultivar(file, components[axis].c_str()), &DBFreeMultivar);
		ASSERT_TRUE(multi) << "Expression component: " << components[axis];
		ASSERT_EQ(multi->nvars, blocks);
		EXPECT_EQ(multi->tensor_rank, DB_VARTYPE_SCALAR);
		EXPECT_STREQ(multi->mmesh_name, "mesh");
		for (int b = 0; b < blocks; ++b) {
			EXPECT_EQ(multi->vartypes[b], DB_QUADVAR);
			EXPECT_EQ(std::string(multi->varnames[b]), "block" + std::to_string(b) + "/" + components[axis]);
			std::unique_ptr<DBquadvar, decltype(&DBFreeQuadvar)> var(DBGetQuadvar(file, multi->varnames[b]), &DBFreeQuadvar);
			ASSERT_TRUE(var);
			EXPECT_EQ(var->nvals, 1);
			EXPECT_EQ(var->centering, DB_ZONECENT);
		}
	}
	EXPECT_EQ(std::string(expressions->defns[index]), "{" + components[0] + "," + components[1] + "," + components[2] + "}");
	// The vector is only an expression: no duplicate vector arrays or objects.
	EXPECT_EQ(DBInqVarExists(file, name.c_str()), 0);
	for (int b = 0; b < blocks; ++b) {
		EXPECT_EQ(DBInqVarExists(file, ("block" + std::to_string(b) + "/" + name).c_str()), 0);
	}
}

void check(std::string const& problem, int dimensions, std::filesystem::path const& root) {
	using std::abs;
	using std::isfinite;

	auto c = octotigerII::test::parseConfig({"--mesh.level=0", "--mesh.cells=4", "--runtime.stopTime=0", "--output.directory=" + root.string()});
	c.verification.gravityReference = "continuum";
	c.hydro.meanMolecularWeight = 0.6;
	if (c.hydroEnabled() && build::massFractions) {
		c.massFractions.enabled = true;
		c.massFractions.species = composition::parseSpecies("fuel:0.7:He;ash:0.3:oxygen;dye:0.2:A=0,Z=0");
	}
	auto patch = initialSnapshot(c, {0, {}});
	// Give every component a different value, including a nonzero analytic
	// error, so component permutation or sign errors cannot pass unnoticed.
	if constexpr (octotigerII::test::hydro) {
		hydro::HydroSystem gas(c.hydro);
		for (auto& state : patch.hydro.values()) {
			auto primitive = gas.reconstructionVariables(state);
			for (int axis = 0; axis < ndim; ++axis) {
				primitive.velocity(axis) = units::Velocity::from_value(0.125 * (axis + 1));
			}
			state = gas.conservedState(primitive);
		}
	}
	if constexpr (octotigerII::test::radiation) {
		for (auto& state : patch.radiation.values()) {
			for (int axis = 0; axis < ndim; ++axis) {
				state.radiativeFlux(axis) = constants::c * state.energy() * (Real(axis + 1) / (10 * ndim));
			}
		}
	}
	if (patch.gravityEnabled) {
		gravity::State state{};
		state.potential() = units::VelocitySquared::from_value(-12);
		for (int axis = 0; axis < ndim; ++axis) {
			state.acceleration(axis) = units::Acceleration::from_value(axis + 1);
		}
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
	ASSERT_TRUE(bool(file)) << "Silo open";
	EXPECT_EQ(DBInqVarExists(file.get(), "MetadataIsTimeVarying"), 0);
	std::unique_ptr<DBmultimesh, decltype(&DBFreeMultimesh)> multi(DBGetMultimesh(file.get(), "mesh"), &DBFreeMultimesh);
	ASSERT_TRUE(multi && multi->nblocks == 1) << "Silo multimesh";
	std::unique_ptr<DBquadmesh, decltype(&DBFreeQuadmesh)> mesh(DBGetQuadmesh(file.get(), "block0/mesh"), &DBFreeQuadmesh);
	ASSERT_TRUE(mesh && mesh->ndims == dimensions) << "Silo dimensional geometry";
	for (int axis = 0; axis < dimensions; ++axis) {
		EXPECT_TRUE(mesh->dims[axis] == 5) << "Silo coordinate count";
		EXPECT_TRUE(mesh->units[axis] && std::string(mesh->units[axis]) == "cm") << "Silo coordinate units";
		auto const* values = static_cast<double const*>(mesh->coords[axis]);
		EXPECT_TRUE(values[0] == units::value(patch.lower[axis]) && values[4] == units::value(patch.lower[axis] + 4.0 * patch.cellWidth)) << "Silo coordinates";
	}
	auto field = [&](char const* name, char const* units, auto expected) {
		std::unique_ptr<DBquadvar, decltype(&DBFreeQuadvar)> var(DBGetQuadvar(file.get(), ("block0/" + std::string(name)).c_str()), &DBFreeQuadvar);
		ASSERT_TRUE(var && var->datatype == DB_DOUBLE && var->centering == DB_ZONECENT && var->nels == int(patch.layout.interiorCellCount()) && var->nvals == 1)
			<< "Silo field metadata";
		EXPECT_TRUE(var->units && std::string(var->units) == units) << "Silo CGS field units";
		EXPECT_TRUE(var->dtime == 0 && var->cycle == 0) << "Silo time/cycle";
		auto const* values = static_cast<double const*>(var->vals[0]);
		std::size_t j = 0;
		patch.layout.forEachInterior([&](mesh::Coordinates const&, std::size_t i) {
			Real const value = units::value(expected(i));
			EXPECT_TRUE(isfinite(values[j]) && abs(values[j] - value) <= 1e-14 * std::max(Real(1), abs(value))) << "Silo field roundtrip";
			++j;
		});
	};
	for (std::size_t s = 0; s < patch.species.size(); ++s) {
		auto const& species = c.massFractions.species[s];
		field(((species.tracer() ? "tracerDensity_" : "partialDensity_") + species.name).c_str(), "g/cm^3",
			[&](std::size_t i) { return patch.species[s].values()[i]; });
		std::string name = (species.tracer() ? "tracer_" : "massFraction_") + species.name;
		std::unique_ptr<DBquadvar, decltype(&DBFreeQuadvar)> fraction(DBGetQuadvar(file.get(), ("block0/" + name).c_str()), &DBFreeQuadvar);
		ASSERT_TRUE(fraction);
		auto values = static_cast<double const*>(fraction->vals[0]);
		for (int i = 0; i < fraction->nels; ++i) EXPECT_NEAR(values[i], species.initialFraction, 1e-14);
	}
	std::set<std::string> expectedExpressions;
	auto vector = [&](char const* name, char const* units, auto expected) {
		checkVectorExpression(file.get(), name, "", 1);
		expectedExpressions.insert(name);
		for (int axis = 0; axis < ndim; ++axis) {
			field((std::string(name) + "XYZ"[axis]).c_str(), units, [&](std::size_t i) { return expected(i, axis); });
		}
	};
	if (patch.hydroEnabled) {
		field("density", "g/cm^3", [&](std::size_t i) { return patch.hydro.values()[i].density(); });
		field("gasEnergy", "erg/cm^3", [&](std::size_t i) { return patch.hydro.values()[i].totalEnergy(); });
#if OCTOTIGERII_HYDRO
		hydro::HydroSystem gas(c.hydro);
		field("internalEnergy", "erg/cm^3", [&](std::size_t i) { return gas.internalEnergy(patch.hydro.values()[i]); });
		field("temperature", "K", [&](std::size_t i) { return gas.temperature(patch.hydro.values()[i]); });
		field("dualEnergy", "g/cm^3", [&](std::size_t i) { return patch.hydro.values()[i].auxiliary(); });
#endif
		vector("momentum", "g/(cm^2 s)", [&](std::size_t i, int axis) { return patch.hydro.values()[i].momentum(axis); });
		vector("velocity", "cm/s", [&](std::size_t i, int axis) { return patch.hydro.values()[i].momentum(axis) / patch.hydro.values()[i].density(); });
	}
	if (patch.radiationEnabled) {
		field("radiationEnergy", "erg/cm^3", [&](std::size_t i) { return patch.radiation.values()[i].energy(); });
		vector("radiationFlux", "erg/(cm^2 s)", [&](std::size_t i, int axis) { return patch.radiation.values()[i].radiativeFlux(axis); });
	}
	if (patch.gravityEnabled) {
		field("potential", "cm^2/s^2", [](std::size_t) { return units::VelocitySquared::from_value(-12); });
		vector("acceleration", "cm/s^2", [](std::size_t, int axis) { return units::Acceleration::from_value(axis + 1); });
	}
	for (auto const& f : verification::sample(patch, c, verification::reference(c, patch.time))) {
		for (auto stem : {"momentum", "velocity", "radiationFlux", "acceleration"}) {
			if (f.name == std::string(stem) + "X") {
				for (auto suffix : {"Exact", "Error"}) {
					checkVectorExpression(file.get(), stem, suffix, 1);
					expectedExpressions.insert(std::string(stem) + suffix);
				}
			}
		}
		std::unique_ptr<DBquadvar, decltype(&DBFreeQuadvar)> exact(DBGetQuadvar(file.get(), ("block0/" + f.name + "Exact").c_str()), &DBFreeQuadvar);
		std::unique_ptr<DBquadvar, decltype(&DBFreeQuadvar)> error(DBGetQuadvar(file.get(), ("block0/" + f.name + "Error").c_str()), &DBFreeQuadvar);
		ASSERT_TRUE(exact && error && exact->nels == int(f.exact.size()) && error->nels == exact->nels) << "Analytic Silo arrays";
		ASSERT_EQ(exact->nvals, 1);
		ASSERT_EQ(error->nvals, 1);
		EXPECT_TRUE(exact->units && f.units == exact->units && error->units && f.units == error->units) << "Analytic Silo CGS units";
		for (std::size_t i = 0; i < f.exact.size(); ++i) {
			EXPECT_TRUE(static_cast<double const*>(exact->vals[0])[i] == f.exact[i]) << "Analytic Silo value";
			EXPECT_TRUE(static_cast<double const*>(error->vals[0])[i] == f.numerical[i] - f.exact[i]) << "Signed Silo error";
		}
		for (auto suffix : {"Exact", "Error"}) {
			std::unique_ptr<DBmultivar, decltype(&DBFreeMultivar)> multiComparison(DBGetMultivar(file.get(), (f.name + suffix).c_str()), &DBFreeMultivar);
			ASSERT_TRUE(multiComparison && multiComparison->nvars == 1) << "Analytic Silo multivar";
			EXPECT_EQ(multiComparison->tensor_rank, DB_VARTYPE_SCALAR);
			EXPECT_STREQ(multiComparison->mmesh_name, "mesh");
		}
	}
	std::unique_ptr<DBdefvars, decltype(&DBFreeDefvars)> expressions(DBGetDefvars(file.get(), "expressions"), &DBFreeDefvars);
	ASSERT_TRUE(expressions);
	EXPECT_EQ(expressions->ndefs, int(expectedExpressions.size()));
	for (auto stem : {"momentum", "velocity", "radiationFlux", "acceleration"}) {
		for (int axis = ndim; axis < 3; ++axis) {
			for (auto suffix : {"", "Exact", "Error"}) {
				auto const name = std::string(stem) + "XYZ"[axis] + suffix;
				EXPECT_EQ(DBInqVarExists(file.get(), name.c_str()), 0) << "Inactive component: " << name;
				EXPECT_EQ(DBInqVarExists(file.get(), ("block0/" + name).c_str()), 0) << "Inactive component: " << name;
			}
		}
	}
	std::cout << problem << ' ' << dimensions << "D Silo roundtrip passed\n";
}

}	 // namespace


TEST(SiloOutput, GeometryFieldsUnitsExactErrorsAndOverwriteRoundTrip) {
	octotigerII::test::TemporaryDirectory directory;
	check(octotigerII::test::problem, ndim, directory.path);
}

TEST(SiloOutput, MixedLevelsContainOnlyLeavesAndReportTheirGeometry) {
	octotigerII::test::TemporaryDirectory directory;
	auto c = octotigerII::test::parseConfig({"--mesh.level=1", "--mesh.cells=4", "--amr.enabled=on", "--amr.maxLevel=2", "--output.directory=" + directory.path.string()});
	std::vector<Snapshot> leaves;
	mesh::BlockLocation const root;
	for (int slot = 0; slot < (1 << ndim); ++slot) {
		auto const child = root.child(slot);
		if (slot == 0) {
			for (int fine = 0; fine < (1 << ndim); ++fine) {
				leaves.push_back(initialSnapshot(c, child.child(fine)));
			}
		} else
			leaves.push_back(initialSnapshot(c, child));
	}
	Output output(c);
	output(leaves, 0, diagnose(leaves, c));
	std::unique_ptr<DBfile, decltype(&DBClose)> file(DBOpen((directory.path / "frame_000000.silo").c_str(), DB_HDF5, DB_READ), &DBClose);
	ASSERT_TRUE(file);
	std::unique_ptr<DBmultimesh, decltype(&DBFreeMultimesh)> multi(DBGetMultimesh(file.get(), "mesh"), &DBFreeMultimesh);
	std::unique_ptr<DBmultivar, decltype(&DBFreeMultivar)> levels(DBGetMultivar(file.get(), "refinementLevel"), &DBFreeMultivar);
	ASSERT_TRUE(multi);
	ASSERT_TRUE(levels);
	EXPECT_EQ(multi->nblocks, int(leaves.size()));
	EXPECT_EQ(levels->nvars, int(leaves.size()));
	for (std::size_t b = 0; b < leaves.size(); ++b) {
		auto const prefix = "block" + std::to_string(b) + "/";
		std::unique_ptr<DBquadmesh, decltype(&DBFreeQuadmesh)> mesh(DBGetQuadmesh(file.get(), (prefix + "mesh").c_str()), &DBFreeQuadmesh);
		std::unique_ptr<DBquadvar, decltype(&DBFreeQuadvar)> level(DBGetQuadvar(file.get(), (prefix + "refinementLevel").c_str()), &DBFreeQuadvar);
		ASSERT_TRUE(mesh);
		ASSERT_TRUE(level);
		EXPECT_EQ(level->nels, int(leaves[b].layout.interiorCellCount()));
		for (int i = 0; i < level->nels; ++i) {
			EXPECT_DOUBLE_EQ(static_cast<double const*>(level->vals[0])[i], leaves[b].location.level);
		}
		for (int d = 0; d < ndim; ++d) {
			ASSERT_EQ(mesh->dims[d], 5);
			auto const* coordinate = static_cast<double const*>(mesh->coords[d]);
			for (int i = 0; i <= 4; ++i) {
				EXPECT_DOUBLE_EQ(coordinate[i], units::value(leaves[b].lower[d] + Real(i) * leaves[b].cellWidth));
			}
		}
	}
}

TEST(SiloOutput, AmrTimeSeriesRefreshesDomainsAndPreservesCoverage) {
	using std::lround;

	octotigerII::test::TemporaryDirectory directory;
	auto c = octotigerII::test::parseConfig({"--mesh.level=1", "--mesh.cells=4", "--amr.enabled=on", "--amr.maxLevel=2", "--output.every=1",
		"--verification.analytic=off", "--output.directory=" + directory.path.string()});
	int constexpr children = 1 << ndim;
	std::array<int, 4> const counts{children, 2 * children - 1, children * children, children};
	Output output(c);
	for (int frame = 0; frame < 4; ++frame) {
		std::vector<Snapshot> leaves;
		mesh::BlockLocation const root;
		for (int slot = 0; slot < children; ++slot) {
			auto const child = root.child(slot);
			if (frame == 2 || (frame == 1 && slot == 0)) {
				for (int fine = 0; fine < children; ++fine) {
					leaves.push_back(initialSnapshot(c, child.child(fine)));
				}
			} else {
				leaves.push_back(initialSnapshot(c, child));
			}
		}
		for (auto& leaf : leaves) {
			leaf.time = units::Time::from_value(0.01 * frame);
		}
		output(leaves, frame, diagnose(leaves, c));
	}

	std::ifstream series(directory.path / "frames.visit");
	for (int frame = 0; frame < 4; ++frame) {
		std::string filename;
		ASSERT_TRUE(bool(std::getline(series, filename)));
		std::unique_ptr<DBfile, decltype(&DBClose)> file(DBOpen((directory.path / filename).c_str(), DB_HDF5, DB_READ), &DBClose);
		ASSERT_TRUE(file);
		// VisIt's Silo reader uses these root flags to invalidate both caches.
		for (char const* name : {"MetadataIsTimeVarying", "ConnectivityIsTimeVarying"}) {
			ASSERT_EQ(DBInqVarExists(file.get(), name), 1) << name;
			ASSERT_EQ(DBGetVarType(file.get(), name), DB_INT);
			ASSERT_EQ(DBGetVarLength(file.get(), name), 1);
			int value = 0;
			ASSERT_EQ(DBReadVar(file.get(), name, &value), 0);
			EXPECT_EQ(value, 1);
		}
		std::unique_ptr<DBmultimesh, decltype(&DBFreeMultimesh)> multi(DBGetMultimesh(file.get(), "mesh"), &DBFreeMultimesh);
		std::unique_ptr<DBmultivar, decltype(&DBFreeMultivar)> levels(DBGetMultivar(file.get(), "refinementLevel"), &DBFreeMultivar);
		ASSERT_TRUE(multi);
		ASSERT_TRUE(levels);
		ASSERT_EQ(multi->nblocks, counts[frame]);
		ASSERT_EQ(levels->nvars, counts[frame]);
		if constexpr (octotigerII::test::hydro) {
			checkVectorExpression(file.get(), "momentum", "", counts[frame]);
			checkVectorExpression(file.get(), "velocity", "", counts[frame]);
		}
		if constexpr (octotigerII::test::radiation) {
			checkVectorExpression(file.get(), "radiationFlux", "", counts[frame]);
		}
		if constexpr (octotigerII::test::gravity) {
			checkVectorExpression(file.get(), "acceleration", "", counts[frame]);
		}
		std::unique_ptr<DBdefvars, decltype(&DBFreeDefvars)> expressions(DBGetDefvars(file.get(), "expressions"), &DBFreeDefvars);
		ASSERT_TRUE(expressions);
		EXPECT_EQ(expressions->ndefs, 2 * int(octotigerII::test::hydro) + int(octotigerII::test::radiation) + int(octotigerII::test::gravity)) << "Analytic expressions disabled";
		// Count coverage on the finest block lattice using only coordinates read
		// from disk: every region must appear once, with no holes or overlaps.
		auto const lattice = mesh::filledCoordinates(4);
		std::vector<int> coverage(mesh::MeshLayout(4).interiorCellCount());
		for (int b = 0; b < multi->nblocks; ++b) {
			std::unique_ptr<DBquadmesh, decltype(&DBFreeQuadmesh)> grid(DBGetQuadmesh(file.get(), multi->meshnames[b]), &DBFreeQuadmesh);
			ASSERT_TRUE(grid);
			ASSERT_EQ(grid->ndims, ndim);
			mesh::Coordinates lower{}, extent{};
			for (int d = 0; d < ndim; ++d) {
				auto const* x = static_cast<double const*>(grid->coords[d]);
				auto toIndex = [&](double value) { return int(lround(4 * (value - units::value(c.mesh.lower)) / units::value(c.mesh.upper - c.mesh.lower))); };
				lower[d] = toIndex(x[0]);
				int const upper = toIndex(x[grid->dims[d] - 1]);
				ASSERT_GE(lower[d], 0);
				ASSERT_LE(upper, 4);
				ASSERT_GT(upper, lower[d]);
				extent[d] = upper - lower[d];
			}
			mesh::forEachCoordinate(extent, [&](auto cell) {
				for (int d = 0; d < ndim; ++d) {
					cell[d] += lower[d];
				}
				++coverage[mesh::linearIndex(cell, lattice)];
			});
		}
		for (int count : coverage) {
			EXPECT_EQ(count, 1) << filename;
		}
	}
}
