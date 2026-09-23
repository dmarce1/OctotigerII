#include "octotigerII/output.hpp"
#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <memory>
#include <silo.h>
#include <sstream>
#include <stdexcept>
#include "octotigerII/profiling.hpp"
#include "octotigerII/verification/analytic.hpp"

namespace octotigerII {

namespace {

	class Variable {
	public:
		std::string name;
		std::string units;
		bool isVector = false;
		std::vector<std::vector<double>> values;
	};

	std::vector<Variable> variables(Snapshot const& b, [[maybe_unused]] Config const& c) {
		std::vector<Variable> result;
		if (c.amr.enabled) result.push_back({"refinementLevel", "", false, {std::vector<double>(b.layout.interiorCellCount(), b.location.level)}});
		auto field = [&](std::string name, std::string units, bool isVector, auto value) {
			int const components = isVector ? ndim : 1;
			Variable v{std::move(name), std::move(units), isVector, std::vector<std::vector<double>>(components)};
			for (int axis = 0; axis < components; ++axis) {
				v.values[axis].reserve(b.layout.interiorCellCount());
				b.layout.forEachInterior([&](mesh::Coordinates const&, std::size_t i) { v.values[axis].push_back(units::value(value(i, axis))); });
			}
			result.push_back(std::move(v));
		};
#if OCTOTIGERII_HYDRO
		{
			field("density", "g/cm^3", false, [&](std::size_t i, int) { return b.hydro.values()[i].density(); });
			field("momentum", "g/(cm^2 s)", true, [&](std::size_t i, int axis) { return b.hydro.values()[i].momentum(axis); });
			field("velocity", "cm/s", true, [&](std::size_t i, int axis) { return b.hydro.values()[i].momentum(axis) / b.hydro.values()[i].density(); });
			field("gasEnergy", "erg/cm^3", false, [&](std::size_t i, int) { return b.hydro.values()[i].totalEnergy(); });
			field("pressure", "dyn/cm^2", false,
				[&](std::size_t i, int) { return hydro::HydroSystem(c.hydro.gamma).reconstructionVariables(b.hydro.values()[i]).pressure(); });
		}
#elif OCTOTIGERII_GRAVITY
		{
			field("density", "g/cm^3", false, [&](std::size_t i, int) { return b.density.values()[i]; });
		}
#endif
		if constexpr (build::radiation) {
			field("radiationEnergy", "erg/cm^3", false, [&](std::size_t i, int) { return b.radiation.values()[i].energy(); });
			field("radiationFlux", "erg/(cm^2 s)", true, [&](std::size_t i, int axis) { return b.radiation.values()[i].radiativeFlux(axis); });
		}
		if constexpr (build::gravity) {
			field("potential", "cm^2/s^2", false, [&](std::size_t i, int) { return b.gravity.values()[i].potential(); });
			field("acceleration", "cm/s^2", true, [&](std::size_t i, int axis) { return b.gravity.values()[i].acceleration(axis); });
		}

		auto const reference = verification::reference(c, b.time);
		auto const comparisons = verification::sample(b, c, reference);
		auto const numericalCount = result.size();
		for (std::size_t v = 0; v < numericalCount; ++v) {
			auto const& numerical = result[v];
			Variable exact{numerical.name + "Exact", numerical.units, numerical.isVector, {}};
			Variable error{numerical.name + "Error", numerical.units, numerical.isVector, {}};
			for (std::size_t axis = 0; axis < numerical.values.size(); ++axis) {
				auto const name = numerical.name + (numerical.isVector ? std::string(1, "XYZ"[axis]) : "");
				auto const match = std::find_if(comparisons.begin(), comparisons.end(), [&](auto const& f) { return f.name == name; });
				if (match == comparisons.end()) break;
				exact.values.push_back(match->exact);
				auto& component = error.values.emplace_back();
				component.reserve(match->exact.size());
				for (std::size_t i = 0; i < match->exact.size(); ++i) {
					component.push_back(match->numerical[i] - match->exact[i]);
				}
			}
			if (exact.values.empty()) continue;
			if (exact.values.size() != numerical.values.size()) throw std::logic_error("Incomplete analytic Silo vector: " + numerical.name);
			result.push_back(std::move(exact));
			result.push_back(std::move(error));
		}
		return result;
	}

	void writeSilo(std::vector<Snapshot> const& patches, Config const& c, std::string const& filename, int cycle, units::Time time) {
		profiling::Region profile("output.silo");
		std::unique_ptr<DBfile, decltype(&DBClose)> file(DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL, "OctotigerII (cgs)", DB_HDF5), &DBClose);
		if (!file) throw std::runtime_error("Cannot create Silo file: " + filename);
		int status = 0;
		if (c.amr.enabled) {
			// VisIt otherwise reuses the first frame's domain list and connectivity
			// after refinement, displaying only a prefix of the new leaf blocks.
			int const enabled = 1;
			int const length = 1;
			status |= DBWrite(file.get(), "MetadataIsTimeVarying", &enabled, &length, 1, DB_INT);
			status |= DBWrite(file.get(), "ConnectivityIsTimeVarying", &enabled, &length, 1, DB_INT);
		}
		std::unique_ptr<DBoptlist, decltype(&DBFreeOptlist)> options(DBMakeOptlist(8), &DBFreeOptlist);
		if (!options) throw std::runtime_error("Cannot allocate Silo options");
		Real outputTime = units::value(time);
		status |= DBAddOption(options.get(), DBOPT_DTIME, &outputTime);
		status |= DBAddOption(options.get(), DBOPT_CYCLE, &cycle);
		char centimeters[] = "cm";
		int const coordinateUnits[] = {DBOPT_XUNITS, DBOPT_YUNITS, DBOPT_ZUNITS};
		for (int axis = 0; axis < ndim; ++axis) {
			status |= DBAddOption(options.get(), coordinateUnits[axis], centimeters);
		}
		std::vector<std::string> meshes;
		auto const fieldList = variables(patches.front(), c);
		std::vector<std::vector<std::string>> names(fieldList.size());
		for (std::size_t b = 0; b < patches.size(); ++b) {
			auto const& patch = patches[b];
			std::string const block = "block" + std::to_string(b);
			status |= DBMkDir(file.get(), block.c_str());
			status |= DBSetDir(file.get(), block.c_str());
			int meshDims[ndim]{}, zoneDims[ndim]{};
			std::array<std::vector<double>, ndim> coords;
			void* coordPointers[ndim]{};
			for (int d = 0; d < ndim; ++d) {
				zoneDims[d] = patch.layout.interiorExtent(d);
				meshDims[d] = zoneDims[d] + 1;
				for (int i = 0; i < meshDims[d]; ++i) {
					coords[d].push_back(units::value(patch.lower[d] + Real(i) * patch.cellWidth));
				}
				coordPointers[d] = coords[d].data();
			}
			status |= DBPutQuadmesh(file.get(), "mesh", nullptr, coordPointers, meshDims, ndim, DB_DOUBLE, DB_COLLINEAR, options.get());
			meshes.push_back(block + "/mesh");
			auto const fields = variables(patch, c);
			for (std::size_t f = 0; f < fields.size(); ++f) {
				status |= DBAddOption(options.get(), DBOPT_UNITS, const_cast<char*>(fields[f].units.c_str()));
				if (fields[f].isVector) {
					std::array<std::string, ndim> componentNames;
					std::array<char const*, ndim> componentPointers{};
					std::array<void const*, ndim> valuePointers{};
					for (int axis = 0; axis < ndim; ++axis) {
						componentNames[axis] = fields[f].name + "_" + "xyz"[axis];
						componentPointers[axis] = componentNames[axis].c_str();
						valuePointers[axis] = fields[f].values[axis].data();
					}
					status |= DBPutQuadvar(file.get(), fields[f].name.c_str(), "mesh", ndim, componentPointers.data(), valuePointers.data(), zoneDims,
						ndim, nullptr, 0, DB_DOUBLE, DB_ZONECENT, options.get());
				} else {
					status |= DBPutQuadvar1(file.get(), fields[f].name.c_str(), "mesh", fields[f].values.front().data(), zoneDims, ndim, nullptr, 0,
						DB_DOUBLE, DB_ZONECENT, options.get());
				}
				names[f].push_back(block + "/" + fields[f].name);
				status |= DBClearOption(options.get(), DBOPT_UNITS);
			}
			status |= DBSetDir(file.get(), "/");
		}
		std::vector<char const*> pointers;
		std::vector<int> types(patches.size(), DB_QUAD_RECT);
		for (auto const& mesh : meshes) {
			pointers.push_back(mesh.c_str());
		}
		status |= DBPutMultimesh(file.get(), "mesh", static_cast<int>(patches.size()), pointers.data(), types.data(), options.get());
		std::fill(types.begin(), types.end(), DB_QUADVAR);
		char meshName[] = "mesh";
		status |= DBAddOption(options.get(), DBOPT_MMESH_NAME, meshName);
		for (std::size_t f = 0; f < names.size(); ++f) {
			pointers.clear();
			for (auto const& name : names[f]) {
				pointers.push_back(name.c_str());
			}
			int tensorRank = fieldList[f].isVector ? DB_VARTYPE_VECTOR : DB_VARTYPE_SCALAR;
			status |= DBAddOption(options.get(), DBOPT_TENSOR_RANK, &tensorRank);
			status |= DBPutMultivar(file.get(), fieldList[f].name.c_str(), static_cast<int>(patches.size()), pointers.data(), types.data(), options.get());
			status |= DBClearOption(options.get(), DBOPT_TENSOR_RANK);
		}
		status |= DBClose(file.release());
		if (status < 0) throw std::runtime_error("Silo write failed: " + filename);
	}

}	 // namespace

Output::Output(Config const& c)
  : config_(c) {
	if (!c.output.enabled) return;
	std::filesystem::create_directories(c.output.directory);
	series_.open(std::filesystem::path(c.output.directory) / "frames.visit");
	series_.exceptions(std::ios::badbit | std::ios::failbit);
}

void Output::operator()(std::vector<Snapshot> const& patches, int step, Diagnostics const& d) {
	if (!config_.output.enabled) return;
	if (step != 0 && step % config_.output.every != 0 && d.time < config_.runtime.stopTime) return;
	std::ostringstream base;
	base << "frame_" << std::setfill('0') << std::setw(6) << frame_++ << ".silo";
	std::string const filename = (std::filesystem::path(config_.output.directory) / base.str()).string();
	writeSilo(patches, config_, filename, step, d.time);
	series_ << base.str() << '\n';
	series_.flush();
	auto const comparison = verification::compare(patches, config_);
	std::ofstream report;
	report.exceptions(std::ios::badbit | std::ios::failbit);
	report.open(std::filesystem::path(config_.output.directory) / "analytic-errors.json");
	comparison.writeJson(report);
}
}	 // namespace octotigerII
