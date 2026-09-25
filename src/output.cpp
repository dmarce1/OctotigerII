#include "octotigerII/output.hpp"
#include <limits>
#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <map>
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
		std::string vectorName;
		std::vector<double> values;
	};

	std::vector<Variable> variables(Snapshot const& b, [[maybe_unused]] Config const& c) {
		std::vector<Variable> result;
		if (c.amr.enabled) result.push_back({"refinementLevel", "", "", std::vector<double>(b.layout.interiorCellCount(), b.location.level)});
		auto field = [&](std::string name, std::string units, bool isVector, auto value) {
			int const components = isVector ? ndim : 1;
			for (int axis = 0; axis < components; ++axis) {
				Variable v{name + (isVector ? std::string(1, "XYZ"[axis]) : ""), units, isVector ? name : "", {}};
				v.values.reserve(b.layout.interiorCellCount());
				b.layout.forEachInterior([&](mesh::Coordinates const&, std::size_t i) { v.values.push_back(units::value(value(i, axis))); });
				result.push_back(std::move(v));
			}
		};
#if OCTOTIGERII_HYDRO
		if (c.hydroEnabled()) {
			field("density", "g/cm^3", false, [&](std::size_t i, int) { return b.hydro.values()[i].density(); });
			for (std::size_t s = 0; s < c.massFractions.species.size(); ++s) if (c.massFractions.enabled) {
				auto const& species = c.massFractions.species[s];
				field((species.tracer() ? "tracerDensity_" : "partialDensity_") + species.name, "g/cm^3", false, [&](std::size_t i, int) { return b.species.at(s).values()[i]; });
				field((species.tracer() ? "tracer_" : "massFraction_") + species.name, "", false,
					[&](std::size_t i, int) { return b.species.at(s).values()[i] / b.hydro.values()[i].density(); });
			}
			field("momentum", "g/(cm^2 s)", true, [&](std::size_t i, int axis) { return b.hydro.values()[i].momentum(axis); });
			field("velocity", "cm/s", true, [&](std::size_t i, int axis) { return b.hydro.values()[i].momentum(axis) / b.hydro.values()[i].density(); });
			field("gasEnergy", "erg/cm^3", false, [&](std::size_t i, int) { return b.hydro.values()[i].totalEnergy(); });
			field("internalEnergy", "erg/cm^3", false, [&](std::size_t i, int) { return hydro::HydroSystem(c.hydro).internalEnergy(b.hydro.values()[i]); });
			field("temperature", "K", false, [&](std::size_t i, int) { return hydro::HydroSystem(c.hydro).temperature(b.hydro.values()[i]); });
			if (c.hydro.dualEnergy.enabled)
				field("dualEnergy", "g/cm^3", false, [&](std::size_t i, int) { return b.hydro.values()[i].auxiliary(); });
			field("pressure", "dyn/cm^2", false,
				[&](std::size_t i, int) { return hydro::HydroSystem(c.hydro).reconstructionVariables(b.hydro.values()[i]).pressure(); });
		}
#endif
#if OCTOTIGERII_GRAVITY
		if (!c.hydroEnabled() && c.gravityEnabled()) {
			field("density", "g/cm^3", false, [&](std::size_t i, int) { return b.density.values()[i]; });
		}
#endif
		if (build::radiation && c.radiationEnabled()) {
			field("radiationEnergy", "erg/cm^3", false, [&](std::size_t i, int) { return b.radiation.values()[i].energy(); });
			field("radiationFlux", "erg/(cm^2 s)", true, [&](std::size_t i, int axis) { return b.radiation.values()[i].radiativeFlux(axis); });
		}
		if (build::gravity && c.gravityEnabled()) {
			field("potential", "cm^2/s^2", false, [&](std::size_t i, int) { return b.gravity.values()[i].potential(); });
			field("acceleration", "cm/s^2", true, [&](std::size_t i, int axis) { return b.gravity.values()[i].acceleration(axis); });
		}

		auto const reference = verification::reference(c, b.time);
		auto const comparisons = verification::sample(b, c, reference);
		auto const numericalCount = result.size();
		for (std::size_t v = 0; v < numericalCount; ++v) {
			auto const& numerical = result[v];
			auto const match = std::find_if(comparisons.begin(), comparisons.end(), [&](auto const& f) { return f.name == numerical.name; });
			if (match == comparisons.end()) continue;
			Variable exact{numerical.name + "Exact", numerical.units, numerical.vectorName.empty() ? "" : numerical.vectorName + "Exact", match->exact};
			Variable error{numerical.name + "Error", numerical.units, numerical.vectorName.empty() ? "" : numerical.vectorName + "Error", {}};
			error.values.reserve(match->exact.size());
			for (std::size_t i = 0; i < match->exact.size(); ++i) {
				error.values.push_back(match->numerical[i] - match->exact[i]);
			}
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
				status |= DBPutQuadvar1(file.get(), fields[f].name.c_str(), "mesh", fields[f].values.data(), zoneDims, ndim, nullptr, 0,
					DB_DOUBLE, DB_ZONECENT, options.get());
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
		int tensorRank = DB_VARTYPE_SCALAR;
		status |= DBAddOption(options.get(), DBOPT_TENSOR_RANK, &tensorRank);
		std::map<std::string, std::vector<std::string>> vectorComponents;
		for (std::size_t f = 0; f < names.size(); ++f) {
			pointers.clear();
			for (auto const& name : names[f]) {
				pointers.push_back(name.c_str());
			}
			status |= DBPutMultivar(file.get(), fieldList[f].name.c_str(), static_cast<int>(patches.size()), pointers.data(), types.data(), options.get());
			if (!fieldList[f].vectorName.empty()) {
				vectorComponents[fieldList[f].vectorName].push_back(fieldList[f].name);
			}
		}
		if (!vectorComponents.empty()) {
			// Define vectors from the root scalar multivars so VisIt can plot both
			// components and vectors without storing a second copy of the arrays.
			std::vector<std::string> definitions;
			std::vector<char const*> expressionNames, expressionValues;
			std::vector<int> expressionTypes(vectorComponents.size(), DB_VARTYPE_VECTOR);
			for (auto const& [name, components] : vectorComponents) {
				if (components.size() != ndim) throw std::logic_error("Incomplete Silo vector expression: " + name);
				std::string definition = "{";
				for (int axis = 0; axis < 3; ++axis) {
					if (axis != 0) definition += ',';
					// VisIt needs three components. Inactive directions are expressions,
					// not stored arrays; explicit zone centering matches the scalars.
					definition += axis < ndim ? components[axis] : "zonal_constant(<mesh>,0)";
				}
				definition += '}';
				definitions.push_back(std::move(definition));
				expressionNames.push_back(name.c_str());
			}
			for (auto const& definition : definitions) {
				expressionValues.push_back(definition.c_str());
			}
			status |= DBPutDefvars(file.get(), "expressions", static_cast<int>(definitions.size()), expressionNames.data(), expressionTypes.data(),
				expressionValues.data(), nullptr);
		}
		status |= DBClose(file.release());
		if (status < 0) throw std::runtime_error("Silo write failed: " + filename);
	}

}	 // namespace

Output::Output(Config const& c)
  : config_(c) {
	std::filesystem::create_directories(c.output.directory);
	conservation_.exceptions(std::ios::badbit | std::ios::failbit);
	conservation_.open(std::filesystem::path(c.output.directory) / "conservation.csv");
	conservation_ << std::scientific << std::setprecision(std::numeric_limits<Real>::max_digits10);
	conservation_ << "step,time_s";
	auto header = [&](std::string const& name) {
		for (auto const* suffix : {"grid", "in", "out", "corrected", "l1", "norm", "drift_scaled"}) conservation_ << ',' << name << '_' << suffix;
	};
	if (c.hydroEnabled() || c.gravityEnabled()) header("mass_g");
	if (c.hydroEnabled()) {
		for (int d = 0; d < ndim; ++d) header(std::string("momentum_") + "xyz"[d] + "_g_cm_s");
		header("gas_energy_erg");
		conservation_ << ",kinetic_energy_erg_grid,thermal_energy_erg_grid";
		if (c.gravityEnabled()) conservation_ << ",potential_energy_erg_grid,potential_energy_erg_in,potential_energy_erg_out,gas_gravity_energy_erg_grid,gas_gravity_energy_erg_corrected,gas_gravity_energy_erg_norm,gas_gravity_energy_drift_scaled,gravity_reciprocity_defect_erg,gravity_regrid_energy_change_erg,gravity_energy_budget_residual_scaled";
	}
	if (c.radiationEnabled()) {
		header("radiation_energy_erg");
		for (int d = 0; d < ndim; ++d) header(std::string("radiation_flux_integral_") + "xyz"[d] + "_erg_cm_s");
	}
	conservation_ << '\n';
	if (!c.output.enabled) return;
	series_.open(std::filesystem::path(c.output.directory) / "frames.visit");
	series_.exceptions(std::ios::badbit | std::ios::failbit);
}

void Output::operator()(std::vector<Snapshot> const& patches, int step, Diagnostics const& d) {
	if (!haveInitial_) { initial_ = d; haveInitial_ = true; }
	conservation_ << step << ',' << units::value(d.time);
	auto write = [&](auto grid, auto inward, auto outward, auto l1, auto initial, auto initialL1) {
		auto const corrected = grid + outward - inward;
		auto const norm = std::max({initialL1, l1, inward + outward});
		Real const drift = norm > decltype(norm){} ? Real((corrected - initial) / norm) : Real(0);
		conservation_ << ',' << units::value(grid) << ',' << units::value(inward) << ',' << units::value(outward)
			<< ',' << units::value(corrected) << ',' << units::value(l1) << ',' << units::value(norm) << ',' << drift;
	};
	auto const& in = d.boundary.inward;
	auto const& out = d.boundary.outward;
	if (config_.hydroEnabled() || config_.gravityEnabled()) write(d.mass, in.mass, out.mass, d.norm.mass, initial_.mass, initial_.norm.mass);
	if (config_.hydroEnabled()) {
		for (int axis = 0; axis < ndim; ++axis) write(d.momentum[axis], in.momentum[axis], out.momentum[axis], d.norm.momentum[axis], initial_.momentum[axis], initial_.norm.momentum[axis]);
		write(d.gasEnergy, in.gasEnergy, out.gasEnergy, d.norm.gasEnergy, initial_.gasEnergy, initial_.norm.gasEnergy);
		conservation_ << ',' << units::value(d.kineticEnergy) << ',' << units::value(d.thermalEnergy);
		if (config_.gravityEnabled()) {
			auto const norm = std::max(initial_.gasGravityNorm, d.gasGravityNorm);
			Real const drift = norm > units::Energy{} ? Real((d.gasGravityEnergy + out.gasEnergy - in.gasEnergy + out.potentialEnergy - in.potentialEnergy - initial_.gasGravityEnergy) / norm) : Real(0);
			conservation_ << ',' << units::value(d.potentialEnergy) << ',' << units::value(in.potentialEnergy) << ',' << units::value(out.potentialEnergy)
				<< ',' << units::value(d.gasGravityEnergy)
				<< ',' << units::value(d.gasGravityEnergy + out.gasEnergy - in.gasEnergy + out.potentialEnergy - in.potentialEnergy)
				<< ',' << units::value(norm) << ',' << drift
				<< ',' << units::value(d.gravityReciprocityDefect) << ',' << units::value(d.gravityRegridEnergyChange)
				<< ',' << (norm > units::Energy{} ? drift - Real((d.gravityReciprocityDefect + d.gravityRegridEnergyChange) / norm) : Real(0));
		}
	}
	if (config_.radiationEnabled()) {
		write(d.radiationEnergy, in.radiationEnergy, out.radiationEnergy, d.norm.radiationEnergy, initial_.radiationEnergy, initial_.norm.radiationEnergy);
		for (int axis = 0; axis < ndim; ++axis) write(d.radiationFlux[axis], in.radiationFlux[axis], out.radiationFlux[axis], d.norm.radiationFlux[axis], initial_.radiationFlux[axis], initial_.norm.radiationFlux[axis]);
	}
	conservation_ << '\n';
	conservation_.flush();
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
