#include "octotigerII/output.hpp"
#include <filesystem>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#ifdef OCTOTIGERII_WITH_SILO
#include <silo.h>
#endif

namespace octotigerII {
namespace {
struct Variable {
	std::string name;
	std::vector<double> values;
};
std::vector<Variable> variables(Snapshot const& b, Config const& c) {
	std::vector<Variable> result;
	auto field = [&](std::string name, auto value) {
		Variable v{std::move(name), {}};
		b.layout.forEachInterior(
			[&](mesh::Coordinates const&, std::size_t i) { v.values.push_back(value(i)); });
		result.push_back(std::move(v));
	};
	if (b.hydroEnabled) {
		std::array<std::string, 5> const names{"density", "momentumX", "momentumY", "momentumZ",
											   "gasEnergy"};
		for (int f = 0; f < 5; ++f)
			field(names[f], [&](std::size_t i) { return b.hydro.values()[i][f]; });
		field("pressure", [&](std::size_t i) {
			return hydro::HydroSystem(c.gamma)
				.reconstructionVariables(b.hydro.values()[i])
				.pressure();
		});
	} else if (b.gravityEnabled)
		field("density", [&](std::size_t i) { return b.density.values()[i]; });
	if (b.radiationEnabled) {
		std::array<std::string, 4> const names{"radiationEnergy", "radiationFluxX",
											   "radiationFluxY", "radiationFluxZ"};
		for (int f = 0; f < 4; ++f)
			field(names[f], [&](std::size_t i) {
				return b.radiation.values()[i][f] * (f == 0 ? 1 : physicalLightSpeed);
			});
	}
	if (b.gravityEnabled) {
		std::array<std::string, 4> const names{"potential", "accelerationX", "accelerationY",
											   "accelerationZ"};
		for (int f = 0; f < 4; ++f)
			field(names[f], [&](std::size_t i) { return b.gravity.values()[i][f]; });
	}
	return result;
}
void writeCsv(std::vector<Snapshot> const& patches, Config const& c, std::string const& filename) {
	std::ofstream file(filename);
	file.exceptions(std::ios::badbit | std::ios::failbit);
	file << "time_s,block,x_cm,y_cm,z_cm,dx_cm";
	for (auto const& v : variables(patches.front(), c))
		file << ',' << v.name;
	file << '\n' << std::scientific << std::setprecision(17);
	for (std::size_t b = 0; b < patches.size(); ++b) {
		auto const& patch = patches[b];
		auto const fields = variables(patch, c);
		std::size_t interior = 0;
		patch.layout.forEachInterior([&](mesh::Coordinates const& cell, std::size_t) {
			auto const point = patch.layout.cellCenter(patch.lower, patch.cellWidth, cell);
			file << patch.time << ',' << b;
			for (Real coordinate : point)
				file << ',' << coordinate;
			file << ',' << patch.cellWidth;
			for (auto const& v : fields)
				file << ',' << v.values[interior];
			file << '\n';
			++interior;
		});
	}
}
#ifdef OCTOTIGERII_WITH_SILO
void writeSilo(std::vector<Snapshot> const& patches, Config const& c, std::string const& filename,
			   int cycle, Real time) {
	std::unique_ptr<DBfile, decltype(&DBClose)> file(
		DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL, "OctotigerII (cgs)", DB_HDF5), &DBClose);
	if (!file)
		throw std::runtime_error("Cannot create Silo file: " + filename);
	std::unique_ptr<DBoptlist, decltype(&DBFreeOptlist)> options(DBMakeOptlist(2), &DBFreeOptlist);
	if (!options)
		throw std::runtime_error("Cannot allocate Silo options");
	int status = DBAddOption(options.get(), DBOPT_DTIME, &time);
	status |= DBAddOption(options.get(), DBOPT_CYCLE, &cycle);
	std::vector<std::string> meshes;
	auto const fieldList = variables(patches.front(), c);
	std::vector<std::vector<std::string>> names(fieldList.size());
	for (std::size_t b = 0; b < patches.size(); ++b) {
		auto const& patch = patches[b];
		std::string const block = "block" + std::to_string(b);
		status |= DBMkDir(file.get(), block.c_str());
		status |= DBSetDir(file.get(), block.c_str());
		int meshDims[3]{1, 1, 1}, zoneDims[3]{1, 1, 1};
		std::array<std::vector<double>, 3> coords;
		void* coordPointers[3]{};
		int const ndim = patch.layout.dimensionCount();
		for (int d = 0; d < ndim; ++d) {
			zoneDims[d] = patch.layout.interiorExtent(d);
			meshDims[d] = zoneDims[d] + 1;
			for (int i = 0; i < meshDims[d]; ++i)
				coords[d].push_back(patch.lower[d] + i * patch.cellWidth);
			coordPointers[d] = coords[d].data();
		}
		status |= DBPutQuadmesh(file.get(), "mesh", nullptr, coordPointers, meshDims, ndim,
								DB_DOUBLE, DB_COLLINEAR, options.get());
		meshes.push_back(block + "/mesh");
		auto const fields = variables(patch, c);
		for (std::size_t f = 0; f < fields.size(); ++f) {
			status |=
				DBPutQuadvar1(file.get(), fields[f].name.c_str(), "mesh", fields[f].values.data(),
							  zoneDims, ndim, nullptr, 0, DB_DOUBLE, DB_ZONECENT, options.get());
			names[f].push_back(block + "/" + fields[f].name);
		}
		status |= DBSetDir(file.get(), "/");
	}
	std::vector<char const*> pointers;
	std::vector<int> types(patches.size(), DB_QUAD_RECT);
	for (auto const& mesh : meshes)
		pointers.push_back(mesh.c_str());
	status |= DBPutMultimesh(file.get(), "mesh", static_cast<int>(patches.size()), pointers.data(),
							 types.data(), options.get());
	std::fill(types.begin(), types.end(), DB_QUADVAR);
	for (std::size_t f = 0; f < names.size(); ++f) {
		pointers.clear();
		for (auto const& name : names[f])
			pointers.push_back(name.c_str());
		status |=
			DBPutMultivar(file.get(), fieldList[f].name.c_str(), static_cast<int>(patches.size()),
						  pointers.data(), types.data(), options.get());
	}
	status |= DBClose(file.release());
	if (status < 0)
		throw std::runtime_error("Silo write failed: " + filename);
}
#endif
} // namespace
Output::Output(Config const& c) : config_(c) {
	if (!c.outputEnabled)
		return;
	std::filesystem::create_directories(c.outputDirectory);
	diagnostics_.open(std::filesystem::path(c.outputDirectory) / "diagnostics.csv");
	diagnostics_.exceptions(std::ios::badbit | std::ios::failbit);
	diagnostics_ << "step,time_s,mass,gas_energy,radiation_energy,min_density,min_pressure,min_"
					"radiation_energy,max_reduced_flux\n";
	diagnostics_ << std::scientific << std::setprecision(17);
	if (c.outputFormat == "silo") {
		series_.open(std::filesystem::path(c.outputDirectory) / "frames.visit");
		series_.exceptions(std::ios::badbit | std::ios::failbit);
	}
}
void Output::operator()(std::vector<Snapshot> const& patches, int step, Diagnostics const& d) {
	if (!config_.outputEnabled)
		return;
	diagnostics_ << step << ',' << d.time << ',' << d.mass << ',' << d.gasEnergy << ','
				 << d.radiationEnergy << ',' << d.minimumDensity << ',' << d.minimumPressure << ','
				 << d.minimumRadiationEnergy << ',' << d.maximumReducedFlux << '\n';
	diagnostics_.flush();
	if (step != 0 && step % config_.outputEvery != 0 && d.time < config_.stopTime)
		return;
	std::ostringstream base;
	base << "frame_" << std::setfill('0') << std::setw(6) << frame_++ << '.'
		 << config_.outputFormat;
	std::string const filename =
		(std::filesystem::path(config_.outputDirectory) / base.str()).string();
	if (config_.outputFormat == "csv")
		writeCsv(patches, config_, filename);
#ifdef OCTOTIGERII_WITH_SILO
	else {
		writeSilo(patches, config_, filename, step, d.time);
		series_ << base.str() << '\n';
		series_.flush();
	}
#else
	else
		throw std::logic_error("Silo support was not compiled");
#endif
}
} // namespace octotigerII
