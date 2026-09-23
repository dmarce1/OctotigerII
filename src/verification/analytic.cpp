#include "octotigerII/verification/analytic.hpp"
#include "octotigerII/profiling.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include "octotigerII/problems.hpp"
#if OCTOTIGERII_GRAVITY
#include "octotigerII/verification/directGravity.hpp"
#endif


namespace octotigerII::verification {


namespace {

	std::string jsonString(std::string const& value) {
		std::string result = "\"";
		char const* hex = "0123456789abcdef";
		for (unsigned char c : value) {
			if (c == '"' || c == '\\') {
				result += '\\';
				result += c;
			} else if (c < 0x20) {
				result += "\\u00";
				result += hex[c / 16];
				result += hex[c % 16];
			} else
				result += c;
		}
		return result + '"';
	}

	void relativeJson(std::ostream& out, Real error, Real norm) {
		if (norm > 0)
			out << error / norm;
		else
			out << "null";
	}

}	 // namespace


Reference reference(Config const& c, units::Time time) {
	Reference result;
	if (c.verification.analytic == "off") {
		result.reason = "Analytic comparison disabled";
		return result;
	}
	if (usesDirectGravity(c)) {
		result.name = "Direct summation of cell-center masses";
		result.reason = "Direct gravity comparison requires all snapshots; no full-grid exact Silo fields are generated";
		return result;
	}
	result = problemReference(c);
	if (time > result.validUntil) result.evaluate = {};
	if (!result.evaluate && (c.verification.analytic == "on" || c.verification.relativeL1Tolerance >= 0))
		throw std::invalid_argument("Analytic comparison unavailable: " + result.reason);
	return result;
}


std::vector<Field> sample(Snapshot const& b, [[maybe_unused]] Config const& c, Reference const& ref) {
	std::vector<Field> result;
	if (!ref.evaluate) return result;
	if (b.time > ref.validUntil) throw std::invalid_argument("Analytic reference sampled after its validity interval");
	std::vector<ExactState> exact;
	b.layout.forEachInterior(
		[&](mesh::Coordinates const& cell, std::size_t) { exact.push_back(ref.evaluate(b.layout.cellCenter(b.lower, b.cellWidth, cell), b.time)); });
	auto field = [&](std::string name, std::string units, auto numerical, auto expected) {
		Field f{std::move(name), std::move(units), {}, {}};
		std::size_t j = 0;
		b.layout.forEachInterior([&](mesh::Coordinates const&, std::size_t i) {
			auto const actual = numerical(i), target = expected(exact[j++]);
			if (!units::finite(actual) || !units::finite(target)) throw std::runtime_error("Nonfinite analytic comparison: " + f.name);
			f.numerical.push_back(units::value(actual));
			f.exact.push_back(units::value(target));
		});
		result.push_back(std::move(f));
	};
#if OCTOTIGERII_HYDRO
	{
		hydro::HydroSystem const gas(c.hydro.gamma);
		field("density", "g/cm^3", [&](std::size_t i) { return b.hydro.values()[i].density(); }, [](auto const& q) { return q.hydro.density(); });
		field(
			"pressure", "dyn/cm^2", [&](std::size_t i) { return gas.reconstructionVariables(b.hydro.values()[i]).pressure(); },
			[](auto const& q) { return q.hydro.pressure(); });
		field(
			"gasEnergy", "erg/cm^3", [&](std::size_t i) { return b.hydro.values()[i].totalEnergy(); },
			[&](auto const& q) { return gas.conservedState(q.hydro).totalEnergy(); });
		for (int axis = 0; axis < ndim; ++axis) {
			field(
				std::string("velocity") + "XYZ"[axis], "cm/s",
				[&](std::size_t i) { return b.hydro.values()[i].momentum(axis) / b.hydro.values()[i].density(); },
				[axis](auto const& q) { return q.hydro.velocity(axis); });
			field(
				std::string("momentum") + "XYZ"[axis], "g/(cm^2 s)", [&](std::size_t i) { return b.hydro.values()[i].momentum(axis); },
				[axis](auto const& q) { return q.hydro.density() * q.hydro.velocity(axis); });
		}
	}
#endif
	if constexpr (build::gravity) {
		field("potential", "cm^2/s^2", [&](std::size_t i) { return b.gravity.values()[i].potential(); }, [](auto const& q) { return q.gravity.potential(); });
		for (int axis = 0; axis < ndim; ++axis)
			field(
				std::string("acceleration") + "XYZ"[axis], "cm/s^2", [&](std::size_t i) { return b.gravity.values()[i].acceleration(axis); },
				[axis](auto const& q) { return q.gravity.acceleration(axis); });
	}
	if constexpr (build::radiation) {
		field(
			"radiationEnergy", "erg/cm^3", [&](std::size_t i) { return b.radiation.values()[i].energy(); }, [](auto const& q) { return q.radiation.energy(); });
		for (int axis = 0; axis < ndim; ++axis)
			field(
				std::string("radiationFlux") + "XYZ"[axis], "erg/(cm^2 s)", [&](std::size_t i) { return b.radiation.values()[i].radiativeFlux(axis); },
				[axis](auto const& q) { return q.radiation.radiativeFlux(axis); });
	}
	return result;
}


Comparison compare(std::vector<Snapshot> const& snapshots, Config const& c) {
	using std::abs;
	using std::sqrt;

	profiling::Region profile("verification.compare");

	if (snapshots.empty()) throw std::invalid_argument("Analytic comparison needs snapshots");

#if OCTOTIGERII_GRAVITY
	if (usesDirectGravity(c)) return compareDirectGravity(snapshots, c);
#endif
	Comparison result;
	result.time = snapshots.front().time;
	result.cellsPerAxis = c.mesh.cells * (1 << c.mesh.level);
	result.multipoleOrder = c.gravity.multipoleOrder;
	result.openingAngle = c.gravity.openingAngle;
	auto const ref = reference(c, result.time);
	result.name = ref.name;
	result.status = ref.evaluate ? "available" : c.verification.analytic == "off" ? "disabled" : "unavailable";
	if (!ref.evaluate) {
		result.reason = ref.reason;
		return result;
	}
	Real measure = 0;
	for (auto const& block : snapshots) {
		if (!units::finite(block.time) || block.time != result.time) throw std::invalid_argument("Analytic comparison needs synchronized finite times");
		Real const volume = units::value(block.layout.cellMeasure(block.cellWidth));
		auto const fields = sample(block, c, ref);
		if (result.fields.empty())
			for (auto const& field : fields)
				result.fields.push_back(ErrorNorm{field.name, field.units});
		result.cells += block.layout.interiorCellCount();
		measure += volume * block.layout.interiorCellCount();
		for (std::size_t f = 0; f < fields.size(); ++f) {
			auto& error = result.fields[f];
			for (std::size_t i = 0; i < fields[f].exact.size(); ++i) {
				Real const e = abs(fields[f].numerical[i] - fields[f].exact[i]), r = abs(fields[f].exact[i]);
				error.l1 += volume * e;
				error.l2 += volume * e * e;
				error.linf = std::max(error.linf, e);
				error.referenceL1 += volume * r;
				error.referenceL2 += volume * r * r;
				error.referenceLinf = std::max(error.referenceLinf, r);
			}
		}
	}
	result.totalCells = result.cells;
	for (auto& error : result.fields) {
		error.l1 /= measure;
		error.l2 = sqrt(error.l2 / measure);
		error.referenceL1 /= measure;
		error.referenceL2 = sqrt(error.referenceL2 / measure);
	}
	return result;
}


void Comparison::enforce(Config const& c) const {
	if (c.verification.relativeL1Tolerance < 0) return;
	if (status != "available") throw std::runtime_error("Analytic accuracy gate requires an available reference");
	for (auto const& f : fields) {
		bool const failed = f.referenceL1 > 0 ? f.l1 / f.referenceL1 > c.verification.relativeL1Tolerance : f.linf > c.verification.absoluteTolerance;
		if (failed) throw std::runtime_error("Analytic accuracy tolerance exceeded: " + f.name);
	}
}


void Comparison::print(std::ostream& out) const {
	out << "Analytic comparison: " << status;
	if (status != "available") {
		out << " (" << reason << ")\n";
		return;
	}
	out << " | " << name << " | cell centers | t=" << units::value(time) << " s | relative L1/L2/Linf\n";
	if (referenceKind == "direct")
		out << "  direct targets=" << cells << "/" << totalCells << " sources=" << sourceCells
			<< (cells < totalCells ? " | sampled norms (Linf is sample maximum)" : " | all targets") << " | mt19937_64 seed=" << randomSeed << '\n';
	for (auto const& f : fields) {
		out << "  " << f.name << " L1=";
		if (f.referenceL1 > 0)
			out << f.l1 / f.referenceL1;
		else
			out << "undefined";
		out << " L2=";
		if (f.referenceL2 > 0)
			out << f.l2 / f.referenceL2;
		else
			out << "undefined";
		out << " Linf=";
		if (f.referenceLinf > 0)
			out << f.linf / f.referenceLinf;
		else
			out << "undefined";
		if (f.referenceL1 == 0 || f.referenceL2 == 0 || f.referenceLinf == 0)
			out << " (zero analytic field; absolute [" << f.units << "] L1=" << f.l1 << " L2=" << f.l2 << " Linf=" << f.linf << ')';
		out << '\n';
		if (referenceKind == "direct" && cells < totalCells) {
			if (f.samplingUncertaintyAvailable)
				out << "    estimated 95% sampling half-width: L1=" << f.relativeL1HalfWidth95 << " L2=" << f.relativeL2HalfWidth95 << '\n';
			if (f.samplingWarning)
				out << "    WARNING: "
					<< (f.samplingUncertaintyAvailable ? "sampling uncertainty exceeds the measured error" : "sampling uncertainty cannot be estimated")
					<< "; increase verification.directSamples\n";
		}
	}
}


void Comparison::writeJson(std::ostream& out) const {
	profiling::Region profile("output.verification_json");
	out << std::scientific << std::setprecision(17);
	out << "{\n  \"schemaVersion\": 3,\n  \"problem\": " << jsonString(build::problem) << ",\n  \"ndim\": " << ndim << ",\n  \"status\": " << jsonString(status)
		<< ",\n  \"reference\": " << jsonString(name) << ",\n  \"reason\": " << jsonString(reason) << ",\n  \"sampling\": " << jsonString(sampling)
		<< ",\n  \"time\": " << units::value(time) << ",\n  \"cells\": " << cells << ",\n  \"cellsPerAxis\": " << cellsPerAxis
		<< ",\n  \"referenceKind\": " << jsonString(referenceKind) << ",\n  \"totalCells\": " << totalCells << ",\n  \"sourceCells\": " << sourceCells
		<< ",\n  \"multipoleOrder\": " << multipoleOrder << ",\n  \"openingAngle\": " << openingAngle << ",\n  \"targetIndices\": [";
	for (std::size_t i = 0; i < targetIndices.size(); ++i)
		out << (i ? "," : "") << targetIndices[i];
	out << "]";
	if (referenceKind == "direct") out << ",\n  \"randomGenerator\": \"mt19937_64\",\n  \"randomSeed\": " << randomSeed;
	out << ",\n  \"fields\": [";
	for (std::size_t i = 0; i < fields.size(); ++i) {
		auto const& f = fields[i];
		out << (i ? "," : "") << "\n    {\"name\": " << jsonString(f.name) << ", \"units\": " << jsonString(f.units) << ", \"L1\": ";
		relativeJson(out, f.l1, f.referenceL1);
		out << ", \"L2\": ";
		relativeJson(out, f.l2, f.referenceL2);
		out << ", \"Linf\": ";
		relativeJson(out, f.linf, f.referenceLinf);
		out << ", \"absoluteL1\": " << f.l1 << ", \"absoluteL2\": " << f.l2 << ", \"absoluteLinf\": " << f.linf << ", \"referenceL1\": " << f.referenceL1
			<< ", \"referenceL2\": " << f.referenceL2 << ", \"referenceLinf\": " << f.referenceLinf;
		if (referenceKind == "direct") {
			out << ", \"samplingWarning\": " << (f.samplingWarning ? "true" : "false") << ", \"relativeL1HalfWidth95\": ";
			if (f.samplingUncertaintyAvailable)
				out << f.relativeL1HalfWidth95;
			else
				out << "null";
			out << ", \"relativeL2HalfWidth95\": ";
			if (f.samplingUncertaintyAvailable)
				out << f.relativeL2HalfWidth95;
			else
				out << "null";
		}
		out << '}';
	}
	out << "\n  ]\n}\n";
}

}	 // namespace octotigerII::verification
