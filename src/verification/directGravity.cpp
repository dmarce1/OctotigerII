#include "octotigerII/verification/directGravity.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <tuple>
#include "octotigerII/gravity/ewald.hpp"
#include "octotigerII/gravity/images.hpp"
#include "octotigerII/profiling.hpp"

namespace octotigerII::verification {

std::vector<std::size_t> directTargets(std::size_t population, std::size_t count, std::uint64_t seed) {
	count = std::min(count, population);
	if (count == population) {
		std::vector<std::size_t> result(count);
		std::iota(result.begin(), result.end(), 0);
		return result;
	}
	std::mt19937_64 generator(seed);
	auto bounded = [&](std::uint64_t bound) {
		// Unsigned wraparound gives 2^64 modulo bound. Reject that many values
		// before taking the remainder, avoiding modulo bias.
		auto const threshold = (std::uint64_t(0) - bound) % bound;
		std::uint64_t value;
		do {
			value = generator();
		} while (value < threshold);
		return value % bound;
	};
	// Floyd's sampling algorithm; memory is proportional to the sample size.
	std::set<std::size_t> selected;
	for (std::size_t j = population - count; j < population; ++j) {
		auto const candidate = static_cast<std::size_t>(bounded(j + 1));
		if (!selected.insert(candidate).second) selected.insert(j);
	}
	return {selected.begin(), selected.end()};
}

ErrorNorm directErrorNorm(Field const& field, std::size_t population) {
	using std::abs;
	using std::isfinite;
	using std::sqrt;

	profiling::Region profile("verification.direct_error_norm");

	auto const count = field.exact.size();
	if (count == 0 || count > population || count != field.numerical.size()) throw std::invalid_argument("Invalid direct norm sample");
	ErrorNorm result{field.name, field.units};
	for (std::size_t i = 0; i < count; ++i) {
		Real const e = abs(field.numerical[i] - field.exact[i]), r = abs(field.exact[i]);
		if (!isfinite(e) || !isfinite(r)) throw std::runtime_error("Nonfinite direct gravity comparison");
		result.l1 += e;
		result.l2 += e * e;
		result.linf = std::max(result.linf, e);
		result.referenceL1 += r;
		result.referenceL2 += r * r;
		result.referenceLinf = std::max(result.referenceLinf, r);
	}
	result.l1 /= count;
	result.l2 = sqrt(result.l2 / count);
	result.referenceL1 /= count;
	result.referenceL2 = sqrt(result.referenceL2 / count);
	if (count == population) {
		result.samplingUncertaintyAvailable = true;
		return result;
	}
	if (count < 2 || result.referenceL1 == 0 || result.referenceL2 == 0) {
		result.samplingWarning = true;
		return result;
	}
	Real const l1 = result.l1 / result.referenceL1, l2 = result.l2 / result.referenceL2;
	Real variance1 = 0, variance2 = 0;
	for (std::size_t i = 0; i < count; ++i) {
		Real const e = abs(field.numerical[i] - field.exact[i]), r = abs(field.exact[i]);
		Real const residual1 = (e - l1 * r) / result.referenceL1;
		Real const residual2 = (e * e - l2 * l2 * r * r) / (result.referenceL2 * result.referenceL2);
		variance1 += residual1 * residual1;
		variance2 += residual2 * residual2;
	}
	// SRS without replacement: Var(mean) = (1-k/N) * s^2/k.
	Real const factor = (1 - Real(count) / Real(population)) / (Real(count) * Real(count - 1));
	result.relativeL1HalfWidth95 = 1.96 * sqrt(factor * variance1);
	result.relativeL2HalfWidth95 = l2 > 0 ? 1.96 * sqrt(factor * variance2) / (2 * l2) : 0;
	result.samplingUncertaintyAvailable = true;
	result.samplingWarning = result.relativeL1HalfWidth95 > l1 || result.relativeL2HalfWidth95 > l2;
	return result;
}

Comparison compareDirectGravity(std::vector<Snapshot> const& snapshots, Config const& c) {
	profiling::Region profile("verification.direct_gravity");
	static_assert(ndim == 3);
	if (snapshots.empty()) throw std::invalid_argument("Direct gravity comparison needs snapshots");
	gravity::ImageGeometry const images(c.mesh.boundary);
	Comparison result;
	result.name = "Direct summation of cell-center masses";
	if (images.active()) result.name = "Direct cell-center image summation with Ewald lattice correction";
	result.referenceKind = "direct";
	result.status = "available";
	result.time = snapshots.front().time;
	int finest = 0;
	std::vector<mesh::BlockLocation> leaves;
	for (auto const& block : snapshots) {
		leaves.push_back(block.location);
		if (block.location.level < 0 || block.location.level > 16 || block.layout.cellsPerActiveDimension() != c.mesh.cells || block.layout.ghostWidth() != 0 ||
			!block.gravityEnabled || block.gravity.values().size() != block.layout.interiorCellCount())
			throw std::invalid_argument("Invalid direct gravity snapshot geometry or fields");
		auto const width = (c.mesh.upper - c.mesh.lower) / Real(c.mesh.cells * (1 << block.location.level));
		if (block.cellWidth != width) throw std::invalid_argument("Inconsistent direct gravity cell width");
		finest = std::max(finest, block.location.level);
	}
	// Reuse mesh validation for complete coverage, bounds, and duplicate or
	// overlapping leaves before treating any cell as an independent mass.
	CartesianTopology const geometry(c, 1, leaves);
	result.cellsPerAxis = c.mesh.cells * (1 << finest);
	result.multipoleOrder = c.gravity.multipoleOrder;
	result.openingAngle = c.gravity.openingAngle;
	result.randomSeed = c.randomSeed;
	int const n = 2 * result.cellsPerAxis;
	auto const h = (c.mesh.upper - c.mesh.lower) / Real(n);
	class Cell {
	public:
		gravity::diagonal::Offset center{};
		units::Mass mass{};
		gravity::State const* field = nullptr;
	};
	std::vector<Cell> cells;
	for (auto const& block : snapshots) {
		if (!units::finite(block.time) || block.time != result.time) throw std::invalid_argument("Direct gravity reference requires synchronized snapshots");
		auto const volume = block.layout.cellMeasure(block.cellWidth);
		block.layout.forEachInterior([&](mesh::Coordinates const& coordinate, std::size_t i) {
			Cell cell;
			for (int d = 0; d < 3; ++d)
				cell.center[d] = (2 * (block.location.coordinates[d] * c.mesh.cells + coordinate[d]) + 1) * (1 << (finest - block.location.level));
			auto const density = block.hydroEnabled ? block.hydro.values()[i].density() : block.density.values()[i];
			cell.mass = density * volume;
			cell.field = &block.gravity.values()[i];
			if (!(cell.mass >= units::Mass{}) || !units::finite(cell.mass)) throw std::invalid_argument("Invalid direct gravity mass");
			cells.push_back(cell);
		});
	}
	// Physical x-fast ordering is independent of Morton block placement and
	// locality count. No covered coarse cell enters the reference sum.
	std::sort(cells.begin(), cells.end(),
		[](auto const& a, auto const& b) { return std::tie(a.center[2], a.center[1], a.center[0]) < std::tie(b.center[2], b.center[1], b.center[0]); });
	result.totalCells = cells.size();
	std::vector<units::Mass> masses;
	std::vector<gravity::State const*> numerical;
	for (auto const& cell : cells) {
		masses.push_back(cell.mass);
		numerical.push_back(cell.field);
	}
	std::vector<std::size_t> sources;
	for (std::size_t i = 0; i < masses.size(); ++i) {
		if (!numerical[i]) throw std::invalid_argument("Missing direct gravity cell");
		if (masses[i] != units::Mass{}) sources.push_back(i);
	}
	result.sourceCells = sources.size();
	std::size_t count = result.totalCells;
	if (c.verification.directSamples > 0)
		count = std::min(count, static_cast<std::size_t>(c.verification.directSamples));
	else if (sources.size() > 0 && count > std::size_t(c.verification.directMaxPairs) / sources.size())
		count = std::min({count, std::size_t(1024), std::max(std::size_t(1), std::size_t(c.verification.directMaxPairs) / sources.size())});
	result.targetIndices = directTargets(result.totalCells, count, c.randomSeed);
	result.cells = count;
	if (count < result.totalCells) result.sampling = "random-cell-center-without-replacement";
	std::array<Field, 4> fields{Field{"potential", "cm^2/s^2", {}, {}}, Field{"accelerationX", "cm/s^2", {}, {}}, Field{"accelerationY", "cm/s^2", {}, {}},
		Field{"accelerationZ", "cm/s^2", {}, {}}};
	auto coordinates = [&](std::size_t id) { return cells.at(id).center; };
	for (auto const target : result.targetIndices) {
		auto const x = coordinates(target);
		std::array<long double, 4> sum{};
		for (auto const source : sources) {
			auto const y = coordinates(source);
			if (images.active()) {
				gravity::diagonal::Coefficients local(5, 0);
				for (auto const& image : images.images()) {
					auto const r = images.centerSeparation(x, y, n, image);
					if (r != gravity::diagonal::Offset{}) gravity::diagonal::addDirect(local, units::value(masses[source]), r, units::value(h));
					if (images.periodic()) gravity::ewald::addDirect(local, units::value(masses[source]), r, images.periods(n), units::value(h));
				}
				sum[0] += units::value(constants::G) * local[0];
				for (int d = 0; d < 3; ++d) {
					std::array<int, 3> e{};
					e[d] = 1;
					sum[d + 1] -= units::value(constants::G) * gravity::diagonal::derivative(local, e[0], e[1], e[2]) / units::value(h);
				}
				continue;
			}
			if (target == source) continue;
			std::array<units::Length, 3> const r{h * Real(x[0] - y[0]), h * Real(x[1] - y[1]), h * Real(x[2] - y[2])};
			auto const distance = units::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
			auto const potential = -constants::G * masses[source] / distance;
			sum[0] += units::value(potential);
			for (int d = 0; d < 3; ++d)
				sum[d + 1] += units::value(potential * r[d] / (distance * distance));
		}
		fields[0].numerical.push_back(units::value(numerical[target]->potential()));
		for (int d = 0; d < 3; ++d)
			fields[d + 1].numerical.push_back(units::value(numerical[target]->acceleration(d)));
		for (int f = 0; f < 4; ++f)
			fields[f].exact.push_back(static_cast<Real>(sum[f]));
	}
	for (auto const& field : fields)
		result.fields.push_back(directErrorNorm(field, result.totalCells));
	return result;
}

}	 // namespace octotigerII::verification
