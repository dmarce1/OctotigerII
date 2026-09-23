#include "octotigerII/verification/directGravity.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
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
	result.cellsPerAxis = c.mesh.cells * (1 << c.mesh.level);
	result.multipoleOrder = c.gravity.multipoleOrder;
	result.openingAngle = c.gravity.openingAngle;
	result.randomSeed = c.randomSeed;
	int const n = result.cellsPerAxis;
	result.totalCells = std::size_t(n) * n * n;
	auto const h = (c.mesh.upper - c.mesh.lower) / Real(n);
	auto const volume = h * h * h;
	// Canonical x-fast global indices make both the draw and source summation
	// independent of block order, task completion order and HPX locality count.
	std::vector<units::Mass> masses(result.totalCells);
	std::vector<gravity::State const*> numerical(result.totalCells, nullptr);
	for (auto const& block : snapshots) {
		if (!units::finite(block.time) || block.time != result.time || block.cellWidth != h || block.location.level != c.mesh.level)
			throw std::invalid_argument("Direct gravity reference requires synchronized uniform snapshots");
		block.layout.forEachInterior([&](mesh::Coordinates cell, std::size_t i) {
			for (int d = 0; d < ndim; ++d) {
				cell[d] += block.location.coordinates[d] * c.mesh.cells;
				if (cell[d] < 0 || cell[d] >= n) throw std::invalid_argument("Invalid direct gravity cell coordinate");
			}
			auto const id = mesh::linearIndex(cell, mesh::filledCoordinates(n));
			if (numerical[id]) throw std::invalid_argument("Duplicate direct gravity cell");
			numerical[id] = &block.gravity.values()[i];
			auto const density = block.hydroEnabled ? block.hydro.values()[i].density() : block.density.values()[i];
			masses[id] = density * volume;
			if (!(masses[id] >= units::Mass{}) || !units::finite(masses[id])) throw std::invalid_argument("Invalid direct gravity mass");
		});
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
	auto coordinates = [n](std::size_t id) { return std::array<int, 3>{int(id % n), int(id / n % n), int(id / (std::size_t(n) * n))}; };
	for (auto const target : result.targetIndices) {
		auto const x = coordinates(target);
		std::array<long double, 4> sum{};
		for (auto const source : sources) {
			auto const y = coordinates(source);
			if (images.active()) {
				gravity::diagonal::Coefficients local(5, 0);
				for (auto const& image : images.images()) {
					auto const r = images.separation(x, y, n, image);
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
