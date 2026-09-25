// Distributed under the Boost Software License, Version 1.0.
#include "octotigerII/gravity/solver.hpp"
#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>
#include "octotigerII/gravity/diagonal/fmm.hpp"
#include "octotigerII/gravity/ewald.hpp"
#include "octotigerII/gravity/images.hpp"
#include "octotigerII/profiling.hpp"

namespace octotigerII::gravity {
static_assert(ndim == 3, "The gravity solver requires a 3D build");

namespace {
	using diagonal::Coefficients;
	using Coordinate = std::array<int, 3>;

	std::size_t index(Coordinate c, int n) {
		return (static_cast<std::size_t>(c[2]) * n + c[1]) * n + c[0];
	}

	Coordinate child(Coordinate c, int slot) {
		return {2 * c[0] + (slot & 1), 2 * c[1] + ((slot >> 1) & 1), 2 * c[2] + ((slot >> 2) & 1)};
	}

	diagonal::Vector childOffset(int slot) {
		return {(slot & 1) ? 0.25 : -0.25, (slot & 2) ? 0.25 : -0.25, (slot & 4) ? 0.25 : -0.25};
	}

	void add(Coefficients& destination, Coefficients const& source) {
		for (std::size_t i = 0; i < destination.size(); ++i)
			destination[i] += source[i];
	}

	class Level {
	public:
		int count;
		Real width;
		std::vector<Coefficients> moments, locals, forceLocals;

		Level(int n, Real h, int coefficients, int forceCoefficients)
		  : count(n)
		  , width(h)
		  , moments(static_cast<std::size_t>(n) * n * n, Coefficients(coefficients))
		  , locals(moments.size(), Coefficients(coefficients))
		  , forceLocals(moments.size(), Coefficients(forceCoefficients)) {}
	};

}	 // namespace

Solution solve(
	std::vector<units::Density> const& density, int n, units::Length cellWidth, int order, Real theta, physics::BoundaryConditions const& boundaries) {
	using std::hypot;
	using std::isfinite;
	using std::sqrt;

	profiling::Region profile("gravity.serial.solve");

	// The diagonal expansion uses dimensionless coefficients, normalized to
	// one gram and one centimeter. Restore dimensions before applying G.
	auto const cm = units::Length::from_value(1);
	auto const gram = units::Mass::from_value(1);
	Real const h = cellWidth / cm;
	if (n < 2 || (n & (n - 1)) || density.size() != static_cast<std::size_t>(n) * n * n || !(h > 0) || !isfinite(h) || order < 1 || order > 10 ||
		!(theta > 0 && theta < 1 / sqrt(3.0)))
		throw std::invalid_argument("Gravity requires a power-of-two cubic mesh, p=1..10, 0<theta<1/sqrt(3)");
	ImageGeometry const images(boundaries);
	std::vector<Level> levels;
	for (int count = 1; count <= n; count *= 2)
		levels.emplace_back(count, h * (n / count), diagonal::coefficientCount(order) + (images.active() ? 1 : 0),
			diagonal::coefficientCount(order + 1) + (images.active() ? 1 : 0));
	auto& leaves = levels.back();
	auto const volume = cellWidth * cellWidth * cellWidth;
	for (std::size_t i = 0; i < density.size(); ++i) {
		if (!(density[i] >= units::Density{}) || !units::finite(density[i]) || !units::finite(density[i] * volume))
			throw std::invalid_argument("Gravity density/mass must be finite and nonnegative");
		leaves.moments[i][0] = density[i] * volume / gram;
	}
	// Upward pass: compact moments about geometric cell centers.
	{
		profiling::Region profile("gravity.serial.m2m");
		for (int level = static_cast<int>(levels.size()) - 2; level >= 0; --level) {
			auto& parent = levels[level];
			auto const& fine = levels[level + 1];
			for (int z = 0; z < parent.count; ++z)
				for (int y = 0; y < parent.count; ++y)
					for (int x = 0; x < parent.count; ++x) {
						Coordinate const c{x, y, z};
						for (int slot = 0; slot < 8; ++slot)
							add(parent.moments[index(c, parent.count)],
								(images.active() ? imageShiftMultipole : diagonal::shiftMultipole)(
									fine.moments[index(child(c, slot), fine.count)], childOffset(slot), 0.5, order));
					}
		}
	}
	Solution result;
	// Symmetric cell-pair traversal. Both cells always have the same width.
	// A pair is processed at its first accepted level; descendants then see
	// only the translated local, so there is no double counting.
	std::function<void(int, Coordinate, Coordinate)> interact;
	interact = [&](int depth, Coordinate a, Coordinate b) {
		auto& level = levels[depth];
		auto const ia = index(a, level.count), ib = index(b, level.count);
		bool const same = a == b;
		diagonal::Offset r{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
		Real const distance = hypot(Real(r[0]), Real(r[1]), Real(r[2]));
		bool const leaf = depth + 1 == static_cast<int>(levels.size());
		if (!same && !leaf && distance * theta > 1) {
			diagonal::getOperator(order, r)->add(level.locals[ia], level.moments[ib], level.width);
			diagonal::getOperator(order, r, true)->add(level.forceLocals[ia], level.moments[ib], level.width);
			for (auto& value : r)
				value = -value;
			diagonal::getOperator(order, r)->add(level.locals[ib], level.moments[ia], level.width);
			diagonal::getOperator(order, r, true)->add(level.forceLocals[ib], level.moments[ia], level.width);
			++result.statistics.multipolePairs;
		} else if (leaf) {
			if (same) return;
			diagonal::addDirect(level.locals[ia], level.moments[ib][0], r, level.width);
			diagonal::addDirect(level.forceLocals[ia], level.moments[ib][0], r, level.width);
			for (auto& value : r)
				value = -value;
			diagonal::addDirect(level.locals[ib], level.moments[ia][0], r, level.width);
			diagonal::addDirect(level.forceLocals[ib], level.moments[ia][0], r, level.width);
			++result.statistics.directPairs;
		} else {
			for (int i = 0; i < 8; ++i)
				for (int j = (same ? i : 0); j < 8; ++j)
					interact(depth + 1, child(a, i), child(b, j));
		}
	};
	{
		profiling::Region interactions("gravity.serial.interactions");
		if (!images.active())
			interact(0, {0, 0, 0}, {0, 0, 0});
		else
			for (unsigned depth = 0; depth < levels.size(); ++depth) {
				auto& level = levels[depth];
				bool const leaf = depth + 1 == levels.size();
				for (int z = 0; z < level.count; ++z)
					for (int y = 0; y < level.count; ++y)
						for (int x = 0; x < level.count; ++x) {
							Coordinate const a{x, y, z};
							auto& local = level.locals[index(a, level.count)];
							auto& forceLocal = level.forceLocals[index(a, level.count)];
							images.interactions(depth, a, leaf, theta, [&](Coordinate b, diagonal::Offset r, unsigned mask, bool correction) {
								auto const& source = level.moments[index(b, level.count)];
								if (source[0] == 0) return;
								if (leaf) {
									if (correction) {
										ewald::addDirect(local, source[0], r, images.periods(level.count), level.width);
										ewald::addDirect(forceLocal, source[0], r, images.periods(level.count), level.width);
									} else {
										diagonal::addDirect(local, source[0], r, level.width);
										diagonal::addDirect(forceLocal, source[0], r, level.width);
									}
									++result.statistics.directPairs;
								} else {
									auto const moment = reflectMultipole(source, mask, order);
									if (correction) {
										ewald::getOperator(order, r, images.periods(level.count))->add(local, moment, level.width);
										ewald::getOperator(order, r, images.periods(level.count), true)->add(forceLocal, moment, level.width);
									} else {
										diagonal::getOperator(order, r)->add(local, moment, level.width);
										diagonal::getOperator(order, r, true)->add(forceLocal, moment, level.width);
									}
									++result.statistics.multipolePairs;
								}
								if (correction) ++result.statistics.ewaldPairs;
								if (mask) ++result.statistics.reflectedPairs;
							});
						}
			}
	}
	// Downward pass: normalized derivatives scale with the level width.
	{
		profiling::Region profile("gravity.serial.l2l");
		for (std::size_t depth = 0; depth + 1 < levels.size(); ++depth) {
			auto const& parent = levels[depth];
			auto& fine = levels[depth + 1];
			for (int z = 0; z < parent.count; ++z)
				for (int y = 0; y < parent.count; ++y)
					for (int x = 0; x < parent.count; ++x) {
						Coordinate const c{x, y, z};
						for (int slot = 0; slot < 8; ++slot) {
							add(fine.locals[index(child(c, slot), fine.count)],
								(images.active() ? imageShiftLocal : diagonal::shiftLocal)(
									parent.locals[index(c, parent.count)], childOffset(slot), 0.5, order));
							add(fine.forceLocals[index(child(c, slot), fine.count)],
								(images.active() ? imageShiftLocal : diagonal::shiftLocal)(
									parent.forceLocals[index(c, parent.count)], childOffset(slot), 0.5, order + 1));
						}
					}
		}
	}
	profiling::Region publish("gravity.serial.l2p");
	result.fields.resize(density.size());
	for (std::size_t i = 0; i < density.size(); ++i) {
		auto& field = result.fields[i];
		auto const& local = leaves.locals[i];
		auto const& forceLocal = leaves.forceLocals[i];
		field.potential() = constants::G * (gram / cm) * local[0];
		field.acceleration(0) = -constants::G * (gram / cm) * diagonal::derivative(forceLocal, 1, 0, 0) / cellWidth;
		field.acceleration(1) = -constants::G * (gram / cm) * diagonal::derivative(forceLocal, 0, 1, 0) / cellWidth;
		field.acceleration(2) = -constants::G * (gram / cm) * diagonal::derivative(forceLocal, 0, 0, 1) / cellWidth;
		if (!finite(field)) throw std::runtime_error("Nonfinite gravity solution");
	}
	return result;
}
}	 // namespace octotigerII::gravity
