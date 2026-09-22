// Distributed under the Boost Software License, Version 1.0.
#include "octotigerII/gravity/solver.hpp"
#include "octotigerII/gravity/diagonal/fmm.hpp"
#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>

namespace octotigerII::gravity {
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
struct Level {
	int count;
	Real width;
	std::vector<Coefficients> moments, locals;
	Level(int n, Real h, int coefficients)
		: count(n), width(h),
		  moments(static_cast<std::size_t>(n) * n * n, Coefficients(coefficients)),
		  locals(moments.size(), Coefficients(coefficients)) {}
};
} // namespace

Solution solve(std::vector<Real> const& density, int n, Real h, int order, Real theta) {
	if (n < 2 || (n & (n - 1)) || density.size() != static_cast<std::size_t>(n) * n * n ||
		!(h > 0) || !std::isfinite(h) || order < 3 || order > 5 ||
		!(theta > 0 && theta < 1 / std::sqrt(3.0)))
		throw std::invalid_argument(
			"Gravity requires a power-of-two cubic mesh, p=3..5, 0<theta<1/sqrt(3)");
	std::vector<Level> levels;
	for (int count = 1; count <= n; count *= 2)
		levels.emplace_back(count, h * (n / count), diagonal::coefficientCount(order));
	auto& leaves = levels.back();
	Real const volume = h * h * h;
	for (std::size_t i = 0; i < density.size(); ++i) {
		if (!(density[i] >= 0) || !std::isfinite(density[i]) || !std::isfinite(density[i] * volume))
			throw std::invalid_argument("Gravity density/mass must be finite and nonnegative");
		leaves.moments[i][0] = density[i] * volume;
	}
	// Upward pass: compact moments about geometric cell centers.
	for (int level = static_cast<int>(levels.size()) - 2; level >= 0; --level) {
		auto& parent = levels[level];
		auto const& fine = levels[level + 1];
		for (int z = 0; z < parent.count; ++z)
			for (int y = 0; y < parent.count; ++y)
				for (int x = 0; x < parent.count; ++x) {
					Coordinate const c{x, y, z};
					for (int slot = 0; slot < 8; ++slot)
						add(parent.moments[index(c, parent.count)],
							diagonal::shiftMultipole(
								fine.moments[index(child(c, slot), fine.count)], childOffset(slot),
								0.5, order));
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
		Real const distance = std::hypot(Real(r[0]), Real(r[1]), Real(r[2]));
		bool const leaf = depth + 1 == static_cast<int>(levels.size());
		if (!same && !leaf && distance * theta > 1) {
			diagonal::getOperator(order, r)->add(level.locals[ia], level.moments[ib], level.width);
			for (auto& value : r)
				value = -value;
			diagonal::getOperator(order, r)->add(level.locals[ib], level.moments[ia], level.width);
			++result.statistics.multipolePairs;
		} else if (leaf) {
			if (same)
				return;
			diagonal::addDirect(level.locals[ia], level.moments[ib][0], r, level.width);
			for (auto& value : r)
				value = -value;
			diagonal::addDirect(level.locals[ib], level.moments[ia][0], r, level.width);
			++result.statistics.directPairs;
		} else {
			for (int i = 0; i < 8; ++i)
				for (int j = (same ? i : 0); j < 8; ++j)
					interact(depth + 1, child(a, i), child(b, j));
		}
	};
	interact(0, {0, 0, 0}, {0, 0, 0});
	// Downward pass: normalized derivatives scale with the level width.
	for (std::size_t depth = 0; depth + 1 < levels.size(); ++depth) {
		auto const& parent = levels[depth];
		auto& fine = levels[depth + 1];
		for (int z = 0; z < parent.count; ++z)
			for (int y = 0; y < parent.count; ++y)
				for (int x = 0; x < parent.count; ++x) {
					Coordinate const c{x, y, z};
					for (int slot = 0; slot < 8; ++slot)
						add(fine.locals[index(child(c, slot), fine.count)],
							diagonal::shiftLocal(parent.locals[index(c, parent.count)],
												 childOffset(slot), 0.5, order));
				}
	}
	result.fields.resize(density.size());
	for (std::size_t i = 0; i < density.size(); ++i) {
		auto& field = result.fields[i];
		auto const& local = leaves.locals[i];
		field.potential() = gravitationalConstant * local[0];
		field.acceleration(0) = -gravitationalConstant * diagonal::derivative(local, 1, 0, 0) / h;
		field.acceleration(1) = -gravitationalConstant * diagonal::derivative(local, 0, 1, 0) / h;
		field.acceleration(2) = -gravitationalConstant * diagonal::derivative(local, 0, 0, 1) / h;
		for (int f = 0; f < 4; ++f)
			if (!std::isfinite(field[f]))
				throw std::runtime_error("Nonfinite gravity solution");
	}
	return result;
}
} // namespace octotigerII::gravity
