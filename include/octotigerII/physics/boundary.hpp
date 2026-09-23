/** @file
 * @brief Per-face physical boundary rules shared by patch and distributed transport.
 */
#pragma once
#include <algorithm>
#include <array>
#include <functional>
#include <stdexcept>
#include <string>
#include "octotigerII/mesh.hpp"


namespace octotigerII::physics {


enum class BoundaryCondition
{
	Outflow,
	Periodic,
	Reflecting,
	Analytic
};


inline BoundaryCondition parseBoundaryCondition(std::string const& value) {
	if (value == "outflow") return BoundaryCondition::Outflow;
	if (value == "periodic") return BoundaryCondition::Periodic;
	if (value == "reflecting") return BoundaryCondition::Reflecting;
	if (value == "analytic") return BoundaryCondition::Analytic;
	throw std::invalid_argument("Expected periodic, reflecting, outflow, or analytic boundary: " + value);
}


/// Donor coordinates and transformations for a single ghost cell.
class BoundaryCell {
public:

	mesh::Coordinates source{};
	unsigned reflectionMask = 0;
	bool analytic = false;
};


/// One independent rule per active face, with paired periodic faces.
class BoundaryConditions {
public:

	std::array<BoundaryCondition, ndim> lower{};
	std::array<BoundaryCondition, ndim> upper{};

	static BoundaryConditions uniform(BoundaryCondition rule) {
		BoundaryConditions result;
		result.lower.fill(rule);
		result.upper.fill(rule);
		return result;
	}

	static BoundaryConditions periodic() {
		return uniform(BoundaryCondition::Periodic);
	}

	bool periodic(int axis) const {
		return lower.at(axis) == BoundaryCondition::Periodic && upper.at(axis) == BoundaryCondition::Periodic;
	}

	bool all(BoundaryCondition rule) const {
		return std::all_of(lower.begin(), lower.end(), [=](auto face) { return face == rule; }) &&
			std::all_of(upper.begin(), upper.end(), [=](auto face) { return face == rule; });
	}

	bool contains(BoundaryCondition rule) const {
		return std::find(lower.begin(), lower.end(), rule) != lower.end() || std::find(upper.begin(), upper.end(), rule) != upper.end();
	}

	void validate() const {
		for (int axis = 0; axis < ndim; ++axis) {
			for (auto face : {lower[axis], upper[axis]}) {
				switch (face) {
				case BoundaryCondition::Outflow:
				case BoundaryCondition::Periodic:
				case BoundaryCondition::Reflecting:
				case BoundaryCondition::Analytic:
					break;
				default:
					throw std::invalid_argument("Invalid physical boundary enum");
				}
			}
			if ((lower[axis] == BoundaryCondition::Periodic) != (upper[axis] == BoundaryCondition::Periodic))
				throw std::invalid_argument(std::string("Periodic boundaries must be paired on the ") + "xyz"[axis] + " axis");
		}
	}

	/// Map a cell in a uniform domain. Call validate() once before mapping a batch.
	/// Analytic takes precedence at edges/corners and uses the original ghost position.
	BoundaryCell map(mesh::Coordinates cell, int count) const {
		if (count < 1) throw std::invalid_argument("Boundary mapping requires a positive cell count");
		BoundaryCell result{cell};
		for (int axis = 0; axis < ndim; ++axis) {
			int& x = result.source[axis];
			while (x < 0 || x >= count) {
				bool const low = x < 0;
				switch (low ? lower[axis] : upper[axis]) {
				case BoundaryCondition::Periodic:
					x = (x % count + count) % count;
					break;
				case BoundaryCondition::Outflow:
					x = std::clamp(x, 0, count - 1);
					break;
				case BoundaryCondition::Reflecting:
					x = low ? -x - 1 : 2 * count - x - 1;
					result.reflectionMask ^= 1u << axis;
					break;
				case BoundaryCondition::Analytic:
					result.analytic = true;
					x = std::clamp(x, 0, count - 1);
					break;
				default:
					throw std::invalid_argument("Invalid physical boundary enum");
				}
			}
		}
		return result;
	}

	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		for (int axis = 0; axis < ndim; ++axis)
			archive & lower[axis] & upper[axis];
	}
};


template <typename State>
using AnalyticBoundary = std::function<State(mesh::PhysicalCoordinates const&, units::Time)>;


/// Apply all crossed reflecting planes; scalars and tangential vectors stay unchanged.
template <typename System>
typename System::State reflectBoundary(typename System::State state, unsigned mask, System const& system) {
	for (int axis = 0; axis < ndim; ++axis)
		if (mask & (1u << axis)) state = system.reflected(state, axis);
	return state;
}


/// Sample prescribed data, rejecting missing functions and invalid physical states.
template <typename System>
typename System::State evaluateBoundary(AnalyticBoundary<typename System::State> const& analytic,
	mesh::PhysicalCoordinates const& position, units::Time time, System const& system) {
	if (!analytic) throw std::invalid_argument("Analytic boundary requires a problem evaluator");
	if (!units::finite(time) || time < units::Time{}) throw std::invalid_argument("Invalid analytic boundary time");
	auto state = analytic(position, time);
	if (!system.admissible(state)) throw std::runtime_error("Analytic boundary produced an inadmissible state");
	return state;
}


} // namespace octotigerII::physics
