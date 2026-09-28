/** @file
 * @brief Rigidly rotating Cartesian mesh geometry with inertial stored vectors.
 */
#pragma once

#include "octotigerII/mesh.hpp"
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace octotigerII::finiteVolume {

/// Logical coordinates rotate about the inertial z axis at constant angular rate.
/// All transformations are orthogonal changes of components, never velocity boosts.
class RotatingFrame {
public:

	explicit RotatingFrame(units::InverseTime omega = {}) : omega_(omega) {
		if (!units::finite(omega_)) throw std::invalid_argument("Grid angular velocity must be finite");
		if constexpr (ndim < 2)
			if (active()) throw std::invalid_argument("Grid rotation requires at least two spatial dimensions");
	}

	bool active() const { return omega_ != units::InverseTime{}; }
	units::InverseTime omega() const { return omega_; }

	template <typename Q>
	std::array<Q, ndim> toInertial(std::array<Q, ndim> value, units::Time time) const {
		return rotate(value, time);
	}

	template <typename Q>
	std::array<Q, ndim> toGrid(std::array<Q, ndim> value, units::Time time) const {
		return rotate(value, -time);
	}

	/// State/flux schemas all place their Cartesian vector in slots 1 through ndim.
	template <typename State>
	State toInertialState(State value, units::Time time) const { return rotateState(value, time); }

	template <typename State>
	State toGridState(State value, units::Time time) const { return rotateState(value, -time); }

	std::array<units::Velocity, ndim> velocityGrid(mesh::PhysicalCoordinates const& position) const {
		std::array<units::Velocity, ndim> result{};
		if constexpr (ndim >= 2) {
			result[0] = -omega_ * position[1];
			result[1] = omega_ * position[0];
		}
		return result;
	}

	units::Velocity normalSpeed(mesh::PhysicalCoordinates const& position, int axis) const {
		return velocityGrid(position).at(axis);
	}

	/// Resolve changing face normals even when material corotates with the mesh.
	units::Time maximumTimestep() const {
		return active() ? Real(0.1) / units::abs(omega_) : units::Time::from_value(std::numeric_limits<Real>::infinity());
	}

private:
	units::InverseTime omega_;

	template <typename Q>
	std::array<Q, ndim> rotate(std::array<Q, ndim> value, units::Time time) const {
		if (!active()) return value;
		if constexpr (ndim >= 2) {
			Real const angle = units::value(omega_ * time);
			Real const cosine = std::cos(angle), sine = std::sin(angle);
			auto const x = value[0], y = value[1];
			value[0] = cosine * x - sine * y;
			value[1] = sine * x + cosine * y;
		}
		return value;
	}

	template <typename State>
	State rotateState(State value, units::Time time) const {
		if (!active()) return value;
		using Q = std::remove_cvref_t<decltype(value.momentum(0))>;
		std::array<Q, ndim> vector{};
		for (int axis = 0; axis < ndim; ++axis) vector[axis] = value.momentum(axis);
		vector = rotate(vector, time);
		for (int axis = 0; axis < ndim; ++axis) value.momentum(axis) = vector[axis];
		return value;
	}
};

} // namespace octotigerII::finiteVolume
