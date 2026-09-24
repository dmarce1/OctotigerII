#pragma once
#include <vector>
#include "octotigerII/config.hpp"

namespace octotigerII::problems {

/// Dimensionless regular Lane-Emden solution, integrated to its first zero.
class LaneEmden {
public:
	class Value {
	public:
		Real theta = 1, mass = 0; // mass = -xi^2 dtheta/dxi
	};
	explicit LaneEmden(Real index);
	Value operator()(Real xi) const;
	Real surface() const { return points_.back().xi; }
	Real surfaceMass() const { return points_.back().value.mass; }
	Real index() const { return index_; }
private:
	class Point {
	public:
		Real xi;
		Value value;
	};
	Real index_;
	std::vector<Point> points_;
};

/// Physical scaling of a vacuum polytrope by surface radius and central density.
class Polytrope {
public:
	class State {
	public:
		units::Density density{};
		units::Pressure pressure{};
		units::Mass enclosedMass{};
		units::VelocitySquared potential{};
	};
	Polytrope(Real index, units::Length radius, units::Density centralDensity);
	State operator()(units::Length radius) const;
	units::Mass mass() const;
	units::Pressure centralPressure() const;
	units::Length scaleLength() const { return scale_; }
private:
	LaneEmden const* solution_;
	units::Length radius_, scale_;
	units::Density centralDensity_;
};

} // namespace octotigerII::problems
