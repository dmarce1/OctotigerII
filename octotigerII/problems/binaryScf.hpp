#pragma once
#include <array>
#include <memory>
#include <iosfwd>
#include <vector>
#include "octotigerII/config.hpp"
#include "octotigerII/hydro/hydroSystem.hpp"

namespace octotigerII::problems {

/// Synchronous, isolated binary constructed on a uniform reference mesh using
/// OctoII's cell-mass gravity. Hydro stores inertial velocities in every frame.
class BinaryScf {
public:
	struct Diagnostics {
		int iterations = 0;
		Real densityResidual = 0, bernoulliResidual = 0, virialResidual = 0;
		Real maxStarBernoulliResidual = 0;
		std::array<Real, 2> mass{}, coreMass{}, center{}, bernoulli{}, centralDensity{};
		Real omega = 0, separation = 0, rotationCenter = 0, l1 = 0, l1Potential = 0;
		Real kinetic = 0, potential = 0, pressureIntegral = 0, orbitalAngularMomentum = 0, spinAngularMomentum = 0;
	};
	struct Cell {
		Real density = 0, pressure = 0;
		std::array<Real, 4> partial{}; // primary core/envelope, donor core/envelope
	};
	explicit BinaryScf(Config const&);
	static std::shared_ptr<BinaryScf const> get(Config const&);
	static void validate(Config const&);
	hydro::PrimitiveState operator()(mesh::PhysicalCoordinates const&) const;
	hydro::PrimitiveState state(Cell const&, mesh::PhysicalCoordinates const&) const;
	Cell sample(mesh::PhysicalCoordinates const&) const;
	/// Conservative overlap of piecewise constant reference cells with a hydro cell.
	Cell average(mesh::PhysicalCoordinates const& center, units::Length width) const;
	Diagnostics const& diagnostics() const { return diagnostics_; }
	units::InverseTime angularVelocity() const;
	Real cellWidth() const { return width_; } // reference, dimensionless
	std::vector<Cell> const& cells() const { return cells_; }
	int resolution() const { return n_; }
	units::Length lengthUnit() const { return length_; }
	units::Density densityUnit() const { return density_; }
	void writeJson(std::ostream&) const;
private:
	Config::BinaryScfOptions parameters_;
	int n_;
	Real width_, lower_, atmosphere_;
	units::Length length_;
	units::Density density_;
	units::Pressure pressure_;
	Diagnostics diagnostics_;
	std::vector<Cell> cells_;
};
} // namespace octotigerII::problems
