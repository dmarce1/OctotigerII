/** @file
 * @brief Standard problem-provided exact solutions and volume-weighted errors.
 */
#pragma once
#include <functional>
#include <iosfwd>
#include <limits>
#include <string>
#include <vector>
#include "octotigerII/subgrid/subgrid.hpp"


namespace octotigerII::verification {


/// Physical CGS values returned by a problem's independent reference solution.
class ExactState {
public:
	hydro::PrimitiveState hydro;
	gravity::State gravity;
	radiation::RadiationSystem::State radiation;
};


/// An empty evaluator explicitly means that no analytic solution is available.
/// References are sampled at cell centers, at the actual snapshot time.
class Reference {
public:
	std::string name;
	std::string reason = "This problem has no implemented analytic solution";
	units::Time validUntil = units::Time::from_value(std::numeric_limits<Real>::infinity());
	std::function<ExactState(mesh::PhysicalCoordinates const&, units::Time)> evaluate;
};


/// CGS scalar arrays at the reporting/Silo boundary; physical solvers stay typed.
class Field {
public:
	std::string name, units;
	std::vector<Real> numerical, exact;
};


class ErrorNorm {
public:
	std::string name, units;
	/// Absolute CGS norms, retained for gates and for zero-reference fields.
	Real l1 = 0, l2 = 0, linf = 0;
	Real referenceL1 = 0, referenceL2 = 0, referenceLinf = 0;
	Real relativeL1HalfWidth95 = 0, relativeL2HalfWidth95 = 0;
	bool samplingUncertaintyAvailable = false, samplingWarning = false;
};


class Comparison {
public:
	std::string name, status, reason;
	units::Time time{};
	std::size_t cells = 0, totalCells = 0, sourceCells = 0;
	std::string referenceKind = "continuum", sampling = "cell-center";
	std::vector<std::size_t> targetIndices;
	std::int64_t randomSeed = 5489;
	int cellsPerAxis = 0, multipoleOrder = 0;
	Real openingAngle = 0;
	std::vector<ErrorNorm> fields;

	/// Fail an explicitly requested accuracy gate. Zero-reference fields use
	/// absoluteTolerance in that field's CGS units; relative norms are undefined.
	void enforce(Config const& config) const;

	void print(std::ostream& out) const;

	void writeJson(std::ostream& out) const;
};


/// Whether gravity uses the discrete snapshot reference (gravity fields only).
inline bool usesDirectGravity(Config const& config) {
	return build::gravity && config.verification.gravityReference == "direct" && config.verification.analytic != "off";
}

/// Apply auto/on/off policy and the reference's time-validity restriction.
Reference reference(Config const& config, units::Time time);

std::vector<Field> sample(Snapshot const& snapshot, Config const& config, Reference const& reference);

Comparison compare(std::vector<Snapshot> const& snapshots, Config const& config);

/// Independent exact solutions used by the corresponding problem hooks.
Reference sodReference(Config const& config);

Reference sphereReference(Config const& config, bool gaussian);

Reference streamingReference(Config const& config);

}	 // namespace octotigerII::verification
