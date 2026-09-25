/** @file
 * @brief Harmonic-equivalent Cartesian expansions and diagonal M2L translations.
 * @ingroup numerics
 */
// Distributed under the Boost Software License, Version 1.0.
#pragma once

#include <array>
#include <complex>
#include <cstddef>
#include <memory>
#include <vector>


namespace octotigerII::gravity::diagonal {
using Vector = std::array<double, 3>;
using Offset = std::array<int, 3>;
using Coefficients = std::vector<double>;
using Field = std::vector<Coefficients>;


// Cartesian polynomial modulo kx^2+ky^2+kz^2. Store only kz exponents 0,1.
// These are harmonic-equivalent moments, NOT a subset of raw moments.
/// Exponent triple in the harmonic-equivalent coefficient basis.
/// Only z=0,1 are stored; higher powers are reduced with kz²=−kx²−ky².
/// @ingroup numerics
class Index {
public:

	int x, y, z;

	int degree() const {
		return x + y + z;
	}
};


/// Return (p+1)² harmonic-equivalent coefficients for a supported order.
int coefficientCount(int p);

/// Return exponent triples in total-degree order for the supported expansion order.
const std::vector<Index>& indices(int p);

/// Form the harmonic-equivalent moments of one normalized point source.
Coefficients point(double mass, const Vector& offset, int p);

// offset=(old center-new center)/new scale; ratio=old scale/new scale.
/// Translate source moments: offset=(old center−new center)/new scale;
/// ratio=old scale/new scale.
Coefficients shiftMultipole(const Coefficients&, const Vector& offset, double ratio, int p);

// Locals store derivatives of the potential with respect to normalized coordinates.
// offset=(new center-old center)/old scale; ratio=new scale/old scale.
/// Translate local derivatives: offset=(new center−old center)/old scale;
/// ratio=new scale/old scale.
Coefficients shiftLocal(const Coefficients&, const Vector& offset, double ratio, int p);

/// Recover a normalized-coordinate potential derivative using harmonic reduction.
double derivative(const Coefficients&, int x, int y, int z);

/// Return potential and three normalized-coordinate derivatives at an offset.
std::array<double, 4> evaluate(const Coefficients&, const Vector& offset, int p);


// Greengard--Rokhlin (1997), section 7: C_XL D C_MX. The plane-wave axis
// is rotated onto the separation. Quadrature is exact for degrees <=2p.
/// Cached diagonal plane-wave multipole-to-local translation.
/// Uses C_XL D C_MX from @ref ref_greengard1997 "Greengard and Rokhlin (1997)".
/// The quadrature axis follows the separation; coefficient storage is Cartesian.
/// @ingroup numerics
class Operator {
public:

	/// forceLocal adds one local degree while retaining source degree p. Its
	/// gradient has degree p on both sides of an interaction, preserving mutual force.
	Operator(int p, Offset separation, bool forceLocal = false);

	// R=(target center-source center)/cell width. Source moments use that width.
	// Outputs potential derivatives with respect to target normalized coordinates.
	/// Accumulate one source expansion into a target local expansion.
	/// separation=(target center−source center)/cellWidth; locals are normalized derivatives.
	void add(Coefficients& local, const Coefficients& moment, double cellWidth) const;

	/// Return the number of plane-wave quadrature samples in this operator.
	std::size_t sampleCount() const {
		return diagonal_.size();
	}

private:

	int sourceCount_, localCount_;
	std::vector<std::complex<double>> toWave_, fromWave_;
	std::vector<double> diagonal_;
};


/// Find or build an immutable cached translation for the order and integer separation.
std::shared_ptr<const Operator> getOperator(int p, Offset separation, bool forceLocal = false);

/// Return the number of cached order/separation translation operators.
std::size_t cachedOperatorCount();

// Exact P2P; self interaction is excluded by the caller.
/// Add the exact non-self point-source potential and derivatives.
/// The caller excludes zero separation; physical G is applied by the typed driver.
void addDirect(Coefficients& local, double sourceMass, Offset separation, double cellWidth);
}	 // namespace octotigerII::gravity::diagonal
