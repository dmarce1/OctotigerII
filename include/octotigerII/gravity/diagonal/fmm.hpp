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
struct Index {
	int x, y, z;
	int degree() const {
		return x + y + z;
	}
};
int coefficientCount(int p);
const std::vector<Index>& indices(int p);
Coefficients point(double mass, const Vector& offset, int p);
// offset=(old center-new center)/new scale; ratio=old scale/new scale.
Coefficients shiftMultipole(const Coefficients&, const Vector& offset, double ratio, int p);
// Locals store derivatives of the potential with respect to normalized coordinates.
// offset=(new center-old center)/old scale; ratio=new scale/old scale.
Coefficients shiftLocal(const Coefficients&, const Vector& offset, double ratio, int p);
double derivative(const Coefficients&, int x, int y, int z);
std::array<double, 4> evaluate(const Coefficients&, const Vector& offset, int p);

// Greengard--Rokhlin (1997), section 7: C_XL D C_MX. The plane-wave axis
// is rotated onto the separation. Quadrature is exact for degrees <=2p.
class Operator {
  public:
	Operator(int p, Offset separation);
	// R=(target center-source center)/cell width. Source moments use that width.
	// Outputs potential derivatives with respect to target normalized coordinates.
	void add(Coefficients& local, const Coefficients& moment, double cellWidth) const;
	std::size_t sampleCount() const {
		return diagonal_.size();
	}

  private:
	int count_;
	std::vector<std::complex<double>> toWave_, fromWave_;
	std::vector<double> diagonal_;
};
std::shared_ptr<const Operator> getOperator(int p, Offset separation);
std::size_t cachedOperatorCount();
// Exact P2P; self interaction is excluded by the caller.
void addDirect(Coefficients& local, double sourceMass, Offset separation, double cellWidth);
} // namespace octotigerII::gravity::diagonal
