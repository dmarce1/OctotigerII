#pragma once
#include <complex>
#include <vector>
#include "octotigerII/units/cgs.hpp"

namespace octotigerII::gravity {
/// Reusable isolated 1/r convolution on a uniform cube, with G=1 and the
/// cell's self potential omitted, exactly as in gravity::solve. Zero padding
/// prevents periodic images. Intended for repeated, potential-only SCF solves.
class IsolatedPotential {
public:
	IsolatedPotential(int cells, Real cellWidth);
	std::vector<Real> operator()(std::vector<Real> const& density) const;
private:
	int n_, padded_;
	std::vector<std::complex<Real>> kernel_;
};
} // namespace octotigerII::gravity
