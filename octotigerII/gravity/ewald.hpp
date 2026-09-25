/** @file
 * Newtonian 1/r Ewald corrections for 1P, 2P and 3P lattices.
 * Distributed under the Boost Software License, Version 1.0.
 */
#pragma once
#include "octotigerII/gravity/diagonal/fmm.hpp"

namespace octotigerII::gravity::ewald {

/// Ordinary derivatives of the correction psi_periodic + 1/r. The nearest
/// image is at the origin; callers must wrap only the periodic coordinates.
/// Units are inverse lengths and their derivatives, with G and mass omitted.
class Derivatives {
public:
	double at(int x, int y, int z) const;
	std::vector<double> values;
	double laplacian = 0;
};

/// Periods=0 denotes an open direction. The optional alpha is in inverse
/// coordinate units; zero selects 2/min(period). quadrature controls the 1P
/// integral and central analytic term, and permits independent convergence tests.
Derivatives derivatives(diagonal::Vector r, diagonal::Vector periods, int degree, double alpha = 0, int quadrature = 64);

/// Dense harmonic M2L: (p+1)^2 source and local entries, hence O(p^4).
/// Carries the additional trace/background terms explicitly in 3P.
class Operator {
public:
	/// forceLocal retains one additional local degree for a mutual order-p force.
	Operator(int order, diagonal::Offset separation, diagonal::Offset periods, bool forceLocal = false);
	void add(diagonal::Coefficients& local, diagonal::Coefficients const& moment, double cellWidth) const;

private:
	int sourceCount_, localCount_;
	std::vector<double> matrix_;
	double laplace_;
};

std::shared_ptr<Operator const> getOperator(int order, diagonal::Offset separation, diagonal::Offset periods, bool forceLocal = false);
void addDirect(diagonal::Coefficients& local, double mass, diagonal::Offset separation, diagonal::Offset periods, double cellWidth);

}	 // namespace octotigerII::gravity::ewald
