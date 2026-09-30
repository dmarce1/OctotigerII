// Distributed under the Boost Software License, Version 1.0.
#pragma once

#include <cmath>
#include "octotigerII/math/Real.hpp"

namespace octotigerII {

// Neumaier summation for finite Real contributions. The correction retains
// roundoff lost by successive additions; all arithmetic remains in Real.
// This improves accuracy but does not make arbitrary summation orders identical.
class CompensatedSum {
public:

	void add(Real value) {
		Real const next = sum_ + value;
		correction_ += std::abs(sum_) >= std::abs(value) ? (sum_ - next) + value : (value - next) + sum_;
		sum_ = next;
	}

	Real value() const { return sum_ + correction_; }

private:
	Real sum_ = 0;
	Real correction_ = 0;
};

} // namespace octotigerII
