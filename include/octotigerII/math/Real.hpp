/** @file
 * @brief Scalar arithmetic type and mathematical constants.
 * @ingroup math
 */
// Distributed under the Boost Software License, Version 1.0.
#pragma once
#include <limits>
#include <numbers>

using Real = double;
inline constexpr Real epsilonR = std::numeric_limits<Real>::epsilon();
inline constexpr Real piR = std::numbers::pi_v<Real>;
