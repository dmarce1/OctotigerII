/** @file
 * @brief Reproducible target sampling and independent isolated P2P verification.
 */
#pragma once
#include "octotigerII/verification/analytic.hpp"


namespace octotigerII::verification {


/// Sorted, unbiased sample without replacement. Uses a specified bounded-integer
/// mapping so results do not depend on the standard library's distributions.
std::vector<std::size_t> directTargets(std::size_t population, std::size_t count, std::uint64_t seed);

/// Relative norms and approximate 95% sampling half-widths (delta method with
/// finite-population correction). Linf is only the observed sample maximum.
ErrorNorm directErrorNorm(Field const& field, std::size_t population);

Comparison compareDirectGravity(std::vector<Snapshot> const& snapshots, Config const& config);

}	 // namespace octotigerII::verification
