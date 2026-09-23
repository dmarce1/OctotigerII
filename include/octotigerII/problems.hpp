/** @file
 * @brief Link-time problem hooks implemented by one selected bin directory.
 * @ingroup runtime
 */
#pragma once
#include "octotigerII/config.hpp"

namespace octotigerII {
class Snapshot;
namespace verification {
	class Reference;
}

/// Every problem declares its exact solution or explicitly returns unavailable.
verification::Reference problemReference(Config const& config);


/// Set problem-local defaults before reading runtime inputs.
void problemDefaults(Config& config);

/// Validate restrictions particular to the selected problem.
void validateProblem(Config const& config);

/// Fill physical initial conditions on one interior snapshot.
void initializeProblem(Snapshot& snapshot, Config const& config);
}	 // namespace octotigerII
