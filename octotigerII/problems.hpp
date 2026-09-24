/** @file
 * @brief Runtime dispatch to the available problem implementations.
 * @ingroup runtime
 */
#pragma once
#include <algorithm>
#include <functional>
#include "octotigerII/config.hpp"

namespace octotigerII {
class Snapshot;
namespace verification {
	class Reference;
	class ExactState;
}

/// Prescribed state at a physical ghost-cell center and stage time, in CGS.
/// An empty function means this problem does not provide analytic boundary data.
using ProblemBoundary = std::function<verification::ExactState(mesh::PhysicalCoordinates const&, units::Time)>;

ProblemBoundary problemBoundary(Config const& config);


/// Every problem declares its exact solution or explicitly returns unavailable.
verification::Reference problemReference(Config const& config);


/// Set problem-local defaults before reading runtime inputs.
void problemDefaults(Config& config);

std::string problemHelp();

/// Validate restrictions particular to the selected problem.
void validateProblem(Config const& config);

/// Keep an unresolved feature visible while constructing the initial hierarchy.
/// Initializers receive the actual cell width through Snapshot::cellWidth.
/// A problem must check that its physical width is resolved on the final mesh.
inline units::Length initialFeatureWidth(units::Length physicalWidth, units::Length cellWidth) {
	return std::max(physicalWidth, cellWidth);
}

/// Fill initial conditions at this snapshot resolution. refinementProbe permits
/// temporary feature widening for startup tagging; final stored data use false.
void initializeProblem(Snapshot& snapshot, Config const& config, bool refinementProbe = false);
}	 // namespace octotigerII
