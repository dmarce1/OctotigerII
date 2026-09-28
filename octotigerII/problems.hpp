/** @file
 * @brief Runtime dispatch to the available problem implementations.
 * @ingroup runtime
 */
#pragma once
#include <algorithm>
#include <functional>
#include <stdexcept>
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

/// Prescribed material and isotropic photon heating in physical grid coordinates,
/// matching the problem initializer. Photon power is physical erg/(cm^3 s): the
/// reduced-speed radiation equation receives (chat/c)*photonPower. It supplies
/// neither gas energy nor momentum directly. Opacity is true absorption in cm^2/g;
/// configured scattering is applied separately by the gray opacity law.
struct RadiationMaterial {
	units::Quantity<2, -1, 0> opacity{};
	units::Quantity<-1, 1, -3> photonPower{};
};
using ProblemRadiationMaterial = std::function<RadiationMaterial(mesh::PhysicalCoordinates const&, units::Time)>;

/// Construct once and share the immutable problem model across cell evaluations.
/// Problems without this optional hook return config.radiation.opacity and zero
/// heating. The analytic model evaluates its state-dependent coefficients
/// separately. Custom material hooks may contain transparent cells with heating.
ProblemRadiationMaterial problemRadiationMaterial(Config const& config);
bool problemHasRadiationMaterial(Config const& config);

inline RadiationMaterial checkedRadiationMaterial(ProblemRadiationMaterial const& material,
	mesh::PhysicalCoordinates const& position, units::Time time) {
	auto const value = material(position, time);
	if (!units::finite(value.opacity) || value.opacity < decltype(value.opacity){} ||
		!units::finite(value.photonPower) || value.photonPower < decltype(value.photonPower){})
		throw std::invalid_argument("Prescribed radiation opacity and photon heating must be finite and nonnegative");
	return value;
}


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
