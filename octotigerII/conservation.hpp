/** @file
 * @brief Volume integrals and cumulative numerical transport at physical boundaries.
 */
#pragma once
#include <array>
#include "octotigerII/buildConfig.hpp"
#include "octotigerII/units/cgs.hpp"

namespace octotigerII {

/// All evolved conserved components, integrated over physical volume.
class ConservedTotals {
public:
	units::Mass mass{};
	units::Energy gasEnergy{}, radiationEnergy{}, potentialEnergy{};
	std::array<units::Momentum, ndim> momentum{};
	std::array<units::Quantity<3, 1, -3>, ndim> radiationFlux{};

	ConservedTotals& operator+=(ConservedTotals const& other) {
		mass += other.mass;
		gasEnergy += other.gasEnergy;
		potentialEnergy += other.potentialEnergy;
		radiationEnergy += other.radiationEnergy;
		for (int d = 0; d < ndim; ++d) {
			momentum[d] += other.momentum[d];
			radiationFlux[d] += other.radiationFlux[d];
		}
		return *this;
	}

	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & mass & gasEnergy & radiationEnergy & potentialEnergy & momentum & radiationFlux;
	}
};

/// Positive/negative parts of outward-oriented transport, per component.
/// Thus grid + outward - inward removes the physical-boundary contribution.
/// Momentum includes pressure traction at walls, not just advected momentum.
class BoundaryTransport {
public:
	ConservedTotals inward, outward;

	BoundaryTransport& operator+=(BoundaryTransport const& other) {
		inward += other.inward;
		outward += other.outward;
		return *this;
	}

	template <typename Archive>
	void serialize(Archive& archive, unsigned) {
		archive & inward & outward;
	}
};

} // namespace octotigerII
