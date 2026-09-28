/** @file
 * @brief M1 transport adapter with physical (E, F) stored quantities.
 * @ingroup numerics
 */
// Distributed under the Boost Software License, Version 1.0.
#pragma once
#include <array>
#include "octotigerII/mesh.hpp"
#include "octotigerII/finiteVolume/solver.hpp"
#include "octotigerII/radiation/m1.hpp"
#include "octotigerII/units/constants.hpp"


namespace octotigerII::radiation {


/// Transport adapter whose state stores E and the physical radiative flux F.
/// E has units erg/cm³; F has units erg/(cm² s). Conversion to F/c occurs only
/// at the calculation interface and always uses physical c. The reduced speed
/// ĉ changes propagation, not the stored flux units or the realizability bound.
/// @ingroup numerics
class RadiationSystem {
public:
	struct InsulatingBoundary {
		bool enabled = false;
		units::Length lower{}, upper{};
	};

	using Method = M1;
	using State = units::ScalarVectorState<units::EnergyDensity, units::EnergyFlux>;
	using Reconstruction = Method::State;
	using Flux = units::ScalarVectorState<units::EnergyFlux, units::EnergyFluxTransport>;

	explicit RadiationSystem(units::Velocity reducedLightSpeed);
	RadiationSystem(units::Velocity reducedLightSpeed, InsulatingBoundary boundary);
	/// True only for a physical domain face, including faces inside a padded
	/// predictor patch. The caller must use an inertial (nonrotating) mesh.
	bool closedEnergyFace(mesh::PhysicalCoordinates const& position, int normal, units::Length cellWidth) const;
	bool closedBoundary() const { return boundary_.enabled; }

	/// Return the configured transport speed ĉ in cm/s.
	units::Velocity reducedLightSpeed() const;

	/// Validate the physical state and form temporary (E, F/c) reconstruction values.
	Reconstruction reconstructionVariables(State const&) const;

	/// Convert a reconstructed (E, F/c) tuple to physical (E, F).
	/// Admissibility is checked separately so slope limiting can test candidate states.
	State conservedState(Reconstruction const&) const;

	/// Return ((ĉ/c)F_n, cĉP_in)-w_n(E,F), retaining physical inertial moments.
	/// The second component transports radiative flux and has units erg/(cm s²).
	Flux physicalFlux(State const&, int normal, units::Velocity faceSpeed = {}) const;

	/// Evaluate HLL with characteristic speeds shifted by the face speed, then scale
	/// its vector flux components by physical c. No radiation rest-frame boost is used.
	Flux riemann(State const&, State const&, int normal, units::Velocity faceSpeed = {}) const;

	/// Reverse only the physical flux component normal to the reflecting surface.
	State reflected(State, int normal) const;

	/// Zero inward normal radiation flux; copy energy and tangential fluxes.
	/// An insulating domain instead reflects normal F for the virtual neighbor;
	/// the face solver separately constrains the accepted energy flux after AP.
	State outflow(State, int normal, bool lower) const;

	/// Return max |lambda-w_n| for M1 normal characteristic speeds computed with ĉ.
	units::Velocity maximumSignalSpeed(State const&, int normal, units::Velocity faceSpeed = {}) const;

	/// Check finite E≥0 and |F|≤cE, with the M1 roundoff tolerance.
	bool admissible(State const&) const;

	/// Exact closed-cone test for conservative interpolation and face limiting.
	/// Do not spend a cell's roundoff allowance while choosing a slope fraction.
	bool admissibleInterpolation(State const&) const;

	/// Repair only roundoff-sized violations in calculation variables; reject larger errors.
	State correctRoundoff(State, State const& updateScale) const;

	/// Use a common face blend and test both adjacent contributions against the M1 cone.
	/// The low-order flux must already be admissible at the supplied timestep.
	/// zeroEnergyFlux preserves the insulating constraint even if the corrected
	/// low-order candidate is rejected and the limiter uses its own fallback.
	Flux limitFlux(State const&, State const&, Flux const&, int normal, units::TimePerLength stepOverCellWidth,
		units::Velocity faceSpeed = {}, Flux const* correctedLowOrder = nullptr, bool zeroEnergyFlux = false) const;

	/// Global-speed Lax-Friedrichs face flux used by the realizability fallback.
	Flux lowOrderFlux(State const&, State const&, int normal, units::Velocity faceSpeed = {}) const;

	/// Form (E, F/c) using physical c, independent of the configured transport speed.
	static Method::State toCalculationState(State const&);

	/// Multiply Q by physical c to form the physical stored flux F.
	static State fromCalculationState(Method::State const&);

	/// Scale the ndim Q-transport fluxes by physical c to transport F.
	static Flux fromCalculationFlux(Method::Flux const&);

	/// Multiply the typed state by velocity to obtain its advective flux.
	static Flux advectiveFlux(State const&, units::Velocity);

	/// Multiply a face-flux difference by Δt/Δx to obtain a typed state increment.
	static State integratedFlux(Flux const&, units::TimePerLength);

	/// Construct and validate the stored state directly from physical energy density and flux.
	static State fromPhysical(units::EnergyDensity energyDensity, std::array<units::EnergyFlux, ndim> const& physicalFlux);

	/// Validate and return the ndim stored physical flux quantities without rescaling.
	static std::array<units::EnergyFlux, ndim> toPhysicalFlux(State const&);

private:

	units::Velocity reducedLightSpeed_;
	InsulatingBoundary boundary_;
};


using Solver = finiteVolume::MusclHancock<RadiationSystem>;
using Fields = mesh::PatchData<RadiationSystem::State>;
}	 // namespace octotigerII::radiation
