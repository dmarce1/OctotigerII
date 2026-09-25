/** @file
 * @brief Ideal-gas Euler states, HLLC fluxes, and positivity safeguards.
 * @ingroup numerics
 */
// Distributed under the Boost Software License, Version 1.0.
#pragma once
#include "octotigerII/mesh.hpp"
#include "octotigerII/config.hpp"
#include "octotigerII/physics/finiteVolume.hpp"
#include "octotigerII/units/cgs.hpp"
#include "octotigerII/units/state.hpp"


namespace octotigerII::hydro {
using ConservedState = units::FluidState<units::Density, units::MomentumDensity, units::EnergyDensity, units::Density>;
using PrimitiveState = units::FluidState<units::Density, units::Velocity, units::Pressure, units::Dimensionless>;
using ConservedFlux = units::FluidState<units::MassFlux, units::MomentumFlux, units::EnergyFlux, units::MassFlux>;


/// Ideal-gas Euler flux adapter with primitive reconstruction.
/// HLLC follows @ref ref_toro1994 "Toro et al. (1994)"; inadmissible or degenerate
/// star states fall back to the HLL flux of @ref ref_harten1983 "Harten et al. (1983)".
/// @ingroup numerics
class HydroSystem {
public:

	using State = ConservedState;
	using Reconstruction = PrimitiveState;
	using Flux = ConservedFlux;

	explicit HydroSystem(Real adiabaticIndex = Real(5) / 3, units::Density densityFloor = units::Density::from_value(1e-14),
		units::Pressure pressureFloor = units::Pressure::from_value(1e-14), DualEnergyOptions dualEnergy = {}, Real meanMolecularWeight = 1);

	explicit HydroSystem(Config::HydroOptions const&);

	/// Selected thermal energy density, using the lower dual-energy threshold.
	units::EnergyDensity internalEnergy(State const&) const;
	units::Temperature temperature(State const&) const;

	/// Encode/decode the normalized entropy density. Positive arguments required.
	units::Density auxiliaryFromInternalEnergy(units::Density, units::EnergyDensity) const;
	units::EnergyDensity internalEnergyFromAuxiliary(State const&) const;

	/// Reset only A, when (E-K)/E exceeds the upper threshold. Never changes E.
	void synchronize(State&) const;

	/// Return the ideal-gas ratio of specific heats.
	Real adiabaticIndex() const;

	/// Convert (rho, rho v, E, A) to (rho, v, P, A/rho); reject inadmissible input.
	PrimitiveState reconstructionVariables(State const&) const;

	/// Initialize/constrain (rho, rho v, E, A) from CGS primitives.
	/// A zero primitive auxiliary requests initialization from rho and pressure;
	/// a positive value preserves independently reconstructed A/rho.
	State conservedState(PrimitiveState const&) const;

	/// Return the Euler flux normal to an axis in [0, ndim).
	Flux physicalFlux(State const&, int normal) const;

	/// Return the HLLC face flux, falling back to HLL for degenerate or inadmissible star states.
	Flux riemann(State const&, State const&, int normal) const;

	/// Reverse normal momentum while preserving density, tangential momenta, and energy.
	State reflected(State, int normal) const;

	/// Zero inward normal momentum; copy density, total energy, and tangential momenta.
	State outflow(State, int normal, bool lower) const;

	/// Return |v_n| plus the adiabatic sound speed.
	units::Velocity maximumSignalSpeed(State const&, int normal) const;

	/// Require finite state values and density/pressure at or above the configured floors.
	bool admissible(State const&) const;

	/// Return the candidate unchanged; hydro does not silently repair a failed state.
	State correctRoundoff(State, State const& updateScale) const;

	/// Blend toward a first-order local Lax–Friedrichs flux when needed.
	/// A common face coefficient preserves conservative flux sharing. Related idea:
	/// @ref ref_hu2013 "Hu et al. (2013)"; this implementation selects it by bisection.
	Flux limitFlux(State const&, State const&, Flux const&, int normal, units::TimePerLength stepOverCellWidth) const;

	/// Multiply a face-flux difference by Δt/Δx to obtain a typed state increment.
	static State integratedFlux(Flux const&, units::TimePerLength);

	/// Multiply the typed state by velocity to obtain its advective flux.
	static Flux advectiveFlux(State const&, units::Velocity);

private:

	DualEnergyOptions dualEnergy_;
	Real meanMolecularWeight_;
	Real adiabaticIndex_;
	units::Density densityFloor_;
	units::Pressure pressureFloor_;

	units::EnergyDensity totalInternalEnergy(State const&) const;

	/// Compute P=(gamma-1)u from the selected internal energy.
	units::Pressure pressure(State const&) const;

	/// Two-wave hydro fallback from @ref ref_harten1983 "Harten et al. (1983)".
	Flux hll(State const&, State const&, int normal) const;
};


using Solver = physics::MusclHancock<HydroSystem>;
using Fields = mesh::PatchData<ConservedState>;
}	 // namespace octotigerII::hydro
