/** @file
 * @brief Euler states, selectable thermodynamic closures, and HLLC fluxes.
 * @ingroup numerics
 */
// Distributed under the Boost Software License, Version 1.0.
#pragma once
#include "octotigerII/mesh.hpp"
#include "octotigerII/config.hpp"
#include "octotigerII/physics/finiteVolume.hpp"
#include "octotigerII/hydro/whiteDwarfEos.hpp"
#include "octotigerII/hydro/helmholtzClosure.hpp"
#include "octotigerII/units/cgs.hpp"
#include "octotigerII/units/state.hpp"


namespace octotigerII::hydro {
using ConservedState = units::FluidState<units::Density, units::MomentumDensity, units::EnergyDensity, units::Density, units::Density, units::Density>;
using PrimitiveState = units::FluidState<units::Density, units::Velocity, units::Pressure, units::Dimensionless, units::Dimensionless, units::Dimensionless>;
using ConservedFlux = units::FluidState<units::MassFlux, units::MomentumFlux, units::EnergyFlux, units::MassFlux, units::MassFlux, units::MassFlux>;


/// Euler flux adapter with primitive reconstruction and selectable EOS.
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

	explicit HydroSystem(Config::HydroOptions const&, bool separateRadiation = false);
	explicit HydroSystem(Config const&);
	bool helmholtz() const { return bool(helmholtz_); }
	bool gasOnly() const { return !helmholtz_ || helmholtz_->gasOnly(); }
	double temperatureFloor() const { return helmholtz_ ? helmholtz_->temperatureFloor() : 0; }
	HelmholtzClosure::Point thermodynamics(State const&, units::EnergyDensity) const;
	State stateFromTemperature(units::Density, units::Temperature, Real abar, Real zbar) const;
	void setComposition(State&, std::vector<units::Density> const&, composition::Options const&) const;
	units::Density auxiliaryFromInternalEnergy(State const&, units::EnergyDensity) const;
	/// Returns the actual energy density added; only call on accepted states.
	units::EnergyDensity applyTemperatureFloor(State&) const;
	static Flux compositionFlux(Flux, State const&, State const&);
	void constrainComposition(State& candidate, State const& donor) const;

	/// Selected internal energy density (thermal-only for the cold white-dwarf model).
	units::EnergyDensity internalEnergy(State const&) const;
	units::Temperature temperature(State const&) const;

	/// Encode/decode the auxiliary: EOS entropy for Helmholtz (which may be signed),
	/// or the positive constant-gamma entropy power for the other closures.
	units::Density auxiliaryFromInternalEnergy(units::Density, units::EnergyDensity) const;
	units::EnergyDensity internalEnergyFromAuxiliary(State const&) const;

	/// Reset only A, when (E-K)/E exceeds the upper threshold. Never changes E.
	void synchronize(State&) const;

	/// Return the configured ideal-gas gamma; not a Helmholtz derivative.
	Real adiabaticIndex() const;

	/// Convert (rho, rho v, E, A) to (rho, v, P, A/rho); reject inadmissible input.
	PrimitiveState reconstructionVariables(State const&) const;

	/// Initialize/constrain (rho, rho v, E, A) from CGS primitives.
	/// A zero primitive auxiliary requests initialization from rho and pressure;
	/// a positive value preserves independently reconstructed A/rho.
	State conservedState(PrimitiveState const&) const;

	/// Return the inertial Euler ALE flux F_n-w_n U, using components along the face axes.
	Flux physicalFlux(State const&, int normal, units::Velocity faceSpeed = {}) const;

	/// Sample HLLC at the moving face, retaining inertial momentum and energy fluxes.
	/// Falls back to HLL for degenerate or inadmissible star states.
	Flux riemann(State const&, State const&, int normal, units::Velocity faceSpeed = {}) const;

	/// Reverse normal momentum while preserving density, tangential momenta, and energy.
	State reflected(State, int normal) const;

	/// Clip inward velocity relative to the face. A moving-face clip preserves internal
	/// energy; the stationary-face overload retains the established boundary rule.
	State outflow(State, int normal, bool lower, units::Velocity faceSpeed = {}) const;

	/// Return |v_n-w_n| plus the adiabatic sound speed.
	units::Velocity maximumSignalSpeed(State const&, int normal, units::Velocity faceSpeed = {}) const;
	units::Velocity adiabaticSoundSpeed(State const&) const;

	/// Require finite state values and density/pressure at or above the configured floors.
	bool admissible(State const&) const;

	/// Return the candidate unchanged; hydro does not silently repair a failed state.
	State correctRoundoff(State, State const& updateScale) const;

	/// Blend toward a first-order local Lax–Friedrichs flux when needed.
	/// A common face coefficient preserves conservative flux sharing. Related idea:
	/// @ref ref_hu2013 "Hu et al. (2013)"; this implementation selects it by bisection.
	Flux limitFlux(State const&, State const&, Flux const&, int normal, units::TimePerLength stepOverCellWidth, units::Velocity faceSpeed = {}) const;

	/// Multiply a face-flux difference by Δt/Δx to obtain a typed state increment.
	static State integratedFlux(Flux const&, units::TimePerLength);

	/// Multiply the typed state by velocity to obtain its advective flux.
	static Flux advectiveFlux(State const&, units::Velocity);

private:

	DualEnergyOptions dualEnergy_;
	Real meanMolecularWeight_;
	Real adiabaticIndex_;
	bool degenerate_ = false;
	WhiteDwarfEos whiteDwarfEos_;
	std::shared_ptr<const HelmholtzClosure> helmholtz_;
	Real defaultAbar_ = 12, defaultZbar_ = 6;
	std::pair<Real, Real> composition(State const&) const;
	HelmholtzClosure::Point minimum(State const&) const;
	units::Density densityFloor_;
	units::Pressure pressureFloor_;

	units::EnergyDensity totalInternalEnergy(State const&) const;
	units::EnergyDensity degenerateEnergy(units::Density) const;
	units::Pressure degeneratePressure(units::Density) const;
	units::Velocity soundSpeed(PrimitiveState const&) const;

	/// Evaluate pressure using the selected thermodynamic closure.
	units::Pressure pressure(State const&) const;

	/// Two-wave hydro fallback from @ref ref_harten1983 "Harten et al. (1983)".
	Flux hll(State const&, State const&, int normal) const;
};


inline HydroSystem::HydroSystem(Config const& config) : HydroSystem(config.hydro, config.radiationEnabled()) {
    if (helmholtz_ && config.massFractions.enabled) {
        auto state = State{};
        auto const rho = units::Density::from_value(1);
        setComposition(state, composition::initialDensities(config.massFractions, rho), config.massFractions);
        auto const c = composition(state);
        defaultAbar_ = c.first; defaultZbar_ = c.second;
    }
}
using Solver = physics::MusclHancock<HydroSystem>;
using Fields = mesh::PatchData<ConservedState>;
}	 // namespace octotigerII::hydro
