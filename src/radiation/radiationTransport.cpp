// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.

#include "octotigerII/radiation/radiationTransport.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>


namespace octotigerII::radiation {

RadiationSystem::RadiationSystem(units::Velocity reducedLightSpeed)
  : RadiationSystem(reducedLightSpeed, InsulatingBoundary{}) {}

RadiationSystem::RadiationSystem(units::Velocity reducedLightSpeed, InsulatingBoundary boundary)
  : reducedLightSpeed_(reducedLightSpeed), boundary_(boundary) {
	if (!(reducedLightSpeed_ > units::Velocity{}) || !units::finite(reducedLightSpeed_)) {
		throw std::invalid_argument("Reduced light speed must be positive and finite");
	}
	if (boundary_.enabled && (!(boundary_.upper > boundary_.lower)
		|| !units::finite(boundary_.lower) || !units::finite(boundary_.upper)))
		throw std::invalid_argument("An insulating radiation boundary needs finite ordered domain bounds");
}

bool RadiationSystem::closedEnergyFace(mesh::PhysicalCoordinates const& position, int normal, units::Length cellWidth) const {
	if (!boundary_.enabled) return false;
	// Account only for coordinate-arithmetic roundoff. The cell-width cap keeps
	// a large coordinate offset from accidentally marking an interior face.
	auto const scale = std::max({units::abs(position.at(normal)), units::abs(boundary_.lower), units::abs(boundary_.upper), cellWidth});
	auto const tolerance = std::min(Real(64) * epsilonR * scale, Real(1e-6) * cellWidth);
	return units::abs(position.at(normal) - boundary_.lower) <= tolerance
		|| units::abs(position.at(normal) - boundary_.upper) <= tolerance;
}

units::Velocity RadiationSystem::reducedLightSpeed() const {
	return reducedLightSpeed_;
}

RadiationSystem::Reconstruction RadiationSystem::reconstructionVariables(State const& state) const {
	auto result = toCalculationState(state);
	Method::checkState(result);
	return result;
}

RadiationSystem::State RadiationSystem::conservedState(Reconstruction const& state) const {
	// A slope trial outside the cone must reach the generic limiter unchanged.
	// Put accepted states near the cone boundary slightly inside it, so the
	// subsequent Hancock arithmetic and grid rotation do not turn roundoff
	// into a first-order fallback. Check the physical round-trip because
	// multiplication by c can move a trial across the tolerance boundary.
	auto const physical = fromCalculationState(state);
	auto const calculation = toCalculationState(physical);
	auto const magnitude = Method::magnitude(calculation);
	constexpr Real interiorMargin = Real(64) * epsilonR;
	if (Method::admissible(calculation) && calculation[0] > units::EnergyDensity{} &&
		magnitude > (Real(1) - interiorMargin) * calculation[0]) {
		auto interior = calculation;
		auto const factor = Real((Real(1) - interiorMargin) * calculation[0] / magnitude);
		for (int i = 1; i < Reconstruction::size(); ++i) interior[i] *= factor;
		return fromCalculationState(interior);
	}
	return physical;
}

RadiationSystem::Flux RadiationSystem::physicalFlux(State const& state, int normal, units::Velocity faceSpeed) const try {
	auto result = fromCalculationFlux(Method::physicalFlux(toCalculationState(state), normal, reducedLightSpeed_).flux);
	if (faceSpeed != units::Velocity{}) result -= advectiveFlux(state, faceSpeed);
	return result;
} catch (std::runtime_error const& error) {
	throw std::runtime_error(std::string("Radiation physical face flux: ") + error.what());
}

RadiationSystem::Flux RadiationSystem::riemann(State const& left, State const& right, int normal, units::Velocity faceSpeed) const {
	if (faceSpeed == units::Velocity{})
		return fromCalculationFlux(Method::hll(toCalculationState(left), toCalculationState(right), normal, reducedLightSpeed_));
	// Only the integration surface moves: E and F retain their inertial meaning.
	auto const ul = Method::canonical(toCalculationState(left)), ur = Method::canonical(toCalculationState(right));
	auto l = Method::physicalFlux(ul, normal, reducedLightSpeed_), r = Method::physicalFlux(ur, normal, reducedLightSpeed_);
	auto const sm = std::min(l.minus, r.minus) - faceSpeed, sp = std::max(l.plus, r.plus) - faceSpeed;
	l.flux -= faceSpeed * ul;
	r.flux -= faceSpeed * ur;
	if (sm >= units::Velocity{}) return fromCalculationFlux(l.flux);
	if (sp <= units::Velocity{}) return fromCalculationFlux(r.flux);
	return fromCalculationFlux(Real(sp / (sp - sm)) * (l.flux - sm * ul) - Real(sm / (sp - sm)) * (r.flux - sp * ur));
}

RadiationSystem::State RadiationSystem::reflected(State state, int normal) const {
	state.radiativeFlux(normal) = -state.radiativeFlux(normal);
	return state;
}

RadiationSystem::State RadiationSystem::outflow(State state, int normal, bool lower) const {
	// The reflected extension supplies admissible virtual neighbors to the face
	// limiter. The accepted face-energy constraint is imposed after AP matching;
	// this ghost transform alone cannot suppress the material advective flux.
	if (boundary_.enabled) return reflected(state, normal);
	auto& flux = state.radiativeFlux(normal);
	if (lower ? flux > units::EnergyFlux{} : flux < units::EnergyFlux{}) flux = {};
	return state;
}

units::Velocity RadiationSystem::maximumSignalSpeed(State const& state, int normal, units::Velocity faceSpeed) const try {
	auto const waves = Method::physicalFlux(toCalculationState(state), normal, reducedLightSpeed_);
	return std::max(units::abs(waves.minus - faceSpeed), units::abs(waves.plus - faceSpeed));
} catch (std::runtime_error const& error) {
	throw std::runtime_error(std::string("Radiation signal speed: ") + error.what());
}

bool RadiationSystem::admissible(State const& state) const {
	return Method::admissible(toCalculationState(state));
}

bool RadiationSystem::admissibleInterpolation(State const& physical) const {
	auto const state = toCalculationState(physical);
	return finite(state) && state[0] >= units::EnergyDensity{} && Method::magnitude(state) <= state[0];
}

RadiationSystem::State RadiationSystem::correctRoundoff(State state, State const& updateScale) const {
	return fromCalculationState(Method::roundoffState(toCalculationState(state), toCalculationState(updateScale)));
}

RadiationSystem::Flux RadiationSystem::limitFlux(
	State const& left, State const& right, Flux const& highOrderFlux, int normal, units::TimePerLength stepOverCellWidth,
	units::Velocity faceSpeed, Flux const* correctedLowOrder, bool zeroEnergyFlux) const {
	if (zeroEnergyFlux && (highOrderFlux.energy() != units::EnergyFlux{} || faceSpeed != units::Velocity{}))
		throw std::invalid_argument("An insulating radiation face requires zero energy flux on a fixed grid");
	if (stepOverCellWidth == units::TimePerLength{}) {
		return highOrderFlux;
	}
	Flux const leftPhysical = physicalFlux(left, normal, faceSpeed);
	Flux const rightPhysical = physicalFlux(right, normal, faceSpeed);
	auto const factor = Real(2 * ndim) * stepOverCellWidth;
	auto validState = [&](State const& physical) {
		// A tolerance based on the brighter neighbor can admit cone overshoot
		// larger than the faint cell's complete update and its roundoff budget.
		// Keep these convex face states inside the cone; the final summed update
		// still has its own scale-aware roundoff repair.
		return admissibleInterpolation(physical);
	};
	auto validFlux = [&](Flux const& flux) {
		return validState(State(left - integratedFlux(flux - leftPhysical, factor))) && validState(State(right + integratedFlux(flux - rightPhysical, factor)));
	};
	auto roundoffValidFlux = [&](Flux const& flux) {
		// A face flux can land just outside the exact cone when an input lies
		// on its boundary. Match the update's local M1 roundoff
		// allowance here without borrowing a tolerance from the other cell.
		return admissible(State(left - integratedFlux(flux - leftPhysical, factor))) &&
			admissible(State(right + integratedFlux(flux - rightPhysical, factor)));
	};
	if (validFlux(highOrderFlux) || roundoffValidFlux(highOrderFlux)) {
		return highOrderFlux;
	}

	Flux lowFlux = lowOrderFlux(left, right, normal, faceSpeed);
	if (zeroEnergyFlux) lowFlux.energy() = {};
	if (correctedLowOrder && (validFlux(*correctedLowOrder) || roundoffValidFlux(*correctedLowOrder))) lowFlux = *correctedLowOrder;
	if (zeroEnergyFlux && lowFlux.energy() != units::EnergyFlux{})
		throw std::invalid_argument("Insulating radiation limiter received a leaking low-order flux");
	if (!validFlux(lowFlux) && !roundoffValidFlux(lowFlux)) {
		throw std::runtime_error("First-order M1 flux violates realizability at this timestep");
	}
	Real low = 0;
	Real high = 1;
	for (int iteration = 0; iteration < 56; ++iteration) {
		Real const fraction = Real(0.5) * (low + high);
		Flux const candidate = Flux(lowFlux + fraction * (highOrderFlux - lowFlux));
		if (validFlux(candidate)) {
			low = fraction;
		} else {
			high = fraction;
		}
	}
	return Flux(lowFlux + low * (highOrderFlux - lowFlux));
}

RadiationSystem::Flux RadiationSystem::lowOrderFlux(State const& left, State const& right, int normal, units::Velocity faceSpeed) const {
	auto const speed = (reducedLightSpeed_ + units::abs(faceSpeed)) * (Real(1) + Method::roundoff);
	return Flux(Real(0.5) * (physicalFlux(left, normal, faceSpeed) + physicalFlux(right, normal, faceSpeed)
		- advectiveFlux(right - left, speed)));
}

RadiationSystem::State RadiationSystem::fromPhysical(units::EnergyDensity energyDensity, std::array<units::EnergyFlux, ndim> const& physicalFlux) {
	State result{};
	result.energy() = energyDensity;
	for (int axis = 0; axis < ndim; ++axis)
		result.radiativeFlux(axis) = physicalFlux[axis];
	Method::checkState(toCalculationState(result));
	return result;
}

std::array<units::EnergyFlux, ndim> RadiationSystem::toPhysicalFlux(State const& state) {
	Method::checkState(toCalculationState(state));
	std::array<units::EnergyFlux, ndim> result{};
	for (int axis = 0; axis < ndim; ++axis)
		result[axis] = state.radiativeFlux(axis);
	return result;
}

RadiationSystem::Method::State RadiationSystem::toCalculationState(State const& state) {
	Method::State result{};
	result[0] = state.energy();
	for (int axis = 0; axis < ndim; ++axis)
		result[axis + 1] = state.radiativeFlux(axis) / constants::c;
	return result;
}

RadiationSystem::State RadiationSystem::fromCalculationState(Method::State const& state) {
	State result{};
	result.energy() = state[0];
	for (int axis = 0; axis < ndim; ++axis)
		result.radiativeFlux(axis) = constants::c * state[axis + 1];
	return result;
}

RadiationSystem::Flux RadiationSystem::fromCalculationFlux(Method::Flux const& flux) {
	Flux result{};
	result.energy() = flux[0];
	for (int axis = 0; axis < ndim; ++axis)
		result.radiativeFlux(axis) = constants::c * flux[axis + 1];
	return result;
}

RadiationSystem::Flux RadiationSystem::advectiveFlux(State const& state, units::Velocity speed) {
	return state * speed;
}

RadiationSystem::State RadiationSystem::integratedFlux(Flux const& flux, units::TimePerLength factor) {
	return flux * factor;
}
}	 // namespace octotigerII::radiation
