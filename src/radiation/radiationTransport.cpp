// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.

#include "octotigerII/radiation/radiationTransport.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>


namespace octotigerII::radiation {

RadiationSystem::RadiationSystem(units::Velocity reducedLightSpeed)
  : reducedLightSpeed_(reducedLightSpeed) {
	if (!(reducedLightSpeed_ > units::Velocity{}) || !units::finite(reducedLightSpeed_)) {
		throw std::invalid_argument("Reduced light speed must be positive and finite");
	}
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
	return fromCalculationState(state);
}

RadiationSystem::Flux RadiationSystem::physicalFlux(State const& state, int normal) const {
	return fromCalculationFlux(Method::physicalFlux(toCalculationState(state), normal, reducedLightSpeed_).flux);
}

RadiationSystem::Flux RadiationSystem::riemann(State const& left, State const& right, int normal) const {
	return fromCalculationFlux(Method::hll(toCalculationState(left), toCalculationState(right), normal, reducedLightSpeed_));
}

RadiationSystem::State RadiationSystem::reflected(State state, int normal) const {
	state.radiativeFlux(normal) = -state.radiativeFlux(normal);
	return state;
}

units::Velocity RadiationSystem::maximumSignalSpeed(State const& state, int normal) const {
	auto const waves = Method::physicalFlux(toCalculationState(state), normal, reducedLightSpeed_);
	return std::max(units::abs(waves.minus), units::abs(waves.plus));
}

bool RadiationSystem::admissible(State const& state) const {
	return Method::admissible(toCalculationState(state));
}

RadiationSystem::State RadiationSystem::correctRoundoff(State state, State const& updateScale) const {
	return fromCalculationState(Method::roundoffState(toCalculationState(state), toCalculationState(updateScale)));
}

RadiationSystem::Flux RadiationSystem::limitFlux(
	State const& left, State const& right, Flux const& highOrderFlux, int normal, units::TimePerLength stepOverCellWidth) const {
	if (stepOverCellWidth == units::TimePerLength{}) {
		return highOrderFlux;
	}
	Flux const leftPhysical = physicalFlux(left, normal);
	Flux const rightPhysical = physicalFlux(right, normal);
	auto const factor = Real(2 * ndim) * stepOverCellWidth;
	auto const tolerance = Method::roundoff * (left.energy() + right.energy());
	auto validState = [&](State const& physical) {
		auto const state = toCalculationState(physical);
		if (Method::admissible(state)) {
			return true;
		}
		if (!units::finite(state[0]) || state[0] < -tolerance) {
			return false;
		}
		units::EnergyDensity magnitude{};
		for (int field = 1; field < State::size(); ++field) {
			if (!units::finite(state[field])) {
				return false;
			}
			magnitude = units::hypot(magnitude, state[field]);
		}
		return magnitude - std::max(units::EnergyDensity{}, state[0]) <= tolerance;
	};
	auto validFlux = [&](Flux const& flux) {
		return validState(State(left - integratedFlux(flux - leftPhysical, factor))) && validState(State(right + integratedFlux(flux - rightPhysical, factor)));
	};
	if (validFlux(highOrderFlux)) {
		return highOrderFlux;
	}

	auto const speed = reducedLightSpeed_ * (Real(1) + Method::roundoff);
	Flux const lowOrderFlux = Flux(Real(0.5) * (leftPhysical + rightPhysical - advectiveFlux(right - left, speed)));
	if (!validFlux(lowOrderFlux)) {
		throw std::runtime_error("First-order M1 flux violates realizability at this timestep");
	}
	Real low = 0;
	Real high = 1;
	for (int iteration = 0; iteration < 56; ++iteration) {
		Real const fraction = Real(0.5) * (low + high);
		Flux const candidate = Flux(lowOrderFlux + fraction * (highOrderFlux - lowOrderFlux));
		if (validFlux(candidate)) {
			low = fraction;
		} else {
			high = fraction;
		}
	}
	return Flux(lowOrderFlux + low * (highOrderFlux - lowOrderFlux));
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
