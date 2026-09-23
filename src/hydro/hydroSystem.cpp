// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.

#include "octotigerII/hydro/hydroSystem.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>


namespace octotigerII::hydro {

HydroSystem::HydroSystem(Real adiabaticIndex, units::Density densityFloor, units::Pressure pressureFloor)
  : adiabaticIndex_(adiabaticIndex)
  , densityFloor_(densityFloor)
  , pressureFloor_(pressureFloor) {
	if (!(adiabaticIndex_ > 1) || !std::isfinite(adiabaticIndex_) || !(densityFloor_ > units::Density{}) || !units::finite(densityFloor_) ||
		!(pressureFloor_ > units::Pressure{}) || !units::finite(pressureFloor_)) {
		throw std::invalid_argument("Hydro EOS must have finite gamma > 1 and finite positive floors");
	}
}

Real HydroSystem::adiabaticIndex() const {
	return adiabaticIndex_;
}

PrimitiveState HydroSystem::reconstructionVariables(ConservedState const& state) const {
	if (!admissible(state)) {
		throw std::runtime_error("Cannot convert an inadmissible hydro state to primitive variables");
	}
	PrimitiveState result;
	result.setDensity(state.density());
	for (int axis = 0; axis < ndim; ++axis) {
		result.setVelocity(axis, state.momentum(axis) / state.density());
	}
	result.setPressure(pressure(state));
	return result;
}

ConservedState HydroSystem::conservedState(PrimitiveState const& state) const {
	if (!(state.density() >= densityFloor_) || !(state.pressure() >= pressureFloor_) || !units::finite(state.density()) || !units::finite(state.pressure())) {
		throw std::invalid_argument("Cannot convert an inadmissible hydro primitive state");
	}
	ConservedState result;
	result.setDensity(state.density());
	units::VelocitySquared speedSquared{};
	for (int axis = 0; axis < ndim; ++axis) {
		if (!units::finite(state.velocity(axis))) {
			throw std::invalid_argument("Hydro primitive velocity must be finite");
		}
		result.setMomentum(axis, state.density() * state.velocity(axis));
		speedSquared += state.velocity(axis) * state.velocity(axis);
	}
	result.setTotalEnergy(state.pressure() / (adiabaticIndex_ - 1) + Real(0.5) * state.density() * speedSquared);
	return result;
}

ConservedFlux HydroSystem::physicalFlux(ConservedState const& state, int normal) const {
	PrimitiveState const primitive = reconstructionVariables(state);
	auto const normalVelocity = primitive.velocity(normal);
	ConservedFlux result;
	result.setMass(state.density() * normalVelocity);
	for (int axis = 0; axis < ndim; ++axis) {
		result.setMomentum(axis, state.momentum(axis) * normalVelocity + (axis == normal ? primitive.pressure() : units::Pressure{}));
	}
	result.setEnergy((state.totalEnergy() + primitive.pressure()) * normalVelocity);
	return result;
}

ConservedFlux HydroSystem::riemann(ConservedState const& left, ConservedState const& right, int normal) const {
	PrimitiveState const leftPrimitive = reconstructionVariables(left);
	PrimitiveState const rightPrimitive = reconstructionVariables(right);
	auto const leftSound = units::sqrt(adiabaticIndex_ * leftPrimitive.pressure() / leftPrimitive.density());
	auto const rightSound = units::sqrt(adiabaticIndex_ * rightPrimitive.pressure() / rightPrimitive.density());
	auto const leftSpeed = std::min(leftPrimitive.velocity(normal) - leftSound, rightPrimitive.velocity(normal) - rightSound);
	auto const rightSpeed = std::max(leftPrimitive.velocity(normal) + leftSound, rightPrimitive.velocity(normal) + rightSound);
	ConservedFlux const leftFlux = physicalFlux(left, normal);
	ConservedFlux const rightFlux = physicalFlux(right, normal);
	if (leftSpeed >= units::Velocity{}) {
		return leftFlux;
	}
	if (rightSpeed <= units::Velocity{}) {
		return rightFlux;
	}

	auto const leftDenominator = leftPrimitive.density() * (leftSpeed - leftPrimitive.velocity(normal));
	auto const rightDenominator = rightPrimitive.density() * (rightSpeed - rightPrimitive.velocity(normal));
	auto const denominator = leftDenominator - rightDenominator;
	auto const denominatorScale = std::max({units::MassFlux::from_value(1), units::abs(leftDenominator), units::abs(rightDenominator)});
	auto const waveScale = std::max({units::Velocity::from_value(1), units::abs(leftSpeed), units::abs(rightSpeed)});
	if (units::abs(denominator) <= 64 * epsilonR * denominatorScale) {
		return hll(left, right, normal);
	}
	auto const contactSpeed = (rightPrimitive.pressure() - leftPrimitive.pressure() +
								  leftPrimitive.density() * leftPrimitive.velocity(normal) * (leftSpeed - leftPrimitive.velocity(normal)) -
								  rightPrimitive.density() * rightPrimitive.velocity(normal) * (rightSpeed - rightPrimitive.velocity(normal))) /
		denominator;
	if (!units::finite(contactSpeed) || contactSpeed < leftSpeed || contactSpeed > rightSpeed) {
		return hll(left, right, normal);
	}

	auto starState = [&](ConservedState const& state, PrimitiveState const& primitive, units::Velocity waveSpeed) {
		auto const waveDifference = waveSpeed - primitive.velocity(normal);
		auto const starDifference = waveSpeed - contactSpeed;
		if (units::abs(starDifference) <= 64 * epsilonR * waveScale || units::abs(waveDifference) <= 64 * epsilonR * waveScale) {
			return ConservedState{};
		}
		ConservedState star;
		star.setDensity(primitive.density() * waveDifference / starDifference);
		for (int axis = 0; axis < ndim; ++axis) {
			auto const velocity = axis == normal ? contactSpeed : primitive.velocity(axis);
			star.setMomentum(axis, star.density() * velocity);
		}
		star.setTotalEnergy(star.density() *
			(state.totalEnergy() / primitive.density() +
				(contactSpeed - primitive.velocity(normal)) * (contactSpeed + primitive.pressure() / (primitive.density() * waveDifference))));
		return star;
	};

	if (contactSpeed >= units::Velocity{}) {
		ConservedState const star = starState(left, leftPrimitive, leftSpeed);
		if (!admissible(star)) {
			return hll(left, right, normal);
		}
		return leftFlux + advectiveFlux(star - left, leftSpeed);
	}
	ConservedState const star = starState(right, rightPrimitive, rightSpeed);
	if (!admissible(star)) {
		return hll(left, right, normal);
	}
	return rightFlux + advectiveFlux(star - right, rightSpeed);
}

ConservedState HydroSystem::reflected(ConservedState state, int normal) const {
	state.setMomentum(normal, -state.momentum(normal));
	return state;
}

ConservedState HydroSystem::outflow(State state, int normal, bool lower) const {
	auto& momentum = state.momentum(normal);
	if (lower ? momentum > units::MomentumDensity{} : momentum < units::MomentumDensity{}) momentum = {};
	return state;
}

units::Velocity HydroSystem::maximumSignalSpeed(ConservedState const& state, int normal) const {
	PrimitiveState const primitive = reconstructionVariables(state);
	return units::abs(primitive.velocity(normal)) + units::sqrt(adiabaticIndex_ * primitive.pressure() / primitive.density());
}

bool HydroSystem::admissible(ConservedState const& state) const {
	if (!finite(state)) return false;
	return state.density() >= densityFloor_ && pressure(state) >= pressureFloor_;
}

ConservedState HydroSystem::correctRoundoff(ConservedState state, State const& updateScale) const {
	(void) updateScale;
	return state;
}

ConservedFlux HydroSystem::limitFlux(
	ConservedState const& left, ConservedState const& right, ConservedFlux const& highOrderFlux, int normal, units::TimePerLength stepOverCellWidth) const {
	if (stepOverCellWidth == units::TimePerLength{}) {
		return highOrderFlux;
	}
	ConservedFlux const leftPhysical = physicalFlux(left, normal);
	ConservedFlux const rightPhysical = physicalFlux(right, normal);
	auto const factor = Real(2 * ndim) * stepOverCellWidth;
	auto validFlux = [&](ConservedFlux const& flux) {
		return admissible(left - integratedFlux(flux - leftPhysical, factor)) && admissible(right + integratedFlux(flux - rightPhysical, factor));
	};
	if (validFlux(highOrderFlux)) {
		return highOrderFlux;
	}
	auto const speed = std::max(maximumSignalSpeed(left, normal), maximumSignalSpeed(right, normal));
	ConservedFlux const lowOrderFlux = Real(0.5) * (leftPhysical + rightPhysical - advectiveFlux(right - left, speed));
	if (!validFlux(lowOrderFlux)) {
		throw std::runtime_error("First-order hydro flux is not positivity preserving at this timestep");
	}
	Real low = 0;
	Real high = 1;
	for (int iteration = 0; iteration < 56; ++iteration) {
		Real const fraction = Real(0.5) * (low + high);
		if (validFlux(lowOrderFlux + fraction * (highOrderFlux - lowOrderFlux))) {
			low = fraction;
		} else {
			high = fraction;
		}
	}
	return lowOrderFlux + low * (highOrderFlux - lowOrderFlux);
}

units::Pressure HydroSystem::pressure(ConservedState const& state) const {
	if (!(state.density() > units::Density{}) || !units::finite(state.density())) {
		return units::Pressure::from_value(-std::numeric_limits<Real>::infinity());
	}
	decltype(units::MomentumDensity{} * units::MomentumDensity{}) momentumSquared{};
	for (int axis = 0; axis < ndim; ++axis) {
		momentumSquared += state.momentum(axis) * state.momentum(axis);
	}
	return (adiabaticIndex_ - 1) * (state.totalEnergy() - Real(0.5) * momentumSquared / state.density());
}

ConservedFlux HydroSystem::hll(ConservedState const& left, ConservedState const& right, int normal) const {
	PrimitiveState const leftPrimitive = reconstructionVariables(left);
	PrimitiveState const rightPrimitive = reconstructionVariables(right);
	auto const leftSound = units::sqrt(adiabaticIndex_ * leftPrimitive.pressure() / leftPrimitive.density());
	auto const rightSound = units::sqrt(adiabaticIndex_ * rightPrimitive.pressure() / rightPrimitive.density());
	auto const leftSpeed = std::min(leftPrimitive.velocity(normal) - leftSound, rightPrimitive.velocity(normal) - rightSound);
	auto const rightSpeed = std::max(leftPrimitive.velocity(normal) + leftSound, rightPrimitive.velocity(normal) + rightSound);
	ConservedFlux const leftFlux = physicalFlux(left, normal);
	ConservedFlux const rightFlux = physicalFlux(right, normal);
	if (leftSpeed >= units::Velocity{}) {
		return leftFlux;
	}
	if (rightSpeed <= units::Velocity{}) {
		return rightFlux;
	}
	return Real(rightSpeed / (rightSpeed - leftSpeed)) * leftFlux - Real(leftSpeed / (rightSpeed - leftSpeed)) * rightFlux +
		advectiveFlux(right - left, leftSpeed * rightSpeed / (rightSpeed - leftSpeed));
}

ConservedFlux HydroSystem::advectiveFlux(ConservedState const& state, units::Velocity speed) {
	ConservedFlux result;
	result.setMass(state.density() * speed);
	for (int axis = 0; axis < ndim; ++axis)
		result.setMomentum(axis, state.momentum(axis) * speed);
	result.setEnergy(state.totalEnergy() * speed);
	return result;
}

ConservedState HydroSystem::integratedFlux(ConservedFlux const& flux, units::TimePerLength factor) {
	ConservedState result;
	result.setDensity(flux.mass() * factor);
	for (int axis = 0; axis < ndim; ++axis)
		result.setMomentum(axis, flux.momentum(axis) * factor);
	result.setTotalEnergy(flux.energy() * factor);
	return result;
}
}	 // namespace octotigerII::hydro
