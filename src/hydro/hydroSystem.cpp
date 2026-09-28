// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.

#include "octotigerII/hydro/hydroSystem.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>


namespace octotigerII::hydro {

namespace {
	// Logarithmic arithmetic avoids intermediate powers overflowing even when
	// their final product is representable. Long double also handles tiny alpha.
	Real positiveExp(long double logarithm) {
		if (!std::isfinite(logarithm) || logarithm < std::log(static_cast<long double>(std::numeric_limits<Real>::min())) ||
			logarithm > std::log(static_cast<long double>(std::numeric_limits<Real>::max())))
			throw std::runtime_error("Dual-energy value outside the positive floating-point range; choose a less extreme exponent");
		return static_cast<Real>(std::exp(logarithm));
	}
}

HydroSystem::HydroSystem(Real adiabaticIndex, units::Density densityFloor, units::Pressure pressureFloor,
	DualEnergyOptions dualEnergy, Real meanMolecularWeight)
  : dualEnergy_(dualEnergy)
  , meanMolecularWeight_(meanMolecularWeight)
  , adiabaticIndex_(adiabaticIndex)
  , densityFloor_(densityFloor)
  , pressureFloor_(pressureFloor) {
	dualEnergy_.validate();
	if (!std::isfinite(meanMolecularWeight_) || !(meanMolecularWeight_ > 0))
		throw std::invalid_argument("Mean molecular weight must be finite and positive");
	if (!(adiabaticIndex_ > 1) || !std::isfinite(adiabaticIndex_) || !(densityFloor_ > units::Density{}) || !units::finite(densityFloor_) ||
		!(pressureFloor_ > units::Pressure{}) || !units::finite(pressureFloor_)) {
		throw std::invalid_argument("Hydro EOS must have finite gamma > 1 and finite positive floors");
	}
}

HydroSystem::HydroSystem(Config::HydroOptions const& options, bool separateRadiation)
  : HydroSystem(options.gamma, units::Density::from_value(1e-14), units::Pressure::from_value(1e-14), options.dualEnergy, options.meanMolecularWeight) {
	if (options.eos != "ideal" && options.eos != "white-dwarf" && options.eos != "helmholtz")
		throw std::invalid_argument("Unsupported hydro EOS");
	if (options.eos == "helmholtz") helmholtz_ = std::make_shared<HelmholtzClosure>(options.helmholtzTable, separateRadiation, options.temperatureFloor);
	degenerate_ = options.eos == "white-dwarf";
	whiteDwarfEos_ = WhiteDwarfEos(options.meanMassPerElectron);
}

std::pair<Real, Real> HydroSystem::composition(State const& state) const {
    if (state.nuclei() == units::Density{} && state.electrons() == units::Density{})
        return {defaultAbar_, defaultZbar_};
    if (!(state.density()>units::Density{}) || !(state.nuclei()>units::Density{}) || !(state.electrons()>units::Density{}))
        throw std::invalid_argument("Invalid Helmholtz composition moments");
    Real const a = Real(state.density() / state.nuclei());
    Real const z = Real(state.electrons() / state.nuclei());
    if (!std::isfinite(a) || !std::isfinite(z) || !(a > 0 && z > 0 && z <= a))
        throw std::invalid_argument("Invalid Helmholtz material composition");
    return {a,z};
}
void HydroSystem::setComposition(State& state, std::vector<units::Density> const& species, composition::Options const& options) const {
    if (!helmholtz_) return;
    if (species.size() != options.species.size()) throw std::invalid_argument("Composition field count mismatch");
    state.density() = {}; state.nuclei() = {}; state.electrons() = {};
    for (std::size_t i=0; i<species.size(); ++i) {
        auto const& s = options.species[i];
        if (s.tracer()) continue;
        if (!units::finite(species[i]) || species[i] < units::Density{}) throw std::invalid_argument("Invalid material partial density");
        state.density() += species[i];
        state.nuclei() += species[i] / s.atomicMass;
        state.electrons() += species[i] * (s.atomicNumber / s.atomicMass);
    }
    (void) composition(state);
}
HelmholtzClosure::Point HydroSystem::minimum(State const& state) const {
    auto const [a,z] = composition(state);
    return helmholtz_->at(units::value(state.density()), helmholtz_->temperatureFloor(), a,z);
}
HelmholtzClosure::Point HydroSystem::thermodynamics(State const& state, units::EnergyDensity energy) const {
    if (!helmholtz_) throw std::logic_error("Helmholtz thermodynamics requested for another EOS");
    auto const [a,z] = composition(state);
    auto const q = helmholtz_->invert(units::value(state.density()), units::value(energy / state.density()), a,z, HelmholtzClosure::Variable::Energy);
    if (!(q.cv > 0 && q.pressure > 0 && q.soundSquared > 0) || !std::isfinite(q.cv) || !std::isfinite(q.soundSquared))
        throw std::runtime_error("Helmholtz state has nonpositive/nonfinite heat capacity, pressure or sound speed squared");
    return q;
}
HydroSystem::State HydroSystem::stateFromTemperature(units::Density rho, units::Temperature temperature, Real a, Real z) const {
    if (!helmholtz_) throw std::logic_error("Temperature initialization requires Helmholtz");
    auto const q = helmholtz_->at(units::value(rho), units::value(temperature), a,z);
    State state;
    state.density() = rho; state.nuclei() = rho/a; state.electrons() = rho*(z/a);
    state.totalEnergy() = rho*units::VelocitySquared::from_value(q.energy);
    if (dualEnergy_.enabled) state.auxiliary() = rho*(q.entropy / 1e8);
    return state;
}
units::Density HydroSystem::auxiliaryFromInternalEnergy(State const& state, units::EnergyDensity u) const {
    if (!helmholtz_) return auxiliaryFromInternalEnergy(state.density(),u);
    return state.density() * (thermodynamics(state,u).entropy / 1e8);
}
units::EnergyDensity HydroSystem::applyTemperatureFloor(State& state) const {
    if (!helmholtz_) return {};
    auto const low = minimum(state);
    auto const floor = state.density()*units::VelocitySquared::from_value(low.energy);
    auto const fromTotal = totalInternalEnergy(state);
    bool const usableAuxiliary = dualEnergy_.enabled && units::finite(state.auxiliary()) &&
        Real(state.auxiliary()/state.density())*1e8 > low.entropy + 1e-12*std::max(Real(1),std::abs(low.entropy));
    // Prefer an independently advected entropy if E-K has lost thermal accuracy.
    if (fromTotal >= floor) {
        if (dualEnergy_.enabled && !usableAuxiliary) state.auxiliary()=auxiliaryFromInternalEnergy(state,fromTotal);
        return {};
    }
    if (usableAuxiliary) return {};
    auto const old = state.totalEnergy();
    state.totalEnergy() += floor-fromTotal;
    if (dualEnergy_.enabled) state.auxiliary() = state.density()*(low.entropy/1e8);
    return state.totalEnergy()-old;
}
void HydroSystem::constrainComposition(State& candidate, State const& donor) const {
    if (!helmholtz_) return;
    candidate.nuclei()=donor.nuclei()*Real(candidate.density()/donor.density());
    candidate.electrons()=donor.electrons()*Real(candidate.density()/donor.density());
}
HydroSystem::Flux HydroSystem::compositionFlux(Flux flux, State const& left, State const& right) {
    auto const& donor = flux.mass() >= units::MassFlux{} ? left : right;
    flux.nuclei() = flux.mass()*Real(donor.nuclei()/donor.density());
    flux.electrons() = flux.mass()*Real(donor.electrons()/donor.density());
    return flux;
}

units::EnergyDensity HydroSystem::degenerateEnergy(units::Density rho) const {
	return degenerate_ ? whiteDwarfEos_.internalEnergy(rho) : units::EnergyDensity{};
}

units::Pressure HydroSystem::degeneratePressure(units::Density rho) const {
	return degenerate_ ? whiteDwarfEos_.pressure(rho) : units::Pressure{};
}

units::Density HydroSystem::auxiliaryFromInternalEnergy(units::Density rho, units::EnergyDensity u) const {
	if (helmholtz_) { State state; state.density()=rho; return auxiliaryFromInternalEnergy(state,u); }
	if (!(rho > units::Density{}) || !(u > units::EnergyDensity{}) || !units::finite(rho) || !units::finite(u))
		throw std::invalid_argument("Dual energy requires finite positive density and internal energy");
	long double const logRho = std::log(static_cast<long double>(units::value(rho)));
	long double const logU = std::log(static_cast<long double>(units::value(u)));
	return units::Density::from_value(positiveExp(logRho + static_cast<long double>(dualEnergy_.exponent) * (logU - adiabaticIndex_ * logRho)));
}

units::EnergyDensity HydroSystem::internalEnergyFromAuxiliary(State const& state) const {
    if (helmholtz_) {
        auto const [a,z] = composition(state);
        auto const q = helmholtz_->invert(units::value(state.density()), Real(state.auxiliary()/state.density())*1e8, a,z, HelmholtzClosure::Variable::Entropy);
        return state.density()*units::VelocitySquared::from_value(q.energy);
    }
	if (!(state.density() > units::Density{}) || !(state.auxiliary() > units::Density{}) ||
		!units::finite(state.density()) || !units::finite(state.auxiliary()))
		throw std::invalid_argument("Dual energy requires finite positive density and auxiliary density");
	long double const logRho = std::log(static_cast<long double>(units::value(state.density())));
	long double const logA = std::log(static_cast<long double>(units::value(state.auxiliary())));
	return units::EnergyDensity::from_value(positiveExp(adiabaticIndex_ * logRho + (logA - logRho) / dualEnergy_.exponent));
}

units::EnergyDensity HydroSystem::internalEnergy(State const& state) const {
	auto const thermal = totalInternalEnergy(state) - degenerateEnergy(state.density());
	auto const resolved = helmholtz_ ? thermal-state.density()*units::VelocitySquared::from_value(minimum(state).energy) : thermal;
	if (!dualEnergy_.enabled || (state.totalEnergy() > units::EnergyDensity{} && resolved > dualEnergy_.pressureThreshold * state.totalEnergy() &&
        (!helmholtz_ || thermal >= state.density()*units::VelocitySquared::from_value(minimum(state).energy))))
		return thermal;
	try {
		return internalEnergyFromAuxiliary(state);
	} catch (std::exception const&) {
		// Trial states in slope/flux limiters must fail admissibility rather
		// than aborting the search for a usable conservative update.
		return units::EnergyDensity::from_value(-std::numeric_limits<Real>::infinity());
	}
}

units::Temperature HydroSystem::temperature(State const& state) const {
	if (!admissible(state)) throw std::runtime_error("Cannot compute temperature of an inadmissible hydro state");
	if (helmholtz_) return units::Temperature::from_value(thermodynamics(state,internalEnergy(state)).temperature);
	return (adiabaticIndex_ - 1) * internalEnergy(state) *
		(meanMolecularWeight_ * constants::atomicMassUnit) / (state.density() * constants::boltzmann);
}

void HydroSystem::synchronize(State& state) const {
    if (helmholtz_) {
        if (!dualEnergy_.enabled) return;
        auto const u = totalInternalEnergy(state);
        if (state.totalEnergy() > units::EnergyDensity{} && u-state.density()*units::VelocitySquared::from_value(minimum(state).energy) > dualEnergy_.syncThreshold*state.totalEnergy() &&
            u >= state.density()*units::VelocitySquared::from_value(minimum(state).energy))
            state.auxiliary()=auxiliaryFromInternalEnergy(state,u);
        return;
    }
	if (!dualEnergy_.enabled) return;
	auto const thermal = totalInternalEnergy(state) - degenerateEnergy(state.density());
	if (state.totalEnergy() > units::EnergyDensity{} && units::finite(thermal) &&
		thermal > dualEnergy_.syncThreshold * state.totalEnergy() &&
		degeneratePressure(state.density()) + (adiabaticIndex_ - 1) * thermal >= pressureFloor_)
		state.auxiliary() = auxiliaryFromInternalEnergy(state.density(), thermal);
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
	result.auxiliary() = state.auxiliary() / state.density();
	result.nuclei() = state.nuclei()/state.density();
	result.electrons() = state.electrons()/state.density();
	return result;
}

ConservedState HydroSystem::conservedState(PrimitiveState const& state) const {
	if (!finite(state) || !(state.density() >= densityFloor_) || !(state.pressure() >= pressureFloor_) || !units::finite(state.density()) || !units::finite(state.pressure())) {
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
    result.nuclei() = state.density()*state.nuclei();
    result.electrons() = state.density()*state.electrons();
    if (helmholtz_) {
        auto const [a,z] = composition(result);
        result.nuclei() = result.density()/a; result.electrons() = result.density()*(z/a);
        auto const q = helmholtz_->invert(units::value(state.density()), units::value(state.pressure()), a,z, HelmholtzClosure::Variable::Pressure, false);
        result.totalEnergy() = result.density()*units::VelocitySquared::from_value(q.energy) + Real(0.5)*state.density()*speedSquared;
        if (dualEnergy_.enabled) result.auxiliary() = state.auxiliary() == units::Dimensionless{} ?
            result.density()*(q.entropy/1e8) : state.density()*state.auxiliary();
        return result;
    }
	auto const thermalPressure = state.pressure() - degeneratePressure(state.density());
	if (!(thermalPressure > units::Pressure{}))
		throw std::invalid_argument("Hydro primitive pressure must exceed cold degeneracy pressure");
	auto const thermal = thermalPressure / (adiabaticIndex_ - 1);
	result.setTotalEnergy(degenerateEnergy(state.density()) + thermal + Real(0.5) * state.density() * speedSquared);
	if (!units::finite(state.auxiliary()) || state.auxiliary() < units::Dimensionless{})
		throw std::invalid_argument("Primitive auxiliary entropy must be finite and nonnegative");
	// Zero is an initialization marker only in primitives, never in an enabled
	// conserved state. Reconstruction retains A/rho independently of pressure.
	if (dualEnergy_.enabled)
		result.auxiliary() = state.auxiliary() == units::Dimensionless{} ?
			auxiliaryFromInternalEnergy(state.density(), thermal) : state.density() * state.auxiliary();
	return result;
}

ConservedFlux HydroSystem::physicalFlux(ConservedState const& state, int normal, units::Velocity faceSpeed) const {
	PrimitiveState const primitive = reconstructionVariables(state);
	auto const normalVelocity = primitive.velocity(normal);
	ConservedFlux result;
	result.setMass(state.density() * normalVelocity);
	result.auxiliary() = state.auxiliary() * normalVelocity;
	result.nuclei() = state.nuclei()*normalVelocity;
	result.electrons() = state.electrons()*normalVelocity;
	for (int axis = 0; axis < ndim; ++axis) {
		result.setMomentum(axis, state.momentum(axis) * normalVelocity + (axis == normal ? primitive.pressure() : units::Pressure{}));
	}
	result.setEnergy((state.totalEnergy() + primitive.pressure()) * normalVelocity);
	if (faceSpeed != units::Velocity{}) result -= advectiveFlux(state, faceSpeed);
	return result;
}

ConservedFlux HydroSystem::riemann(ConservedState const& left, ConservedState const& right, int normal, units::Velocity faceSpeed) const {
	if (faceSpeed != units::Velocity{}) {
		// Solve in the face's translating normal frame, then transform the flux
		// back to inertial momentum and energy. The pressure work is p*v_n.
		auto boost = [&](ConservedState value) {
			value.totalEnergy() += -faceSpeed * value.momentum(normal) + Real(0.5) * value.density() * faceSpeed * faceSpeed;
			value.momentum(normal) -= value.density() * faceSpeed;
			return value;
		};
		auto result = riemann(boost(left), boost(right), normal);
		result.energy() += faceSpeed * result.momentum(normal) + Real(0.5) * faceSpeed * faceSpeed * result.mass();
		result.momentum(normal) += faceSpeed * result.mass();
		return result;
	}
	PrimitiveState const leftPrimitive = reconstructionVariables(left);
	PrimitiveState const rightPrimitive = reconstructionVariables(right);
	auto const leftSound = soundSpeed(leftPrimitive);
	auto const rightSound = soundSpeed(rightPrimitive);
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
		star.auxiliary() = state.auxiliary() * (star.density() / state.density());
		star.nuclei() = state.nuclei()*(star.density()/state.density());
		star.electrons() = state.electrons()*(star.density()/state.density());
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

ConservedState HydroSystem::outflow(State state, int normal, bool lower, units::Velocity faceSpeed) const {
	auto& momentum = state.momentum(normal);
	auto const meshMomentum = state.density() * faceSpeed;
	if (lower ? momentum > meshMomentum : momentum < meshMomentum) {
		if (faceSpeed != units::Velocity{})
			state.totalEnergy() += Real(0.5) * (meshMomentum * meshMomentum - momentum * momentum) / state.density();
		momentum = meshMomentum;
	}
	return state;
}

units::Velocity HydroSystem::maximumSignalSpeed(ConservedState const& state, int normal, units::Velocity faceSpeed) const {
	PrimitiveState const primitive = reconstructionVariables(state);
	return units::abs(primitive.velocity(normal) - faceSpeed) + soundSpeed(primitive);
}

units::Velocity HydroSystem::adiabaticSoundSpeed(ConservedState const& state) const {
	return soundSpeed(reconstructionVariables(state));
}

units::Velocity HydroSystem::soundSpeed(PrimitiveState const& primitive) const {
    if (helmholtz_) {
        State state; state.density()=primitive.density();
        state.nuclei()=primitive.density()*primitive.nuclei(); state.electrons()=primitive.density()*primitive.electrons();
        auto const [a,z] = composition(state);
        auto const q = helmholtz_->invert(units::value(primitive.density()),units::value(primitive.pressure()),a,z,HelmholtzClosure::Variable::Pressure);
        if (!(q.soundSquared>0) || !std::isfinite(q.soundSquared)) throw std::runtime_error("Invalid Helmholtz sound speed");
        return units::Velocity::from_value(std::sqrt(q.soundSquared));
    }
	auto const thermalPressure = primitive.pressure() - degeneratePressure(primitive.density());
	auto const coldDerivative = degenerate_ ? whiteDwarfEos_.pressureDerivative(primitive.density()) : units::VelocitySquared{};
	return units::sqrt(coldDerivative + adiabaticIndex_*thermalPressure/primitive.density());
}

bool HydroSystem::admissible(ConservedState const& state) const {
    if (helmholtz_) {
        if (!finite(state) || !(state.density()>=densityFloor_)) return false;
        try { (void)thermodynamics(state,internalEnergy(state)); return true; }
        catch (std::exception const&) { return false; }
    }
	if (!finite(state)) return false;
	if (!(state.density() >= densityFloor_)) return false;
	// Conservative gravity work can undershoot E in a dilute atmosphere.
	// With dual energy, thermodynamics comes from the positive entropy auxiliary
	// in that regime; do not clip E and destroy the gas-plus-gravity balance.
	if (!dualEnergy_.enabled && !(state.totalEnergy() > units::EnergyDensity{})) return false;
	if (dualEnergy_.enabled && !(state.auxiliary() > units::Density{})) return false;
	if (!(internalEnergy(state) > units::EnergyDensity{})) return false;
	auto const p = pressure(state);
	return units::finite(p) && p >= pressureFloor_;
}

ConservedState HydroSystem::correctRoundoff(ConservedState state, State const& updateScale) const {
	(void) updateScale;
	return state;
}

ConservedFlux HydroSystem::limitFlux(
	ConservedState const& left, ConservedState const& right, ConservedFlux const& highOrderFlux, int normal, units::TimePerLength stepOverCellWidth, units::Velocity faceSpeed) const {
	if (stepOverCellWidth == units::TimePerLength{}) {
		return highOrderFlux;
	}
	ConservedFlux const leftPhysical = physicalFlux(left, normal, faceSpeed);
	ConservedFlux const rightPhysical = physicalFlux(right, normal, faceSpeed);
	auto const factor = Real(2 * ndim) * stepOverCellWidth;
	auto validFlux = [&](ConservedFlux const& flux) {
		auto const consistent = compositionFlux(flux,left,right);
		return admissible(left - integratedFlux(consistent - leftPhysical, factor)) && admissible(right + integratedFlux(consistent - rightPhysical, factor));
	};
	if (validFlux(highOrderFlux)) {
		return highOrderFlux;
	}
	auto const speed = std::max(maximumSignalSpeed(left, normal, faceSpeed), maximumSignalSpeed(right, normal, faceSpeed));
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
	if (helmholtz_) return units::Pressure::from_value(thermodynamics(state,internalEnergy(state)).pressure);
	return degeneratePressure(state.density()) + (adiabaticIndex_ - 1) * internalEnergy(state);
}

units::EnergyDensity HydroSystem::totalInternalEnergy(ConservedState const& state) const {
	if (!(state.density() > units::Density{}) || !units::finite(state.density())) {
		return units::Pressure::from_value(-std::numeric_limits<Real>::infinity());
	}
	decltype(units::MomentumDensity{} * units::MomentumDensity{}) momentumSquared{};
	for (int axis = 0; axis < ndim; ++axis) {
		momentumSquared += state.momentum(axis) * state.momentum(axis);
	}
	return state.totalEnergy() - Real(0.5) * momentumSquared / state.density();
}

ConservedFlux HydroSystem::hll(ConservedState const& left, ConservedState const& right, int normal) const {
	PrimitiveState const leftPrimitive = reconstructionVariables(left);
	PrimitiveState const rightPrimitive = reconstructionVariables(right);
	auto const leftSound = soundSpeed(leftPrimitive);
	auto const rightSound = soundSpeed(rightPrimitive);
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
	result.auxiliary() = state.auxiliary() * speed;
	result.nuclei() = state.nuclei()*speed;
	result.electrons() = state.electrons()*speed;
	return result;
}

ConservedState HydroSystem::integratedFlux(ConservedFlux const& flux, units::TimePerLength factor) {
	ConservedState result;
	result.setDensity(flux.mass() * factor);
	for (int axis = 0; axis < ndim; ++axis)
		result.setMomentum(axis, flux.momentum(axis) * factor);
	result.setTotalEnergy(flux.energy() * factor);
	result.auxiliary() = flux.auxiliary() * factor;
	result.nuclei() = flux.nuclei()*factor;
	result.electrons() = flux.electrons()*factor;
	return result;
}
}	 // namespace octotigerII::hydro
