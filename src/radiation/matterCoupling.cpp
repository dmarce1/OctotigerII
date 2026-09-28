// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#include "octotigerII/radiation/matterCoupling.hpp"
#include "octotigerII/math/FpeGuard.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace octotigerII::radiation {
namespace {
using Number = long double;
constexpr int count = ndim + 1;
using Vector = std::array<Number, count>;
using Matrix = std::array<Vector, count>;
constexpr Number gamma = 0.2928932188134524755991556378951509607L;

Number norm(Vector const& x) {
	Number value = 0;
	for (auto v : x) {
		if (!std::isfinite(v)) return std::numeric_limits<Number>::infinity();
		value = std::max(value, std::abs(v));
	}
	return value;
}

// All numerical variables below are dimensionless (E,F/c)/scale. Conversion
// from/to physical CGS quantities is confined to the public interface.
struct Problem {
	Number rho, finalRho, c, chat, scale, thermal, finalThermal, temperaturePerEnergy, opticalInterval, referenceOpacity;
	OpacityLaw opacities;
	std::array<Number, ndim> velocity, finalVelocity;
	Vector initial, drive;
    hydro::HydroSystem const* system = nullptr;
    hydro::ConservedState material, finalMaterial;
    hydro::ConservedState materialAt(Number time) const {
        auto q = material;
        q.density() = units::Density::from_value(Real(density(time)));
        q.nuclei() = material.nuclei() + Real(time)*(finalMaterial.nuclei()-material.nuclei());
        q.electrons() = material.electrons() + Real(time)*(finalMaterial.electrons()-material.electrons());
        return q;
    }
    Number tolerance() const { return system->helmholtz() ? 2e-12L : 256 * std::numeric_limits<Number>::epsilon(); }

	Number density(Number time) const { return rho + time * (finalRho - rho); }
	Number baselineVelocity(int d, Number time) const {
		return ((1 - time) * rho * velocity[d] + time * finalRho * finalVelocity[d]) / density(time);
	}
	Number thermalEnergy(Vector const& x, Number time) const {
		Number const currentRho = density(time);
		Number result = thermal + time * (finalThermal - thermal)
			- (c / chat) * scale * (x[0] - initial[0] - time * drive[0]);
		for (int d = 0; d < ndim; ++d) {
			Number const dv = finalVelocity[d] - velocity[d];
			// Stable convex kinetic-energy remainder for linearly driven mass
			// and momentum; equivalent to linearly interpolating the thermal
			// discrepancy from the conservative total energy at the endpoints.
			result += time * (1 - time) * rho * finalRho / (2 * currentRho) * dv * dv;
			Number const impulse = -scale * (x[d + 1] - initial[d + 1] - time * drive[d + 1]) / chat;
			// Difference of kinetic energies without subtracting large squares.
			result -= impulse * (baselineVelocity(d, time) + impulse / (2 * currentRho));
		}
		return result;
	}

	bool admissible(Vector const& x, Number time) const {
		for (auto v : x) { if (!std::isfinite(v)) return false; }
		if (!(x[0] >= 0)) return false;
        if (system->helmholtz()) {
            try { (void)system->thermodynamics(materialAt(time),units::EnergyDensity::from_value(Real(thermalEnergy(x,time)))); }
            catch (std::exception const&) { return false; }
        } else if (!(thermalEnergy(x,time)>0)) return false;
		Number magnitude = 0;
		for (int d = 1; d < count; ++d) { magnitude = std::hypot(magnitude, x[d]); }
		return magnitude <= x[0] * (1 + 64 * std::numeric_limits<Real>::epsilon());
	}

	Vector rate(Vector const& x, Number time, Matrix* derivative = nullptr) const {
		std::array<Number, ndim> beta{}, f{};
		Number f2 = 0, betaQ = 0, betaF = 0;
		for (int d = 0; d < ndim; ++d) {
			Number const impulse = -scale * (x[d + 1] - initial[d + 1] - time * drive[d + 1]) / chat;
			beta[d] = (baselineVelocity(d, time) + impulse / density(time)) / c;
			if (x[0] > 0) f[d] = x[d + 1] / x[0];
			f2 += f[d] * f[d];
			betaQ += beta[d] * x[d + 1];
			betaF += beta[d] * f[d];
		}
		Number const root = std::sqrt(4 - 3 * f2);
		Number const isotropic = x[0] * (1 - f2) / (root + 1);
		Number const directed = x[0] * 3 / (root + 2);
		Number const internal = thermalEnergy(x, time);
		Number temperature = internal * temperaturePerEnergy * rho / density(time);
        Number temperatureDerivative = temperaturePerEnergy * rho / density(time);
        if (system->helmholtz()) {
            auto const q = system->thermodynamics(materialAt(time), units::EnergyDensity::from_value(Real(internal)));
            temperature = q.temperature;
            temperatureDerivative = 1 / (density(time)*q.cv);
            // Clamped trial states emit at the floor. The accepted floor pass
            // accounts for any energy added; the radiation exchange stays paired.
            if (q.temperature == system->temperatureFloor() && internal < density(time)*q.energy) temperatureDerivative = 0;
        }
		Number const temperature2 = temperature * temperature;
		Number const emission = Number(units::value(constants::radiation)) * temperature2 * temperature2 / scale;
		auto const state = materialAt(time);
		Real const ye = state.electrons() > units::Density{} ? Real(state.electrons() / state.density()) : -1;
		auto const coefficients = opacities.evaluate(Real(density(time)), Real(temperature), ye);
		Number const absorption = units::value(coefficients.planckAbsorption) / referenceOpacity;
		Number const diffusionExtinction = units::value(coefficients.fluxExtinction) / referenceOpacity;
		Number const streamingCorrection = std::max(Number(0), absorption - diffusionExtinction);
		// The Rosseland flux mean applies in diffusion. At a free-streaming
		// cone, using it with a much larger Planck energy mean would remove
		// energy faster than momentum and drive |F| beyond cE. Recover the
		// Planck flux scale continuously as the reduced flux approaches one.
		Number const extinction = diffusionExtinction + streamingCorrection * f2;
		Number const temperatureSlope = opacities.ionizedGas ? -3.5L / temperature : 0;
		Number const absorptionSlope = temperatureSlope * absorption;
		Number const diffusionSlope = temperatureSlope * Number(units::value(coefficients.rosselandAbsorption)) / referenceOpacity;
		Vector result{};
		result[0] = absorption * (emission - x[0]) + (2 * absorption - extinction) * betaQ;
		// Mixed-frame first-order force with kappa_E approximated by the
		// Planck absorption mean and kappa_F by the total flux extinction.
		for (int d = 0; d < ndim; ++d) {
			result[d + 1] = -extinction * x[d + 1]
				+ ((extinction - absorption) * x[0] + absorption * emission + extinction * isotropic) * beta[d]
				+ extinction * directed * f[d] * betaF;
		}
		for (auto& value : result) { value *= density(time) / rho; }
		if (derivative) {
			// Analytic derivatives keep Newton directions tangent to a cold
			// streaming beam. A finite-difference closure derivative can point
			// outside the cone even when the exact source preserves it.
			Number const energyConversion = scale * c / chat;
			Number const betaConversion = scale / (chat * density(time) * c);
			for (int column = 0; column < count; ++column) {
				Number const dEnergy = column == 0 ? 1 : 0;
				Number const dInternal = column == 0 ? -energyConversion : energyConversion * beta[column - 1];
				Number const dEmission = 4 * emission * temperatureDerivative * dInternal / temperature;
				Number const dAbsorption = absorptionSlope * temperatureDerivative * dInternal;
				Number const dDiffusion = diffusionSlope * temperatureDerivative * dInternal;
				std::array<Number, ndim> dBeta{}, dFluxFactor{};
				Number dFluxSquared = 0, dBetaFluxFactor = 0, dBetaFlux = 0;
				for (int d = 0; d < ndim; ++d) {
					Number const dFlux = column == d + 1 ? 1 : 0;
					dBeta[d] = -betaConversion * dFlux;
					if (x[0] > 0) dFluxFactor[d] = (dFlux - f[d] * dEnergy) / x[0];
					dFluxSquared += 2 * f[d] * dFluxFactor[d];
					dBetaFluxFactor += dBeta[d] * f[d] + beta[d] * dFluxFactor[d];
					dBetaFlux += dBeta[d] * x[d + 1] + beta[d] * dFlux;
				}
				Number const dExtinction = dDiffusion + (streamingCorrection > 0
					? (dAbsorption - dDiffusion) * f2 + streamingCorrection * dFluxSquared : 0);
				Number const dRoot = -1.5L * dFluxSquared / root;
				Number const dIsotropic = (dEnergy * (1 - f2) - x[0] * dFluxSquared - isotropic * dRoot) / (root + 1);
				Number const dDirected = (3 * dEnergy - directed * dRoot) / (root + 2);
				(*derivative)[0][column] = dAbsorption * (emission - x[0])
					+ absorption * (dEmission - dEnergy)
					+ (2 * dAbsorption - dExtinction) * betaQ
					+ (2 * absorption - extinction) * dBetaFlux;
				for (int d = 0; d < ndim; ++d) {
					Number const dFlux = column == d + 1 ? 1 : 0;
					Number const carried = (extinction - absorption) * x[0] + absorption * emission + extinction * isotropic;
					Number const dCarried = (dExtinction - dAbsorption) * x[0]
						+ (extinction - absorption) * dEnergy + dAbsorption * emission
						+ absorption * dEmission + dExtinction * isotropic + extinction * dIsotropic;
					(*derivative)[d + 1][column] = -dExtinction * x[d + 1] - extinction * dFlux
						+ dCarried * beta[d] + carried * dBeta[d]
						+ (dExtinction * directed + extinction * dDirected) * f[d] * betaF
						+ extinction * directed * (dFluxFactor[d] * betaF + f[d] * dBetaFluxFactor);
				}
				for (int row = 0; row < count; ++row) { (*derivative)[row][column] *= density(time) / rho; }
			}
		}
		return result;
	}
};

bool solveLinear(Matrix a, Vector b, Vector& x) {
	for (int k = 0; k < count; ++k) {
		int pivot = k;
		for (int i = k + 1; i < count; ++i) { if (std::abs(a[i][k]) > std::abs(a[pivot][k])) pivot = i; }
		if (!(std::abs(a[pivot][k]) > std::numeric_limits<Number>::min()) || !std::isfinite(a[pivot][k])) return false;
		std::swap(a[pivot], a[k]);
		std::swap(b[pivot], b[k]);
		for (int i = k + 1; i < count; ++i) {
			Number const factor = a[i][k] / a[k][k];
			for (int j = k + 1; j < count; ++j) { a[i][j] -= factor * a[k][j]; }
			b[i] -= factor * b[k];
		}
	}
	for (int i = count - 1; i >= 0; --i) {
		Number rhs = b[i];
		for (int j = i + 1; j < count; ++j) { rhs -= a[i][j] * x[j]; }
		x[i] = rhs / a[i][i];
		if (!std::isfinite(x[i])) return false;
	}
	return true;
}

bool implicitStage(Problem const& problem, Vector const& base, Number duration, Number time, Vector& state) {
	Number const step = duration * problem.opticalInterval;
	Number const inverse = 1 / (1 + step), weight = step / (1 + step);
	auto residual = [&](Vector const& value, Vector* source = nullptr) {
		auto const rhs = problem.rate(value, time);
		if (source) *source = rhs;
		Vector result{};
		for (int j = 0; j < count; ++j) { result[j] = inverse * (value[j] - base[j] - duration * problem.drive[j]) - weight * rhs[j]; }
		return result;
	};
	// Holding radiation at its old value after a large incoming transport drive
	// assigns that energy to trial gas heat through conservation. The enormous
	// trial emission can then stall cone-constrained line searches. Try including
	// transport in the starting point, but preserve an already-balanced guess
	// with a smaller residual, including its accurately represented tiny fluxes.
	// This changes only the initial iterate, never the implicit stage equation.
	auto transported = base;
	for (int j = 0; j < count; ++j) { transported[j] += duration * problem.drive[j]; }
	if (problem.admissible(transported, time) && (!problem.admissible(state, time)
		|| norm(residual(transported)) < norm(residual(state)))) state = transported;
	if (!problem.admissible(state,time)) return false;
	for (int iteration = 0; iteration < 64; ++iteration) {
		Vector source{};
		auto const r = residual(state, &source);
		Number const error = norm(r);
		if (!std::isfinite(error)) return false;
		// Scale by the terms in the radiation equations, not the gas energy
		// used to normalize variables. An absolute dimensionless tolerance would
		// erase a weak radiation field in a cold, gas-dominated cell. Include
		// the separate source magnitudes so stiff equilibrium cancellation does
		// not demand precision beyond the represented radiation state.
		Number const residualScale = inverse * (norm(state) + norm(base) + duration * norm(problem.drive))
			+ weight * (problem.density(time) / problem.rho * norm(state) + norm(source));
		if (error <= problem.tolerance() * residualScale) return problem.admissible(state, time);
		Matrix jacobian{};
		problem.rate(state, time, &jacobian);
		for (int i = 0; i < count; ++i) { for (int j = 0; j < count; ++j) {
			jacobian[i][j] = (i == j ? inverse : 0) - weight * jacobian[i][j];
		}
		}
		Vector rhs{}, update{};
		for (int j = 0; j < count; ++j) { rhs[j] = -r[j]; }
		if (!solveLinear(jacobian, rhs, update)) return false;
		// In a radiation-dominated equilibrium, reconstructing gas heat from
		// conserved energy makes the residual sensitive to tiny radiation
		// roundoff. Also accept a resolved Newton correction; its scale is the
		// radiation state itself, so this introduces no gas-relative floor.
		if (norm(update) <= problem.tolerance() * norm(state)
			&& problem.admissible(state, time)) return true;
		bool accepted = false;
		for (Number fraction = 1; fraction >= std::ldexp(Number(1), -40); fraction /= 2) {
			auto candidate = state;
			for (int j = 0; j < count; ++j) { candidate[j] += fraction * update[j]; }
			if (!problem.admissible(candidate, time)) continue;
			Number const nextError = norm(residual(candidate));
			if (std::isfinite(nextError) && nextError < error) { state = candidate; accepted = true; break; }
		}
		if (!accepted) return false;
	}
	return false;
}

Vector advance(Problem const& problem, Vector const& old, Number time, Number interval, int depth, unsigned& attempts) {
	if (++attempts > 4096) throw std::runtime_error("Radiation matter exchange exceeded its source subdivision budget");
	auto stage = old;
	if (implicitStage(problem, old, gamma * interval, time + gamma * interval, stage)) {
		Vector base{};
		// Stage-one identity avoids multiplying a cancellation-sized residual
		// by a very large opacity in a stiff second stage.
		for (int j = 0; j < count; ++j) { base[j] = old[j] + ((1 - gamma) / gamma) * (stage[j] - old[j]); }
		auto result = stage;
		if (implicitStage(problem, base, gamma * interval, time + interval, result)) return result;
	}
	if (depth == 24 || !(interval / 2 > 0)) throw std::runtime_error("Radiation matter exchange has no converged admissible source step");
	auto const middle = advance(problem, old, time, interval / 2, depth + 1, attempts);
	return advance(problem, middle, time + interval / 2, interval / 2, depth + 1, attempts);
}
} // namespace

void couple(hydro::ConservedState& gas, RadiationSystem::State& radiation,
	hydro::HydroSystem const& system, Opacity opacity, Real ratio, units::Time interval) {
	coupleForcedWithOpacityLaw(gas, radiation, {}, {}, system, OpacityLaw{false, opacity}, ratio, interval);
}

void coupleWithOpacityLaw(hydro::ConservedState& gas, RadiationSystem::State& radiation,
	hydro::HydroSystem const& system, OpacityLaw const& opacities, Real ratio, units::Time interval) {
	coupleForcedWithOpacityLaw(gas, radiation, {}, {}, system, opacities, ratio, interval);
}

void coupleForced(hydro::ConservedState& gas, RadiationSystem::State& radiation,
	hydro::ConservedState const& gasIncrement, RadiationSystem::State const& radiationIncrement,
	hydro::HydroSystem const& system, Opacity opacity, Real ratio, units::Time interval) {
	coupleForcedWithOpacityLaw(gas, radiation, gasIncrement, radiationIncrement, system, OpacityLaw{false, opacity}, ratio, interval);
}

void coupleForcedWithOpacityLaw(hydro::ConservedState& gas, RadiationSystem::State& radiation,
	hydro::ConservedState const& gasIncrement, RadiationSystem::State const& radiationIncrement,
	hydro::HydroSystem const& system, OpacityLaw const& opacities, Real ratio, units::Time interval) {
	if (!units::finite(opacities.constantAbsorption) || opacities.constantAbsorption < Opacity{}
		|| !units::finite(opacities.constantScattering) || opacities.constantScattering < Opacity{}
		|| !std::isfinite(ratio) || !(ratio > 0 && ratio <= 1)
		|| !units::finite(interval) || interval < units::Time{}) throw std::invalid_argument("Invalid radiation matter-exchange parameters");
	if (interval == units::Time{}) return;
	RadiationSystem transport(ratio * constants::c);
	if (!system.admissible(gas) || !transport.admissible(radiation)) throw std::invalid_argument("Inadmissible radiation matter-exchange input");
	auto const drivenGas = hydro::ConservedState(gas + gasIncrement);
	auto const drivenRadiation = RadiationSystem::State(radiation + radiationIncrement);
	if (!finite(gasIncrement) || !finite(radiationIncrement) || !(drivenGas.density() > units::Density{}))
		throw std::invalid_argument("Invalid radiation matter-exchange transport drive");
	if (!opacities.ionizedGas && opacities.constantAbsorption == Opacity{} && opacities.constantScattering == Opacity{}) {
		if (!system.admissible(drivenGas) || !transport.admissible(drivenRadiation)) throw std::runtime_error("Inadmissible uncoupled transport drive");
		gas = drivenGas; radiation = drivenRadiation; return;
	}
	FpeGuard guard;
	Problem problem{};
	if (!system.gasOnly()) throw std::invalid_argument("Radiation coupling requires the gas-only Helmholtz closure");
	problem.system = &system; problem.material = gas; problem.finalMaterial = drivenGas;
	problem.opacities = opacities;
	problem.rho = units::value(gas.density());
	problem.finalRho = units::value(drivenGas.density());
	problem.c = units::value(constants::c);
	problem.chat = ratio * problem.c;
	problem.thermal = units::value(system.internalEnergy(gas));
	problem.finalThermal = units::value(system.internalEnergy(drivenGas));
	if (!std::isfinite(problem.finalThermal)) throw std::runtime_error("Invalid driven gas thermal energy");
	problem.temperaturePerEnergy = units::value(system.temperature(gas)) / problem.thermal;
	problem.scale = units::value(radiation.energy()) + ratio * problem.thermal;
	if (!(problem.scale > 0) || !std::isfinite(problem.scale)) throw std::runtime_error("Invalid radiation matter-exchange energy scale");
	problem.initial[0] = units::value(radiation.energy()) / problem.scale;
	problem.drive[0] = units::value(radiationIncrement.energy()) / problem.scale;
	for (int d = 0; d < ndim; ++d) {
		problem.velocity[d] = units::value(gas.momentum(d)) / problem.rho;
		problem.finalVelocity[d] = units::value(drivenGas.momentum(d)) / problem.finalRho;
		problem.initial[d + 1] = units::value(radiation.radiativeFlux(d)) / (problem.c * problem.scale);
		problem.drive[d + 1] = units::value(radiationIncrement.radiativeFlux(d)) / (problem.c * problem.scale);
	}
	auto const initialOpacity = opacities.evaluate(gas, system);
	problem.referenceOpacity = std::max(Number(units::value(initialOpacity.planckAbsorption)),
		Number(units::value(initialOpacity.fluxExtinction)));
	if (!(problem.referenceOpacity > 0)) throw std::runtime_error("Analytic radiation opacity vanished at the input state");
	problem.opticalInterval = Number(units::value(interval)) * problem.chat * problem.rho * problem.referenceOpacity;
	if (!std::isfinite(problem.opticalInterval)) throw std::runtime_error("Unrepresentable radiation matter-exchange interval");
	unsigned attempts = 0;
	Vector result;
	try {
		result = advance(problem, problem.initial, 0, 1, 0, attempts);
	} catch (std::runtime_error const& error) {
		std::ostringstream details;
		details.precision(17);
		details << error.what() << ": dt=" << units::value(interval) << " absorption="
			<< units::value(initialOpacity.planckAbsorption) << " extinction=" << units::value(initialOpacity.fluxExtinction) << " ratio=" << ratio;
		auto append = [&](char const* name, auto const& state) {
			details << ' ' << name << "=[";
			bool first = true;
			state.forEach([&](auto, auto value) { if (!first) details << ','; first = false; details << units::value(value); });
			details << ']';
		};
		append("gas", gas); append("radiation", radiation);
		append("gasDrive", gasIncrement); append("radiationDrive", radiationIncrement);
		throw std::runtime_error(details.str());
	}
	auto nextGas = drivenGas;
	auto nextRadiation = radiation;
	nextRadiation.energy() = units::EnergyDensity::from_value(Real(result[0] * problem.scale));
	// Use the actual representable radiation changes, not separately rounded
	// source evaluations. Opposite gas increments then cancel to roundoff.
	auto const energyChange = nextRadiation.energy() - radiation.energy() - radiationIncrement.energy();
	nextGas.totalEnergy() -= energyChange / ratio;
	Number thermal = problem.finalThermal - Number(units::value(energyChange)) / ratio;
	for (int d = 0; d < ndim; ++d) {
		nextRadiation.radiativeFlux(d) = units::EnergyFlux::from_value(Real(result[d + 1] * problem.scale * problem.c));
		auto const impulse = -(nextRadiation.radiativeFlux(d) - radiation.radiativeFlux(d) - radiationIncrement.radiativeFlux(d)) / (ratio * constants::c * constants::c);
		nextGas.momentum(d) += impulse;
		Number const dp = units::value(impulse);
		thermal -= dp * (problem.finalVelocity[d] + dp / (2 * problem.finalRho));
	}
	if ((!system.helmholtz() && !(thermal > 0)) || !std::isfinite(thermal)) throw std::runtime_error("Radiation matter exchange exhausted gas thermal energy");
	nextGas.auxiliary() = system.auxiliaryFromInternalEnergy(nextGas, units::EnergyDensity::from_value(Real(thermal)));
	if (!system.admissible(nextGas) || !transport.admissible(nextRadiation)) throw std::runtime_error("Radiation matter exchange produced an inadmissible state");
	gas = nextGas;
	radiation = nextRadiation;
}
} // namespace octotigerII::radiation
