#include "testSupport.hpp"
#include "octotigerII/radiation/matterCoupling.hpp"
#include <cmath>
#include <iostream>
#include <limits>

using namespace octotigerII;
using namespace octotigerII::radiation;

namespace {
constexpr Real density = 1e-7, temperature = 1e5, opacityValue = 0.4;

hydro::ConservedState gasState(hydro::HydroSystem const& system, Real rho = density, Real t = temperature, Real velocity = 0) {
	hydro::PrimitiveState primitive{};
	primitive.density() = units::Density::from_value(rho);
	primitive.pressure() = primitive.density() * constants::boltzmann * units::Temperature::from_value(t) / constants::atomicMassUnit;
	primitive.velocity(0) = units::Velocity::from_value(velocity);
	return system.conservedState(primitive);
}
RadiationSystem::State radiationState(Real factor = 1) {
	RadiationSystem::State state{};
	state.energy() = units::EnergyDensity::from_value(factor * units::value(constants::radiation) * std::pow(temperature, 4));
	return state;
}
units::Time interval(Real opticalTime, Real ratio = 1) {
	return units::Time::from_value(opticalTime / (ratio * units::value(constants::c) * density * opacityValue));
}
void invariants(hydro::ConservedState const& beforeGas, RadiationSystem::State const& beforeRad,
	hydro::ConservedState const& afterGas, RadiationSystem::State const& afterRad, Real ratio) {
	auto const energy = beforeGas.totalEnergy() + beforeRad.energy() / ratio;
	EXPECT_NEAR(units::value(afterGas.totalEnergy() + afterRad.energy() / ratio - energy), 0,
		32 * epsilonR * units::value(energy));
	for (int d = 0; d < ndim; ++d) {
		auto const old = beforeGas.momentum(d) + beforeRad.radiativeFlux(d) / (ratio * constants::c * constants::c);
		auto const next = afterGas.momentum(d) + afterRad.radiativeFlux(d) / (ratio * constants::c * constants::c);
		auto const scale = units::abs(beforeGas.momentum(d)) + units::abs(afterGas.momentum(d))
			+ (units::abs(beforeRad.radiativeFlux(d)) + units::abs(afterRad.radiativeFlux(d))) / (ratio * constants::c * constants::c);
		EXPECT_NEAR(units::value(next - old), 0, 32 * epsilonR * units::value(scale));
	}
}

TEST(RadiationCoupling, ThermalEquilibriumRemainsStationaryForStiffSteps) {
	hydro::HydroSystem system;
	for (Real ratio : {1.0, 0.125}) { for (Real depth : {1e-6, 1.0, 1e9}) {
		auto gas = gasState(system);
		auto rad = radiationState();
		auto const beforeGas = gas;
		auto const beforeRad = rad;
		couple(gas, rad, system, Opacity::from_value(opacityValue), ratio, interval(depth, ratio));
		test::expectStateNear(gas, beforeGas, 2e-13);
		test::expectStateNear(rad, beforeRad, 2e-13);
		invariants(beforeGas, beforeRad, gas, rad, ratio);
	}
	}
}

TEST(RadiationCoupling, AnalyticMeansSeparateThermalAbsorptionAndFluxExtinction) {
	OpacityLaw law;
	law.ionizedGas = true;
	law.hydrogenFraction = 0.7;
	law.metalFraction = 0.02;
	auto const means = law.evaluate(1e-7, 1e5);
	Real const rosseland = 3.68e22 * 1.7 * 0.98 * 1e-7 * std::pow(1e5, -3.5);
	EXPECT_NEAR(units::value(means.planckAbsorption), 37 * rosseland, 1e-13 * 37 * rosseland);
	EXPECT_NEAR(units::value(means.fluxExtinction), rosseland + 0.34, 1e-13);
	EXPECT_NEAR(units::value(means.rosselandAbsorption), rosseland, 1e-13 * rosseland);
	EXPECT_NEAR(units::value(means.scattering), 0.34, 1e-13);
	EXPECT_NEAR(units::value(law.evaluate(2e-7, 1e5).planckAbsorption),
		2 * units::value(means.planckAbsorption), 1e-13 * units::value(means.planckAbsorption));
	EXPECT_NEAR(units::value(law.evaluate(1e-7, 2e5).planckAbsorption),
		std::pow(2, -3.5) * units::value(means.planckAbsorption), 1e-13 * units::value(means.planckAbsorption));
	EXPECT_NEAR(units::value(law.evaluate(1e-7, 1e5, 0.5).fluxExtinction), rosseland + 0.2, 1e-13);
}

TEST(RadiationCoupling, ElasticScatteringTransfersMomentumWithoutStationaryHeating) {
	hydro::HydroSystem system;
	auto gas = gasState(system);
	auto rad = radiationState(0.5);
	rad.radiativeFlux(0) = 0.1 * constants::c * rad.energy();
	auto const beforeGas = gas;
	auto const beforeRad = rad;
	OpacityLaw law;
	law.constantScattering = Opacity::from_value(opacityValue);
	coupleWithOpacityLaw(gas, rad, system, law, 1, interval(0.2));
	EXPECT_NEAR(units::value(rad.energy() / beforeRad.energy()), 1, 1e-9);
	EXPECT_NEAR(units::value(rad.radiativeFlux(0) / beforeRad.radiativeFlux(0)), std::exp(-0.2), 5e-4);
	EXPECT_GT(gas.momentum(0), beforeGas.momentum(0));
	invariants(beforeGas, beforeRad, gas, rad, 1);
}

TEST(RadiationCoupling, AnalyticOpacityPreservesStreamingConeWithLargePlanckMean) {
	hydro::HydroSystem system;
	auto gas = gasState(system);
	auto rad = radiationState(100);
	rad.radiativeFlux(0) = 0.9999 * constants::c * rad.energy();
	auto const beforeGas = gas;
	auto const beforeRad = rad;
	OpacityLaw law;
	law.ionizedGas = true;
	auto const means = law.evaluate(gas, system);
	ASSERT_GT(means.planckAbsorption, means.fluxExtinction);
	auto const step = units::Time::from_value(0.01 /
		(units::value(constants::c) * density * units::value(means.planckAbsorption)));
	ASSERT_NO_THROW(coupleWithOpacityLaw(gas, rad, system, law, 1, step));
	EXPECT_TRUE(RadiationSystem(constants::c).admissible(rad));
	invariants(beforeGas, beforeRad, gas, rad, 1);
}

TEST(RadiationCoupling, LorentzBoostedLteRemainsStationaryForStiffSteps) {
	hydro::HydroSystem system;
	for (Real ratio : {1.0, 0.125}) { for (Real beta : {0.001, 0.05}) { for (Real depth : {1e-6, 1.0, 1e9}) {
		SCOPED_TRACE(ratio);
		SCOPED_TRACE(beta);
		SCOPED_TRACE(depth);
		auto gas = gasState(system, density, temperature, beta * units::value(constants::c));
		auto const t = units::value(system.temperature(gas));
		auto const comoving = units::EnergyDensity::from_value(units::value(constants::radiation) * std::pow(t, 4));
		Real const boostSquared = 1 / (1 - beta * beta);
		RadiationSystem::State rad;
		// Independent Lorentz transformation of isotropic comoving radiation;
		// it satisfies the retained first-order source equations exactly under M1.
		rad.energy() = boostSquared * (1 + beta * beta / 3) * comoving;
		rad.radiativeFlux(0) = (4 * beta * boostSquared / 3) * constants::c * comoving;
		auto const beforeGas = gas;
		auto const beforeRad = rad;
		couple(gas, rad, system, Opacity::from_value(opacityValue), ratio, interval(depth, ratio));
		test::expectStateNear(gas, beforeGas, 3e-12);
		test::expectStateNear(rad, beforeRad, 3e-12);
		invariants(beforeGas, beforeRad, gas, rad, ratio);
	}
	}
	}
}

TEST(RadiationCoupling, BoostedTransverseDiffusionHasNoSpuriousThermalSource) {
	if constexpr (ndim < 2) GTEST_SKIP() << "A transverse diffusion flux needs at least two dimensions";
	hydro::HydroSystem system;
	// At increasing opacity the diffusive flux scales inversely with optical
	// depth, while its balancing pressure-force impulse remains finite. Keeping
	// depth*f0 bounded also keeps the uncoupled gas drive thermally admissible.
	std::array<std::array<Real, 2>, 8> const cases{{
		{1, 1e-7}, {1, 1e-3}, {1e6, 1e-7}, {1e6, 1e-3},
		{1e9, 1e-9}, {1e9, 1e-7}, {1e12, 1e-12}, {1e12, 1e-9}}};
	Real largestThermalError = 0, largestRadiationError = 0, largestFluxError = 0, largestDiffusionDrift = 0;
	for (bool oblique : {false, true}) { for (Real beta : {0.001, 0.05}) { for (auto const& parameters : cases) {
		Real const depth = parameters[0], diffusionFraction = parameters[1];
		SCOPED_TRACE(oblique);
		SCOPED_TRACE(beta);
		SCOPED_TRACE(diffusionFraction);
		SCOPED_TRACE(depth);
		std::array<Real, ndim> velocityDirection{}, diffusionDirection{};
		velocityDirection[0] = oblique ? std::sqrt(Real(0.5)) : 1;
		velocityDirection[1] = oblique ? std::sqrt(Real(0.5)) : 0;
		diffusionDirection[0] = -velocityDirection[1];
		diffusionDirection[1] = velocityDirection[0];
		auto primitive = system.reconstructionVariables(gasState(system));
		for (int d = 0; d < ndim; ++d) { primitive.velocity(d) = (beta * velocityDirection[d]) * constants::c; }
		auto gas = system.conservedState(primitive);
		auto const thermal = system.internalEnergy(gas);
		auto const t = units::value(system.temperature(gas));
		auto const comovingEnergy = units::EnergyDensity::from_value(units::value(constants::radiation) * std::pow(t, 4));
		auto const comovingFlux = diffusionFraction * constants::c * comovingEnergy;
		Real const gammaSquared = 1 / (1 - beta * beta), gamma = std::sqrt(gammaSquared);
		Real const f2 = diffusionFraction * diffusionFraction;
		// Boost an independently specified comoving M1 field whose flux is
		// perpendicular to the material velocity, as for a uniformly rotating
		// axisymmetric star. A is its pressure transverse to the comoving flux.
		auto const transversePressure = ((1 - f2) / (1 + std::sqrt(4 - 3 * f2))) * comovingEnergy;
		RadiationSystem::State rad;
		rad.energy() = gammaSquared * (comovingEnergy + beta * beta * transversePressure);
		for (int d = 0; d < ndim; ++d) {
			rad.radiativeFlux(d) = gammaSquared * (comovingEnergy + transversePressure) * primitive.velocity(d)
				+ (gamma * diffusionDirection[d]) * comovingFlux;
		}
		auto const beforeGas = gas;
		auto const beforeRad = rad;
		hydro::ConservedState gasDrive{};
		RadiationSystem::State radDrive{};
		// The implemented gray source has Q=0 and G=chi*F0/(gamma*c).
		// Prescribe the opposing pressure/gravity drives independently, so the
		// complete forced ODE is stationary even for a stiff optical interval.
		auto uncoupledThermal = thermal;
		for (int d = 0; d < ndim; ++d) {
			radDrive.radiativeFlux(d) = (depth * diffusionDirection[d] / gamma) * comovingFlux;
			gasDrive.momentum(d) = -radDrive.radiativeFlux(d) / (constants::c * constants::c);
			// v.dot(delta p)=0 analytically; its two terms cancel exactly for
			// both selected directions. Only the quadratic kinetic term remains.
			uncoupledThermal -= gasDrive.momentum(d) * gasDrive.momentum(d) / (2.0 * gas.density());
		}
		ASSERT_GT(uncoupledThermal, units::EnergyDensity{});
		gasDrive.auxiliary() = system.auxiliaryFromInternalEnergy(gas.density(), uncoupledThermal) - gas.auxiliary();
		coupleForced(gas, rad, gasDrive, radDrive, system, Opacity::from_value(opacityValue), 1, interval(depth));
		test::expectStateNear(gas, beforeGas, 3e-11);
		EXPECT_NEAR(Real(system.internalEnergyFromAuxiliary(gas) / thermal), 1, 3e-11);
		test::expectStateNear(rad, beforeRad, 3e-9);
		invariants(hydro::ConservedState(beforeGas + gasDrive), RadiationSystem::State(beforeRad + radDrive), gas, rad, 1);
		if (depth == 1e12) {
			largestThermalError = std::max(largestThermalError, std::abs(Real(system.internalEnergyFromAuxiliary(gas) / thermal) - 1));
			largestRadiationError = std::max(largestRadiationError, std::abs(Real(rad.energy() / beforeRad.energy()) - 1));
			long double transverseDrift = 0;
			for (int d = 0; d < ndim; ++d) {
				largestFluxError = std::max(largestFluxError, std::abs(Real((rad.radiativeFlux(d) - beforeRad.radiativeFlux(d)) / (constants::c * comovingEnergy))));
				transverseDrift += static_cast<long double>(diffusionDirection[d]) * units::value(rad.radiativeFlux(d) - beforeRad.radiativeFlux(d));
			}
			largestDiffusionDrift = std::max(largestDiffusionDrift, Real(std::abs(transverseDrift) / units::value(gamma * comovingFlux)));
		}
	}
	}
	}
	// Stiffness must not amplify the equilibrium rounding errors. The flux
	// bound is on the natural cB scale: an arbitrarily tiny diffusion component
	// embedded in oblique advective flux cannot have uniform relative accuracy.
	EXPECT_LE(largestThermalError, 128 * epsilonR);
	EXPECT_LE(largestRadiationError, 128 * epsilonR);
	EXPECT_LE(largestFluxError, 128 * epsilonR);
	std::cout << "[ stiff moving LTE ] optical_interval=1e12 max_relative_thermal_drift=" << largestThermalError
		<< " max_relative_Er_drift=" << largestRadiationError << " max_flux_drift_over_cB=" << largestFluxError
		<< " max_relative_transverse_diffusion_drift=" << largestDiffusionDrift << '\n';
}

TEST(RadiationCoupling, StiffThermalRelaxationRecoversIndependentEquilibriumRoot) {
	hydro::HydroSystem system;
	for (Real ratio : {1.0, 0.125}) { for (Real radiationFactor : {0.0, 0.01, 100.0}) {
		auto gas = gasState(system);
		auto rad = radiationState(radiationFactor);
		auto const beforeGas = gas;
		auto const beforeRad = rad;
		Real const total = units::value(gas.totalEnergy() + rad.energy() / ratio);
		Real const heatCapacity = density * units::value(constants::boltzmann / constants::atomicMassUnit) / (system.adiabaticIndex() - 1);
		Real lo = 0, hi = total / heatCapacity;
		for (int i = 0; i < 100; ++i) {
			Real const mid = (lo + hi) / 2;
			if (heatCapacity * mid + units::value(constants::radiation) * std::pow(mid, 4) / ratio > total) hi = mid;
			else lo = mid;
		}
		couple(gas, rad, system, Opacity::from_value(opacityValue), ratio, interval(1e9, ratio));
		EXPECT_NEAR(units::value(system.temperature(gas)), (lo + hi) / 2, 2e-7 * (lo + hi) / 2);
		EXPECT_TRUE(system.admissible(gas));
		EXPECT_TRUE(RadiationSystem(ratio * constants::c).admissible(rad));
		invariants(beforeGas, beforeRad, gas, rad, ratio);
	}
	}
}

TEST(RadiationCoupling, GasDominatedCellDoesNotImposeArtificialRadiationFloor) {
	hydro::HydroSystem system;
	for (Real ratio : {1.0, 0.125}) {
		auto gas = gasState(system, 1, 1);
		RadiationSystem::State rad;
		auto const initial = gas;
		auto const t = units::value(system.temperature(gas));
		auto const equilibrium = units::EnergyDensity::from_value(units::value(constants::radiation) * std::pow(t, 4));
		auto const dt = units::Time::from_value(1e9 / (ratio * units::value(constants::c) * 1e-10));
		couple(gas, rad, system, Opacity::from_value(1e-10), ratio, dt);
		// Er/u is ~6e-23 here. Gas temperature is unchanged to representable
		// accuracy, but radiation must still relax to its nonzero LTE value.
		EXPECT_NEAR(Real(rad.energy() / equilibrium), 1, 1e-7);
		invariants(initial, {}, gas, rad, ratio);
	}
}

TEST(RadiationCoupling, ThermalRelaxationConvergesAtSecondOrder) {
	hydro::HydroSystem system;
	for (Real ratio : {1.0, 0.125}) {
		auto const initialGas = gasState(system);
		auto const initialRad = radiationState(0.1);
		Real const total = units::value(initialGas.totalEnergy() + initialRad.energy() / ratio);
		Real const temperatureCoefficient = temperature / units::value(system.internalEnergy(initialGas));
		auto rhs = [&](Real e) { return units::value(constants::radiation) * std::pow(temperatureCoefficient * (total - e / ratio), 4) - e; };
		Real reference = units::value(initialRad.energy());
		Real const h = 0.5 / 8192;
		for (int i = 0; i < 8192; ++i) {
			Real const k1 = rhs(reference), k2 = rhs(reference + h * k1 / 2), k3 = rhs(reference + h * k2 / 2), k4 = rhs(reference + h * k3);
			reference += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
		}
		Real previous = 0;
		for (int steps : {16, 32, 64}) {
			auto gas = initialGas;
			auto rad = initialRad;
			for (int n = 0; n < steps; ++n) { couple(gas, rad, system, Opacity::from_value(opacityValue), ratio, interval(0.5 / steps, ratio)); }
			Real const error = std::abs(units::value(rad.energy()) - reference);
			if (previous > 0) { EXPECT_GT(previous / error, 3.7); }
			previous = error;
			invariants(initialGas, initialRad, gas, rad, ratio);
		}
	}
}

TEST(RadiationCoupling, MomentumRelaxationMatchesIndependentLinearSolution) {
	hydro::HydroSystem system;
	for (Real ratio : {1.0, 0.125}) {
		auto gas = gasState(system);
		auto rad = radiationState();
		rad.radiativeFlux(0) = 1e-5 * constants::c * rad.energy();
		auto const beforeGas = gas;
		auto const beforeRad = rad;
		Real const feedback = 4 * units::value(rad.energy()) / (3 * density * ratio * std::pow(units::value(constants::c), 2));
		Real const expectedFraction = (feedback + std::exp(-(1 + feedback))) / (1 + feedback);
		for (int n = 0; n < 256; ++n) { couple(gas, rad, system, Opacity::from_value(opacityValue), ratio, interval(1.0 / 256, ratio)); }
		EXPECT_NEAR(Real(rad.radiativeFlux(0) / beforeRad.radiativeFlux(0)), expectedFraction, 3e-7);
		EXPECT_GT(gas.momentum(0), units::MomentumDensity{});
		invariants(beforeGas, beforeRad, gas, rad, ratio);
	}
}

TEST(RadiationCoupling, MovingGasIncludesWorkAndFullM1Tensor) {
	hydro::HydroSystem system;
	auto gas = gasState(system, density, temperature, 1e6);
	auto rad = radiationState();
	rad.radiativeFlux(0) = 0.7 * constants::c * rad.energy();
	auto const beforeGas = gas;
	auto const beforeRad = rad;
	Real const ratio = 0.25;
	auto const dt = interval(1e-6, ratio);
	auto const pressure = M1::physicalFlux(RadiationSystem::toCalculationState(rad), 0, constants::c).flux[1] / constants::c;
	auto const velocity = gas.momentum(0) / gas.density();
	auto const sigma = gas.density() * Opacity::from_value(opacityValue);
	auto const blackbodyEnergy = units::EnergyDensity::from_value(units::value(constants::radiation) * std::pow(units::value(system.temperature(gas)), 4));
	auto const force = sigma * (rad.radiativeFlux(0) - (blackbodyEnergy + pressure) * velocity) / constants::c;
	auto const work = -sigma * velocity * rad.radiativeFlux(0) / constants::c;
	couple(gas, rad, system, Opacity::from_value(opacityValue), ratio, dt);
	EXPECT_NEAR(Real((gas.momentum(0) - beforeGas.momentum(0)) / (dt * force)), 1, 3e-6);
	EXPECT_NEAR(Real((gas.totalEnergy() - beforeGas.totalEnergy()) / (dt * work)), 1, 3e-5);
	invariants(beforeGas, beforeRad, gas, rad, ratio);
}

TEST(RadiationCoupling, ColdMovingPureBeamPreservesConeAndAbsorbsAtFirstOrderRate) {
	hydro::HydroSystem system;
	for (Real ratio : {1.0, 0.125}) { for (Real beta : {0.001, 0.01}) { for (Real depth : {1e-4, 0.1, 10.0}) {
		SCOPED_TRACE(ratio);
		SCOPED_TRACE(beta);
		SCOPED_TRACE(depth);
		auto gas = gasState(system, 1, 1, beta * units::value(constants::c));
		RadiationSystem::State rad;
		rad.energy() = units::EnergyDensity::from_value(1e8);
		rad.radiativeFlux(0) = constants::c * rad.energy();
		auto const beforeGas = gas;
		auto const beforeRad = rad;
		auto const beforeThermal = system.internalEnergy(gas);
		auto const opacity = Opacity::from_value(1e-10);
		auto const dt = units::Time::from_value(depth / (ratio * units::value(constants::c) * units::value(gas.density() * opacity)));
		couple(gas, rad, system, opacity, ratio, dt);
		EXPECT_TRUE(system.admissible(gas));
		EXPECT_TRUE(RadiationSystem(ratio * constants::c).admissible(rad));
		EXPECT_LT(rad.energy(), beforeRad.energy());
		EXPECT_GT(system.internalEnergyFromAuxiliary(gas), beforeThermal);
		if (depth == 1e-4) {
			// Emission is negligible and radiation inertia is <1e-12 of gas
			// inertia. Both beam moments decay initially at chi*chat*(1-beta),
			// whereas replacing aT^4 by Er gives mismatched attenuation rates.
			Real const expected = 1 - beta;
			EXPECT_NEAR(Real((beforeRad.energy() - rad.energy()) / beforeRad.energy()) / depth, expected, 6e-5);
			EXPECT_NEAR(Real((beforeRad.radiativeFlux(0) - rad.radiativeFlux(0)) / beforeRad.radiativeFlux(0)) / depth, expected, 6e-5);
		}
		invariants(beforeGas, beforeRad, gas, rad, ratio);
	}
	}
	}
}

TEST(RadiationCoupling, HighKineticDualEnergyRetainsSourceHeating) {
	hydro::HydroSystem system;
	auto gas = gasState(system, 1.0, 1.0, 1e8);
	auto rad = radiationState(1000);
	auto const beforeGas = gas;
	auto const beforeRad = rad;
	auto const beforeThermal = system.internalEnergy(gas);
	couple(gas, rad, system, Opacity::from_value(1e-10), 1, units::Time::from_value(1e-4));
	auto const expected = beforeThermal - (rad.energy() - beforeRad.energy());
	units::EnergyDensity kineticChange{};
	for (int d = 0; d < ndim; ++d) {
		auto const impulse = gas.momentum(d) - beforeGas.momentum(d);
		kineticChange += impulse * (beforeGas.momentum(d) + impulse / 2.0) / gas.density();
	}
	EXPECT_GT(system.internalEnergyFromAuxiliary(gas), beforeThermal);
	EXPECT_NEAR(units::value(system.internalEnergyFromAuxiliary(gas) - expected + kineticChange), 0, 2e-7 * units::value(expected));
	invariants(beforeGas, beforeRad, gas, rad, 1);
}

TEST(RadiationCoupling, JointThermalAndMomentumSolvePreservesTheRadiationCone) {
	hydro::HydroSystem system;
	for (Real ratio : {1.0, 0.125}) { for (Real depth : {0.1, 10.0, 1e9}) { for (Real factor : {0.001, 1000.0}) {
		auto gas = gasState(system);
		auto rad = radiationState(factor);
		rad.radiativeFlux(0) = constants::c * rad.energy();
		auto const beforeGas = gas;
		auto const beforeRad = rad;
		couple(gas, rad, system, Opacity::from_value(opacityValue), ratio, interval(depth, ratio));
		EXPECT_TRUE(system.admissible(gas));
		EXPECT_TRUE(RadiationSystem(ratio * constants::c).admissible(rad));
		invariants(beforeGas, beforeRad, gas, rad, ratio);
	}
	}
	}
}

TEST(RadiationCoupling, ForcedStiffDiffusionBalancePreservesFlux) {
	hydro::HydroSystem system;
	for (Real ratio : {1.0, 0.125}) { for (Real depth : {1.0, 1e4, 1e8}) {
		auto gas = gasState(system);
		auto rad = radiationState();
		rad.radiativeFlux(0) = 1e-8 * constants::c * rad.energy();
		auto const beforeGas = gas;
		auto const beforeRad = rad;
		hydro::ConservedState gasDrive{};
		RadiationSystem::State radDrive{};
		radDrive.radiativeFlux(0) = depth * rad.radiativeFlux(0);
		gasDrive.momentum(0) = -radDrive.radiativeFlux(0) / (ratio * constants::c * constants::c);
		// Keep the uncoupled baseline's conservative thermal defect zero.
		auto drivenGas = hydro::ConservedState(gas + gasDrive);
		drivenGas.auxiliary() = system.auxiliaryFromInternalEnergy(drivenGas.density(),
			drivenGas.totalEnergy() - drivenGas.momentum(0) * drivenGas.momentum(0) / (2.0 * drivenGas.density()));
		gasDrive.auxiliary() = drivenGas.auxiliary() - gas.auxiliary();
		coupleForced(gas, rad, gasDrive, radDrive, system, Opacity::from_value(opacityValue), ratio, interval(depth, ratio));
		test::expectStateNear(rad, beforeRad, 2e-8);
		EXPECT_NEAR(units::value(gas.momentum(0)), 0, 1e-9 * std::max(Real(1), units::value(units::abs(gasDrive.momentum(0)))));
		invariants(hydro::ConservedState(beforeGas + gasDrive), RadiationSystem::State(beforeRad + radDrive), gas, rad, ratio);
	}
	}
}

TEST(RadiationCoupling, ForcedDensityAndThermalDrivingConvergeAtSecondOrder) {
	hydro::HydroSystem system;
	Real const ratio = 0.25;
	auto const initialGas = gasState(system);
	auto const initialRad = radiationState(0.1);
	hydro::ConservedState gasDrive{};
	RadiationSystem::State radDrive{};
	gasDrive.density() = 0.2 * initialGas.density();
	gasDrive.totalEnergy() = 0.1 * initialGas.totalEnergy();
	radDrive.energy() = 0.2 * initialRad.energy();
	Real const gasEnergy = units::value(initialGas.totalEnergy()), radEnergy = units::value(initialRad.energy());
	Real const gasDelta = units::value(gasDrive.totalEnergy()), radDelta = units::value(radDrive.energy());
	Real const temperatureCoefficient = temperature / gasEnergy;
	auto rhs = [&](Real time, Real energy) {
		Real const rhoRatio = 1 + 0.2 * time;
		Real const thermal = gasEnergy + time * gasDelta - (energy - radEnergy - time * radDelta) / ratio;
		return radDelta + 0.5 * rhoRatio * (units::value(constants::radiation) * std::pow(temperatureCoefficient * thermal / rhoRatio, 4) - energy);
	};
	Real reference = radEnergy;
	Real const h = 1.0 / 8192;
	for (int i = 0; i < 8192; ++i) {
		Real const time = i * h;
		Real const k1 = rhs(time, reference), k2 = rhs(time + h / 2, reference + h * k1 / 2);
		Real const k3 = rhs(time + h / 2, reference + h * k2 / 2), k4 = rhs(time + h, reference + h * k3);
		reference += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
	}
	Real previous = 0;
	for (int steps : {16, 32, 64}) {
		auto gas = initialGas;
		auto rad = initialRad;
		for (int n = 0; n < steps; ++n) {
			coupleForced(gas, rad, gasDrive / Real(steps), radDrive / Real(steps), system,
				Opacity::from_value(opacityValue), ratio, interval(0.5 / steps, ratio));
		}
		Real const error = std::abs(units::value(rad.energy()) - reference);
		if (previous > 0) { EXPECT_GT(previous / error, 3.7); }
		previous = error;
		EXPECT_NEAR(Real(gas.density() / initialGas.density()), 1.2, 3e-14);
		invariants(hydro::ConservedState(initialGas + gasDrive), RadiationSystem::State(initialRad + radDrive), gas, rad, ratio);
	}
}

TEST(RadiationCoupling, ZeroOpacityAppliesOnlyTransportAndInvalidInputIsAtomic) {
	hydro::HydroSystem system;
	auto gas = gasState(system);
	auto rad = radiationState();
	auto const beforeGas = gas;
	auto const beforeRad = rad;
	hydro::ConservedState dg{};
	RadiationSystem::State dr{};
	dg.totalEnergy() = 0.1 * gas.totalEnergy();
	dr.energy() = 0.1 * rad.energy();
	coupleForced(gas, rad, dg, dr, system, {}, 0.25, interval(1));
	test::expectStateNear(gas, hydro::ConservedState(beforeGas + dg), 0);
	test::expectStateNear(rad, RadiationSystem::State(beforeRad + dr), 0);
	auto const unchangedGas = gas;
	auto const unchangedRad = rad;
	EXPECT_THROW(couple(gas, rad, system, Opacity::from_value(-1), 1, interval(1)), std::invalid_argument);
	EXPECT_THROW(couple(gas, rad, system, Opacity::from_value(1), 0, interval(1)), std::invalid_argument);
	test::expectStateNear(gas, unchangedGas, 0);
	test::expectStateNear(rad, unchangedRad, 0);
}

TEST(RadiationCoupling, UnrealizableThinForcingIsRejectedWithoutClipping) {
	hydro::HydroSystem system;
	auto gas = gasState(system, 0.01, 1402.743);
	RadiationSystem::State rad;
	rad.energy() = units::EnergyDensity::from_value(0.01465342324909228);
	auto const beforeGas = gas;
	auto const beforeRad = rad;
	hydro::ConservedState gasDrive;
	gasDrive.density() = units::Density::from_value(0.03308338909895038);
	gasDrive.momentum(0) = units::MomentumDensity::from_value(-649373.99947350135);
	gasDrive.totalEnergy() = units::EnergyDensity::from_value(36856256158355.766);
	gasDrive.auxiliary() = 3.308338909895038 * gas.auxiliary();
	RadiationSystem::State radDrive;
	radDrive.energy() = units::EnergyDensity::from_value(24495335094.028572);
	radDrive.radiativeFlux(0) = units::EnergyFlux::from_value(-1.8392121911067523e21);
	// This reproduced an unlimited opacity-transition predictor: |delta F|/c
	// exceeds delta Er by a factor 2.5 while the optical interval is ~1e-5.
	// A local source cannot repair the nonphysical driving; shared transport
	// fluxes must be limited before they supply the predictor.
	EXPECT_THROW(coupleForced(gas, rad, gasDrive, radDrive, system, Opacity::from_value(1e-10), 1,
		units::Time::from_value(9.6291660077323533e-05)), std::runtime_error);
	gas.forEach([&](auto field, auto value) { EXPECT_EQ(value, beforeGas.template get<field>()); });
	rad.forEach([&](auto field, auto value) { EXPECT_EQ(value, beforeRad.template get<field>()); });
}

TEST(RadiationCoupling, BrightTransportIntoColdAtmosphereRemainsConservative) {
	if constexpr (ndim != 3) GTEST_SKIP() << "Exact three-dimensional rotating-star failure reproduction";
	hydro::HydroSystem system(Real(5) / 3, units::Density::from_value(1e-14),
		units::Pressure::from_value(1e-14), {}, 0.6);
	hydro::ConservedState gas, gasDrive;
	RadiationSystem::State rad, radDrive;
	gas.density() = units::Density::from_value(1e-12);
	gas.totalEnergy() = units::EnergyDensity::from_value(0.91280058491899252);
	gas.auxiliary() = units::Density::from_value(91280058.491899446);
	rad.energy() = units::EnergyDensity::from_value(2.813571076209457);
	gasDrive.density() = units::Density::from_value(4.674985305102448e-13);
	gasDrive.momentum(0) = units::MomentumDensity::from_value(-1.4642902632026057e-06);
	gasDrive.momentum(1) = units::MomentumDensity::from_value(-4.2661144358431261e-06);
	gasDrive.momentum(2) = units::MomentumDensity::from_value(5.6336284020378444e-11);
	gasDrive.totalEnergy() = units::EnergyDensity::from_value(-190.82539523408113);
	gasDrive.auxiliary() = units::Density::from_value(150201.99726428092);
	radDrive.energy() = units::EnergyDensity::from_value(165010.74860785846);
	radDrive.radiativeFlux(0) = units::EnergyFlux::from_value(-4946914226080577.0);
	radDrive.radiativeFlux(1) = units::EnergyFlux::from_value(-2020169074622.3882);
	radDrive.radiativeFlux(2) = units::EnergyFlux::from_value(-13611896605.153202);
	auto const drivenGas = hydro::ConservedState(gas + gasDrive);
	auto const drivenRad = RadiationSystem::State(rad + radDrive);
	ASSERT_TRUE(system.admissible(drivenGas));
	ASSERT_TRUE(RadiationSystem(constants::c).admissible(drivenRad));
	// This accepted atmosphere drive is optically thin (chi*c*dt=0.001656),
	// with a bright nearly streaming incoming beam and positive selected gas
	// heat. An old-radiation Newton guess falsely assigns the incoming energy
	// to trial gas heat, exhausting subdivision despite an admissible solution.
	ASSERT_NO_THROW(coupleForced(gas, rad, gasDrive, radDrive, system,
		Opacity::from_value(0.34), 1, units::Time::from_value(0.16247832095792569)));
	EXPECT_TRUE(system.admissible(gas));
	EXPECT_TRUE(RadiationSystem(constants::c).admissible(rad));
	EXPECT_GT(system.internalEnergyFromAuxiliary(gas), system.internalEnergy(drivenGas));
	EXPECT_LT(rad.energy(), drivenRad.energy());
	invariants(drivenGas, drivenRad, gas, rad, 1);
}
} // namespace
