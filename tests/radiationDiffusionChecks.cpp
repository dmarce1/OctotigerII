#include "testSupport.hpp"
#include <gtest/gtest.h>
#include <iostream>
#include "octotigerII/radiation/diffusionFlux.hpp"
#include "octotigerII/buildConfig.hpp"
#if OCTOTIGERII_HYDRO
#include "octotigerII/radiation/coupledPatch.hpp"
#endif

using namespace octotigerII;
using namespace octotigerII::radiation;

namespace {
RadiationSystem::State state(Real energy, Real reducedFlux = 0) {
	RadiationSystem::State result{};
	result.energy() = units::EnergyDensity::from_value(energy);
	result.radiativeFlux(0) = reducedFlux * constants::c * result.energy();
	return result;
}

TEST(RadiationDiffusion, ZeroExtinctionPreservesEveryThinFluxComponentExactly) {
	RadiationSystem const system(Real(0.2) * constants::c);
	auto const left = state(1, .4), right = state(1.1, -.2);
	auto const thin = system.riemann(left, right, 0, Real(.03) * constants::c);
	auto const corrected = diffusionCorrectedFlux(system, thin, left, right, left, right,
		{}, {}, {}, {}, 0, units::Length::from_value(1), Real(.03) * constants::c);
	test::expectStateNear(corrected, thin, 0);
}

TEST(RadiationDiffusion, MovingSourceEquilibriumAndAleAdvectionAreNotAttenuated) {
	// The faster case checks the algebraic manifold's subluminal domain only;
	// it is not a claim that the nonrelativistic evolution is accurate there.
	for (Real velocityScale : {Real(.03), Real(.6)})
	for (Real ratio : {Real(1), Real(.1)}) for (Real tau : {Real(.01), Real(10), Real(1e8)}) {
		RadiationSystem const system(ratio * constants::c);
		MaterialVelocity velocity{};
		Real beta2 = 0;
		for (int d = 0; d < ndim; ++d) {
			Real const beta = velocityScale / (d + 1);
			velocity[d] = beta * constants::c;
			beta2 += beta * beta;
		}
		auto const equilibrium = materialEquilibriumMoments(units::EnergyDensity::from_value(2), velocity);
		auto const emission = equilibrium.energy() * (Real(3) * (1 - beta2) / (3 + beta2));
		auto thermalBalance = equilibrium.energy() - emission;
		for (int d = 0; d < ndim; ++d)
			thermalBalance -= velocity[d] * equilibrium.radiativeFlux(d) / (constants::c * constants::c);
		EXPECT_NEAR(units::value(thermalBalance), 0, 2e-14);
		for (int normal = 0; normal < ndim; ++normal) {
			auto const w = Real(.012) * constants::c;
			auto const thin = system.physicalFlux(equilibrium, normal, w);
			auto const corrected = diffusionCorrectedFlux(system, thin, equilibrium, equilibrium,
				equilibrium, equilibrium, Extinction::from_value(tau), Extinction::from_value(tau),
				velocity, velocity, normal, units::Length::from_value(1), w);
			test::expectStateNear(corrected, thin, 2e-14);
			// Independent momentum-source check: F_n=(aT^4 delta_nj+P_nj)v_j.
			auto const pressure = M1::physicalFlux(RadiationSystem::toCalculationState(equilibrium), normal, constants::c);
			auto sourceEquilibrium = emission * velocity[normal];
			for (int d = 0; d < ndim; ++d) sourceEquilibrium += pressure.flux[d + 1] * (velocity[d] / constants::c);
			EXPECT_NEAR(units::value(sourceEquilibrium / constants::c),
				units::value(equilibrium.radiativeFlux(normal) / constants::c), 2e-14);
		}
	}
}

TEST(RadiationDiffusion, OpaquePressureIgnoresTransientHancockAnisotropy) {
	RadiationSystem const system(constants::c);
	auto const transient = state(1, .7);
	auto const thin = system.physicalFlux(transient, 0);
	auto const corrected = diffusionCorrectedFlux(system, thin, transient, transient, transient, transient,
		Extinction::from_value(1e8), Extinction::from_value(1e8), {}, {}, 0, units::Length::from_value(1));
	EXPECT_NEAR(units::value(corrected.radiativeFlux(0) / (constants::c * constants::c)), Real(1) / 3, 2e-15);
	EXPECT_LT(std::abs(units::value(corrected.energy() / constants::c)), 1e-15);
}

TEST(RadiationDiffusion, ThinLimitPerturbationIsQuadraticInCellWidth) {
	RadiationSystem const system(constants::c);
	auto const face = state(1, .1);
	auto const thin = system.physicalFlux(face, 0);
	Real previous = 0;
	for (Real dx : {Real(.02), Real(.01), Real(.005)}) {
		auto const corrected = diffusionCorrectedFlux(system, thin, state(1 - .05 * dx), state(1 + .05 * dx), face, face,
			Extinction::from_value(1), Extinction::from_value(1), {}, {}, 0, units::Length::from_value(dx));
		Real const error = std::abs(units::value((corrected.energy() - thin.energy()) / constants::c));
		if (previous > 0) EXPECT_NEAR(previous / error, 4, .002);
		previous = error;
	}
}

TEST(RadiationDiffusion, FixedMeshFourierDiffusionCoefficientHasCorrectOpaqueLimit) {
	// Exercise the production reconstruction and the final shared-face hook.
	// This is a spatial-operator test; source/transport time integration is a
	// separate requirement. Zero initial F deliberately makes standard HLL's
	// energy transport entirely numerical, exposing its spurious diffusion.
	constexpr int cells = 8;
	auto const dx = units::Length::from_value(Real(1) / cells);
	Real const wave = Real(2) * piR / cells;
	for (Real ratio : {Real(1), Real(.1)}) {
		RadiationSystem const system(ratio * constants::c);
		Solver const solver(system);
		mesh::PatchData<RadiationSystem::State> patch(mesh::MeshLayout(cells, 2), dx);
		patch.layout().forEachInterior([&](auto const& cell, auto) {
			patch.atInterior(cell) = state(1 + .05 * std::cos(wave * (cell[0] + .5)));
		});
		physics::fillGhostCells(patch, physics::BoundaryConditions::periodic(), system);
		Real previousError = 1;
		for (Real tau : {Real(10), Real(100), Real(1000), Real(10000)}) {
			Solver::Workspace workspace;
			auto correct = [&](auto const& flux, auto const& centerLeft, auto const& centerRight,
				auto const& faceLeft, auto const& faceRight, auto const&, auto const&, int normal,
				units::Length width, units::Velocity speed) {
				return diffusionCorrectedFlux(system, flux, centerLeft, centerRight, faceLeft, faceRight,
					Extinction::from_value(tau / units::value(dx)), Extinction::from_value(tau / units::value(dx)),
					{}, {}, normal, width, speed);
			};
			solver.advanceInto(patch, {}, workspace, [](auto const&, auto const&) {},
				[](auto const& u, auto const&) { return u; }, false, correct);
			// The same correction must survive the production shared-face limiter
			// at a nonzero light-crossing CFL, not just at an instantaneous flux.
			Solver::Workspace limited;
			auto const dt = (Real(.1) / ndim) * dx / system.reducedLightSpeed();
			solver.advanceInto(patch, dt, limited,
				[&](auto const&, auto const& u) { EXPECT_TRUE(system.admissible(u)); },
				[](auto const& u, auto const&) { return u; }, false, correct);
			for (int axis = 0; axis < ndim; ++axis)
				for (std::size_t face = 0; face < workspace.fluxes[axis].size(); ++face)
					test::expectStateNear(limited.fluxes[axis][face], workspace.fluxes[axis][face], 0);
			Real projection = 0, norm = 0, energyRate = 0;
			patch.layout().forEachInterior([&](auto const& cell, auto) {
				auto upper = cell; ++upper[0];
				auto const delta = workspace.fluxes[0][patch.layout().faceIndex(0, cell)].energy()
					- workspace.fluxes[0][patch.layout().faceIndex(0, upper)].energy();
				Real const mode = std::cos(wave * (cell[0] + .5));
				Real const derivative = units::value(delta / system.reducedLightSpeed()) / units::value(dx);
				projection += derivative * mode;
				norm += .05 * mode * mode;
				energyRate += derivative;
			});
			Real const measured = -projection / norm;
			Real const exactDiscrete = Real(4) * std::pow(std::sin(wave / 2), 2) / (Real(3) * tau * units::value(dx));
			Real const error = std::abs(measured / exactDiscrete - 1);
			std::cout << "diffusion spatial tau=" << tau << " ratio=" << ratio << " Dmeasured/Dexact=" << measured / exactDiscrete << '\n';
			EXPECT_LT(error, previousError / 5);
			EXPECT_LT(std::abs(energyRate), 1e-14);
			previousError = error;
		}
		EXPECT_LT(previousError, 1e-4);
	}
}

TEST(RadiationDiffusion, InstantaneousProbeLimitsSharpOpticalTransitionOverPredictorInterval) {
	// A bright, opaque core beside a faint transparent atmosphere exposes a
	// pressure kick larger than the AP energy flux can keep realizable. The
	// source is weak in the atmosphere, so it cannot repair that predictor.
	constexpr int cells = 8;
	auto const dx = units::Length::from_value(1);
	RadiationSystem const system(constants::c);
	Solver const solver(system);
	Fields patch(mesh::MeshLayout(cells, 2), dx);
	patch.layout().forEachInterior([&](auto const& cell, auto) {
		patch.atInterior(cell) = state(cell[0] < cells / 2 ? 1 : 1e-14);
	});
	physics::fillGhostCells(patch, physics::BoundaryConditions::periodic(), system);
	auto const predictorInterval = (Real(.3) / ndim) * dx / constants::c;
	auto correct = [&](auto const& flux, auto const& centerLeft, auto const& centerRight,
		auto const& faceLeft, auto const& faceRight, auto const&, auto const&, int normal,
		units::Length width, units::Velocity speed) {
		auto extinction = [](auto const& u) {
			return Extinction::from_value(units::value(u.energy()) > .5 ? 100 : 1e-10);
		};
		return diffusionCorrectedFlux(system, flux, centerLeft, centerRight, faceLeft, faceRight,
			extinction(centerLeft), extinction(centerRight), {}, {}, normal, width, speed);
	};
	auto unchanged = [](auto const& u, auto const&) { return u; };
	auto checkZeroUpdate = [&](auto const& cell, auto const& u) {
		test::expectStateNear(u, patch.atInterior(cell), 0);
	};
	Solver::Workspace unbounded, limited;
	solver.advanceInto(patch, {}, unbounded, checkZeroUpdate, unchanged, false, correct);
	solver.advanceInto(patch, {}, limited, checkZeroUpdate, unchanged, false, correct, predictorInterval);
	int unboundedFailures = 0;
	patch.layout().forEachInterior([&](auto const& cell, auto) {
		auto update = [&](auto const& workspace) {
			RadiationSystem::Flux divergence{};
			for (int axis = 0; axis < ndim; ++axis) {
				auto upper = cell; ++upper[axis];
				divergence += workspace.fluxes[axis][patch.layout().faceIndex(axis, cell)]
					- workspace.fluxes[axis][patch.layout().faceIndex(axis, upper)];
			}
			auto const increment = RadiationSystem::integratedFlux(divergence, predictorInterval / dx);
			return RadiationSystem::State(patch.atInterior(cell) + increment);
		};
		if (!system.admissible(update(unbounded))) ++unboundedFailures;
		auto const limitedState = update(limited);
		auto const increment = RadiationSystem::State(limitedState - patch.atInterior(cell));
		EXPECT_TRUE(system.admissible(system.correctRoundoff(limitedState, componentAbs(increment))));
	});
	EXPECT_GT(unboundedFailures, 0);
}

#if OCTOTIGERII_HYDRO && OCTOTIGERII_NDIM == 1
TEST(RadiationDiffusion, TransitionOpticalDepthsRemainAdmissibleWithVaryingDensityAndTemperature) {
	constexpr int cells = 16;
	auto const dx = units::Length::from_value(Real(1) / cells);
	hydro::HydroSystem const gasSystem(Real(5) / 3, units::Density::from_value(1e-30), units::Pressure::from_value(1e-30));
	for (Real ratio : {Real(1), Real(.1)}) for (Real tau : {Real(.03), Real(.3), Real(1), Real(3)}) {
		RadiationSystem const radSystem(ratio * constants::c);
		auto const opacity = Opacity::from_value(tau / (1e-17 * units::value(dx)));
		hydro::Fields gas(mesh::MeshLayout(cells, 4), dx);
		Fields rad(mesh::MeshLayout(cells, 4), dx);
		gas.layout().forEachInterior([&](auto const& cell, auto) {
			Real const phase = Real(2) * piR * (cell[0] + .5) / cells;
			Real const t = 100 * (1 + .2 * std::cos(phase));
			hydro::PrimitiveState primitive{};
			primitive.density() = units::Density::from_value(1e-17 * (1 + .4 * std::sin(phase)));
			primitive.pressure() = primitive.density() * constants::boltzmann * units::Temperature::from_value(t) / constants::atomicMassUnit;
			primitive.velocity(0) = Real(.001 * std::sin(phase)) * radSystem.reducedLightSpeed();
			gas.atInterior(cell) = gasSystem.conservedState(primitive);
			MaterialVelocity velocity{};
			velocity[0] = primitive.velocity(0);
			Real const beta = primitive.velocity(0) / constants::c;
			Real const lorentzEnergyFactor = (3 + beta * beta) / (3 * (1 - beta * beta));
			rad.atInterior(cell) = materialEquilibriumMoments(
				units::EnergyDensity::from_value(lorentzEnergyFactor * units::value(constants::radiation) * std::pow(t, 4)), velocity);
		});
		auto total = [&]() {
			long double sum = 0;
			gas.layout().forEachInterior([&](auto const& cell, auto) {
				sum += units::value(gas.atInterior(cell).totalEnergy()) + units::value(rad.atInterior(cell).energy()) / ratio;
			});
			return sum;
		};
		auto const before = total();
		auto const dt = Real(.3) * dx / radSystem.reducedLightSpeed();
		for (int step = 0; step < 24; ++step) {
			physics::fillGhostCells(gas, physics::BoundaryConditions::periodic(), gasSystem);
			physics::fillGhostCells(rad, physics::BoundaryConditions::periodic(), radSystem);
			auto nextGas = gas;
			auto nextRad = rad;
			CoupledPatchWorkspace workspace;
			advanceCoupledPatch(gas, rad, gasSystem, radSystem, opacity, dt, workspace,
				[&](auto const& cell, auto const& g, auto const& r) {
					EXPECT_TRUE(gasSystem.admissible(g));
					EXPECT_TRUE(radSystem.admissible(r));
					nextGas.atInterior(cell) = g;
					nextRad.atInterior(cell) = r;
				});
			gas = std::move(nextGas);
			rad = std::move(nextRad);
		}
		EXPECT_NEAR(Real((total() - before) / before), 0, 5e-13);
	}
}

TEST(RadiationDiffusion, CoupledFourierModeDecaysAtThePhysicalEquilibriumDiffusionRate) {
	constexpr int cells = 16;
	constexpr Real rho = 1e-17, temperature = 100, perturbation = 1e-5;
	auto const dx = units::Length::from_value(Real(1) / cells);
	hydro::HydroSystem const gasSystem(Real(5) / 3, units::Density::from_value(1e-30), units::Pressure::from_value(1e-30));
	Real const equilibriumE = units::value(constants::radiation) * std::pow(temperature, 4);
	Real const equilibriumU = rho * units::value(constants::boltzmann / constants::atomicMassUnit) * temperature / (gasSystem.adiabaticIndex() - 1);
	Real const wave = Real(2) * piR / cells;
	for (Real ratio : {Real(1), Real(.1)}) for (Real tau : {Real(10), Real(100), Real(1000)}) {
		RadiationSystem const radiationSystem(ratio * constants::c);
		hydro::Fields gas(mesh::MeshLayout(cells, 4), dx);
		Fields radiation(mesh::MeshLayout(cells, 4), dx);
		Real const extinction = tau / units::value(dx);
		auto const opacity = Opacity::from_value(extinction / rho);
		gas.layout().forEachInterior([&](auto const& cell, auto) {
			Real const phase = wave * (cell[0] + .5);
			Real const t = temperature * (1 + perturbation * std::cos(phase));
			hydro::PrimitiveState primitive{};
			primitive.density() = units::Density::from_value(rho);
			primitive.pressure() = primitive.density() * constants::boltzmann * units::Temperature::from_value(t) / constants::atomicMassUnit;
			gas.atInterior(cell) = gasSystem.conservedState(primitive);
			auto& rad = radiation.atInterior(cell);
			rad.energy() = units::EnergyDensity::from_value(units::value(constants::radiation) * std::pow(t, 4));
			Real const derivative = -4 * equilibriumE * perturbation * std::sin(phase) * std::sin(wave) / units::value(dx);
			rad.radiativeFlux(0) = units::EnergyFlux::from_value(-units::value(constants::c) * derivative / (3 * extinction));
		});
		auto fill = [&](auto& g, auto& r) {
			physics::fillGhostCells(g, physics::BoundaryConditions::periodic(), gasSystem);
			physics::fillGhostCells(r, physics::BoundaryConditions::periodic(), radiationSystem);
		};
		auto mode = [&]() {
			long double amplitude = 0;
			gas.layout().forEachInterior([&](auto const& cell, auto) {
				long double const energy = units::value(gas.atInterior(cell).totalEnergy()) + units::value(radiation.atInterior(cell).energy()) / ratio;
				amplitude += (energy - equilibriumU - equilibriumE / ratio) * std::cos(wave * (cell[0] + .5));
			});
			return amplitude;
		};
		auto total = [&]() {
			long double sum = 0;
			gas.layout().forEachInterior([&](auto const& cell, auto) {
				sum += units::value(gas.atInterior(cell).totalEnergy()) + units::value(radiation.atInterior(cell).energy()) / ratio;
			});
			return sum;
		};
		auto const amplitude0 = mode(), total0 = total();
		auto const dt = Real(.3) * dx / radiationSystem.reducedLightSpeed();
		constexpr int steps = 8;
		for (int step = 0; step < steps; ++step) {
			fill(gas, radiation);
			auto nextGas = gas;
			auto nextRadiation = radiation;
			CoupledPatchWorkspace workspace;
			advanceCoupledPatch(gas, radiation, gasSystem, radiationSystem, opacity, dt, workspace,
				[&](auto const& cell, auto const& g, auto const& r) {
					nextGas.atInterior(cell) = g;
					nextRadiation.atInterior(cell) = r;
				});
			gas = std::move(nextGas);
			radiation = std::move(nextRadiation);
		}
		Real const measured = -std::log(Real(mode() / amplitude0)) / units::value(Real(steps) * dt);
		Real const diffusivity = units::value(radiationSystem.reducedLightSpeed()) / (3 * extinction)
			* (4 * equilibriumE) / (4 * equilibriumE + ratio * equilibriumU);
		Real const expected = diffusivity * 4 * std::pow(std::sin(wave / 2), 2) / std::pow(units::value(dx), 2);
		std::cout << "coupled Fourier tau=" << tau << " ratio=" << ratio << " decay/exact=" << measured / expected << '\n';
		EXPECT_NEAR(measured / expected, 1, tau == 10 ? .025 : .01);
		EXPECT_NEAR(Real((total() - total0) / total0), 0, 5e-13);
	}
}
#endif

} // namespace
