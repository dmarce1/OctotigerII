#include "testSupport.hpp"
#include <gtest/gtest.h>
#include <iostream>
#include "octotigerII/hydro/hydroSystem.hpp"
#include "octotigerII/radiation/radiationTransport.hpp"
#include "octotigerII/finiteVolume/frame.hpp"

using namespace octotigerII;

namespace {

template <typename System>
void uniformOnRotatingMesh(System const& system, typename System::State const& uniform) {
	if constexpr (ndim < 2) return;
	mesh::PhysicalCoordinates lower{};
	lower.fill(units::Length::from_value(-1));
	mesh::PatchData<typename System::State> patch(mesh::MeshLayout(6, 2), units::Length::from_value(1.0 / 3), lower);
	std::fill(patch.values().begin(), patch.values().end(), uniform);
	patch.timeState().time = units::Time::from_value(0.37);
	finiteVolume::RotatingFrame const frame(units::InverseTime::from_value(0.8));
	finiteVolume::MusclHancock<System> solver(system, frame);
	for (int step = 0; step < 8; ++step) {
		auto const dt = solver.stableTimestep(patch, 0.25);
		EXPECT_LE(units::value(frame.omega() * dt), 0.1);
		solver.advanceWithBoundaryUpdater(patch, dt, [&](auto& p, auto) {
			mesh::forEachCoordinate(p.layout().extents(), [&](auto const& cell) {
				if (!p.layout().isInterior(cell)) p.atStorage(cell) = uniform;
			});
		});
		patch.layout().forEachInterior([&](auto const& cell, auto) {
			test::expectStateNear(patch.atInterior(cell), uniform, 2e-12);
			EXPECT_TRUE(system.admissible(patch.atInterior(cell)));
		});
	}
}

#if OCTOTIGERII_HYDRO
hydro::ConservedState gasState(hydro::HydroSystem const& gas, Real density = 1, Real velocity = 0.3) {
	hydro::PrimitiveState primitive;
	primitive.density() = units::Density::from_value(density);
	primitive.pressure() = units::Pressure::from_value(2);
	primitive.velocity(0) = units::Velocity::from_value(velocity);
	if constexpr (ndim >= 2) primitive.velocity(1) = units::Velocity::from_value(-0.2);
	return gas.conservedState(primitive);
}

TEST(RotatingTransport, ZeroAngularVelocityMatchesExistingUpdateExactly) {
	hydro::HydroSystem const gas(1.4);
	mesh::PatchData<hydro::ConservedState> patch(mesh::MeshLayout(6, 2), units::Length::from_value(0.25));
	mesh::forEachCoordinate(patch.layout().extents(), [&](auto const& cell) {
		patch.atStorage(cell) = gasState(gas, 1 + 0.03 * cell[0], 0.2 + 0.01 * cell[0]);
	});
	auto const dt = units::Time::from_value(0.002);
	hydro::Solver::Workspace legacyWork, zeroWork;
	auto legacy = patch, zero = patch;
	hydro::Solver(gas).advanceInto(patch, dt, legacyWork, [&](auto const& cell, auto const& value) { legacy.atInterior(cell) = value; });
	hydro::Solver(gas, finiteVolume::RotatingFrame{}, units::Time::from_value(42)).advanceInto(
		patch, dt, zeroWork, [&](auto const& cell, auto const& value) { zero.atInterior(cell) = value; });
	patch.layout().forEachInterior([&](auto const& cell, auto) {
		legacy.atInterior(cell).forEach([&](auto field, auto const& value) { EXPECT_EQ(value, zero.atInterior(cell).template get<field>()); });
	});
	for (int axis = 0; axis < ndim; ++axis)
		for (std::size_t face = 0; face < legacyWork.fluxes[axis].size(); ++face)
			legacyWork.fluxes[axis][face].forEach([&](auto field, auto const& value) {
				EXPECT_EQ(value, zeroWork.fluxes[axis][face].template get<field>());
			});
}

TEST(RotatingTransport, MovingHydroFaceRetainsInertialPressureWork) {
	hydro::HydroSystem const gas(1.4);
	auto const state = gasState(gas, 1.3, 3);
	auto const speed = units::Velocity::from_value(3);
	auto const physical = gas.physicalFlux(state, 0, speed);
	auto const numerical = gas.riemann(state, state, 0, speed);
	EXPECT_NEAR(units::value(physical.mass()), 0, 1e-14);
	EXPECT_NEAR(units::value(physical.momentum(0)), 2, 1e-14);
	EXPECT_NEAR(units::value(physical.energy()), 6, 1e-13);
	test::expectStateNear(numerical, physical, 1e-12);
}

TEST(RotatingTransport, MovingHydroFaceSamplesTheCorrectSideOfContact) {
	hydro::HydroSystem const gas(1.4);
	auto const left = gasState(gas, 1, 0), right = gasState(gas, 2, 0);
	auto const speed = units::Velocity::from_value(0.7);
	auto const flux = gas.riemann(left, right, 0, speed);
	test::expectStateNear(flux, gas.physicalFlux(right, 0, speed), 1e-12);
}

TEST(RotatingTransport, UniformInertialHydroSurvivesChangingNormals) {
	hydro::HydroSystem const gas(1.4);
	uniformOnRotatingMesh(gas, gasState(gas));
}

TEST(RotatingTransport, CorotatingGasUsesRelativeCharacteristicSpeeds) {
	if constexpr (ndim < 2) return;
	hydro::HydroSystem const gas(1.4);
	finiteVolume::RotatingFrame const frame(units::InverseTime::from_value(1));
	mesh::PhysicalCoordinates lower{};
	lower.fill(units::Length::from_value(100));
	mesh::PatchData<hydro::ConservedState> patch(mesh::MeshLayout(8, 2), units::Length::from_value(0.125), lower);
	patch.layout().forEachInterior([&](auto const& cell, auto) {
		auto primitive = gas.reconstructionVariables(gasState(gas));
		auto const position = patch.layout().cellCenter(lower, patch.cellWidth(), cell);
		auto const velocity = frame.velocityGrid(position);
		for (int axis = 0; axis < ndim; ++axis) primitive.velocity(axis) = velocity[axis];
		patch.atInterior(cell) = gas.conservedState(primitive);
	});
	auto const stationary = hydro::Solver(gas).stableTimestep(patch, 0.25);
	auto const rotating = hydro::Solver(gas, frame).stableTimestep(patch, 0.25);
	EXPECT_GT(rotating, Real(10) * stationary);
	EXPECT_LE(rotating, frame.maximumTimestep());
}

TEST(RotatingTransport, AdvectedEntropyPulseHasSecondOrderTemporalConvergence) {
	if constexpr (ndim < 2) return;
	// A compact smooth entropy perturbation has an exact inertial solution:
	// constant pressure and velocity, with density translated by v*t. Its
	// asymmetric shape also turns relative to the logical mesh. All support
	// stays well inside the box, so prescribed ghost values are uniform.
	hydro::HydroSystem const gas(1.4);
	finiteVolume::RotatingFrame const frame(units::InverseTime::from_value(1.1));
	using Patch = mesh::PatchData<hydro::ConservedState>;
	mesh::PhysicalCoordinates lower{};
	lower.fill(units::Length::from_value(-1));
	constexpr int cells = 12;
	auto const width = units::Length::from_value(2.0 / cells);
	auto const stop = units::Time::from_value(0.04);
	std::array<Real, 3> const velocity{0.35, 0.15, -0.05}, center{0.12, -0.1, 0}, radius{0.55, 0.7, 0.65};
	auto exact = [&](mesh::Coordinates const& cell, units::Time time) {
		Real density = 0;
		// Tensor-product two-point Gauss quadrature gives the same cell-average
		// initial/reference interpretation on the rigidly rotating volumes.
		for (unsigned corner = 0; corner < (1u << ndim); ++corner) {
			mesh::PhysicalCoordinates xi{};
			for (int d = 0; d < ndim; ++d)
				xi[d] = lower[d] + (cell[d] + 0.5 + ((corner & (1u << d)) ? 1 : -1) / std::sqrt(12.0)) * width;
			auto const x = frame.toInertial(xi, time);
			Real r2 = 0;
			for (int d = 0; d < ndim; ++d) {
				Real const offset = (units::value(x[d]) - velocity[d] * units::value(time) - center[d]) / radius[d];
				r2 += offset * offset;
			}
			density += 1 + (r2 < 1 ? 0.2 * std::exp(-r2 / (1 - r2)) : 0);
		}
		hydro::PrimitiveState primitive;
		primitive.density() = units::Density::from_value(density / Real(1u << ndim));
		primitive.pressure() = units::Pressure::from_value(1);
		for (int d = 0; d < ndim; ++d) primitive.velocity(d) = units::Velocity::from_value(velocity[d]);
		return gas.conservedState(primitive);
	};
	auto fill = [&](Patch& patch, units::Time time) {
		mesh::forEachCoordinate(patch.layout().extents(), [&](auto cell) {
			if (patch.layout().isInterior(cell)) return;
			auto logical = cell;
			for (int d = 0; d < ndim; ++d) logical[d] -= patch.layout().ghostWidth();
			patch.atStorage(cell) = exact(logical, time);
		});
	};
	auto integrate = [&](int steps) {
		Patch patch(mesh::MeshLayout(cells, 2), width, lower);
		patch.layout().forEachInterior([&](auto const& cell, auto) { patch.atInterior(cell) = exact(cell, {}); });
		auto const dt = stop / Real(steps);
		for (int step = 0; step < steps; ++step) {
			auto const time = Real(step) * dt;
			fill(patch, time);
			hydro::Solver const solver(gas, frame, time);
			hydro::Solver::Workspace first, second;
			solver.advanceInto(patch, {}, first, [](auto const&, auto const&) {});
			auto midpoint = patch;
			patch.layout().forEachInterior([&](auto const& cell, auto) {
				hydro::ConservedFlux divergence;
				for (int d = 0; d < ndim; ++d) {
					auto upper = cell;
					++upper[d];
					divergence += first.fluxes[d][patch.layout().faceIndex(d, cell)] - first.fluxes[d][patch.layout().faceIndex(d, upper)];
				}
				midpoint.atInterior(cell) += (Real(0.5) * dt / width) * divergence;
			});
			fill(midpoint, time + Real(0.5) * dt);
			auto next = patch;
			solver.advanceInto(patch, dt, second, [&](auto const& cell, auto const& value) { next.atInterior(cell) = value; },
				[&](auto const&, auto const& storage) { return midpoint.atStorage(storage); }, false);
			patch = std::move(next);
		}
		return patch;
	};
	auto const reference = integrate(128);
	std::array<Real, 3> errors{};
	int run = 0;
	for (int steps : {4, 8, 16}) {
		auto const actual = integrate(steps);
		actual.layout().forEachInterior([&](auto const& cell, auto) {
			errors[run] += std::abs(units::value(actual.atInterior(cell).density() - reference.atInterior(cell).density()));
		});
		++run;
	}
	Real analyticError = 0, evolvedSignal = 0;
	reference.layout().forEachInterior([&](auto const& cell, auto) {
		analyticError += std::abs(units::value(reference.atInterior(cell).density() - exact(cell, stop).density()));
		evolvedSignal += std::abs(units::value(exact(cell, stop).density() - exact(cell, {}).density()));
	});
	std::cout << "Rotating entropy pulse temporal density errors: " << errors[0] << ", " << errors[1] << ", " << errors[2]
		<< "; orders " << std::log2(errors[0] / errors[1]) << ", " << std::log2(errors[1] / errors[2])
		<< "; analytic error / evolved signal " << analyticError / evolvedSignal << '\n';
	EXPECT_GT(errors[0], 1e-9);
	EXPECT_LT(errors[1], 0.3 * errors[0]);
	EXPECT_LT(errors[2], 0.3 * errors[1]);
	EXPECT_LT(analyticError, 0.5 * evolvedSignal);
}
#endif

#if OCTOTIGERII_RADIATION
TEST(RotatingTransport, RadiationAleSubtractsStoredPhysicalMoments) {
	radiation::RadiationSystem const radiation(Real(0.1) * constants::c);
	radiation::RadiationSystem::State state;
	state.energy() = units::EnergyDensity::from_value(2);
	state.radiativeFlux(0) = Real(0.3) * constants::c * state.energy();
	auto const speed = Real(0.04) * constants::c;
	auto const expected = radiation.physicalFlux(state, 0) - radiation.advectiveFlux(state, speed);
	test::expectStateNear(radiation.physicalFlux(state, 0, speed), expected, 1e-12);
	test::expectStateNear(radiation.riemann(state, state, 0, speed), expected, 1e-12);
}

TEST(RotatingTransport, RadiationFaceFasterThanReducedLightUsesRightState) {
	auto const chat = Real(0.1) * constants::c;
	radiation::RadiationSystem const radiation(chat);
	radiation::RadiationSystem::State left, right;
	left.energy() = units::EnergyDensity::from_value(1);
	right.energy() = units::EnergyDensity::from_value(2);
	left.radiativeFlux(0) = Real(0.2) * constants::c * left.energy();
	right.radiativeFlux(0) = -Real(0.1) * constants::c * right.energy();
	auto const speed = Real(1.4) * chat;
	test::expectStateNear(radiation.riemann(left, right, 0, speed), radiation.physicalFlux(right, 0, speed), 1e-12);
}

TEST(RotatingTransport, UniformInertialRadiationSurvivesChangingNormals) {
	// Small reduced speed makes the face rotate appreciably during the steps.
	radiation::RadiationSystem const radiation(units::Velocity::from_value(1));
	radiation::RadiationSystem::State state;
	state.energy() = units::EnergyDensity::from_value(2);
	state.radiativeFlux(0) = Real(0.3) * constants::c * state.energy();
	if constexpr (ndim >= 2) state.radiativeFlux(1) = Real(0.2) * constants::c * state.energy();
	uniformOnRotatingMesh(radiation, state);
}
#endif

} // namespace
