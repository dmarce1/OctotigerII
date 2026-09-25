#include "testSupport.hpp"
#include <gtest/gtest.h>
#include <iostream>
#include <limits>
#include "octotigerII/buildConfig.hpp"
#include "octotigerII/hydro/hydroSystem.hpp"
#include "octotigerII/physics/finiteVolume.hpp"
#include "octotigerII/radiation/radiationTransport.hpp"

using namespace octotigerII;

namespace {

// Only instantiate the physics selected by this executable's problem manifest.
#if OCTOTIGERII_HYDRO && OCTOTIGERII_RADIATION
using Systems = ::testing::Types<hydro::HydroSystem, radiation::RadiationSystem>;
#elif OCTOTIGERII_HYDRO
using Systems = ::testing::Types<hydro::HydroSystem>;
#else
using Systems = ::testing::Types<radiation::RadiationSystem>;
#endif


template <typename System>
class FiniteVolume : public ::testing::Test {
protected:

	System system = [] {
		if constexpr (std::is_same_v<System, hydro::HydroSystem>) return System(1.4);
		else return System(0.25 * constants::c);
	}();

	typename System::State state(Real scale = 1) const {
		if constexpr (std::is_same_v<System, hydro::HydroSystem>) {
			hydro::PrimitiveState p;
			p.density() = units::Density::from_value(scale);
			p.pressure() = units::Pressure::from_value(2);
			for (int d = 0; d < ndim; ++d) p.velocity(d) = units::Velocity::from_value(0.2 * (d + 1));
			return system.conservedState(p);
		} else {
			typename System::State u;
			u.energy() = units::EnergyDensity::from_value(scale);
			for (int d = 0; d < ndim; ++d) u.radiativeFlux(d) = (0.2 / ndim) * constants::c * u.energy();
			return u;
		}
	}

	mesh::PatchData<typename System::State> patch(int ghosts = 2) const {
		mesh::PatchData<typename System::State> result(mesh::MeshLayout(8, ghosts), units::Length::from_value(0.25));
		result.layout().forEachInterior([&](auto const& cell, auto) { result.atInterior(cell) = state(); });
		return result;
	}
};


TYPED_TEST_SUITE(FiniteVolume, Systems);


TYPED_TEST(FiniteVolume, UniformStateIsPreservedForEveryLimiter) {
	for (auto limiter : {physics::Limiter::Minmod, physics::Limiter::VanLeer, physics::Limiter::MinmodTheta}) {
		auto patch = this->patch();
		physics::MusclHancock<TypeParam> solver(this->system, limiter, 1.3);
		auto const dt = solver.stableTimestep(patch, 0.3);
		ASSERT_GT(dt, units::Time{});
		for (int i = 0; i < 4; ++i) {
			auto result = solver.advance(patch, dt, physics::BoundaryConditions::periodic());
			EXPECT_NEAR(Real(result.timeInterval.duration() / dt), 1, 8 * epsilonR);
			for (int d = 0; d < ndim; ++d) EXPECT_EQ(result.faceFluxes[d].size(), patch.layout().faceCount(d));
		}
		patch.layout().forEachInterior([&](auto const& cell, auto) { test::expectStateNear(patch.atInterior(cell), this->state()); });
		EXPECT_EQ(patch.timeState().step, 4u);
		EXPECT_EQ(patch.timeState().time, 4.0 * dt);
	}
}


TYPED_TEST(FiniteVolume, PredictorTransformationUsesCellCentersAndPreservesUpdateBase) {
	auto patch = this->patch();
	physics::fillGhostCells(patch, physics::BoundaryConditions::periodic(), this->system);
	auto const original = patch.values();
	auto transform = [](typename TypeParam::State state, mesh::Coordinates const& storage) {
		Real phase = 0;
		for (int axis = 0; axis < ndim; ++axis) phase += (axis + 1) * (storage[axis] - 1.5) / 8;
		return state * (1 + 0.05 * std::sin(2 * piR * phase));
	};
	auto predictorPatch = patch;
	mesh::forEachCoordinate(patch.layout().extents(), [&](auto const& storage) {
		predictorPatch.atStorage(storage) = transform(patch.atStorage(storage), storage);
	});
	physics::MusclHancock<TypeParam> solver(this->system);
	auto const dt = 0.01 * solver.stableTimestep(predictorPatch, 0.3);
	typename physics::MusclHancock<TypeParam>::Workspace actualWorkspace, referenceWorkspace;
	auto actual = patch.values();
	auto reference = predictorPatch.values();
	solver.advanceInto(patch, dt, actualWorkspace,
		[&](auto const& cell, auto const& state) { actual[patch.layout().index(patch.layout().storageCoordinates(cell))] = state; }, transform);
	solver.advanceInto(predictorPatch, dt, referenceWorkspace,
		[&](auto const& cell, auto const& state) { reference[patch.layout().index(patch.layout().storageCoordinates(cell))] = state; });
	for (int axis = 0; axis < ndim; ++axis)
		for (std::size_t face = 0; face < actualWorkspace.fluxes[axis].size(); ++face)
			test::expectStateNear(actualWorkspace.fluxes[axis][face], referenceWorkspace.fluxes[axis][face]);
	patch.layout().forEachInterior([&](auto const& cell, auto index) {
		auto const expected = reference[index] - predictorPatch.atInterior(cell) + patch.atInterior(cell);
		test::expectStateNear(actual[index], expected);
	});
	for (std::size_t index = 0; index < original.size(); ++index)
		test::expectStateNear(patch.values()[index], original[index], 0);
}

#if OCTOTIGERII_HYDRO
TEST(FiniteVolumeMidpoint, NumericalFluxPredictorHasSecondOrderAtFixedMesh) {
	// A Hancock spatial predictor is not the derivative of the numerical
	// semidiscrete flux operator. The explicit-midpoint path must use that
	// operator and must not apply a second Hancock time prediction afterward.
	hydro::HydroSystem const gas(1.4);
	physics::MusclHancock<hydro::HydroSystem> const solver(gas);
	using Patch = mesh::PatchData<hydro::ConservedState>;
	auto integrate = [&](int steps) {
		Patch patch(mesh::MeshLayout(8, 2), units::Length::from_value(0.25));
		patch.layout().forEachInterior([&](auto const& cell, auto) {
			hydro::PrimitiveState primitive;
			primitive.density() = units::Density::from_value(1 + 0.1 * std::sin(2 * piR * (cell[0] + 0.5) / 8));
			primitive.pressure() = units::Pressure::from_value(2);
			primitive.velocity(0) = units::Velocity::from_value(0.2);
			patch.atInterior(cell) = gas.conservedState(primitive);
		});
		auto const dt = units::Time::from_value(0.04 / steps);
		for (int step = 0; step < steps; ++step) {
			physics::fillGhostCells(patch, physics::BoundaryConditions::periodic(), gas);
			physics::MusclHancock<hydro::HydroSystem>::Workspace first, second;
			solver.advanceInto(patch, {}, first, [](auto const&, auto const&) {});
			auto midpoint = patch;
			patch.layout().forEachInterior([&](auto const& cell, auto) {
				hydro::ConservedFlux divergence;
				for (int d = 0; d < ndim; ++d) {
					auto upper = cell; ++upper[d];
					divergence += first.fluxes[d][patch.layout().faceIndex(d, cell)] - first.fluxes[d][patch.layout().faceIndex(d, upper)];
				}
				midpoint.atInterior(cell) += (0.5 * dt / patch.cellWidth()) * divergence;
			});
			physics::fillGhostCells(midpoint, physics::BoundaryConditions::periodic(), gas);
			auto next = patch;
			solver.advanceInto(patch, dt, second,
				[&](auto const& cell, auto const& value) { next.atInterior(cell) = value; },
				[&](auto const&, auto const& cell) { return midpoint.atStorage(cell); }, false);
			patch = std::move(next);
		}
		return patch;
	};
	auto const reference = integrate(128);
	std::array<Real, 3> errors{};
	int run = 0;
	for (int steps : {2, 4, 8}) {
		auto const actual = integrate(steps);
		reference.layout().forEachInterior([&](auto const& cell, auto) {
			errors[run] += std::abs(units::value(actual.atInterior(cell).density() - reference.atInterior(cell).density()));
		});
		++run;
	}
	std::cout << "Fixed-mesh explicit-midpoint density errors: " << errors[0] << ", " << errors[1] << ", " << errors[2]
		<< "; orders " << std::log2(errors[0] / errors[1]) << ", " << std::log2(errors[1] / errors[2]) << '\n';
	EXPECT_GT(errors[0], 1e-10);
	EXPECT_LT(errors[1], 0.28 * errors[0]);
	EXPECT_LT(errors[2], 0.28 * errors[1]);
}
#endif


TYPED_TEST(FiniteVolume, NonuniformPeriodicEvolutionConservesEveryComponent) {
	for (auto limiter : {physics::Limiter::Minmod, physics::Limiter::VanLeer, physics::Limiter::MinmodTheta}) {
		auto patch = this->patch();
		patch.layout().forEachInterior([&](auto const& cell, auto) {
			Real phase = 0;
			for (int d = 0; d < ndim; ++d) phase += (d + 1) * (cell[d] + 0.5) / 8;
			patch.atInterior(cell) = this->state(1 + 0.1 * std::sin(2 * piR * phase));
		});
		typename TypeParam::State before{}, norm{}, after{};
		patch.layout().forEachInterior([&](auto const& cell, auto) {
			before += patch.atInterior(cell);
			norm += componentAbs(patch.atInterior(cell));
		});
		physics::MusclHancock<TypeParam> solver(this->system, limiter, 1.3);
		for (int i = 0; i < 4; ++i) solver.advance(patch, solver.stableTimestep(patch, 0.3), physics::BoundaryConditions::periodic());
		patch.layout().forEachInterior([&](auto const& cell, auto) {
			EXPECT_TRUE(this->system.admissible(patch.atInterior(cell)));
			after += patch.atInterior(cell);
		});
		before.forEach([&](auto f, auto q) {
			// Entropy is explicitly synchronized; Euler conserved fields still
			// obey flux conservation. Separate dual-energy tests disable sync.
			if constexpr (std::is_same_v<TypeParam, hydro::HydroSystem> && int(f) == ndim + 2) return;
			EXPECT_LE(units::abs(after.template get<f>() - q), 2e-12 * norm.template get<f>()) << "component " << int(f);
		});
	}
}


TYPED_TEST(FiniteVolume, GhostRulesCoverFacesEdgesAndCorners) {
	for (auto boundary : {physics::BoundaryCondition::Periodic, physics::BoundaryCondition::Outflow, physics::BoundaryCondition::Inflow,
			 physics::BoundaryCondition::Reflecting}) {
		auto patch = this->patch();
		patch.layout().forEachInterior([&](auto const& cell, auto) {
			patch.atInterior(cell) = this->state(1 + mesh::linearIndex(cell, mesh::filledCoordinates(8)));
		});
		physics::BoundaryConditions boundaries;
		boundaries.lower.fill(boundary);
		boundaries.upper.fill(boundary);
		physics::fillGhostCells(patch, boundaries, this->system);
		mesh::forEachCoordinate(patch.layout().extents(), [&](auto const& storage) {
			mesh::Coordinates source{};
			std::array<bool, ndim> reflect{};
			for (int d = 0; d < ndim; ++d) {
				int const x = storage[d] - 2;
				reflect[d] = x < 0 || x >= 8;
				if (boundary == physics::BoundaryCondition::Periodic)
					source[d] = (x + 8) % 8;
				else if (boundary == physics::BoundaryCondition::Outflow || boundary == physics::BoundaryCondition::Inflow)
					source[d] = std::clamp(x, 0, 7);
				else
					source[d] = x < 0 ? -1 - x : x >= 8 ? 15 - x : x;
			}
			auto expected = patch.atInterior(source);
			if (boundary == physics::BoundaryCondition::Reflecting)
				for (int d = 0; d < ndim; ++d)
					if (reflect[d]) expected = this->system.reflected(expected, d);
			if (boundary == physics::BoundaryCondition::Outflow)
				for (int d = 0; d < ndim; ++d)
					if (storage[d] < 2) {
						if constexpr (std::is_same_v<TypeParam, hydro::HydroSystem>)
							expected.momentum(d) = {};
						else
							expected.radiativeFlux(d) = {};
					}
			test::expectStateNear(patch.atStorage(storage), expected, 0);
		});
	}
}


TYPED_TEST(FiniteVolume, CflScalesWithCellWidthAndRejectsInvalidInputs) {
	auto patch = this->patch();
	physics::MusclHancock<TypeParam> solver(this->system);
	auto const dt = solver.stableTimestep(patch, 0.3);
	mesh::PatchData<typename TypeParam::State> wider(patch.layout(), 2.0 * patch.cellWidth());
	wider.values() = patch.values();
	EXPECT_EQ(solver.stableTimestep(wider, 0.3), 2.0 * dt);
	EXPECT_EQ(solver.stableTimestep(patch, 0.15), 0.5 * dt);
	for (Real cfl : {Real(0), Real(-1), Real(0.51), std::numeric_limits<Real>::quiet_NaN()})
		EXPECT_THROW(solver.stableTimestep(patch, cfl), std::invalid_argument);
	for (Real invalid : {Real(0), Real(-1), std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		EXPECT_THROW(solver.advance(patch, units::Time::from_value(invalid), physics::BoundaryConditions::periodic()), std::invalid_argument);
		EXPECT_EQ(patch.timeState().step, 0u);
		patch.layout().forEachInterior([&](auto const& cell, auto) { test::expectStateNear(patch.atInterior(cell), this->state(), 0); });
	}
	auto shallow = this->patch(1);
	EXPECT_THROW(solver.advance(shallow, dt, physics::BoundaryConditions::periodic()), std::invalid_argument);
	EXPECT_THROW((physics::MusclHancock<TypeParam>(this->system, physics::Limiter::MinmodTheta, 0.9)), std::invalid_argument);
	EXPECT_THROW((physics::MusclHancock<TypeParam>(this->system, physics::Limiter::MinmodTheta, 2.1)), std::invalid_argument);
}


TYPED_TEST(FiniteVolume, BoundaryUpdaterReceivesBeginningAndEndTimes) {
	auto patch = this->patch();
	physics::MusclHancock<TypeParam> solver(this->system);
	auto const dt = solver.stableTimestep(patch, 0.3);
	patch.timeState().time = units::Time::from_value(0.125);
	auto const start = patch.timeState().time;
	std::vector<units::Time> requested;
	auto result = solver.advanceWithBoundaryUpdater(patch, dt, [&](auto& data, auto time) {
		requested.push_back(time);
		physics::fillGhostCells(data, physics::BoundaryConditions::periodic(), this->system);
	});
	ASSERT_EQ(requested.size(), 2u);
	EXPECT_EQ(requested[0], start);
	EXPECT_EQ(requested[1], start + dt);
	EXPECT_EQ(result.timeInterval.begin, start);
	EXPECT_EQ(result.timeInterval.end, start + dt);
}

} // namespace
