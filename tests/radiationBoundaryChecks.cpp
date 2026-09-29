#include "testSupport.hpp"
#include "octotigerII/radiation/diffusionFlux.hpp"
#include "octotigerII/runtime.hpp"
#include "octotigerII/simulation.hpp"
#include <cmath>
#include <set>

using namespace octotigerII;
namespace {
using System = radiation::RadiationSystem;
using Flux = System::Flux;
using State = System::State;

Config configuration(bool adaptive = false) {
	auto config = parseConfig({"--problem.name=radiation-matter", "--mesh.cells=4", "--mesh.level=1",
		"--mesh.periodic=off", "--radiation.closedBoundary=on", "--radiation.lightSpeedRatio=1",
		"--radiation.opacity=10000", "--runtime.stopTime=0", "--output.enabled=off"});
	config.amr.enabled = adaptive;
	config.amr.minLevel = 1;
	config.amr.maxLevel = 2;
	config.amr.bufferCells = 0;
	config.amr.shadowTolerance = 0;
	config.validate();
	return config;
}

auto apCorrection(System const& system, radiation::MaterialVelocity const& velocity) {
	return [&, velocity](Flux const& flux, State const& left, State const& right, State const& faceLeft,
		State const& faceRight, auto const&, auto const&, int axis, units::Length width, units::Velocity speed) {
		return radiation::diffusionCorrectedFlux(system, flux, left, right, faceLeft, faceRight,
			radiation::Extinction::from_value(1e5), radiation::Extinction::from_value(1e5),
			velocity, velocity, axis, width, speed);
	};
}
}

TEST(RadiationBoundary, OptionIsExplicitAndRejectsUnsupportedGeometry) {
	auto config = configuration();
	EXPECT_TRUE(config.radiation.closedBoundary);
	EXPECT_TRUE(config.mesh.boundary.all(finiteVolume::BoundaryCondition::Outflow));
	EXPECT_EQ(config.frame.omega, units::InverseTime{});
	auto baseline = parseConfig({"--problem.name=radiation-matter", "--output.enabled=off"});
	EXPECT_FALSE(baseline.radiation.closedBoundary);
	EXPECT_THROW(parseConfig({"--problem.name=sod", "--radiation.closedBoundary=on"}), std::invalid_argument);
	EXPECT_THROW(parseConfig({"--problem.name=radiation-matter", "--radiation.closedBoundary=on"}), std::invalid_argument);
	EXPECT_THROW(parseConfig({"--problem.name=radiation-matter", "--radiation.closedBoundary=yes"}), std::invalid_argument);
	if constexpr (ndim >= 2) {
		config.frame.omega = units::InverseTime::from_value(1);
		EXPECT_THROW(config.validate(), std::invalid_argument);
	}
}

TEST(RadiationBoundary, InsulationOverridesApAdvectionAndPreservesPressureFluxInProbes) {
	auto const width = units::Length::from_value(1);
	System const closed(constants::c, {true, {}, 4.0 * width});
	System const open(constants::c);
	radiation::Fields patch(mesh::MeshLayout(4, 2), width);
	State state;
	state.energy() = units::EnergyDensity::from_value(1);
	std::fill(patch.values().begin(), patch.values().end(), state);
	radiation::MaterialVelocity velocity{};
	velocity[0] = Real(.01) * constants::c;
	auto const dt = Real(1e-5) * width / constants::c;
	auto unchanged = [](State input, auto const&) { return input; };
	auto discard = [](auto const&, State const&) {};
	radiation::Solver::Workspace closedWork, openWork;
	radiation::Solver(closed).advanceInto(patch, {}, closedWork, discard, unchanged, false, apCorrection(closed, velocity), dt);
	radiation::Solver(open).advanceInto(patch, {}, openWork, discard, unchanged, false, apCorrection(open, velocity), dt);
	for (int axis = 0; axis < ndim; ++axis) {
		mesh::forEachCoordinate(patch.layout().faceExtents(axis), [&](auto const& face) {
			auto const index = patch.layout().faceIndex(axis, face);
			auto const& flux = closedWork.fluxes[axis][index];
			bool const boundary = face[axis] == 0 || face[axis] == 4;
			if (boundary) EXPECT_EQ(flux.energy(), units::EnergyFlux{});
			else EXPECT_EQ(flux.energy(), openWork.fluxes[axis][index].energy());
			for (int d = 0; d < ndim; ++d) {
				EXPECT_EQ(flux.radiativeFlux(d), openWork.fluxes[axis][index].radiativeFlux(d));
			}
			if (axis == 0) EXPECT_GT(openWork.fluxes[axis][index].energy(), units::EnergyFlux{});
		});
	}
	// Predictor patches extend beyond the physical box. Constraint detection
	// must use physical face coordinates rather than their patch-edge indices.
	mesh::PhysicalCoordinates lower{};
	for (auto& x : lower) { x = -2.0 * width; }
	radiation::Fields padded(mesh::MeshLayout(8, 2), width, lower);
	std::fill(padded.values().begin(), padded.values().end(), state);
	radiation::Solver(closed).advanceInto(padded, {}, closedWork, discard, unchanged, false, apCorrection(closed, velocity), dt);
	mesh::Coordinates face{};
	for (int index : {0, 2, 6, 8}) {
		face[0] = index;
		auto const flux = closedWork.fluxes[0][padded.layout().faceIndex(0, face)].energy();
		if (index == 2 || index == 6) EXPECT_EQ(flux, units::EnergyFlux{});
		else EXPECT_GT(flux, units::EnergyFlux{});
	}
}

TEST(RadiationBoundary, RealizabilityFallbackCannotReintroduceEnergyLeakage) {
	System system(constants::c);
	State left, right;
	left.energy() = right.energy() = units::EnergyDensity::from_value(1);
	left.radiativeFlux(0) = Real(.1) * constants::c * left.energy();
	right.radiativeFlux(0) = -left.radiativeFlux(0);
	Flux high = system.riemann(left, right, 0);
	high.energy() = {};
	high.radiativeFlux(0) *= 1e4;
	Flux low = system.lowOrderFlux(left, right, 0);
	low.energy() = {};
	auto const fraction = units::TimePerLength::from_value(.01 / units::value(constants::c));
	auto const limited = system.limitFlux(left, right, high, 0, fraction, {}, &low, true);
	EXPECT_EQ(limited.energy(), units::EnergyFlux{});
	EXPECT_LT(units::abs(limited.radiativeFlux(0)), units::abs(high.radiativeFlux(0)));
	// Force rejection of the supplied low-order momentum candidate as well.
	auto const fallback = system.limitFlux(left, right, high, 0, fraction, {}, &high, true);
	EXPECT_EQ(fallback.energy(), units::EnergyFlux{});
	EXPECT_LT(units::abs(fallback.radiativeFlux(0)), units::abs(high.radiativeFlux(0)));
	EXPECT_GT(fallback.radiativeFlux(0), units::EnergyFluxTransport{});
}

TEST(RadiationBoundary, CoupledRuntimeClosesEnergyAndMomentumLedgersWithFreeGasBoundaries) {
	for (bool adaptive : {false, true}) {
		auto config = configuration(adaptive);
		bool refine = adaptive;
		refinement::Criteria criteria{[&](refinement::CellView const& cell) {
			return refine && cell.center[0] < config.mesh.lower + Real(.25) * (config.mesh.upper - config.mesh.lower)
				&& cell.level < 2 ? Real(2) : Real(0);
		}};
		Runtime runtime(config, criteria);
		if (adaptive) {
			std::set<int> levels;
			for (auto const& block : runtime.snapshots()) { levels.insert(block.location.level); }
			ASSERT_EQ(levels.size(), 2u);
		}
		auto const before = diagnose(runtime.snapshots(), config);
		for (int step = 0; step < 3; ++step) { runtime.advanceCoupled(Real(.2) * runtime.stableTimestep()); }
		if (adaptive) {
			// Exercise the shadow predictor and regridding boundary path as well.
			refine = false;
			runtime.regrid(units::Time{}, true);
			refine = true;
			runtime.regrid(units::Time{}, true);
			runtime.advanceCoupled(Real(.1) * runtime.stableTimestep());
		}
		auto const after = diagnose(runtime.snapshots(), config);
		auto const boundary = runtime.boundaryTransport();
		EXPECT_EQ(boundary.outward.radiationEnergy, units::Energy{});
		EXPECT_EQ(boundary.inward.radiationEnergy, units::Energy{});
		EXPECT_EQ(runtime.radiationSourceEnergy(), units::Energy{});
		EXPECT_GT(boundary.outward.mass, units::Mass{});
		auto const gasEnergy = boundary.outward.gasEnergy - boundary.inward.gasEnergy;
		EXPECT_NEAR(Real((after.physicalTotalEnergy + gasEnergy - before.physicalTotalEnergy) / before.physicalTotalEnergyNorm), 0, 2e-12);
		EXPECT_NEAR(Real((after.mass + boundary.outward.mass - boundary.inward.mass - before.mass) / before.mass), 0, 2e-12);
		for (int d = 0; d < ndim; ++d) {
				auto const momentum = boundary.outward.momentum[d] - boundary.inward.momentum[d]
					+ (boundary.outward.radiationFlux[d] - boundary.inward.radiationFlux[d]) / (constants::c * constants::c);
				// A transverse component can be exactly zero in every cell while
				// opposing wall-pressure impulses are large. Scale their cancellation
				// by the separately accepted fluxes, not their roundoff-sized net.
				auto const transportedNorm = boundary.outward.momentum[d] + boundary.inward.momentum[d]
					+ (boundary.outward.radiationFlux[d] + boundary.inward.radiationFlux[d]) / (constants::c * constants::c);
				auto const scale = std::max({before.physicalTotalMomentumNorm[d], after.physicalTotalMomentumNorm[d], transportedNorm});
				if (scale > units::Momentum{})
					EXPECT_NEAR(Real((after.physicalTotalMomentum[d] + momentum - before.physicalTotalMomentum[d]) / scale), 0, 3e-12)
						<< "adaptive=" << adaptive << " axis=" << d;
		}
	}
}
