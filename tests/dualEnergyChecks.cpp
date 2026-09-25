#include "testSupport.hpp"
#include "octotigerII/hydro/hydroSystem.hpp"
#include "octotigerII/amr/hierarchy.hpp"
#include "octotigerII/simulation.hpp"
#include <limits>

using namespace octotigerII;
using namespace octotigerII::hydro;

namespace {
PrimitiveState primitive(Real rho, Real pressure, Real speed = 0) {
	PrimitiveState p;
	p.density() = units::Density::from_value(rho);
	p.pressure() = units::Pressure::from_value(pressure);
	p.velocity(0) = units::Velocity::from_value(speed);
	return p;
}
void relative(Real a, Real b, Real tolerance = 3e-13) {
	EXPECT_NEAR(a, b, tolerance * std::abs(b));
}

TEST(DualEnergy, ArbitraryNonzeroExponentAndAdiabaticCompression) {
	for (Real alpha : {Real(-2), Real(-0.5), Real(0.25), Real(1), Real(1 / 1.4), Real(2)}) {
		Config::HydroOptions options;
		options.dualEnergy.exponent = alpha;
		HydroSystem gas(options);
		for (Real rho : {Real(1e-8), Real(2), Real(1e8)}) {
			auto u = gas.conservedState(primitive(rho, 3.0));
			Real const thermal = 3 / (options.gamma - 1);
			relative(units::value(u.auxiliary()), rho * std::pow(thermal / std::pow(rho, options.gamma), alpha));
			relative(units::value(gas.internalEnergyFromAuxiliary(u)), thermal);
			// Fixed A/rho under compression gives u proportional to rho^gamma.
			u.density() *= 3;
			u.auxiliary() *= 3;
			relative(units::value(gas.internalEnergyFromAuxiliary(u)), thermal * std::pow(3, options.gamma));
		}
	}
}

TEST(DualEnergy, PressureTemperatureAndSynchronizationUseIndependentThresholds) {
	Config::HydroOptions options;
	options.meanMolecularWeight = 0.6;
	HydroSystem gas(options);
	for (Real thermal : {Real(0.5), Real(2), Real(50), Real(101)}) {
		auto u = gas.conservedState(primitive(1, 0.1, std::sqrt(2 * (1000 - thermal))));
		u.totalEnergy() = units::EnergyDensity::from_value(1000);
		u.auxiliary() = gas.auxiliaryFromInternalEnergy(u.density(), units::EnergyDensity::from_value(0.25));
		auto const original = u;
		Real const expected = thermal > 1 ? thermal : 0.25;
		relative(units::value(gas.internalEnergy(u)), expected);
		relative(units::value(gas.reconstructionVariables(u).pressure()), (options.gamma - 1) * expected);
		relative(units::value(gas.temperature(u)), (options.gamma - 1) * expected * 0.6 * units::value(constants::atomicMassUnit) / units::value(constants::boltzmann));
		// Reconstructing from pressure must not overwrite the independent A.
		auto const reconstructed = gas.conservedState(gas.reconstructionVariables(u));
		EXPECT_EQ(reconstructed.auxiliary(), original.auxiliary());
		gas.synchronize(u);
		EXPECT_EQ(u.totalEnergy(), original.totalEnergy());
		EXPECT_EQ(u.density(), original.density());
		for (int d = 0; d < ndim; ++d) EXPECT_EQ(u.momentum(d), original.momentum(d));
		if (thermal > 100) relative(units::value(gas.internalEnergyFromAuxiliary(u)), thermal);
		else EXPECT_EQ(u.auxiliary(), original.auxiliary());
	}
}

TEST(DualEnergy, ThresholdEqualityUsesAuxiliaryAndDoesNotSynchronize) {
	Config::HydroOptions options;
	options.dualEnergy.pressureThreshold = 0.5;
	options.dualEnergy.syncThreshold = 0.75;
	HydroSystem gas(options);
	auto u = gas.conservedState(primitive(1, 1, 32));
	u.totalEnergy() = units::EnergyDensity::from_value(1024); // K=512 exactly
	relative(units::value(gas.internalEnergy(u)), 2.5);
	options.dualEnergy.pressureThreshold = 0.125;
	options.dualEnergy.syncThreshold = 0.5;
	gas = HydroSystem(options);
	auto const original = u;
	gas.synchronize(u);
	EXPECT_EQ(u.auxiliary(), original.auxiliary());
	u.totalEnergy() += units::EnergyDensity::from_value(1);
	gas.synchronize(u);
	relative(units::value(gas.internalEnergyFromAuxiliary(u)), 513);
}

TEST(DualEnergy, NonpositiveConservativeEnergyUsesAuxiliaryWithoutClipping) {
	HydroSystem gas(1.4);
	auto u = gas.conservedState(primitive(1, 1, 2));
	for (Real energy : {Real(0), Real(-1e6)}) {
		u.totalEnergy() = units::EnergyDensity::from_value(energy);
		auto const before = u;
		EXPECT_TRUE(gas.admissible(u));
		relative(units::value(gas.internalEnergy(u)), 2.5);
		gas.synchronize(u);
		test::expectStateNear(u, before, 0);
		EXPECT_GT(gas.temperature(u), units::Temperature{});
	}
}

TEST(DualEnergy, NegativeOrCancelledTotalThermalEnergyStillHasPositivePressure) {
	HydroSystem gas(1.4);
	auto u = gas.conservedState(primitive(1, 1, 1e10));
	ASSERT_EQ(units::value(u.totalEnergy()), 5e19); // Thermal part was rounded away.
	for (Real offset : {Real(0), Real(-16384)}) {
		u.totalEnergy() = units::EnergyDensity::from_value(5e19 + offset);
		EXPECT_TRUE(gas.admissible(u));
		relative(units::value(gas.reconstructionVariables(u).pressure()), 1);
		auto const a = u.auxiliary();
		gas.synchronize(u);
		EXPECT_EQ(u.auxiliary(), a);
	}
}

TEST(DualEnergy, DisabledModeUsesOnlyTotalEnergy) {
	Config::HydroOptions options;
	options.dualEnergy.enabled = false;
	HydroSystem gas(options);
	auto u = gas.conservedState(primitive(1, 2, 1));
	EXPECT_EQ(u.auxiliary(), units::Density{});
	relative(units::value(gas.internalEnergy(u)), 5);
	auto const before = u;
	gas.synchronize(u);
	test::expectStateNear(u, before, 0);
	u.totalEnergy() = units::EnergyDensity::from_value(0.5);
	u.auxiliary() = units::Density::from_value(100);
	EXPECT_FALSE(gas.admissible(u));
}

TEST(DualEnergy, InvalidAuxiliaryAndUnrepresentablePowersAreRejected) {
	HydroSystem gas(1.4);
	for (Real a : {Real(0), Real(-1), std::numeric_limits<Real>::infinity()}) {
		auto u = gas.conservedState(primitive(1, 1, 1e10));
		u.auxiliary() = units::Density::from_value(a);
		EXPECT_FALSE(gas.admissible(u));
	}
	Config::HydroOptions options;
	options.dualEnergy.exponent = 1e300;
	EXPECT_THROW(HydroSystem(options).conservedState(primitive(1, 10)), std::runtime_error);
	for (Real alpha : {Real(0), std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		options.dualEnergy.exponent = alpha;
		EXPECT_THROW(HydroSystem{options}, std::invalid_argument);
	}
}

TEST(DualEnergy, AuxiliaryFluxIsAdvectedThroughHllcContactAndReflection) {
	HydroSystem gas(1.4);
	auto l = gas.conservedState(primitive(1, 1, 0.2));
	auto r = gas.conservedState(primitive(2, 1, 0.2));
	l.auxiliary() *= 2;
	r.auxiliary() *= 3;
	auto const f = gas.riemann(l, r, 0);
	relative(units::value(f.auxiliary()), 0.2 * units::value(l.auxiliary()));
	EXPECT_EQ(gas.reflected(l, 0).auxiliary(), l.auxiliary());
	test::expectStateNear(HydroSystem::integratedFlux(HydroSystem::advectiveFlux(l, units::Velocity::from_value(2)), units::TimePerLength::from_value(0.5)), l);
}

TEST(DualEnergy, ColdPeriodicTransportConservesAuxiliaryWithoutSynchronization) {
	for (Real alpha : {Real(-0.5), Real(1), Real(2)}) {
		Config::HydroOptions options;
		options.dualEnergy.exponent = alpha;
		HydroSystem gas(options);
		Solver solver(gas);
		Fields patch(mesh::MeshLayout(8, 2), units::Length::from_value(0.125));
		ConservedState before{}, after{};
		patch.layout().forEachInterior([&](auto const& cell, auto) {
			auto u = gas.conservedState(primitive(1 + 0.01 * std::sin(2 * piR * (cell[0] + 0.5) / 8), 1, 1e10));
			patch.atInterior(cell) = u;
			before += u;
		});
		for (int step = 0; step < 5; ++step)
			solver.advance(patch, solver.stableTimestep(patch, 0.3), physics::BoundaryConditions::periodic());
		patch.layout().forEachInterior([&](auto const& cell, auto) {
			auto const u = patch.atInterior(cell);
			EXPECT_TRUE(gas.admissible(u));
			EXPECT_NEAR(units::value(gas.reconstructionVariables(u).pressure()), 1, 0.002);
			after += u;
		});
		test::expectStateNear(after, before, 2e-12);
	}
}

TEST(DualEnergy, IndependentAuxiliarySurvivesAdvectionBetweenThresholds) {
	HydroSystem gas(1.4);
	Solver solver(gas);
	Fields patch(mesh::MeshLayout(8, 2), units::Length::from_value(0.125));
	units::Density before{}, after{};
	patch.layout().forEachInterior([&](auto const& cell, auto) {
		auto u = gas.conservedState(primitive(1 + 0.01 * std::sin(2 * piR * (cell[0] + 0.5) / 8), 2, 50));
		u.auxiliary() *= 0.5;
		patch.atInterior(cell) = u;
		before += u.auxiliary();
	});
	for (int step = 0; step < 4; ++step)
		solver.advance(patch, solver.stableTimestep(patch, 0.3), physics::BoundaryConditions::periodic());
	patch.layout().forEachInterior([&](auto const& cell, auto) {
		auto const u = patch.atInterior(cell);
		after += u.auxiliary();
		EXPECT_NEAR(units::value(gas.reconstructionVariables(u).pressure()), 2, 1e-8);
		EXPECT_LT(gas.internalEnergyFromAuxiliary(u), 0.6 * gas.internalEnergy(u));
	});
	relative(units::value(after), units::value(before), 2e-12);
}

TEST(DualEnergy, RefineAndRestrictPreserveIndependentAuxiliary) {
	auto c = octotigerII::parseConfig({"--problem.name=sod", "--mesh.cells=4", "--mesh.level=0", "--output.enabled=off"});
	c.hydro.dualEnergy.exponent = -0.5;
	HydroSystem gas(c.hydro);
	auto parent = initialSnapshot(c, {});
	ConservedState before{}, after{};
	parent.layout.forEachInterior([&](auto const& cell, std::size_t i) {
		auto u = gas.conservedState(primitive(1 + 0.03 * cell[0], 1, 1e8));
		u.auxiliary() *= 0.7;
		parent.hydro.values()[i] = u;
		before += u;
	});
	amr::Hierarchy hierarchy(c, {parent});
	std::vector<Snapshot> children;
	for (int child = 0; child < (1 << ndim); ++child) children.push_back(hierarchy.transfer(mesh::BlockLocation{}.child(child)));
	for (auto const& b : children) b.layout.forEachInterior([&](auto const&, std::size_t i) {
		EXPECT_TRUE(gas.admissible(b.hydro.values()[i]));
		after += b.hydro.values()[i] / Real(1 << ndim);
	});
	test::expectStateNear(after, before, 3e-12);
	auto restored = amr::Hierarchy(c, children).transfer({});
	parent.layout.forEachInterior([&](auto const&, std::size_t i) {
		test::expectStateNear(restored.hydro.values()[i], parent.hydro.values()[i], 3e-12);
	});
}

TEST(DualEnergy, WarmOwningPatchSynchronizesAfterEveryStep) {
	HydroSystem gas(1.4);
	Solver solver(gas);
	Fields patch(mesh::MeshLayout(4, 2), units::Length::from_value(0.25));
	auto u = gas.conservedState(primitive(1, 1));
	u.auxiliary() *= 0.5;
	patch.layout().forEachInterior([&](auto const& cell, auto) { patch.atInterior(cell) = u; });
	solver.advance(patch, solver.stableTimestep(patch, 0.3), physics::BoundaryConditions::periodic());
	patch.layout().forEachInterior([&](auto const& cell, auto) {
		relative(units::value(gas.internalEnergyFromAuxiliary(patch.atInterior(cell))), 2.5);
	});
}

TEST(DualEnergy, RuntimeSynchronizesAfterRefluxAtEachTimestep) {
	for (bool adaptive : {false, true}) {
		auto c = octotigerII::parseConfig({"--problem.name=sod", "--mesh.cells=4", "--mesh.level=1", "--output.enabled=off"});
		c.hydro.dualEnergy.exponent = -0.5;
		c.amr.enabled = adaptive;
		c.amr.maxLevel = 2;
		c.amr.shadowTolerance = 0;
		c.amr.bufferCells = 0;
		Runtime runtime(c, {[](refinement::CellView const& v) { return v.center[0] < units::Length::from_value(0.25) && v.level < 2 ? Real(2) : Real(0); }});
		if (adaptive) {
			auto const blocks = runtime.snapshots();
			bool coarse = false, fine = false;
			for (auto const& b : blocks) { coarse |= b.location.level == 1; fine |= b.location.level == 2; }
			ASSERT_TRUE(coarse && fine);
		}
		HydroSystem gas(c.hydro);
		for (int step = 0; step < 3; ++step) {
			runtime.advance(runtime.stableTimestep());
			for (auto const& block : runtime.snapshots()) block.layout.forEachInterior([&](auto const&, std::size_t i) {
				auto const& u = block.hydro.values()[i];
				units::EnergyDensity kinetic{};
				for (int d = 0; d < ndim; ++d) kinetic += 0.5 * u.momentum(d) * u.momentum(d) / u.density();
				auto const thermal = u.totalEnergy() - kinetic;
				if (thermal > c.hydro.dualEnergy.syncThreshold * u.totalEnergy())
					relative(units::value(gas.internalEnergyFromAuxiliary(u)), units::value(thermal));
			});
		}
	}
}
} // namespace
