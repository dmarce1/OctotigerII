#include <gtest/gtest.h>
#include <set>
#include "octotigerII/amr/hierarchy.hpp"
#include "octotigerII/runtime.hpp"
#include "octotigerII/simulation.hpp"
#include "octotigerII/subgrid/view.hpp"
#include "testSupport.hpp"
#if OCTOTIGERII_GRAVITY
#include "octotigerII/verification/directGravity.hpp"
#endif

using namespace octotigerII;
namespace {
Config configuration() {
	auto c = parseConfig({"--mesh.cells=4", "--mesh.level=1", "--output.enabled=off"});
	c.mesh.lower = units::Length::from_value(-0.5);
	c.mesh.upper = units::Length::from_value(0.5);
	c.amr.enabled = true;
	c.amr.maxLevel = 2;
	c.amr.shadowTolerance = 0;
	c.amr.bufferCells = 0;
	c.amr.regridEvery = 1;
	c.mesh.boundary = physics::BoundaryConditions::periodic();
	return c;
}
refinement::Criterion corner(Config const& c, bool* enabled = nullptr) {
	return [c, enabled](refinement::CellView const& cell) {
		if (enabled && !*enabled) return Real(0);
		for (int d = 0; d < ndim; ++d)
			if (cell.center[d] < c.mesh.lower + 0.12 * (c.mesh.upper - c.mesh.lower) || cell.center[d] >= c.mesh.lower + 0.26 * (c.mesh.upper - c.mesh.lower))
				return Real(0);
		return cell.level < c.amr.maxLevel ? Real(2) : Real(0.5);
	};
}
std::vector<Snapshot> uniform(Config const& c) {
	std::vector<Snapshot> result;
	mesh::forEachCoordinate(mesh::filledCoordinates(1 << c.mesh.level), [&](auto const& x) { result.push_back(initialSnapshot(c, {c.mesh.level, x})); });
	return result;
}
void compareTotals(Diagnostics const& a, Diagnostics const& b) {
	EXPECT_NEAR(units::value(a.mass), units::value(b.mass), 5e-13 * std::max(Real(1), units::value(a.mass)));
	EXPECT_NEAR(units::value(a.gasEnergy), units::value(b.gasEnergy), 5e-13 * std::max(Real(1), units::value(a.gasEnergy)));
	EXPECT_NEAR(units::value(a.radiationEnergy), units::value(b.radiationEnergy), 5e-13 * std::max(Real(1), units::value(a.radiationEnergy)));
	for (int d = 0; d < ndim; ++d)
		EXPECT_NEAR(units::value(a.momentum[d]), units::value(b.momentum[d]), 5e-13 * std::max(Real(1), units::value(a.mass)));
}

TEST(Amr, OptionsValidateAndRoundTrip) {
	auto c = parseConfig({"--amr.enabled=on", "--amr.minLevel=0", "--amr.maxLevel=3", "--amr.regridEvery=7", "--amr.shadowTolerance=0.03",
		"--amr.signalBuffer=1.5", "--amr.bufferCells=2"});
	EXPECT_TRUE(c.amr.enabled);
	EXPECT_EQ(c.amr.maxLevel, 3);
	EXPECT_EQ(c.amr.regridEvery, 7);
	EXPECT_DOUBLE_EQ(c.amr.shadowTolerance, 0.03);
	for (auto const* option : {"--amr.regridEvery=0", "--amr.maxLevel=17", "--amr.signalBuffer=0.5", "--amr.coarsenFactor=1", "--amr.shadowFloor=0"})
		EXPECT_THROW(parseConfig({option}), std::invalid_argument);
}

TEST(Amr, MortonDestinationAndConservativeRefineCoarsen) {
	auto c = configuration();
	auto const initial = diagnose(uniform(c), c);
	bool enabled = true;
	Runtime runtime(c, {corner(c, &enabled)});
	auto refined = runtime.snapshots();
	ASSERT_GT(refined.size(), std::size_t(1 << ndim));
	std::set<int> levels;
	std::uint64_t previous = 0;
	for (auto const& block : refined) {
		levels.insert(block.location.level);
		EXPECT_GE(mesh::mortonKey(block.location), previous);
		previous = mesh::mortonKey(block.location);
	}
	EXPECT_EQ(levels.size(), 2u);
	EXPECT_GT(runtime.shadowCellCount(), refined.size() * refined.front().layout.interiorCellCount());
	compareTotals(initial, diagnose(refined, c));
	enabled = false;
	EXPECT_TRUE(runtime.regrid(units::Time{}, true));
	EXPECT_EQ(runtime.size(), std::size_t(1 << ndim));
	compareTotals(initial, diagnose(runtime.snapshots(), c));
}

TEST(Amr, FailedCriterionPreservesPublishedMeshAndState) {
	auto c = configuration();
	bool fail = false;
	Runtime runtime(c, {[&](refinement::CellView const&) -> Real {
		if (fail) throw std::runtime_error("criterion failure");
		return 0;
	}});
	auto const before = diagnose(runtime.snapshots(), c);
	auto const generation = runtime.generation();
	auto const count = runtime.size();
	fail = true;
	EXPECT_THROW(runtime.regrid(units::Time{}, true), std::runtime_error);
	EXPECT_EQ(runtime.generation(), generation);
	EXPECT_EQ(runtime.size(), count);
	compareTotals(before, diagnose(runtime.snapshots(), c));
}

TEST(Amr, DeepRefinementBalancesFacesEdgesAndPeriodicCorners) {
	using std::abs;
	auto c = configuration();
	c.amr.maxLevel = 4;
	auto leaves = uniform(c);
	amr::Hierarchy const source(c, leaves);
	// A tag at the domain corner forces balancing across periodic faces,
	// edges, and corners, through more than one level of refinement.
	refinement::Criteria criteria{[&](refinement::CellView const& cell) {
		for (int d = 0; d < ndim; ++d)
			if (cell.center[d] > c.mesh.lower + 0.08 * (c.mesh.upper - c.mesh.lower)) return Real(0);
		return Real(2);
	}};
	for (int pass = 0; pass < c.amr.maxLevel; ++pass) {
		amr::Hierarchy const current(c, leaves);
		auto const next = amr::selectMesh(c, leaves, current, {}, criteria, false);
		if (!next.changed) break;
		leaves.clear();
		for (auto const& location : next.leaves)
			leaves.push_back(source.transfer(location));
	}
	bool finest = false, coarse = false;
	for (auto const& a : leaves) {
		finest = finest || a.location.level == c.amr.maxLevel;
		coarse = coarse || a.location.level < c.amr.maxLevel - 1;
		for (auto const& b : leaves) {
			bool touching = true;
			for (int d = 0; d < ndim; ++d) {
				bool adjacent = false;
				Real const al = Real(a.location.coordinates[d]) / (1 << a.location.level);
				Real const au = Real(a.location.coordinates[d] + 1) / (1 << a.location.level);
				Real const bl = Real(b.location.coordinates[d]) / (1 << b.location.level);
				Real const bu = Real(b.location.coordinates[d] + 1) / (1 << b.location.level);
				for (int shift = -1; shift <= 1; ++shift)
					adjacent = adjacent || (al <= bu + shift && bl + shift <= au);
				touching = touching && adjacent;
			}
			if (touching) EXPECT_LE(abs(a.location.level - b.location.level), 1);
		}
	}
	EXPECT_TRUE(finest);
	EXPECT_TRUE(coarse);
}

TEST(Amr, MassThresholdRespectsFinestLevel) {
	if constexpr (!build::hydro && !build::gravity) {
		GTEST_SKIP();
	} else {
		auto c = configuration();
		auto leaves = uniform(c);
		units::Mass maximum{};
		for (auto const& block : leaves)
			block.layout.forEachInterior([&](auto const&, std::size_t i) {
				auto const density = build::hydro ? block.hydro.values()[i].density() : block.density.values()[i];
				maximum = std::max(maximum, density * block.layout.cellMeasure(block.cellWidth));
			});
		c.amr.maxCellMass = maximum / Real((1 << ndim) * 2);
		Runtime runtime(c);
		bool capped = false;
		for (auto const& block : runtime.snapshots()) {
			EXPECT_LE(block.location.level, c.amr.maxLevel);
			block.layout.forEachInterior([&](auto const&, std::size_t i) {
				auto const density = build::hydro ? block.hydro.values()[i].density() : block.density.values()[i];
				if (density * block.layout.cellMeasure(block.cellWidth) > c.amr.maxCellMass) {
					EXPECT_EQ(block.location.level, c.amr.maxLevel);
					capped = true;
				}
			});
		}
		EXPECT_TRUE(capped);
	}
}

TEST(Amr, SignalBufferWrapsPeriodicBoundaryAndGrowsWithHorizon) {
	auto c = configuration();
	if constexpr (!build::hydro && !build::radiation) {
		GTEST_SKIP();
	} else {
		auto leaves = uniform(c);
		amr::Hierarchy hierarchy(c, leaves);
		refinement::Criteria criteria{corner(c)};
		auto const small = amr::selectMesh(c, leaves, hierarchy, {}, criteria, false);
		auto const view = hierarchy.sample(leaves.front().location, mesh::Coordinates{});
		auto speed = view.signalSpeed[0];
		ASSERT_GT(speed, units::Velocity{});
		auto const large = amr::selectMesh(c, leaves, hierarchy, 0.7 * (c.mesh.upper - c.mesh.lower) / speed, criteria, false);
		EXPECT_GT(large.leaves.size(), small.leaves.size());
		bool opposite = false;
		for (auto const& leaf : large.leaves)
			if (leaf.level == 2 && leaf.coordinates[0] == 3) opposite = true;
		EXPECT_TRUE(opposite);
	}
}

TEST(Amr, MixedLevelTransportConservesPeriodicIntegrals) {
	if constexpr (!build::hydro && !build::radiation) {
		GTEST_SKIP();
	} else {
		auto c = configuration();
		Runtime runtime(c, {corner(c)});
		auto const initial = diagnose(runtime.snapshots(), c);
		for (int i = 0; i < 4; ++i)
			runtime.advance(0.3 * runtime.stableTimestep());
		compareTotals(initial, diagnose(runtime.snapshots(), c));
	}
}

TEST(Amr, RegridCadenceAndTravelBudget) {
	if constexpr (!build::hydro && !build::radiation) {
		GTEST_SKIP();
	} else {
		auto c = configuration();
		c.amr.regridEvery = 3;
		int calls = 0;
		Runtime runtime(c, {[&](refinement::CellView const&) {
			++calls;
			return Real(0);
		}});
		calls = 0;
		for (int step = 0; step < 2; ++step) {
			auto const dt = 0.1 * runtime.stableTimestep();
			runtime.advance(dt);
			EXPECT_FALSE(runtime.regrid(dt));
			EXPECT_EQ(calls, 0);
		}
		auto const dt = 0.1 * runtime.stableTimestep();
		runtime.advance(dt);
		runtime.regrid(dt);
		EXPECT_GT(calls, 0);
		calls = 0;
		runtime.regrid(10.0 * runtime.stableTimestep());
		EXPECT_GT(calls, 0);
	}
}

TEST(Amr, EvolvedShadowDiffersFromFreshRestriction) {
	if constexpr (!build::hydro && !build::radiation) {
		GTEST_SKIP();
	} else {
		auto c = configuration();
		auto initial = uniform(c);
		// A smooth translating density/radiation profile gives distinct coarse
		// and fine truncation errors without relying on a discontinuity.
		for (auto& block : initial)
			block.layout.forEachInterior([&](auto const& cell, std::size_t i) {
				using std::sin;
				auto const x = Real((block.layout.cellCenter(block.lower, block.cellWidth, cell)[0] - c.mesh.lower) / (c.mesh.upper - c.mesh.lower));
				Real const wave = 1 + 0.2 * sin(2 * piR * x);
				if constexpr (build::hydro) {
					hydro::PrimitiveState q;
					q.density() = units::Density::from_value(wave);
					q.pressure() = units::Pressure::from_value(1);
					q.velocity(0) = units::Velocity::from_value(0.5);
					block.hydro.values()[i] = hydro::HydroSystem(c.hydro.gamma).conservedState(q);
				}
				if constexpr (build::radiation) {
					block.radiation.values()[i].energy() = units::EnergyDensity::from_value(wave);
					for (int d = 0; d < ndim; ++d)
						block.radiation.values()[i].radiativeFlux(d) = {};
					block.radiation.values()[i].radiativeFlux(0) = constants::c * block.radiation.values()[i].energy();
				}
			});
		amr::Hierarchy evolved(c, initial);
		auto const dt = units::Time::from_value(build::radiation ? 1e-13 : 0.001);
		auto const old = evolved.shadow({evolved.cellLevel(c.mesh.level), mesh::filledCoordinates(0)});
		evolved.advance(dt);
		auto const after = evolved.shadow({evolved.cellLevel(c.mesh.level), mesh::filledCoordinates(0)});
		Real difference = 0;
		if constexpr (build::hydro) difference += units::value(units::abs(after.hydro.density() - old.hydro.density()));
		if constexpr (build::radiation) difference += units::value(units::abs(after.radiation.energy() - old.radiation.energy()));
		EXPECT_GT(difference, 1e-10);
		evolved.refreshLeaves(initial);
		auto const preserved = evolved.shadow({evolved.cellLevel(c.mesh.level), mesh::filledCoordinates(0)});
		test::expectStateNear(preserved.hydro, after.hydro);
		test::expectStateNear(preserved.radiation, after.radiation);
	}
}

#if OCTOTIGERII_GRAVITY
TEST(Amr, AdaptiveFmmMatchesDirectAcrossMixedLevelsAndImages) {
	for (int boundary = 0; boundary < 7; ++boundary) {
		// HPX may resume this test on a different OS thread after an action.
		// Keep case context on assertions instead of GoogleTest's TLS trace stack.
		auto c = configuration();
		c.gravity.multipoleOrder = 4;
		c.gravity.openingAngle = 0.3;
		c.verification.directSamples = 16;
		c.mesh.boundary = {};
		int const periodicAxes = boundary <= 3 ? boundary : (boundary == 4 ? 0 : boundary - 4);
		for (int d = 0; d < periodicAxes; ++d)
			c.mesh.boundary.lower[d] = c.mesh.boundary.upper[d] = physics::BoundaryCondition::Periodic;
		if (boundary >= 4) c.mesh.boundary.lower[2] = physics::BoundaryCondition::Reflecting;
		if (boundary == 6) c.mesh.boundary.upper[2] = physics::BoundaryCondition::Reflecting;
		Runtime runtime(c, {corner(c)});
		std::set<int> levels;
		for (auto const& block : runtime.snapshots())
			levels.insert(block.location.level);
		ASSERT_EQ(levels.size(), 2u) << "boundary " << boundary;
		auto const stats = runtime.solveGravity();
		EXPECT_GT(stats.multipolePairs, 0u) << "boundary " << boundary;
		if (periodicAxes > 0 || boundary == 6) EXPECT_GT(stats.ewaldPairs, 0u) << "boundary " << boundary;
		if (boundary >= 4) EXPECT_GT(stats.reflectedPairs, 0u) << "boundary " << boundary;
		auto const result = verification::compareDirectGravity(runtime.snapshots(), c);
		for (auto const& field : result.fields)
			EXPECT_LT(field.l2, 0.003 * field.referenceL2 + 1e-13) << "boundary " << boundary << ", " << field.name;
	}
}
#endif
}	 // namespace
