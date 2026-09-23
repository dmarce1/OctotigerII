#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include "octotigerII/problems.hpp"
#include "octotigerII/simulation.hpp"
#include "testSupport.hpp"

using namespace octotigerII;

namespace {

std::vector<hydro::ConservedState> flatten(std::vector<Snapshot> const& snapshots, int count) {
	std::vector<hydro::ConservedState> result(mesh::MeshLayout(count).interiorCellCount());
	for (auto const& block : snapshots)
		block.layout.forEachInterior([&](auto cell, auto i) {
			for (int d = 0; d < ndim; ++d)
				cell[d] += block.location.coordinates[d] * block.layout.cellsPerActiveDimension();
			result[mesh::linearIndex(cell, mesh::filledCoordinates(count))] = block.hydro.values()[i];
		});
	return result;
}

TEST(RayleighTaylor, DefaultsAndPhysicalParameterValidation) {
	auto c = parseConfig({});
	EXPECT_TRUE(build::hydro);
	EXPECT_FALSE(build::gravity);
	for (int d = 0; d < ndim - 1; ++d)
		EXPECT_TRUE(c.mesh.boundary.periodic(d));
	EXPECT_EQ(c.mesh.boundary.lower[ndim - 1], physics::BoundaryCondition::Reflecting);
	EXPECT_EQ(c.mesh.boundary.upper[ndim - 1], physics::BoundaryCondition::Reflecting);
	EXPECT_LT(c.mesh.lower, units::Length{});
	EXPECT_LT(c.hydro.acceleration[ndim - 1], units::Acceleration{});
	for (auto argument : {"--rayleighTaylor.densityLower=0", "--rayleighTaylor.densityUpper=0.5", "--rayleighTaylor.interfacePressure=0.01",
			 "--rayleighTaylor.perturbation=-1", "--rayleighTaylor.perturbation=nan"})
		EXPECT_THROW(parseConfig({argument}), std::invalid_argument);
	auto const key = std::string("--hydro.acceleration.") + "xyz"[ndim - 1];
	EXPECT_THROW(parseConfig({key + "=0.1"}), std::invalid_argument);
	EXPECT_THROW(parseConfig({key + "=nan"}), std::invalid_argument);
	EXPECT_NO_THROW(parseConfig({key + "=0", "--rayleighTaylor.perturbation=0"}));
}

TEST(RayleighTaylor, InitialLayersHaveHydrostaticPressureAndZeroMeanSeed) {
	auto c = parseConfig({"--mesh.cells=8", "--mesh.level=0", "--mesh.lower=-1", "--mesh.upper=1", "--rayleighTaylor.densityLower=1.25",
		"--rayleighTaylor.densityUpper=2.75", "--rayleighTaylor.interfacePressure=3.5", "--rayleighTaylor.perturbation=0.02"});
	auto block = initialSnapshot(c, {});
	hydro::HydroSystem gas(c.hydro.gamma);
	Real meanVelocity = 0, seedSquared = 0;
	block.layout.forEachInterior([&](auto cell, auto i) {
		auto const p = gas.reconstructionVariables(block.hydro.values()[i]);
		EXPECT_EQ(p.density(), cell[ndim - 1] < 4 ? c.rayleighTaylor.densityLower : c.rayleighTaylor.densityUpper);
		for (int d = 0; d < ndim - 1; ++d)
			EXPECT_EQ(p.velocity(d), units::Velocity{});
		meanVelocity += units::value(p.velocity(ndim - 1));
		seedSquared += units::value(p.velocity(ndim - 1)) * units::value(p.velocity(ndim - 1));
		if (cell[ndim - 1] != 3 && cell[ndim - 1] != 7) {
			++cell[ndim - 1];
			auto const next = gas.reconstructionVariables(block.hydro.atInterior(cell));
			auto const gradient = (next.pressure() - p.pressure()) / block.cellWidth;
			EXPECT_NEAR(units::value(gradient / p.density()), units::value(c.hydro.acceleration[ndim - 1]), 2e-14);
		}
	});
	EXPECT_NEAR(meanVelocity, 0, 1e-13);
	EXPECT_GT(seedSquared, 0);
}

TEST(RayleighTaylor, ExternalKickPreservesInternalEnergyAndRestrictsTimestep) {
	auto c = parseConfig({"--mesh.cells=4", "--mesh.level=1", "--output.enabled=off"});
	Runtime runtime(c);
	auto const before = flatten(runtime.snapshots(), 8);
	auto const dt = units::Time::from_value(0.02);
	hydro::HydroSystem gas(c.hydro.gamma);
	runtime.kickGravity(dt);
	auto const after = flatten(runtime.snapshots(), 8);
	for (std::size_t i = 0; i < before.size(); ++i) {
		auto p = gas.reconstructionVariables(before[i]), q = gas.reconstructionVariables(after[i]);
		EXPECT_EQ(q.density(), p.density());
		EXPECT_NEAR(units::value(q.pressure()), units::value(p.pressure()), 2e-14);
		for (int d = 0; d < ndim; ++d)
			EXPECT_NEAR(units::value(q.velocity(d) - p.velocity(d)), units::value(dt * c.hydro.acceleration[d]), 2e-16);
	}
	EXPECT_EQ(runtime.snapshots().front().time, units::Time{});
	c.hydro.acceleration[ndim - 1] = units::Acceleration::from_value(-1000);
	c.rayleighTaylor.interfacePressure = units::Pressure::from_value(2000);
	Runtime forced(c);
	auto unforcedConfig = c;
	unforcedConfig.hydro.acceleration.fill({});
	Runtime unforced(unforcedConfig);
	EXPECT_LT(forced.stableTimestep(), unforced.stableTimestep());
	EXPECT_LE(forced.stableTimestep(), 0.2 * units::sqrt(units::Length::from_value(0.125) / units::Acceleration::from_value(1000)));
}

TEST(RayleighTaylor, EvolvingGravityAndWallsAgreeAcrossPartitions) {
	auto fine = parseConfig({"--mesh.cells=4", "--mesh.level=1", "--output.enabled=off", "--runtime.stopTime=0.05"});
	auto single = fine;
	single.mesh.cells = 8;
	single.mesh.level = 0;
	auto const tiled = run(fine), whole = run(single);
	auto actual = flatten(tiled.snapshots, 8), expected = flatten(whole.snapshots, 8);
	ASSERT_EQ(tiled.steps, whole.steps);
	EXPECT_GT(tiled.steps, 1);
	for (std::size_t i = 0; i < actual.size(); ++i)
		test::expectStateNear(actual[i], expected[i], 2e-12);
	EXPECT_NEAR(Real(tiled.final.mass / tiled.initial.mass), 1, 2e-13);
	EXPECT_GT(tiled.final.minimumPressure, units::Pressure{});
}

TEST(RayleighTaylor, PerturbationGrowsUnderGravity) {
	using std::cos;

	// Allow the initial compressible adjustment to pass before measuring growth.
	auto c = parseConfig({"--mesh.cells=16", "--mesh.level=0", "--output.enabled=off", "--runtime.stopTime=5"});
	auto amplitude = [&](std::vector<Snapshot> const& snapshots) {
		Real projection = 0, norm = 0;
		for (auto const& b : snapshots)
			b.layout.forEachInterior([&](auto cell, auto i) {
				auto const point = b.layout.cellCenter(b.lower, b.cellWidth, cell);
				Real mode = 1;
				for (int d = 0; d < ndim - 1; ++d)
					mode *= cos(2 * piR * units::value(point[d]));
				projection += mode * units::value(b.hydro.values()[i].momentum(ndim - 1) / b.hydro.values()[i].density());
				norm += mode * mode;
			});
		return projection / norm;
	};
	auto const initial = amplitude({initialSnapshot(c, {})});
	auto const unstable = run(c);
	EXPECT_GT(amplitude(unstable.snapshots), 1.5 * initial);
	EXPECT_NEAR(Real(unstable.final.mass / unstable.initial.mass), 1, 3e-12);
	EXPECT_GT(unstable.final.minimumDensity, units::Density{});
	EXPECT_GT(unstable.final.minimumPressure, units::Pressure{});
	c.hydro.acceleration.fill({});
	auto const control = run(c);
	EXPECT_GT(amplitude(unstable.snapshots), 1.5 * amplitude(control.snapshots));
}

}	 // namespace
