#include <gtest/gtest.h>
#include <hpx/include/serialization.hpp>
#include <iostream>
#include <stdexcept>
#include "octotigerII/subgrid/fluxPacket.hpp"
#include "octotigerII/subgrid/subgrid.hpp"
#include "octotigerII/subgrid/topology.hpp"
using namespace octotigerII;


namespace {

template <typename State>
void same(State const& a, State const& b) {
	a.forEach([&](auto f, auto const& q) { EXPECT_TRUE(q == b.template get<f>()); });
}

template <typename State>
void samePatch(mesh::PatchData<State> const& a, mesh::PatchData<State> const& b) {
	EXPECT_TRUE(a.cellWidth() == b.cellWidth() && a.lower() == b.lower());
	EXPECT_TRUE(a.timeState().time == b.timeState().time && a.timeState().stepSize == b.timeState().stepSize);
	ASSERT_EQ(a.values().size(), b.values().size());
	for (std::size_t i = 0; i < a.values().size(); ++i)
		same(a.values()[i], b.values()[i]);
}

void check() {
	auto c = parseConfig({"--mesh.level=0", "--output.enabled=off"});
	c.verification.analytic = "on";
	c.randomSeed = 987654321;
	c.verification.directSamples = 23;
	c.verification.directMaxPairs = 12345;
	c.verification.gravityReference = "continuum";
	c.verification.relativeL1Tolerance = 0.05;
	c.verification.absoluteTolerance = 2e-12;
	c.rayleighTaylor.perturbation = units::Velocity::from_value(0.031);
	auto before = initialSnapshot(c, {0, {}});
	// Round-trip all option types even in builds without hydro.
	c.hydro.acceleration[ndim - 1] = units::Acceleration::from_value(-0.17);
	before.time = units::Time::from_value(0.125);
	before.hydro = hydro::Fields(before.layout, before.cellWidth, before.lower);
	before.gravity = gravity::Fields(before.layout, before.cellWidth, before.lower);
	before.hydro.timeState().completeStep(before.time);
	before.gravity.timeState().completeStep(before.time);
	// Include every stored field type even though no coupled problem exists yet.
	before.radiation = radiation::Fields(before.layout, before.cellWidth, before.lower);
	before.radiation.timeState().completeStep(before.time);
	before.density = mesh::PatchData<units::Density>(before.layout, before.cellWidth, before.lower);
	before.layout.forEachInterior([&](mesh::Coordinates const&, std::size_t i) {
		before.radiation.values()[i].energy() = units::EnergyDensity::from_value(3);
		before.radiation.values()[i].radiativeFlux(0) = units::EnergyFlux::from_value(1);
		before.gravity.values()[i].potential() = units::VelocitySquared::from_value(-12);
		before.gravity.values()[i].acceleration(0) = units::Acceleration::from_value(2);
		before.hydro.values()[i].density() = units::Density::from_value(1);
		before.hydro.values()[i].totalEnergy() = units::EnergyDensity::from_value(2.5);
		before.density.values()[i] = units::Density::from_value(17);
	});
	FieldFluxPacket<hydro::ConservedFlux> gasFlux;
	gasFlux.interval = {before.time, before.time + units::Time::from_value(0.1)};
	hydro::ConservedFlux gasFace{};
	gasFace.mass() = units::MassFlux::from_value(4);
	gasFace.momentum(0) = units::MomentumFlux::from_value(7);
	gasFace.energy() = units::EnergyFlux::from_value(9);
	gasFlux.fluxes = {{gasFace}};
	FieldFluxPacket<radiation::RadiationSystem::Flux> radFlux;
	radiation::RadiationSystem::Flux radFace{};
	radFace.energy() = units::EnergyFlux::from_value(11);
	radFace.radiativeFlux(0) = units::EnergyFluxTransport::from_value(13);
	radFlux.fluxes = {{radFace}};
	std::vector<char> buffer;
	{
		hpx::serialization::output_archive archive(buffer);
		archive & c & before & gasFlux & radFlux;
	}
	Config restored;
	Snapshot after;
	decltype(gasFlux) restoredGas;
	decltype(radFlux) restoredRad;
	{
		hpx::serialization::input_archive archive(buffer);
		archive & restored & after & restoredGas & restoredRad;
	}
	EXPECT_TRUE(restored.mesh.lower == c.mesh.lower && restored.mesh.upper == c.mesh.upper && restored.runtime.stopTime == c.runtime.stopTime);
	EXPECT_TRUE(restored.randomSeed == c.randomSeed && restored.verification.gravityReference == c.verification.gravityReference &&
		restored.verification.directSamples == c.verification.directSamples && restored.verification.directMaxPairs == c.verification.directMaxPairs);
	EXPECT_TRUE(restored.verification.analytic == c.verification.analytic && restored.verification.relativeL1Tolerance == c.verification.relativeL1Tolerance &&
		restored.verification.absoluteTolerance == c.verification.absoluteTolerance);
	EXPECT_TRUE(before.time == after.time && before.cellWidth == after.cellWidth && before.lower == after.lower);
	EXPECT_EQ(restored.hydro.acceleration, c.hydro.acceleration);
	EXPECT_EQ(restored.rayleighTaylor.densityLower, c.rayleighTaylor.densityLower);
	EXPECT_EQ(restored.rayleighTaylor.densityUpper, c.rayleighTaylor.densityUpper);
	EXPECT_EQ(restored.rayleighTaylor.interfacePressure, c.rayleighTaylor.interfacePressure);
	EXPECT_EQ(restored.rayleighTaylor.perturbation, c.rayleighTaylor.perturbation);
	samePatch(before.hydro, after.hydro);
	samePatch(before.radiation, after.radiation);
	samePatch(before.gravity, after.gravity);
	EXPECT_TRUE(before.density.values() == after.density.values());
	EXPECT_TRUE(restoredGas.interval.begin == gasFlux.interval.begin && restoredGas.interval.end == gasFlux.interval.end);
	same(gasFlux.fluxes[0][0], restoredGas.fluxes[0][0]);
	same(radFlux.fluxes[0][0], restoredRad.fluxes[0][0]);
}

}	 // namespace


TEST(Serialization, TypedConfigSnapshotsAndFluxPacketsRoundTrip) {
	check();
}


TEST(Serialization, PerFaceBoundariesAndHaloPlansRoundTrip) {
	auto c = parseConfig({"--mesh.periodic=off", "--mesh.cells=4", "--mesh.level=1"});
	CartesianTopology topology(c, 2);
	// Serialize all boundary kinds even in builds whose gravity policy disallows them.
	c.mesh.boundary.lower[0] = physics::BoundaryCondition::Analytic;
	c.mesh.boundary.upper[0] = physics::BoundaryCondition::Reflecting;
	if (ndim > 1) c.mesh.boundary.lower[1] = c.mesh.boundary.upper[1] = physics::BoundaryCondition::Periodic;
	auto plan = makeHaloPlan(c, topology.blocks(), 0);
	std::vector<char> buffer;
	{
		hpx::serialization::output_archive archive(buffer);
		archive & c & plan;
	}
	Config restored;
	HaloPlan restoredPlan;
	{
		hpx::serialization::input_archive archive(buffer);
		archive & restored & restoredPlan;
	}
	EXPECT_EQ(restored.mesh.boundary.lower, c.mesh.boundary.lower);
	EXPECT_EQ(restored.mesh.boundary.upper, c.mesh.boundary.upper);
	EXPECT_EQ(restoredPlan.ghostCount, plan.ghostCount);
	EXPECT_EQ(restoredPlan.ghostIndices, plan.ghostIndices);
	EXPECT_EQ(restoredPlan.reflectionMasks, plan.reflectionMasks);
	EXPECT_EQ(restoredPlan.outflowLowerMasks, plan.outflowLowerMasks);
	EXPECT_EQ(restoredPlan.outflowUpperMasks, plan.outflowUpperMasks);
	ASSERT_EQ(restoredPlan.analyticGhosts.size(), plan.analyticGhosts.size());
	for (std::size_t i = 0; i < plan.analyticGhosts.size(); ++i) {
		EXPECT_EQ(restoredPlan.analyticGhosts[i].destination, plan.analyticGhosts[i].destination);
		EXPECT_EQ(restoredPlan.analyticGhosts[i].position, plan.analyticGhosts[i].position);
	}
	ASSERT_EQ(restoredPlan.reads.size(), plan.reads.size());
	for (std::size_t i = 0; i < plan.reads.size(); ++i) {
		auto const& actual = restoredPlan.reads[i];
		auto const& expected = plan.reads[i];
		EXPECT_EQ(actual.range.partition, expected.range.partition);
		EXPECT_EQ(actual.range.offset, expected.range.offset);
		EXPECT_EQ(actual.range.count, expected.range.count);
		ASSERT_EQ(actual.copies.size(), expected.copies.size());
		for (std::size_t j = 0; j < expected.copies.size(); ++j) {
			EXPECT_EQ(actual.copies[j].source, expected.copies[j].source);
			EXPECT_EQ(actual.copies[j].destination, expected.copies[j].destination);
		}
	}
}

TEST(Serialization, InflowAndDirectionalOutflowMasksRoundTrip) {
	auto c = parseConfig({"--mesh.periodic=off", "--mesh.cells=4", "--mesh.level=0"});
	c.mesh.boundary.upper.fill(physics::BoundaryCondition::Inflow);
	CartesianTopology topology(c, 1);
	auto const plan = makeHaloPlan(c, topology.blocks(), 0);
	std::vector<char> buffer;
	{
		hpx::serialization::output_archive archive(buffer);
		archive & c & plan;
	}
	Config restored;
	HaloPlan restoredPlan;
	{
		hpx::serialization::input_archive archive(buffer);
		archive & restored & restoredPlan;
	}
	EXPECT_EQ(restored.mesh.boundary.lower, c.mesh.boundary.lower);
	EXPECT_EQ(restored.mesh.boundary.upper, c.mesh.boundary.upper);
	EXPECT_EQ(restoredPlan.outflowLowerMasks, plan.outflowLowerMasks);
	EXPECT_EQ(restoredPlan.outflowUpperMasks, plan.outflowUpperMasks);
	EXPECT_NE(std::count_if(plan.outflowLowerMasks.begin(), plan.outflowLowerMasks.end(), [](auto mask) { return mask != 0; }), 0);
}
