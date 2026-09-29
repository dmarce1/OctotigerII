#include "testSupport.hpp"
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
	auto c = test::parseConfig({"--mesh.level=0", "--output.enabled=off"});
	c.verification.analytic = "on";
	c.randomSeed = 987654321;
	c.timestep.refinement = false;
	c.gravity.timeIntegration = "conventional";
	c.gravity.energyTreatment = "naive";
	c.gravity.conserveRegridEnergy = false;
	c.hydro.dualEnergy = {true, -0.5, 0.003, 0.2};
	c.hydro.meanMolecularWeight = 0.6;
	c.radiation.lightSpeedRatio = 0.25;
	c.radiation.diagnosticLength = 7.5e8;
	c.radiation.initialEnergyRatio = 0.3;
	if (c.hydroEnabled() && build::massFractions) {
		c.massFractions.enabled = true;
		c.massFractions.species = composition::parseSpecies("gas:1:He=70%,O=30%;dye:2:A=0,Z=0");
	}
	c.amr.refineDensity = units::Density::from_value(0.03);
	c.star.radius = units::Length::from_value(8e8);
	c.star.center[0] = units::Length::from_value(2e8);
	c.scf.primaryMass = units::Mass::from_value(3e33);
	c.scf.separation = units::Length::from_value(7e10);
	c.scf.referenceWidth = 2.2;
	c.scf.massRatio = .6;
	c.scf.coreIndex = {3, 2.7}; c.scf.envelopeIndex = {1.5, 1.2};
	c.scf.interfaceFraction = {.12,.4}; c.scf.densityJump = {2,1.1}; c.scf.fill = {.95,1};
	c.scf.cells = 128; c.scf.history = 3; c.scf.maxIterations = 777;
	c.scf.tolerance = 3e-6; c.scf.virialTolerance = .02; c.scf.relaxation = .3; c.scf.atmosphereFraction = 2e-10;
	c.scf.evolveOrbits = .25; c.scf.framesPerOrbit = 100;
	c.radiatingStar.centralGasFraction = 0.7;
	c.radiatingStar.rotationFraction = 0.15;
	c.radiatingStar.opticalDepthScale = 321;
	c.radiatingStar.opacityCutoffFraction = 0.021;
	c.radiatingStar.radialCells = 192;
	c.radiatingStar.angularPoints = 24;
	c.radiatingStar.multipoles = 10;
	c.radiatingStar.structureTolerance = 3e-9;
	c.verification.directSamples = 23;
	c.verification.directMaxPairs = 12345;
	c.verification.gravityReference = "continuum";
	c.verification.relativeL1Tolerance = 0.05;
	c.verification.absoluteTolerance = 2e-12;
	c.rayleighTaylor.perturbation = units::Velocity::from_value(0.031);
	auto before = initialSnapshot(c, {0, {}});
	// Round-trip all option types even in builds without hydro.
	c.radiation.enabled = true;
	c.radiation.closedBoundary = true;
	c.radiation.opacity = 0.4;
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
		before.hydro.values()[i].auxiliary() = units::Density::from_value(0.7);
		before.density.values()[i] = units::Density::from_value(17);
	});
	FieldFluxPacket<hydro::ConservedFlux> gasFlux;
	gasFlux.interval = {before.time, before.time + units::Time::from_value(0.1)};
	hydro::ConservedFlux gasFace{};
	gasFace.mass() = units::MassFlux::from_value(4);
	gasFace.momentum(0) = units::MomentumFlux::from_value(7);
	gasFace.energy() = units::EnergyFlux::from_value(9);
	gasFace.auxiliary() = units::MassFlux::from_value(3);
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
	EXPECT_EQ(restored.timestep.refinement, c.timestep.refinement);
	EXPECT_EQ(restored.scf.primaryMass, c.scf.primaryMass);
	EXPECT_EQ(restored.scf.separation, c.scf.separation);
	EXPECT_EQ(restored.scf.referenceWidth, c.scf.referenceWidth);
	EXPECT_EQ(restored.scf.massRatio, c.scf.massRatio);
	EXPECT_EQ(restored.scf.coreIndex, c.scf.coreIndex);
	EXPECT_EQ(restored.scf.envelopeIndex, c.scf.envelopeIndex);
	EXPECT_EQ(restored.scf.interfaceFraction, c.scf.interfaceFraction);
	EXPECT_EQ(restored.scf.densityJump, c.scf.densityJump);
	EXPECT_EQ(restored.scf.fill, c.scf.fill);
	EXPECT_EQ(restored.scf.cells, c.scf.cells);
	EXPECT_EQ(restored.scf.history, c.scf.history);
	EXPECT_EQ(restored.scf.maxIterations, c.scf.maxIterations);
	EXPECT_EQ(restored.scf.tolerance, c.scf.tolerance);
	EXPECT_EQ(restored.scf.virialTolerance, c.scf.virialTolerance);
	EXPECT_EQ(restored.scf.relaxation, c.scf.relaxation);
	EXPECT_EQ(restored.scf.atmosphereFraction, c.scf.atmosphereFraction);
	EXPECT_EQ(restored.scf.evolveOrbits, c.scf.evolveOrbits);
	EXPECT_EQ(restored.scf.framesPerOrbit, c.scf.framesPerOrbit);
	EXPECT_TRUE(restored.mesh.lower == c.mesh.lower && restored.mesh.upper == c.mesh.upper && restored.runtime.stopTime == c.runtime.stopTime);
	EXPECT_TRUE(restored.randomSeed == c.randomSeed && restored.verification.gravityReference == c.verification.gravityReference &&
		restored.verification.directSamples == c.verification.directSamples && restored.verification.directMaxPairs == c.verification.directMaxPairs);
	EXPECT_TRUE(restored.verification.analytic == c.verification.analytic && restored.verification.relativeL1Tolerance == c.verification.relativeL1Tolerance &&
		restored.verification.absoluteTolerance == c.verification.absoluteTolerance);
	EXPECT_TRUE(before.time == after.time && before.cellWidth == after.cellWidth && before.lower == after.lower);
	EXPECT_EQ(restored.hydro.acceleration, c.hydro.acceleration);
	EXPECT_EQ(restored.gravity.timeIntegration, c.gravity.timeIntegration);
	EXPECT_EQ(restored.gravity.energyTreatment, c.gravity.energyTreatment);
	EXPECT_EQ(restored.gravity.conserveRegridEnergy, c.gravity.conserveRegridEnergy);
	EXPECT_EQ(restored.hydro.dualEnergy.enabled, c.hydro.dualEnergy.enabled);
	EXPECT_EQ(restored.hydro.dualEnergy.exponent, c.hydro.dualEnergy.exponent);
	EXPECT_EQ(restored.hydro.dualEnergy.pressureThreshold, c.hydro.dualEnergy.pressureThreshold);
	EXPECT_EQ(restored.hydro.dualEnergy.syncThreshold, c.hydro.dualEnergy.syncThreshold);
	EXPECT_EQ(restored.hydro.meanMolecularWeight, c.hydro.meanMolecularWeight);
	EXPECT_EQ(restored.radiation.enabled, c.radiation.enabled);
	EXPECT_EQ(restored.radiation.closedBoundary, c.radiation.closedBoundary);
	EXPECT_EQ(restored.radiation.opacity, c.radiation.opacity);
	EXPECT_EQ(restored.radiation.lightSpeedRatio, c.radiation.lightSpeedRatio);
	EXPECT_EQ(restored.radiation.diagnosticLength, c.radiation.diagnosticLength);
	EXPECT_EQ(restored.radiation.initialEnergyRatio, c.radiation.initialEnergyRatio);
	EXPECT_EQ(restored.amr.refineDensity, c.amr.refineDensity);
	EXPECT_EQ(restored.star.radius, c.star.radius);
	EXPECT_EQ(restored.star.center, c.star.center);
	EXPECT_EQ(restored.star.centralDensity, c.star.centralDensity);
	EXPECT_EQ(restored.star.polytropicIndex, c.star.polytropicIndex);
	EXPECT_EQ(restored.star.atmosphereFraction, c.star.atmosphereFraction);
	EXPECT_EQ(restored.radiatingStar.centralGasFraction, c.radiatingStar.centralGasFraction);
	EXPECT_EQ(restored.radiatingStar.rotationFraction, c.radiatingStar.rotationFraction);
	EXPECT_EQ(restored.radiatingStar.opticalDepthScale, c.radiatingStar.opticalDepthScale);
	EXPECT_EQ(restored.radiatingStar.opacityCutoffFraction, c.radiatingStar.opacityCutoffFraction);
	EXPECT_EQ(restored.radiatingStar.radialCells, c.radiatingStar.radialCells);
	EXPECT_EQ(restored.radiatingStar.angularPoints, c.radiatingStar.angularPoints);
	EXPECT_EQ(restored.radiatingStar.multipoles, c.radiatingStar.multipoles);
	EXPECT_EQ(restored.radiatingStar.structureTolerance, c.radiatingStar.structureTolerance);
	EXPECT_EQ(restored.rayleighTaylor.densityLower, c.rayleighTaylor.densityLower);
	EXPECT_EQ(restored.rayleighTaylor.densityUpper, c.rayleighTaylor.densityUpper);
	EXPECT_EQ(restored.rayleighTaylor.interfacePressure, c.rayleighTaylor.interfacePressure);
	EXPECT_EQ(restored.rayleighTaylor.perturbation, c.rayleighTaylor.perturbation);
	EXPECT_EQ(restored.massFractions.enabled, c.massFractions.enabled);
	ASSERT_EQ(restored.massFractions.species.size(), c.massFractions.species.size());
	ASSERT_EQ(before.species.size(), after.species.size());
	for (std::size_t s = 0; s < c.massFractions.species.size(); ++s) {
		auto const& a = c.massFractions.species[s]; auto const& b = restored.massFractions.species[s];
		EXPECT_EQ(a.name, b.name); EXPECT_EQ(a.atomicMass, b.atomicMass); EXPECT_EQ(a.atomicNumber, b.atomicNumber);
		EXPECT_EQ(a.initialFraction, b.initialFraction);
		ASSERT_EQ(a.mixture.size(), b.mixture.size());
		for (std::size_t i = 0; i < a.mixture.size(); ++i) {
			EXPECT_EQ(a.mixture[i].atomicNumber, b.mixture[i].atomicNumber);
			EXPECT_EQ(a.mixture[i].massFraction, b.mixture[i].massFraction);
		}
		EXPECT_EQ(before.species[s].values(), after.species[s].values());
	}
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
	auto c = test::parseConfig({"--mesh.periodic=off", "--mesh.cells=4", "--mesh.level=1"});
	CartesianTopology topology(c, 2);
	// Serialize all boundary kinds even in builds whose gravity policy disallows them.
	c.mesh.boundary.lower[0] = finiteVolume::BoundaryCondition::Analytic;
	c.mesh.boundary.upper[0] = finiteVolume::BoundaryCondition::Reflecting;
	if (ndim > 1) c.mesh.boundary.lower[1] = c.mesh.boundary.upper[1] = finiteVolume::BoundaryCondition::Periodic;
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
	auto c = test::parseConfig({"--mesh.periodic=off", "--mesh.cells=4", "--mesh.level=0"});
	c.mesh.boundary.upper.fill(finiteVolume::BoundaryCondition::Inflow);
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
