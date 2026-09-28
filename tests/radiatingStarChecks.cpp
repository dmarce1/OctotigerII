#include "testSupport.hpp"
#include "octotigerII/problems.hpp"
#include "octotigerII/problems/radiatingStarStructure.hpp"
#include "octotigerII/radiation/m1.hpp"
#include "octotigerII/radiation/opacity.hpp"
#include "octotigerII/runtime.hpp"
#include "octotigerII/simulation.hpp"
#include "octotigerII/verification/analytic.hpp"
#include <cmath>
#include <iostream>

using namespace octotigerII;
namespace {
Real value(auto q) { return units::value(q); }
Config configuration(int level = 1) {
	return test::parseConfig({"--mesh.cells=4", "--mesh.level=" + std::to_string(level), "--amr.enabled=off",
		"--output.enabled=off", "--runtime.stopTime=0"});
}
problems::RadiatingStarStructure structure(Config const& c) {
	problems::RadiatingStarStructure::Parameters p;
	p.eos = {c.star.polytropicIndex, c.star.centralDensity,
		c.radiatingStar.centralGasFraction, c.hydro.meanMolecularWeight};
	p.spinFractionOfSphericalBreakup = c.radiatingStar.rotationFraction;
	p.radialCells = c.radiatingStar.radialCells;
	p.angularPoints = c.radiatingStar.angularPoints;
	p.maxMultipole = c.radiatingStar.multipoles;
	p.tolerance = c.radiatingStar.structureTolerance;
	return problems::RadiatingStarStructure(p);
}
}

TEST(RadiatingStar, DefaultsUseFullLightSpeedUniformOpacityAndNoPhotonSource) {
	auto const c = configuration();
	EXPECT_TRUE(c.hydroEnabled()); EXPECT_TRUE(c.radiationEnabled()); EXPECT_TRUE(c.gravityEnabled());
	EXPECT_FALSE(c.radiation.enabled);
	EXPECT_EQ(c.radiation.lightSpeedRatio, 1);
	EXPECT_TRUE(c.radiation.closedBoundary);
	EXPECT_EQ(c.radiation.opacity, .34);
	EXPECT_EQ(c.frame.omega, units::InverseTime{});
	EXPECT_EQ(c.star.polytropicIndex, 3.5);
	EXPECT_EQ(c.star.atmosphereFraction, 1e-12);
	EXPECT_EQ(c.radiatingStar.rotationFraction, .2);
	EXPECT_EQ(c.mesh.lower, -1.2 * c.star.radius);
	EXPECT_EQ(c.mesh.upper, 1.2 * c.star.radius);
	EXPECT_FALSE(problemHasRadiationMaterial(c));
	auto const material = problemRadiationMaterial(c);
	for (Real fraction : {Real(0), Real(.7), Real(1.1)}) {
		mesh::PhysicalCoordinates position{}; position[0] = fraction * c.star.radius;
		for (Real t : {Real(0), Real(100)}) {
			auto const m = material(position, units::Time::from_value(t));
			EXPECT_EQ(value(m.opacity), .34);
			EXPECT_EQ(value(m.photonPower), 0);
		}
	}
}

TEST(RadiatingStar, AnalyticOpacitySelectsStateDependentModelWithoutDefaultConstant) {
	auto const c = test::parseConfig({"--radiation.opacityModel=ionized-gas", "--mesh.cells=4",
		"--amr.enabled=off", "--output.enabled=off", "--runtime.stopTime=0"});
	EXPECT_EQ(c.radiation.opacity, 0);
	EXPECT_EQ(c.radiation.opacityModel, "ionized-gas");
	EXPECT_TRUE(radiation::radiationCouplingEnabled(c));
	EXPECT_THROW(test::parseConfig({"--radiation.opacityModel=ionized-gas", "--radiation.opacity=0.34"}), std::invalid_argument);
	auto const fixed = configuration();
	mesh::PhysicalCoordinates x{};
	x[0] = 0.3 * fixed.star.radius;
	auto const fixedState = problemBoundary(fixed)(x, {});
	auto const analyticState = problemBoundary(c)(x, {});
	auto const gasTemperature = analyticState.hydro.pressure() * c.hydro.meanMolecularWeight * constants::atomicMassUnit
		/ (analyticState.hydro.density() * constants::boltzmann);
	auto const extinction = radiation::opacityLaw(c, {}).evaluate(value(analyticState.hydro.density()),
		value(gasTemperature)).fluxExtinction;
	ASSERT_NE(fixedState.radiation.radiativeFlux(0), units::EnergyFlux{});
	EXPECT_NEAR(value(analyticState.radiation.radiativeFlux(0) / fixedState.radiation.radiativeFlux(0)),
		0.34 / value(extinction), 1e-10);
}

TEST(RadiatingStar, BoostedDiffusionMomentsHaveConsistentThermalAndMechanicalForces) {
	auto const c = configuration(); auto const star = structure(c); auto const evaluate = problemBoundary(c);
	radiation::RadiationSystem radiation(constants::c);
	for (Real fraction : {Real(.1), Real(.3), Real(.6), Real(.98), Real(1.05)}) {
		mesh::PhysicalCoordinates x{};
		x[0] = .8 * fraction * star.equatorialRadius();
		x[1] = .6 * fraction * star.equatorialRadius();
		auto const q = evaluate(x, {});
		EXPECT_TRUE(radiation.admissible(q.radiation));
		auto const later = evaluate(x, units::Time::from_value(10));
		test::expectStateNear(q.radiation, later.radiation, 0);
		auto const T = q.hydro.pressure() * c.hydro.meanMolecularWeight * constants::atomicMassUnit
			/ (q.hydro.density() * constants::boltzmann);
		Real const B = value(constants::radiation * boost::units::pow<4>(T));
		Real const E = value(q.radiation.energy()), light = value(constants::c);
		Real vF = 0;
		for (int d = 0; d < ndim; ++d) { vF += value(q.hydro.velocity(d)) * value(q.radiation.radiativeFlux(d)); }
		EXPECT_NEAR((E - B - vF / (light * light)) / E, 0, 3e-14);
		if (fraction < .9) {
			EXPECT_NEAR(Real(q.hydro.velocity(0) / (-star.angularVelocity() * x[1])), 1, 3e-14);
			EXPECT_NEAR(Real(q.hydro.velocity(1) / (star.angularVelocity() * x[0])), 1, 3e-14);
		}
	}
	// Independent Cartesian pressure differences plus the actual M1 mixed-frame
	// force verify radiation support, rather than checking the SCF identity alone.
	mesh::PhysicalCoordinates x{};
	x[0] = .23 * star.equatorialRadius(); x[1] = .17 * star.equatorialRadius(); x[2] = .19 * star.equatorialRadius();
	auto const q = evaluate(x, {});
	auto const calc = radiation::RadiationSystem::toCalculationState(q.radiation);
	auto const T = q.hydro.pressure() * c.hydro.meanMolecularWeight * constants::atomicMassUnit
		/ (q.hydro.density() * constants::boltzmann);
	Real const B = value(constants::radiation * boost::units::pow<4>(T));
	Real const rho = value(q.hydro.density()), light = value(constants::c), omega = value(star.angularVelocity());
	auto const spacing = 1e-5 * star.equatorialRadius();
	for (int axis = 0; axis < ndim; ++axis) {
		auto const pressure = radiation::M1::physicalFlux(calc, axis, constants::c).flux;
		Real Pv = 0;
		for (int d = 0; d < ndim; ++d) { Pv += value(pressure[d+1] / constants::c) * value(q.hydro.velocity(d)); }
		Real const force = c.radiation.opacity * rho / light
			* (value(q.radiation.radiativeFlux(axis)) - Pv - B * value(q.hydro.velocity(axis)));
		auto low = x, high = x; low[axis] -= spacing; high[axis] += spacing;
		Real const derivative = value(evaluate(high, {}).hydro.pressure() - evaluate(low, {}).hydro.pressure()) / (2 * value(spacing));
		Real const gravity = rho * value(q.gravity.acceleration(axis));
		Real const centrifugal = axis < 2 ? rho * omega * omega * value(x[axis]) : 0;
		Real const scale = std::max({std::abs(derivative), std::abs(gravity), std::abs(force)});
		EXPECT_NEAR((derivative - gravity - force - centrifugal) / scale, 0, 2e-6);
	}
}

TEST(RadiatingStar, RejectsReducedLightSpeedLeakingWallsAndObsoleteOpacityProfiles) {
	EXPECT_THROW(test::parseConfig({"--radiation.lightSpeedRatio=.1"}), std::invalid_argument);
	EXPECT_THROW(test::parseConfig({"--radiation.closedBoundary=off"}), std::invalid_argument);
	EXPECT_THROW(test::parseConfig({"--frame.omega=1e-5"}), std::invalid_argument);
	EXPECT_THROW(test::parseConfig({"--radiatingStar.opticalDepthScale=240"}), std::invalid_argument);
	EXPECT_THROW(test::parseConfig({"--star.radius=1e12"}), std::invalid_argument);
}

TEST(RadiatingStar, ShortEvolutionHasNoPhotonHeaterOrEscapingRadiationEnergy) {
	units::Time interval{};
	for (int resolution = 0; resolution < 2; ++resolution) {
		auto const c = configuration(resolution + 1); Runtime runtime(c); runtime.solveGravity();
		auto const initial = runtime.snapshots(); auto const before = diagnose(initial, c);
		if (resolution == 0) interval = 2.0 * runtime.stableTimestep();
		units::Time elapsed{}; int steps = 0;
		while (elapsed < interval) {
			auto const dt = std::min(runtime.stableTimestep(), interval - elapsed);
			runtime.advanceCoupled(dt); elapsed += dt; ++steps;
		}
		auto const snapshots = runtime.snapshots(); auto const after = diagnose(snapshots, c);
		auto const boundary = runtime.boundaryTransport();
		EXPECT_EQ(boundary.outward.radiationEnergy, units::Energy{});
		EXPECT_EQ(boundary.inward.radiationEnergy, units::Energy{});
		EXPECT_EQ(runtime.radiationSourceEnergy(), units::Energy{});
		auto const outward = boundary.outward.gasEnergy - boundary.inward.gasEnergy
			+ boundary.outward.potentialEnergy - boundary.inward.potentialEnergy;
		Real const energyError = Real((after.physicalTotalEnergy + outward - before.physicalTotalEnergy) / before.physicalTotalEnergyNorm);
		Real const massError = Real((after.mass + boundary.outward.mass - boundary.inward.mass - before.mass) / before.mass);
		EXPECT_NEAR(energyError, 0, 3e-12);
		EXPECT_NEAR(massError, 0, 3e-12);
		EXPECT_GT(after.minimumDensity, units::Density{});
		EXPECT_GT(after.minimumPressure, units::Pressure{});
		EXPECT_GE(after.minimumRadiationEnergy, units::EnergyDensity{});
		EXPECT_LE(after.maximumReducedFlux, 1 + radiation::M1::roundoff);
		long double densityChange = 0, densityNorm = 0, meridionalKinetic = 0;
		long double negativeGasEnergy = 0, negativeGasEnergyMass = 0;
		ASSERT_EQ(initial.size(), snapshots.size());
		for (std::size_t b = 0; b < initial.size(); ++b) {
			auto const& block = snapshots[b];
			Real const volume = std::pow(value(block.cellWidth), ndim);
			block.layout.forEachInterior([&](auto const& cell, std::size_t i) {
				auto const& q = block.hydro.values()[i];
				Real const rho = value(q.density());
				if (q.totalEnergy() < units::EnergyDensity{}) {
					negativeGasEnergy -= value(q.totalEnergy()) * volume;
					negativeGasEnergyMass += rho * volume;
				}
				densityChange += std::abs(rho - value(initial[b].hydro.values()[i].density())) * volume;
				densityNorm += value(initial[b].hydro.values()[i].density()) * volume;
				auto const x = block.layout.cellCenter(block.lower, block.cellWidth, cell);
				Real const dx = value(x[0] - c.star.center[0]), dy = value(x[1] - c.star.center[1]);
				Real const cylindricalRadius = std::hypot(dx, dy);
				Real const radialVelocity = cylindricalRadius > 0
					? (dx * value(q.momentum(0)) + dy * value(q.momentum(1))) / (cylindricalRadius * rho) : 0;
				Real const verticalVelocity = value(q.momentum(2)) / rho;
				meridionalKinetic += rho * (radialVelocity * radialVelocity + verticalVelocity * verticalVelocity) * volume;
			});
		}
		Real const meridionalRmsMach = std::sqrt(meridionalKinetic
			/ (c.hydro.gamma * (c.hydro.gamma - 1) * value(after.thermalEnergy)));
		Real const radiationDrift = Real((after.radiationEnergy - before.radiationEnergy) / before.radiationEnergy);
		std::cout << "opaque rotating star cells=" << (8 << resolution) << " steps=" << steps << " time=" << value(elapsed)
			<< " densityL1=" << densityChange / densityNorm << " meridionalRmsMach=" << meridionalRmsMach
			<< " radiationDrift=" << radiationDrift << " energyBudget=" << energyError << " massBudget=" << massError
			<< " negativeGasEnergyNorm=" << negativeGasEnergy / value(before.physicalTotalEnergyNorm)
			<< " negativeGasEnergyMassFraction=" << negativeGasEnergyMass / value(before.mass)
			<< " maximumFlux=" << after.maximumReducedFlux << '\n';
	}
}
