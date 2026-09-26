#include "testSupport.hpp"
#include "octotigerII/gravity/rotationWork.hpp"
#include "octotigerII/runtime.hpp"
#include "octotigerII/simulation.hpp"
#include <cmath>
#include <iostream>
#include <set>
#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/runtime_distributed/find_all_localities.hpp>
#endif

using namespace octotigerII;

namespace {
std::vector<storage::Locality> owners() {
#ifdef OCTOTIGERII_WITH_HPX
	return hpx::find_all_localities();
#else
	return {0};
#endif
}

std::vector<mesh::BlockLocation> leaves(bool adaptive) {
	std::vector<mesh::BlockLocation> result;
	for (int slot = 0; slot < 8; ++slot) {
		auto child = mesh::BlockLocation{}.child(slot);
		if (adaptive && slot == 0) for (int sub = 0; sub < 8; ++sub) result.push_back(child.child(sub));
		else result.push_back(child);
	}
	return result;
}

void checkPairWork(bool adaptive, Real opening) {
	auto c = test::parseConfig({"--mesh.cells=4", "--mesh.level=1", "--output.enabled=off"});
	c.amr.enabled = adaptive;
	c.gravity.multipoleOrder = 5; c.gravity.openingAngle = opening;
	auto locality = owners();
	CartesianTopology topology(c, locality.size(), leaves(adaptive));
	FieldRepository fields(c, topology.storageLayout(), locality);
	storage::PartitionSet store(locality);
	gravity::RotationWorkWorkspace scratch(topology.storageLayout(), store);
	storage::Field<units::Density> density(topology.storageLayout(), store, 1, "rotationTest.density");
	storage::ColumnFields<gravity::State> acceleration(topology.storageLayout(), store, "rotationTest.acceleration", false, 2);
	storage::Field<units::VelocitySquared> potential(topology.storageLayout(), store, 2, "rotationTest.work");
	gravity::FieldSolver solver(c, topology.blocks(), fields.directory(), locality);
	std::size_t count = 0;
	for (auto const& b : topology.blocks()) {
		auto rho = density.handle().output(b.interior, 0);
		for (std::size_t i = 0; i < b.interior.count; ++i)
			rho.data()[i] = units::Density::from_value(2 + Real(((count + i) * 17 + 31) % 127) / 53);
		density.handle().commit(b.interior, 0, rho); count += b.interior.count;
	}
	for (unsigned subset = 0; subset < 2; ++subset) {
		gravity::FieldSolveRequest request;
		request.density = density.handle(); request.output = acceleration.handle(); request.outputBank = subset;
		if (subset) for (auto const& b : topology.blocks()) {
			request.sources.push_back({0, 0, b.id % 2 ? Real(1) : Real(0), 0});
			request.targets.push_back(b.id % 2);
			if (!(b.id % 2)) {
				auto g = acceleration.handle().output(b.interior, subset);
				for (std::size_t i = 0; i < b.interior.count; ++i) g.put(i, {});
				acceleration.handle().commit(b.interior, subset, g);
			}
		}
		solver.solve(request);
		gravity::rotationWorkPotential(solver, request, topology.blocks(), density.handle(), acceleration.handle(), subset,
			scratch, potential.handle(), subset, c.mesh.upper - c.mesh.lower);
	}
	std::array<long double, 3> sum{}, norm{}, error{}, raw{};
	for (auto const& b : topology.blocks()) {
		auto rho = density.handle().read(b.interior, 0).get();
		auto full = potential.handle().read(b.interior, 0).get(), fast = potential.handle().read(b.interior, 1).get();
		auto g = acceleration.handle().read(b.interior, 0).get(), gf = acceleration.handle().read(b.interior, 1).get();
		b.layout.forEachInterior([&](auto const& cell, std::size_t i) {
			auto const x = b.lower[0] + (Real(cell[0]) + 0.5) * b.cellWidth;
			auto const y = b.lower[1] + (Real(cell[1]) + 0.5) * b.cellWidth;
			std::array values{full.data()[i], fast.data()[i], full.data()[i] - fast.data()[i]};
			std::array<gravity::State, 3> forces{g.at(i), gf.at(i), g.at(i) - gf.at(i)};
			auto const mass = rho.data()[i] * b.layout.cellMeasure(b.cellWidth);
			for (int kind = 0; kind < 3; ++kind) {
				long double const actual = units::value(mass * values[kind]);
				long double const physical = units::value(mass * (x * forces[kind].acceleration(1) - y * forces[kind].acceleration(0)));
				sum[kind] += actual; raw[kind] += physical; norm[kind] += std::abs(physical); error[kind] += std::abs(actual - physical);
			}
		});
	}
	for (int kind = 0; kind < 3; ++kind) {
		ASSERT_GT(norm[kind], 0);
		std::cout << "rotation work adaptive=" << adaptive << " opening=" << opening << " component=" << kind
			<< " raw residual=" << raw[kind] / norm[kind] << " balanced residual=" << sum[kind] / norm[kind]
			<< " local correction=" << error[kind] / norm[kind] << '\n';
		EXPECT_LT(std::abs(sum[kind]) / norm[kind], 4e-14L);
		if (opening < 0.1) EXPECT_LT(error[kind] / norm[kind], 4e-14L);
	}
}

Config configuration(std::string const& method, Real omega) {
	auto c = parseConfig({"--problem.name=polytrope", "--mesh.cells=4", "--mesh.level=1",
		"--mesh.lower=-2e8", "--mesh.upper=2e8", "--amr.enabled=on", "--amr.minLevel=1", "--amr.maxLevel=2",
		"--amr.refineDensity=0", "--amr.maxCellMass=0", "--amr.shadowTolerance=0", "--amr.bufferCells=0",
		"--output.enabled=off", "--verification.analytic=off", "--runtime.stopTime=0"});
	c.frame.omega = units::InverseTime::from_value(omega);
	c.star.center[0] = units::Length::from_value(3.5e7);
	c.star.center[1] = units::Length::from_value(-2.0e7);
	c.gravity.timeIntegration = method == "global" ? "hierarchical" : method;
	c.timestep.refinement = method != "global";
	return c;
}

refinement::Criterion criterion(bool& refine) {
	return [&refine](refinement::CellView const& cell) {
		bool lower = true;
		for (auto x : cell.center) lower = lower && x < units::Length{};
		return refine && lower && cell.level < 2 ? Real(2) : Real(0);
	};
}

void expectBudget(Runtime const& runtime, Config const& c, Diagnostics const& before) {
	auto const after = diagnose(runtime.snapshots(), c);
	auto const b = runtime.boundaryTransport();
	EXPECT_NEAR(Real((after.mass + b.outward.mass - b.inward.mass - before.mass) / before.mass), 0, 8e-13);
	auto const energy = after.gasGravityEnergy + b.outward.gasEnergy - b.inward.gasEnergy
		+ b.outward.potentialEnergy - b.inward.potentialEnergy;
	EXPECT_NEAR(Real((energy - before.gasGravityEnergy) / before.gasGravityNorm), 0, 1e-12);
}

std::array<Real, 3> difference(std::vector<Snapshot> const& actual, std::vector<Snapshot> const& reference) {
	std::array<long double, 3> delta{}, norm{};
	if (actual.size() != reference.size()) throw std::logic_error("Rotation convergence changed the leaf mesh");
	for (std::size_t b = 0; b < actual.size(); ++b) {
		auto const volume = units::value(actual[b].layout.cellMeasure(actual[b].cellWidth));
		for (std::size_t i = 0; i < actual[b].hydro.values().size(); ++i) {
			auto const& a = actual[b].hydro.values()[i]; auto const& r = reference[b].hydro.values()[i];
			delta[0] += volume * std::abs(units::value(a.density() - r.density())); norm[0] += volume * std::abs(units::value(r.density()));
			delta[2] += volume * std::abs(units::value(a.totalEnergy() - r.totalEnergy())); norm[2] += volume * std::abs(units::value(r.totalEnergy()));
			for (int d = 0; d < ndim; ++d) {
				delta[1] += volume * std::abs(units::value(a.momentum(d) - r.momentum(d)));
				norm[1] += volume * std::abs(units::value(r.momentum(d)));
			}
		}
	}
	return {Real(delta[0] / norm[0]), Real(delta[1] / norm[1]), Real(delta[2] / norm[2])};
}
}

TEST(RotatingGravityWork, DirectAndMultipoleFullNestedAndRungWorkBalance) {
	for (bool adaptive : {false, true}) for (Real opening : {Real(0.01), Real(0.5)}) checkPairWork(adaptive, opening);
}

TEST(RotatingGravityIntegration, GlobalAndRefinedModesConserveWithEitherRotationSign) {
	for (auto const* mode : {"global", "hierarchical", "conventional"}) for (Real omega : {Real(-0.1), Real(0.1)}) {
		std::cout << "rotating gravity mode=" << mode << " omega=" << omega << '\n';
		auto c = configuration(mode, omega); bool refine = true;
		Runtime runtime(c, {criterion(refine)}); runtime.solveGravity();
		auto const before = diagnose(runtime.snapshots(), c);
		auto const dt = 0.2 * runtime.stableTimestep();
		for (int n = 0; n < 3; ++n) { runtime.advanceGravity(dt); expectBudget(runtime, c, before); }
		if (c.timestep.refinement) {
			auto steps = runtime.statistics().levelSteps;
			ASSERT_GT(steps.size(), 2u); EXPECT_EQ(steps[1], 3u); EXPECT_GE(steps[2], 6u);
		}
		if (c.gravity.timeIntegration == "hierarchical") {
			auto after = diagnose(runtime.snapshots(), c); auto b = runtime.boundaryTransport();
			for (int d = 0; d < ndim; ++d) {
				auto norm = after.norm.momentum[d] + b.outward.momentum[d] + b.inward.momentum[d];
				if (norm > units::Momentum{}) EXPECT_NEAR(Real((after.momentum[d] + b.outward.momentum[d] - b.inward.momentum[d] - before.momentum[d]) / norm), 0, 1e-12);
			}
		}
	}
}

TEST(RotatingGravityIntegration, RefineAndCoarsenPreserveEnergy) {
	for (auto const* mode : {"hierarchical", "conventional"}) {
		auto c = configuration(mode, 0.1); bool refine = true;
		Runtime runtime(c, {criterion(refine)}); runtime.solveGravity();
		auto const before = diagnose(runtime.snapshots(), c);
		runtime.advanceGravity(0.1 * runtime.stableTimestep());
		refine = false; ASSERT_TRUE(runtime.regrid({}, true)); runtime.solveGravity(); expectBudget(runtime, c, before);
		refine = true; ASSERT_TRUE(runtime.regrid({}, true)); runtime.solveGravity(); expectBudget(runtime, c, before);
		runtime.advanceGravity(0.1 * runtime.stableTimestep()); expectBudget(runtime, c, before);
	}
}

TEST(RotatingGravityIntegration, NaiveWorkRemainsIndependent) {
	for (auto const* mode : {"global", "hierarchical", "conventional"}) {
		auto c = configuration(mode, -0.1); c.gravity.energyTreatment = "naive"; bool refine = true;
		Runtime runtime(c, {criterion(refine)}); runtime.solveGravity();
		runtime.advanceGravity(0.1 * runtime.stableTimestep());
		for (auto const& b : runtime.snapshots()) for (auto const& state : b.hydro.values())
			EXPECT_TRUE(hydro::HydroSystem(c.hydro).admissible(state));
	}
}

TEST(RotatingGravityIntegration, TemporalConvergence) {
	for (auto const* mode : {"hierarchical", "conventional"}) {
		auto c = configuration(mode, 0.2); bool refine = true;
		Runtime initial(c, {criterion(refine)}); initial.solveGravity();
		auto const duration = 0.6 * initial.stableTimestep();
		auto integrate = [&](int steps) {
			Runtime runtime(c, {criterion(refine)}); runtime.solveGravity();
			for (int n = 0; n < steps; ++n) runtime.advanceGravity(duration / Real(steps));
			return runtime.snapshots();
		};
		auto const reference = integrate(32);
		std::array<std::array<Real, 3>, 3> errors;
		int n = 0; for (int steps : {2, 4, 8}) errors[n++] = difference(integrate(steps), reference);
		for (int component = 0; component < 3; ++component) {
			std::cout << "rotating " << mode << " temporal component=" << component << " errors=" << errors[0][component]
				<< ',' << errors[1][component] << ',' << errors[2][component] << " orders=" << std::log2(errors[0][component] / errors[1][component])
				<< ',' << std::log2(errors[1][component] / errors[2][component]) << '\n';
			EXPECT_GT(errors[0][component], 1e-13);
			EXPECT_LT(errors[1][component], 0.35 * errors[0][component]);
			EXPECT_LT(errors[2][component], 0.35 * errors[1][component]);
		}
	}
}
