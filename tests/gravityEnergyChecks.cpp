#include "testSupport.hpp"
#include "octotigerII/runtime.hpp"
#include "octotigerII/amr/hierarchy.hpp"
#include "octotigerII/gravity/fieldSolver.hpp"
#include "octotigerII/subgrid/topology.hpp"
#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/runtime_distributed/find_all_localities.hpp>
#endif
#include "octotigerII/simulation.hpp"
using namespace octotigerII;
namespace {
void exercise(bool species, bool adaptive, bool periodic) {
	auto c = parseConfig({"--problem.name=collapse", "--mesh.cells=4", "--mesh.level=1", "--output.enabled=off", "--verification.analytic=off"});
	if (periodic) c.mesh.boundary = physics::BoundaryConditions::periodic();
	c.massFractions.enabled = species;
	if (species) c.massFractions.species = composition::parseSpecies("gas:0.7:He;oxygen:0.3:O;dye:0.4:A=0,Z=0");
	c.amr.enabled = adaptive; c.amr.maxLevel = 2; c.amr.shadowTolerance = 0; c.amr.bufferCells = 0;
	refinement::Criteria criteria{[](refinement::CellView const& cell) { return cell.center[0] < units::Length::from_value(-6e8) && cell.level < 2 ? Real(2) : Real(0); }};
	Runtime runtime(c, criteria);
	runtime.solveGravity();
	auto initial = diagnose(runtime.snapshots(), c);
	units::Energy solverResidual{};
	for (int step = 0; step < 3; ++step) {
		auto old = runtime.snapshots();
		auto dt = std::min(0.1 * runtime.stableTimestep(), units::Time::from_value(0.1));
		runtime.beginGravityEnergy();
		runtime.kickGravity(dt / 2.0);
		runtime.advance(dt);
		runtime.solveGravity();
		runtime.kickGravity(dt / 2.0);
		auto momentum = diagnose(runtime.snapshots(), c).momentum;
		runtime.finishGravityEnergy(dt);
		auto current = runtime.snapshots();
		auto final = diagnose(current, c);
		// Endpoint work is exactly conservative for a reciprocal potential operator.
		// Measure the existing FMM operator defect independently of the gas update.
		for (std::size_t b = 0; b < old.size(); ++b)
			for (std::size_t i = 0; i < old[b].hydro.values().size(); ++i)
				solverResidual += 0.5 * old[b].layout.cellMeasure(old[b].cellWidth) *
					(old[b].hydro.values()[i].density() * current[b].gravity.values()[i].potential() -
					 current[b].hydro.values()[i].density() * old[b].gravity.values()[i].potential());
		for (int d = 0; d < ndim; ++d) EXPECT_EQ(momentum[d], final.momentum[d]);
		auto b = runtime.boundaryTransport();
		auto energy = final.gasGravityEnergy + b.outward.gasEnergy - b.inward.gasEnergy + b.outward.potentialEnergy - b.inward.potentialEnergy;
		Real const drift = Real((energy - initial.gasGravityEnergy) / initial.gasGravityNorm);
		EXPECT_NEAR(drift, Real(solverResidual / initial.gasGravityNorm), 4e-13) << "step=" << step << " species=" << species << " AMR=" << adaptive << " periodic=" << periodic;
		EXPECT_LT(std::abs(drift), 4e-13);
	}
}
}
TEST(GravityEnergy, IsolatedUniformWithAndWithoutSpecies) { exercise(false, false, false); exercise(true, false, false); }
TEST(GravityEnergy, IsolatedAdaptiveWithAndWithoutSpecies) { exercise(true, true, false); exercise(false, true, false); }
TEST(GravityEnergy, PeriodicUniformWithSpecies) { exercise(true, false, true); }

TEST(GravityEnergy, PeriodicAdaptiveWithSpecies) { exercise(true, true, true); }

TEST(GravityEnergy, RegridRefineAndCoarsenIndependentOfTimeTreatment) {
	for (auto const* treatment : {"mullen", "naive"}) for (bool conserve : {false, true}) for (bool species : {false, true}) {
		auto c = parseConfig({"--problem.name=collapse", "--mesh.cells=4", "--mesh.level=1",
			"--output.enabled=off", "--verification.analytic=off", "--amr.enabled=on", "--amr.minLevel=1",
			"--amr.maxLevel=2", "--amr.shadowTolerance=0", "--amr.bufferCells=0", "--amr.signalBuffer=1"});
		c.gravity.energyTreatment = treatment;
		c.gravity.conserveRegridEnergy = conserve;
		c.massFractions.enabled = species;
		if (species) c.massFractions.species = composition::parseSpecies("gas:0.7:He;oxygen:0.3:O;dye:0.4:A=0,Z=0");
		bool refine = false;
		Runtime runtime(c, { [&refine](refinement::CellView const& cell) {
			return refine && cell.center[0] < units::Length::from_value(-2e8) && cell.level < 2 ? Real(2) : Real(0);
		} });
		runtime.solveGravity();
		for (bool refining : {true, false}) {
			amr::Hierarchy const source(c, runtime.snapshots());
			auto const old = diagnose(runtime.snapshots(), c);
			auto const oldSize = runtime.size();
			refine = refining;
			ASSERT_TRUE(runtime.regrid(units::Time{}, true));
			if (refining) EXPECT_GT(runtime.size(), oldSize);
			else EXPECT_LT(runtime.size(), oldSize);
			if (conserve) {
				EXPECT_THROW(runtime.stableTimestep(), std::logic_error);
				EXPECT_THROW(runtime.advance(units::Time::from_value(1e-4)), std::logic_error);
				EXPECT_THROW(runtime.regrid(units::Time{}, true), std::logic_error);
			}
			runtime.solveGravity();
			auto const now = diagnose(runtime.snapshots(), c);
			if (conserve) for (auto const& block : runtime.snapshots()) {
				auto const expected = source.transferGasGravityEnergy(block.location);
				for (std::size_t i = 0; i < expected.size(); ++i) {
					auto const& gas = block.hydro.values()[i];
					auto const potential = 0.5 * gas.density() * block.gravity.values()[i].potential();
					auto const recovered = gas.totalEnergy() + potential;
					auto const scale = units::abs(gas.totalEnergy()) + units::abs(potential);
					EXPECT_NEAR(Real((recovered - expected[i]) / scale), 0, 8e-16);
				}
			}
			EXPECT_NEAR(Real((now.mass - old.mass) / old.mass), 0, 2e-14);
			auto const drift = Real((now.gasGravityEnergy - old.gasGravityEnergy) / old.gasGravityNorm);
			if (conserve) EXPECT_NEAR(drift, 0, 2e-14) << treatment << " refining=" << refining;
			else {
				EXPECT_NEAR(Real((now.gasEnergy - old.gasEnergy) / old.gasGravityNorm), 0, 2e-14);
				EXPECT_GT(std::abs(drift), 1e-8);
			}
			// No double subtraction on a repeated field solve.
			runtime.solveGravity();
			auto const repeated = diagnose(runtime.snapshots(), c);
			EXPECT_EQ(now.gasEnergy, repeated.gasEnergy);
		}
	}
}

TEST(GravityEnergy, NaiveMatchesUnbracketedLegacyKicks) {
	auto c = parseConfig({"--problem.name=collapse", "--mesh.cells=4", "--mesh.level=1",
		"--gravity.energyTreatment=naive", "--output.enabled=off", "--verification.analytic=off"});
	Runtime legacy(c), selected(c);
	legacy.solveGravity(); selected.solveGravity();
	auto const dt = units::Time::from_value(0.01);
	for (int step = 0; step < 3; ++step) {
		selected.beginGravityEnergy();
		for (auto* runtime : {&legacy, &selected}) {
			runtime->kickGravity(dt / 2.0);
			runtime->advance(dt);
			runtime->solveGravity();
			runtime->kickGravity(dt / 2.0);
		}
		selected.finishGravityEnergy(dt);
		auto const a = legacy.snapshots(), b = selected.snapshots();
		for (std::size_t block = 0; block < a.size(); ++block)
			for (std::size_t i = 0; i < a[block].hydro.values().size(); ++i)
				a[block].hydro.values()[i].forEach([&](auto f, auto value) { EXPECT_EQ(value, b[block].hydro.values()[i].template get<f>()); });
	}
	EXPECT_GT(selected.boundaryTransport().outward.potentialEnergy + selected.boundaryTransport().inward.potentialEnergy, units::Energy{});
}

TEST(GravityEnergy, AdaptivePotentialOperatorIsReciprocal) {
	// Independent bilinear test using two unrelated mass distributions. No gas
	// evolution or energy correction can conceal an asymmetric gravity operator.
	for (int order : {2, 5}) {
		auto c = parseConfig({"--problem.name=collapse", "--mesh.cells=4", "--mesh.level=1", "--amr.enabled=on",
			"--output.enabled=off", "--verification.analytic=off"});
		c.gravity.multipoleOrder = order;
		std::vector<mesh::BlockLocation> leaves;
		mesh::BlockLocation const root{};
		for (int slot = 0; slot < 8; ++slot) {
			auto const child = root.child(slot);
			if (slot == 0) for (int sub = 0; sub < 8; ++sub) leaves.push_back(child.child(sub));
			else leaves.push_back(child);
		}
#ifdef OCTOTIGERII_WITH_HPX
		auto owners = hpx::find_all_localities();
#else
		std::vector<storage::Locality> owners{0};
#endif
		CartesianTopology topology(c, owners.size(), leaves);
		FieldRepository repository(c, topology.storageLayout(), owners);
		auto const& fields = repository.directory();
		gravity::FieldSolver solver(c, topology.blocks(), fields, owners);
		std::array<std::vector<long double>, 2> masses, potentials;
		for (int pattern = 0; pattern < 2; ++pattern) {
			std::size_t index = 0;
			for (auto const& block : topology.blocks()) {
				auto gas = fields.hydro.output(block.interior, 0);
				auto gravity = fields.gravity.output(block.interior, 0);
				for (std::size_t i = 0; i < block.interior.count; ++i, ++index) {
					hydro::PrimitiveState primitive{};
					primitive.density() = units::Density::from_value(1 + Real((index * (pattern ? 43 : 17) + 31) % 127) / 53);
					primitive.pressure() = units::Pressure::from_value(1e12);
					gas.put(i, hydro::HydroSystem(c.hydro).conservedState(primitive));
					gravity.put(i, {});
					masses[pattern].push_back(units::value(primitive.density() * block.layout.cellMeasure(block.cellWidth)));
				}
				fields.hydro.commit(block.interior, 0, gas);
				fields.gravity.commit(block.interior, 0, gravity);
			}
			solver.solve(0);
			for (auto const& block : topology.blocks()) {
				auto const field = fields.gravity.read(block.interior, 1).get();
				for (std::size_t i = 0; i < block.interior.count; ++i) potentials[pattern].push_back(units::value(field.at(i).potential()));
			}
		}
		long double ab = 0, ba = 0;
		for (std::size_t i = 0; i < masses[0].size(); ++i) {
			ab += masses[0][i] * potentials[1][i];
			ba += masses[1][i] * potentials[0][i];
		}
		EXPECT_LT(std::abs((ab - ba) / ab), 3e-14L) << "order=" << order;
	}
}
