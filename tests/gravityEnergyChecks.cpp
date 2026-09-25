#include "testSupport.hpp"
#include "octotigerII/runtime.hpp"
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
		EXPECT_LT(std::abs(drift), adaptive ? 1e-9 : 4e-13);
	}
}
}
TEST(GravityEnergy, IsolatedUniformWithAndWithoutSpecies) { exercise(false, false, false); exercise(true, false, false); }
TEST(GravityEnergy, IsolatedAdaptiveWithAndWithoutSpecies) { exercise(true, true, false); exercise(false, true, false); }
TEST(GravityEnergy, PeriodicUniformWithSpecies) { exercise(true, false, true); }

TEST(GravityEnergy, PeriodicAdaptiveWithSpecies) { exercise(true, true, true); }
