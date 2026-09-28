#include "testSupport.hpp"
#include "octotigerII/problems.hpp"
#include "octotigerII/radiation/coupledPatch.hpp"
#include "octotigerII/runtime.hpp"
#include "octotigerII/simulation.hpp"
#include "octotigerII/output.hpp"
#include "octotigerII/radiation/couplingDiagnostics.hpp"
#include <fstream>
#include <map>
#include <sstream>
#include <iostream>
#include <set>

using namespace octotigerII;

TEST(PhotonHeating, ExistingProblemsKeepConstantOpacityAndNoPhotonSource) {
	auto config = parseConfig({"--problem.name=radiation-matter", "--radiation.opacity=1.25"});
	EXPECT_FALSE(problemHasRadiationMaterial(config));
	auto const material = problemRadiationMaterial(config);
	mesh::PhysicalCoordinates position{};
	for (int d = 0; d < ndim; ++d) position[d] = units::Length::from_value(Real(d + 1) * 1e8);
	auto const value = checkedRadiationMaterial(material, position, units::Time::from_value(4));
	EXPECT_EQ(units::value(value.opacity), 1.25);
	EXPECT_EQ(units::value(value.photonPower), 0);
}

TEST(PhotonHeating, PrescribedPhotonDriveHasCorrectReducedSpeedEnergyAndNoDirectMomentum) {
	using PowerDensity = units::Quantity<-1, 1, -3>;
	hydro::HydroSystem gasSystem;
	hydro::PrimitiveState primitive;
	primitive.density() = units::Density::from_value(1e-7);
	auto const temperature = units::Temperature::from_value(1e5);
	primitive.pressure() = primitive.density() * constants::boltzmann * temperature / constants::atomicMassUnit;
	auto const initialGas = gasSystem.conservedState(primitive);
	radiation::RadiationSystem::State initialRadiation;
	initialRadiation.energy() = constants::radiation * boost::units::pow<4>(temperature);
	mesh::MeshLayout layout(4, 4);
	auto const time = units::Time::from_value(0.02);
	auto const interval = units::Time::from_value(0.01);
	auto const scale = units::Time::from_value(0.04);
	auto const power = PowerDensity::from_value(1e8);
	for (Real ratio : {Real(1), Real(0.2)}) for (Real opacity : {Real(0), Real(1)}) {
		radiation::RadiationSystem radSystem(ratio * constants::c);
		hydro::Fields gas(layout, units::Length::from_value(1e11));
		radiation::Fields rad(layout, gas.cellWidth());
		gas.values().assign(layout.cellCount(), initialGas);
		rad.values().assign(layout.cellCount(), initialRadiation);
		ProblemRadiationMaterial const material = [=](auto const&, units::Time at) {
			return RadiationMaterial{radiation::Opacity::from_value(opacity), power * (1 + at / scale)};
		};
		auto const physicalInjection = interval * power * (1 + (time + interval / 2.0) / scale);
		auto const initialEnergy = initialGas.totalEnergy() + initialRadiation.energy() / ratio;
		radiation::CoupledPatchWorkspace workspace;
		radiation::advanceCoupledPatch(gas, rad, gasSystem, radSystem, {}, interval, workspace,
			[&](auto const&, hydro::ConservedState const& updatedGas, radiation::RadiationSystem::State const& updatedRad) {
				auto const residual = updatedGas.totalEnergy() + updatedRad.energy() / ratio - initialEnergy - physicalInjection;
				EXPECT_NEAR(units::value(residual), 0, 64 * epsilonR * units::value(initialEnergy + physicalInjection));
				for (int d = 0; d < ndim; ++d) {
					EXPECT_EQ(units::value(updatedGas.momentum(d)), 0);
					EXPECT_EQ(units::value(updatedRad.radiativeFlux(d)), 0);
				}
				if (opacity == 0) {
					test::expectStateNear(updatedGas, initialGas);
					EXPECT_NEAR(units::value(updatedRad.energy() - initialRadiation.energy() - ratio * physicalInjection),
						0, 32 * epsilonR * units::value(updatedRad.energy()));
				} else {
					EXPECT_GT(updatedGas.totalEnergy(), initialGas.totalEnergy());
					EXPECT_TRUE(gasSystem.admissible(updatedGas));
					EXPECT_TRUE(radSystem.admissible(updatedRad));
				}
			}, finiteVolume::RotatingFrame{}, time, {}, material);
	}
}

TEST(PhotonHeating, NegativePrescribedHeatingOrOpacityIsRejected) {
	for (bool negativeOpacity : {false, true}) {
		ProblemRadiationMaterial material = [=](auto const&, auto) {
			return RadiationMaterial{units::Quantity<2, -1, 0>::from_value(negativeOpacity ? -1 : 0),
				units::Quantity<-1, 1, -3>::from_value(negativeOpacity ? 0 : -1)};
		};
		EXPECT_THROW(checkedRadiationMaterial(material, {}, {}), std::invalid_argument);
	}
}

TEST(PhotonHeating, CsvSeparatesPrescribedSourceFromBoundaryTransport) {
	test::TemporaryDirectory directory;
	auto config = parseConfig({"--problem.name=radiation-matter", "--mesh.cells=4", "--mesh.level=0", "--output.enabled=off"});
	config.radiation.lightSpeedRatio = 0.2;
	config.output.directory = directory.path.string();
	auto block = initialSnapshot(config, {});
	auto const before = diagnose({block}, config);
	{
		Output output(config);
		output({block}, 0, before);
		for (auto& rad : block.radiation.values()) rad.energy() *= 1.2;
		block.time = units::Time::from_value(1);
		auto after = diagnose({block}, config);
		after.radiationSourceEnergy = 0.2 * before.radiationEnergy;
		output({block}, 1, after);
	}
	std::ifstream input(directory.path / "conservation.csv");
	std::string header, row;
	ASSERT_TRUE(std::getline(input, header));
	ASSERT_TRUE(std::getline(input, row));
	ASSERT_TRUE(std::getline(input, row));
	std::istringstream names(header), values(row);
	std::map<std::string, Real> columns;
	std::string name, value;
	while (std::getline(names, name, ',')) {
		ASSERT_TRUE(std::getline(values, value, ','));
		ASSERT_TRUE(columns.emplace(name, std::stod(value)).second);
	}
	auto const injection = units::value(0.2 * before.radiationEnergy);
	EXPECT_NEAR(columns.at("radiation_source_energy_erg"), injection, 4 * epsilonR * injection);
	EXPECT_NEAR(columns.at("rsla_source_energy_erg"), injection / config.radiation.lightSpeedRatio,
		8 * epsilonR * injection / config.radiation.lightSpeedRatio);
	for (auto const* field : {"radiation_energy_erg", "physical_total_energy_erg", "rsla_total_energy_erg"}) {
		EXPECT_EQ(columns.at(std::string(field) + "_in"), 0);
		EXPECT_EQ(columns.at(std::string(field) + "_out"), 0);
		EXPECT_NEAR(columns.at(std::string(field) + "_drift_scaled"), 0, 3e-14);
	}
}

TEST(PhotonHeating, OpticalDepthDiagnosticUsesPrescribedLocalOpacity) {
	auto config = parseConfig({"--problem.name=radiation-matter", "--mesh.cells=4", "--mesh.level=0"});
	auto const block = initialSnapshot(config, {});
	auto const& gas = block.hydro.values().front();
	auto const& rad = block.radiation.values().front();
	auto const transparent = radiation::couplingDiagnostics(gas, rad, block.cellWidth, config, {});
	auto const opaque = radiation::couplingDiagnostics(gas, rad, block.cellWidth, config,
		radiation::Opacity::from_value(2 * config.radiation.opacity));
	EXPECT_EQ(transparent.cellOpticalDepth, 0);
	EXPECT_EQ(transparent.trappingParameter, 0);
	auto const expected = gas.density() * radiation::Opacity::from_value(2 * config.radiation.opacity) * block.cellWidth;
	EXPECT_NEAR(opaque.cellOpticalDepth, units::value(expected), 4 * epsilonR * opaque.cellOpticalDepth);
}

namespace {
Config sourceConfiguration(Real ratio, Real opacity, bool adaptive = false) {
	auto config = parseConfig({"--problem.name=photon-source", "--mesh.cells=4", "--output.enabled=off"});
	config.radiation.lightSpeedRatio = ratio;
	config.radiation.opacity = opacity;
	config.mesh.level = adaptive ? 1 : 0;
	config.amr.enabled = adaptive;
	config.amr.minLevel = config.mesh.level;
	config.amr.maxLevel = 2;
	config.amr.bufferCells = 0;
	config.amr.shadowTolerance = 0;
	config.runtime.stopTime = {};
	config.validate();
	return config;
}

units::Volume volumeOf(std::vector<Snapshot> const& blocks) {
	units::Volume result{};
	for (auto const& block : blocks) result += Real(block.layout.interiorCellCount()) * block.layout.cellMeasure(block.cellWidth);
	return result;
}

void sourceBudget(Runtime const& runtime, Diagnostics const& before, units::Volume volume, Config const& config) {
	auto const snapshots = runtime.snapshots();
	auto const after = diagnose(snapshots, config);
	auto const expected = config.radiation.lightSpeedRatio * units::Quantity<-1, 1, -3>::from_value(1) * after.time * volume;
	auto const injected = runtime.radiationSourceEnergy();
	EXPECT_NEAR(Real(injected / expected), 1, 3e-13);
	auto const boundary = runtime.boundaryTransport();
	auto const boundaryEnergy = boundary.outward.gasEnergy - boundary.inward.gasEnergy
		+ (boundary.outward.radiationEnergy - boundary.inward.radiationEnergy) / config.radiation.lightSpeedRatio;
	Real const error = Real((after.rslaTotalEnergy + boundaryEnergy - injected / config.radiation.lightSpeedRatio - before.rslaTotalEnergy)
		/ before.rslaTotalEnergyNorm);
	EXPECT_NEAR(error, 0, 3e-12);
	std::cout << "photon-source energy residual ratio=" << config.radiation.lightSpeedRatio
		<< " opacity=" << config.radiation.opacity << " residual=" << error << '\n';
}
}

TEST(PhotonHeating, TransparentRuntimeInjectsAtZeroOpacityAndRunReportsItsLedger) {
	for (Real ratio : {Real(1), Real(0.2)}) {
		auto config = sourceConfiguration(ratio, 0);
		ASSERT_TRUE(problemHasRadiationMaterial(config));
		Runtime runtime(config);
		auto const initial = runtime.snapshots();
		auto const before = diagnose(initial, config);
		auto const volume = volumeOf(initial);
		auto const step = 0.2 * runtime.stableTimestep();
		// The explicit local-exchange API does not advance time or add photons.
		runtime.coupleRadiation(step);
		EXPECT_EQ(runtime.radiationSourceEnergy(), units::Energy{});
		for (int i = 0; i < 3; ++i) runtime.advanceCoupled(step);
		sourceBudget(runtime, before, volume, config);
		auto const actual = runtime.snapshots();
		auto const injection = ratio * units::Quantity<-1, 1, -3>::from_value(1) * (3.0 * step);
		for (std::size_t cell = 0; cell < actual[0].hydro.values().size(); ++cell) {
			test::expectStateNear(actual[0].hydro.values()[cell], initial[0].hydro.values()[cell]);
			EXPECT_NEAR(units::value(actual[0].radiation.values()[cell].energy() - initial[0].radiation.values()[cell].energy() - injection),
				0, 32 * epsilonR * units::value(actual[0].radiation.values()[cell].energy()));
			for (int d = 0; d < ndim; ++d) EXPECT_EQ(actual[0].radiation.values()[cell].radiativeFlux(d), units::EnergyFlux{});
		}
		config.runtime.stopTime = step;
		auto const result = run(config);
		EXPECT_NEAR(Real(result.final.radiationSourceEnergy / (ratio * units::Quantity<-1, 1, -3>::from_value(1) * step * volume)), 1, 2e-13);
	}
}

TEST(PhotonHeating, AcceptedSourceLedgerSurvivesSubcyclingAndRegridding) {
	for (auto const [ratio, opacity] : {std::pair<Real, Real>{1, 1}, {0.2, 1}, {0.2, 0}}) {
		auto config = sourceConfiguration(ratio, opacity, true);
		bool refine = true;
		refinement::Criteria criteria{[&](refinement::CellView const& cell) {
			return refine && cell.center[0] < config.mesh.lower + 0.25 * (config.mesh.upper - config.mesh.lower)
				&& cell.level < 2 ? Real(2) : Real(0);
		}};
		Runtime runtime(config, criteria);
		auto const initial = runtime.snapshots();
		std::set<int> levels;
		for (auto const& block : initial) levels.insert(block.location.level);
		ASSERT_EQ(levels, (std::set<int>{1, 2}));
		auto const before = diagnose(initial, config);
		auto const volume = volumeOf(initial);
		runtime.advanceCoupled(0.2 * runtime.stableTimestep());
		sourceBudget(runtime, before, volume, config);
		auto const steps = runtime.statistics().levelSteps;
		ASSERT_GT(steps.size(), 2u);
		EXPECT_GE(steps[2], 2 * steps[1]);
		auto const ledger = runtime.radiationSourceEnergy();
		refine = false;
		runtime.regrid({}, true);
		EXPECT_EQ(runtime.radiationSourceEnergy(), ledger);
		sourceBudget(runtime, before, volume, config);
		refine = true;
		runtime.regrid({}, true);
		EXPECT_EQ(runtime.radiationSourceEnergy(), ledger);
		runtime.advanceCoupled(0.2 * runtime.stableTimestep());
		sourceBudget(runtime, before, volume, config);
	}
}

TEST(PhotonHeating, FailedNumericalIntervalRestoresPhotonLedgerAndState) {
	auto config = sourceConfiguration(0.2, 1);
	config.hydro.acceleration[0] = units::Acceleration::from_value(100);
	Runtime runtime(config), reference(config);
	auto const step = 0.1 * runtime.stableTimestep();
	runtime.advanceCoupled(step);
	reference.advanceCoupled(step);
	auto const before = runtime.snapshots();
	auto const ledger = runtime.radiationSourceEnergy();
	auto const generation = runtime.generation();
	EXPECT_GT(ledger, units::Energy{});
	// A positive finite interval reaches numerical stages but grossly exceeds
	// the acceleration/transport bounds. No provisional photons may survive it.
	EXPECT_THROW(runtime.advanceCoupled(1e10 * step), std::exception);
	EXPECT_EQ(runtime.radiationSourceEnergy(), ledger);
	EXPECT_EQ(runtime.generation(), generation);
	auto same = [](auto const& expected, auto const& actual) {
		ASSERT_EQ(expected.size(), actual.size());
		for (std::size_t b = 0; b < expected.size(); ++b) {
			EXPECT_EQ(expected[b].time, actual[b].time);
			for (std::size_t i = 0; i < expected[b].hydro.values().size(); ++i) {
				expected[b].hydro.values()[i].forEach([&](auto f, auto value) { EXPECT_EQ(value, actual[b].hydro.values()[i].template get<f>()); });
				expected[b].radiation.values()[i].forEach([&](auto f, auto value) { EXPECT_EQ(value, actual[b].radiation.values()[i].template get<f>()); });
			}
		}
	};
	same(before, runtime.snapshots());
	runtime.advanceCoupled(step);
	reference.advanceCoupled(step);
	same(reference.snapshots(), runtime.snapshots());
	EXPECT_EQ(runtime.radiationSourceEnergy(), reference.radiationSourceEnergy());
}
