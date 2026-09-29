#include <octotigerII/helmholtz/helmholtz.hpp>
#include "octotigerII/hydro/hydroSystem.hpp"
#include "octotigerII/radiation/matterCoupling.hpp"
#include "octotigerII/runtime.hpp"
#include "octotigerII/simulation.hpp"
#include "testSupport.hpp"
using namespace octotigerII;
namespace {
Config config(bool radiation = false) {
    auto c = parseConfig({"--problem.name=helmholtz-advection", "--output.enabled=off", "--verification.analytic=off"});
    c.hydro.eos = "helmholtz";
    c.radiation.enabled = radiation;
    c.radiation.opacity = radiation ? 1 : 0;
    c.massFractions.enabled = true;
    c.massFractions.species = composition::parseSpecies("carbon:0.5:A=12,Z=6;oxygen:0.5:A=16,Z=8");
    return c;
}
}    // namespace
TEST(HelmholtzHydro, InversionsSoundSpeedAndPhotonExclusion) {
    EXPECT_EQ(Config::HydroOptions{}.eos, "ideal");
    EXPECT_FALSE(hydro::HydroSystem(Config::HydroOptions{}).helmholtz());
    helmholtz::Eos reference;
    for (bool gasOnly : {false, true}) {
        Config::HydroOptions options;
        options.eos = "helmholtz";
        hydro::HydroSystem gas(options, gasOnly);
        for (double rho : {1e-5, 1e3, 1e6}) {
            for (double t : {1e5, 1e7, 1e9}) {
                SCOPED_TRACE(rho);
                SCOPED_TRACE(t);
                SCOPED_TRACE(gasOnly);
                auto state = gas.stateFromTemperature(units::Density::from_value(rho), units::Temperature::from_value(t), 12, 6);
                auto q = reference.evaluate(rho, t, 12, 6);
                EXPECT_NEAR(units::value(gas.temperature(state)) / t, 1, 2e-8);
                auto primitive = gas.reconstructionVariables(state);
                double const p = gasOnly ? q.pgas : q.ptot;
                EXPECT_NEAR(units::value(primitive.pressure()) / p, 1, 2e-10);
                EXPECT_NEAR(std::pow(units::value(gas.adiabaticSoundSpeed(state)), 2) / ((gasOnly ? q.gam1_gas : q.gam1) * p / rho), 1, 2e-8);
                auto rebuilt = gas.conservedState(primitive);
                EXPECT_NEAR(units::value(rebuilt.totalEnergy() / state.totalEnergy()), 1, 2e-9);
                auto fromEntropy = gas.internalEnergyFromAuxiliary(state);
                EXPECT_NEAR(units::value(fromEntropy / state.totalEnergy()), 1, 2e-8);
            }
        }
    }
}
TEST(HelmholtzHydro, AdiabaticSoundSpeedMatchesEntropyDerivative) {
    Config::HydroOptions options;
    options.eos = "helmholtz";
    hydro::HydroSystem gas(options, true);
    auto state = gas.stateFromTemperature(units::Density::from_value(1e5), units::Temperature::from_value(1e8), 12, 6);
    auto pressure = [&](double factor) {
        auto next = state;
        next.density() *= factor;
        next.nuclei() *= factor;
        next.electrons() *= factor;
        next.auxiliary() *= factor;
        next.totalEnergy() = gas.internalEnergyFromAuxiliary(next);
        return units::value(gas.reconstructionVariables(next).pressure());
    };
    double const slope = (pressure(1.0001) - pressure(0.9999)) / 20;
    EXPECT_NEAR(slope / std::pow(units::value(gas.adiabaticSoundSpeed(state)), 2), 1, 2e-5);
}
TEST(HelmholtzHydro, FloorAndDualEnergyRecovery) {
    Config::HydroOptions options;
    options.eos = "helmholtz";
    hydro::HydroSystem gas(options, true);
    auto cold = gas.stateFromTemperature(units::Density::from_value(1), units::Temperature::from_value(1000), 12, 6);
    auto state = cold;
    state.totalEnergy() *= 0.9;
    state.auxiliary() -= units::Density::from_value(1);
    auto const before = state.totalEnergy();
    auto const correction = gas.applyTemperatureFloor(state);
    EXPECT_GT(units::value(correction), 0);
    EXPECT_EQ(state.totalEnergy() - before, correction);
    EXPECT_NEAR(units::value(gas.temperature(state)), 1000, 1e-6);
    EXPECT_EQ(gas.applyTemperatureFloor(state), units::EnergyDensity{});
    auto warm = gas.stateFromTemperature(units::Density::from_value(1), units::Temperature::from_value(1e6), 12, 6);
    warm.totalEnergy() = before;
    EXPECT_EQ(gas.applyTemperatureFloor(warm), units::EnergyDensity{});
    EXPECT_NEAR(units::value(gas.temperature(warm)) / 1e6, 1, 2e-9);
}
TEST(HelmholtzHydro, SpeciesDetermineCompositionAndTracersAreExcluded) {
    auto c = config();
    c.massFractions.species = composition::parseSpecies("helium:0.3:A=4,Z=2;iron:0.7:A=56,Z=26;tag:0.2:A=0,Z=0");
    hydro::HydroSystem gas(c);
    auto state = gas.stateFromTemperature(units::Density::from_value(1e5), units::Temperature::from_value(1e8), 12, 6);
    std::vector<units::Density> species{units::Density::from_value(3e4), units::Density::from_value(7e4), units::Density::from_value(1e9)};
    gas.setComposition(state, species, c.massFractions);
    EXPECT_DOUBLE_EQ(units::value(state.density()), 1e5);
    EXPECT_NEAR(units::value(state.nuclei()), 3e4 / 4 + 7e4 / 56, 1e-10);
    EXPECT_NEAR(units::value(state.electrons()), 3e4 * 2 / 4 + 7e4 * 26 / 56, 1e-10);
}
TEST(HelmholtzHydro, RadiationLteAndRelaxationConserveEnergy) {
    auto c = config(true);
    hydro::HydroSystem gas(c);
    for (double ratio : {1.0, 0.1}) {
        for (double factor : {1.0, 0.5, 2.0}) {
            auto state = gas.stateFromTemperature(units::Density::from_value(1e-5), units::Temperature::from_value(1e6), 12, 6);
            radiation::RadiationSystem::State rad;
            rad.energy() = units::EnergyDensity::from_value(factor * units::value(constants::radiation) * 1e24);
            auto const energy = state.totalEnergy() + rad.energy() / ratio;
            radiation::couple(state, rad, gas, radiation::Opacity::from_value(1), ratio, units::Time::from_value(1e-3));
            EXPECT_NEAR(units::value(state.totalEnergy() + rad.energy() / ratio - energy), 0, 1e-13 * units::value(energy));
            EXPECT_TRUE(gas.admissible(state));
            if (factor == 1) EXPECT_NEAR(units::value(gas.temperature(state)) / 1e6, 1, 1e-9);
        }
    }
    Config::HydroOptions options;
    options.eos = "helmholtz";
    hydro::HydroSystem photons(options, false);
    auto state = photons.stateFromTemperature(units::Density::from_value(1), units::Temperature::from_value(1e6), 12, 6);
    radiation::RadiationSystem::State rad;
    rad.energy() = units::EnergyDensity::from_value(1);
    EXPECT_THROW(radiation::couple(state, rad, photons, radiation::Opacity::from_value(1), 1, units::Time::from_value(1e-6)), std::invalid_argument);
}
TEST(HelmholtzHydro, CoupledRuntimeUsesSpeciesAndClosesBudget) {
    auto c = config(true);
    c.mesh.cells = 4;
    c.validate();
    Runtime runtime(c);
    auto before = diagnose(runtime.snapshots(), c);
    runtime.advanceCoupled(std::min(runtime.stableTimestep(), units::Time::from_value(1e-5)));
    auto snapshots = runtime.snapshots();
    auto after = diagnose(snapshots, c);
    EXPECT_NEAR(units::value(after.rslaTotalEnergy - before.rslaTotalEnergy - runtime.eosFloorEnergy()), 0, 2e-11 * units::value(before.rslaTotalEnergy));
    for (auto const& block : snapshots) {
        block.layout.forEachInterior([&](auto const&, std::size_t i) {
            auto const& state = block.hydro.values()[i];
            EXPECT_NEAR(units::value(state.nuclei()), units::value(block.species[0].values()[i] / 12.0 + block.species[1].values()[i] / 16.0), 1e-12);
        });
    }
}

TEST(HelmholtzHydro, AmrSubcyclingAndRegridPreserveSpeciesAndEnergy) {
    auto c = config(true);
    c.mesh.cells = 4;
    c.mesh.level = 1;
    c.amr.enabled = true;
    c.amr.maxLevel = 2;
    c.amr.bufferCells = 0;
    c.amr.signalBuffer = 1;
    c.amr.shadowTolerance = 0;
    bool refine = true;
    refinement::Criteria criteria{[&](refinement::CellView const& cell) {
        return refine && cell.center[0] < c.mesh.lower + 0.3 * (c.mesh.upper - c.mesh.lower) && cell.level < 2 ? Real(2) : Real(0);
    }};
    Runtime runtime(c, criteria);
    auto totals = [&]() {
        std::array<units::Mass, 2> result{};
        for (auto const& block : runtime.snapshots()) {
            block.layout.forEachInterior([&](auto const&, std::size_t i) {
                for (int s = 0; s < 2; ++s) {
                    result[s] += block.layout.cellMeasure(block.cellWidth) * block.species[s].values()[i];
                }
                auto const& state = block.hydro.values()[i];
                EXPECT_NEAR(units::value(state.nuclei()), units::value(block.species[0].values()[i] / 12.0 + block.species[1].values()[i] / 16.0), 1e-10);
            });
        }
        return result;
    };
    auto initial = diagnose(runtime.snapshots(), c);
    auto mass = totals();
    runtime.advanceCoupled(std::min(runtime.stableTimestep(), units::Time::from_value(1e-5)));
    auto final = diagnose(runtime.snapshots(), c);
    auto after = totals();
    for (int s = 0; s < 2; ++s) {
        EXPECT_NEAR(Real(after[s] / mass[s]), 1, 2e-12);
    }
    EXPECT_NEAR(Real((final.rslaTotalEnergy - runtime.eosFloorEnergy()) / initial.rslaTotalEnergy), 1, 2e-11);
    refine = false;
    runtime.regrid(units::Time{}, true);
    after = totals();
    for (int s = 0; s < 2; ++s) {
        EXPECT_NEAR(Real(after[s] / mass[s]), 1, 2e-12);
    }
    EXPECT_THROW(runtime.advanceCoupled(units::Time::from_value(-1)), std::invalid_argument);
    EXPECT_EQ(runtime.eosFloorEnergy(), units::Energy{});
}
TEST(HelmholtzHydro, AcceptedPatchFloorLedgerMatchesEnergyChange) {
    Config::HydroOptions options;
    options.eos = "helmholtz";
    hydro::HydroSystem gas(options, true);
    mesh::MeshLayout layout(4, 2);
    hydro::Fields patch(layout, units::Length::from_value(1));
    auto state = gas.stateFromTemperature(units::Density::from_value(1), units::Temperature::from_value(1000), 12, 6);
    state.totalEnergy() *= 0.9;
    state.auxiliary() -= units::Density::from_value(1);
    layout.forEachInterior([&](auto const& cell, std::size_t) { patch.atInterior(cell) = state; });
    auto const step = hydro::Solver(gas).advance(patch, units::Time::from_value(1e-12), finiteVolume::BoundaryConditions::periodic());
    units::Energy change{};
    layout.forEachInterior([&](auto const& cell, std::size_t) {
        change += layout.cellMeasure(patch.cellWidth()) * (patch.atInterior(cell).totalEnergy() - state.totalEnergy());
        EXPECT_NEAR(units::value(gas.temperature(patch.atInterior(cell))), 1000, 1e-6);
    });
    EXPECT_GT(units::value(change), 0);
    EXPECT_EQ(step.eosFloorEnergy, change);
    EXPECT_EQ(step.eosFloorCells, layout.interiorCellCount());
}

TEST(HelmholtzHydro, StiffDrivenRadiationAtMixedCompositionConverges) {
    if constexpr (ndim != 1) GTEST_SKIP() << "Recorded 1D source state";
    Config::HydroOptions options;
    options.eos = "helmholtz";
    hydro::HydroSystem system(options, true);
    hydro::ConservedState gas, drive;
    gas.density() = units::Density::from_value(1000.0000012694189);
    gas.momentum(0) = units::MomentumDensity::from_value(999998593.91741753);
    gas.totalEnergy() = units::EnergyDensity::from_value(1.4698427502633623e18);
    gas.auxiliary() = units::Density::from_value(3550.6790691007191);
    gas.nuclei() = units::Density::from_value(176.82297814391799);
    gas.electrons() = units::Density::from_value(488.74199722341154);
    drive.density() = units::Density::from_value(1.5874126687747615e-7);
    drive.momentum(0) = units::MomentumDensity::from_value(-1.6677356038061912);
    drive.totalEnergy() = units::EnergyDensity::from_value(395278179.45271981);
    drive.auxiliary() = units::Density::from_value(5.6353788834110471e-7);
    drive.nuclei() = units::Density::from_value(-2.560296502752957e-8);
    drive.electrons() = units::Density::from_value(6.9326282399959017e-8);
    radiation::RadiationSystem::State rad, radDrive;
    rad.energy() = units::EnergyDensity::from_value(574473957651106);
    rad.radiativeFlux(0) = units::EnergyFlux::from_value(7.6596418053497094e20);
    radDrive.energy() = units::EnergyDensity::from_value(69099.793877394928);
    radDrive.radiativeFlux(0) = units::EnergyFlux::from_value(-7.8225224237395521e22);
    auto const energy = gas.totalEnergy() + drive.totalEnergy() + rad.energy() + radDrive.energy();
    radiation::coupleForced(gas, rad, drive, radDrive, system, radiation::Opacity::from_value(1), 1, units::Time::from_value(0.00014443192705813527));
    EXPECT_NEAR(Real((gas.totalEnergy() + rad.energy()) / energy), 1, 2e-14);
}
