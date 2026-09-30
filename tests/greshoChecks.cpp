#include "testSupport.hpp"
#include "octotigerII/problems.hpp"
#include "octotigerII/simulation.hpp"
#include "octotigerII/verification/analytic.hpp"
#include <cmath>
#include <iostream>
#include <map>

using namespace octotigerII;

TEST(Gresho, EquilibriumAndContinuity) {
    auto c = test::parseConfig({});
    EXPECT_TRUE(c.hydroEnabled());
    EXPECT_FALSE(c.gravityEnabled());
    EXPECT_FALSE(c.radiationEnabled());
    auto ref = problemReference(c);
    auto sample = [&](Real r) {
        mesh::PhysicalCoordinates x{};
        x[0] = units::Length::from_value(r);
        return ref.evaluate(x, units::Time::from_value(1)).hydro;
    };
    EXPECT_DOUBLE_EQ(units::value(sample(0).pressure()), 5);
    EXPECT_NEAR(units::value(sample(0.2).velocity(1)), 1, 1e-15);
    EXPECT_NEAR(units::value(sample(0.4).pressure()), 3 + 4 * std::log(2.), 1e-14);
    for (Real r : {0.05, 0.15, 0.25, 0.35, 0.45}) {
        Real const eps = 1e-6;
        auto q = sample(r);
        Real const dpdr = units::value(sample(r + eps).pressure() - sample(r - eps).pressure()) / (2 * eps);
        Real const v = units::value(q.velocity(1));
        EXPECT_NEAR(dpdr, units::value(q.density()) * v * v / r, 2e-8);
        EXPECT_EQ(q.velocity(0), units::Velocity{});
        if constexpr (ndim == 3) EXPECT_EQ(q.velocity(2), units::Velocity{});
    }
    for (Real r : {0.2, 0.4}) {
        EXPECT_NEAR(units::value(sample(r - 1e-10).pressure()), units::value(sample(r + 1e-10).pressure()), 2e-9);
        EXPECT_NEAR(units::value(sample(r - 1e-10).velocity(1)), units::value(sample(r + 1e-10).velocity(1)), 2e-9);
    }
    for (auto argument : {"--amr.refineSpeed=-1", "--amr.refineSpeed=nan", "--radiation.enabled=on", "--hydro.acceleration.x=1"})
        EXPECT_THROW(test::parseConfig({argument}), std::invalid_argument);
}

TEST(Gresho, MixedMeshSubcyclesAndConserves) {
    auto c = test::parseConfig({"--mesh.cells=8", "--mesh.level=2", "--amr.enabled=on", "--amr.minLevel=2",
        "--amr.maxLevel=3", "--amr.refineSpeed=0.9", "--amr.bufferCells=0", "--amr.regridEvery=1", "--output.enabled=off"});
    Runtime runtime(c);
    EXPECT_EQ(runtime.shadowCellCount(), 0u);
    std::map<int, int> levels;
    for (auto const& b : runtime.snapshots()) ++levels[b.location.level];
    for (auto [level, count] : levels) std::cout << "Gresho initial level=" << level << " blocks=" << count << '\n';
    ASSERT_EQ(levels.size(), 2u);
    ASSERT_GT(levels[2], 0);
    ASSERT_GT(levels[3], 0);
    auto before = diagnose(runtime.snapshots(), c);
    auto const dt = runtime.stableTimestep();
    runtime.advance(dt);
    auto counts = runtime.statistics().levelSteps;
    ASSERT_EQ(counts.size(), 4u);
    EXPECT_EQ(counts[2], 1u);
    EXPECT_GE(counts[3], 2u);
    std::cout << "Gresho accepted level-2 steps=" << counts[2] << " level-3 steps=" << counts[3] << '\n';
    runtime.regrid(runtime.stableTimestep(), true);
    auto after = diagnose(runtime.snapshots(), c);
    EXPECT_NEAR(units::value(after.mass - before.mass), 0, 3e-12);
    EXPECT_NEAR(units::value(after.gasEnergy - before.gasEnergy), 0, 3e-11);
    for (int d = 0; d < ndim; ++d) EXPECT_NEAR(units::value(after.momentum[d] - before.momentum[d]), 0, 3e-12);
    EXPECT_GT(after.minimumDensity, units::Density{});
    EXPECT_GT(after.minimumPressure, units::Pressure{});
    for (auto const& b : runtime.snapshots()) EXPECT_EQ(b.time, dt);
}

#if OCTOTIGERII_NDIM == 2
TEST(Gresho, VelocityErrorDecreasesWithResolution) {
    Real errors[2]{};
    for (int i = 0; i < 2; ++i) {
        auto c = test::parseConfig({"--mesh.cells=" + std::to_string(8 << i), "--mesh.level=1",
            "--runtime.stopTime=0.05", "--output.enabled=off"});
        auto result = run(c);
        auto comparison = verification::compare(result.snapshots, c);
        for (auto const& f : comparison.fields)
            if (f.name.find("velocity") != std::string::npos) errors[i] += f.l1;
        std::cout << "Gresho resolution=" << (16 << i) << " velocity L1 sum=" << errors[i] << '\n';
    }
    EXPECT_GT(errors[0], 0);
    EXPECT_GT(errors[1], 0);
    EXPECT_LT(errors[1], errors[0]);
}
#endif
