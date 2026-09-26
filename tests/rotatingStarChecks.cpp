#include "testSupport.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <fstream>
#include <numbers>
#include "octotigerII/problems.hpp"
#include "octotigerII/problems/rotatingStar.hpp"
#include "octotigerII/runtime.hpp"
using namespace octotigerII;

TEST(RotatingStar, OriginalTableAndOblateShapeArePreserved) {
	auto const centerCell = problems::RotatingStar::profile(.005, .005);
	EXPECT_DOUBLE_EQ(centerCell.density, 1.0000118570968204);
	EXPECT_DOUBLE_EQ(centerCell.internalEnergy, .25219165291024365);
	auto const interior = problems::RotatingStar::profile(.205, .205);
	EXPECT_NEAR(interior.density, .5249722914563402, 2e-15);
	EXPECT_NEAR(interior.internalEnergy, .08615545127846742, 2e-15);
	EXPECT_GT(problems::RotatingStar::profile(.705, .005).density, .02);
	EXPECT_EQ(problems::RotatingStar::profile(.005, .705).density, 0);
	EXPECT_EQ(problems::RotatingStar::profile(2, 0).density, 0);
	EXPECT_EQ(problems::RotatingStar::profile(0, 2).density, 0);
	EXPECT_DOUBLE_EQ(problems::RotatingStar::profile(-.205, -.205).density, interior.density);
}

TEST(RotatingStar, DensityPressureVelocityAndKineticEnergyScaleInCgs) {
	auto const length = units::Length::from_value(1e9);
	auto const density = units::Density::from_value(1e6);
	problems::RotatingStar const star(length, density);
	EXPECT_NEAR(units::value(star.angularVelocity()), .13319147161856015, 1e-16);
	mesh::PhysicalCoordinates x{};
	x[0] = .305 * length;
	x[1] = .405 * length;
	x[2] = .105 * length;
	auto const primitive = star(x);
	EXPECT_NEAR(Real(primitive.velocity(0) / (-star.angularVelocity()*x[1])), 1, 1e-14);
	EXPECT_NEAR(Real(primitive.velocity(1) / (star.angularVelocity()*x[0])), 1, 1e-14);
	auto const raw = problems::RotatingStar::profile(std::hypot(.305, .405), .105);
	hydro::HydroSystem const gas(problems::RotatingStar::gamma);
	auto const q = gas.conservedState(primitive);
	auto const kinetic = (q.momentum(0)*q.momentum(0)+q.momentum(1)*q.momentum(1))/(2.0*q.density());
	EXPECT_NEAR(Real((q.energy()-kinetic)/(raw.internalEnergy*star.energyUnit())), 1, 2e-15);
	problems::RotatingStar const scaled(2.0*length, 4.0*density);
	for (auto& coordinate : x) coordinate *= 2.0;
	auto const other = scaled(x);
	EXPECT_NEAR(Real(other.density()/primitive.density()), 4, 1e-14);
	EXPECT_NEAR(Real(other.pressure()/primitive.pressure()), 64, 1e-13);
	EXPECT_NEAR(Real(other.velocity(1)/primitive.velocity(1)), 4, 1e-14);
}

TEST(RotatingStar, VacuumAtmosphereUsesOriginalFloorsAndIsInitiallyInertialStationary) {
	auto const length = units::Length::from_value(1e9);
	auto const density = units::Density::from_value(1e6);
	problems::RotatingStar const star(length, density);
	mesh::PhysicalCoordinates x{};
	x[0] = 1.5 * length;
	auto const value = star(x);
	EXPECT_EQ(value.density(), 1e-10*density);
	EXPECT_EQ(value.pressure(), (problems::RotatingStar::gamma-1)*1e-10*star.energyUnit());
	for (int d = 0; d < 3; ++d) EXPECT_EQ(value.velocity(d), units::Velocity{});
}

TEST(RotatingStar, ScfEnthalpyBalancesIndependentRingPotentialAndCentrifugalPotential) {
	// Independent cylindrical volume quadrature of the original density table.
	// This checks that the imported profile is the oblate SCF equilibrium, not a
	// spherical profile with an arbitrary spin pasted onto it. Finite table and
	// quadrature resolution set the tolerance, rather than floating-point error.
	auto potential = [](Real r, Real z) {
		constexpr int azimuths = 96;
		Real phi = 0;
		for (int i = 0; i < 100; ++i) for (int k = 0; k < 100; ++k) {
			Real const sourceR = (i+.5)/100, sourceZ = (k+.5)/100;
			Real const rho = problems::RotatingStar::profile(sourceR, sourceZ).density;
			if (rho == 0) continue;
			Real const mass = rho*sourceR*1e-4*2*std::numbers::pi_v<Real>/azimuths;
			for (int a = 0; a < azimuths; ++a) {
				Real const angle = (a+.5)*2*std::numbers::pi_v<Real>/azimuths;
				Real const transverse = r*r+sourceR*sourceR-2*r*sourceR*std::cos(angle);
				phi -= mass/std::sqrt(transverse+(z-sourceZ)*(z-sourceZ));
				phi -= mass/std::sqrt(transverse+(z+sourceZ)*(z+sourceZ));
			}
		}
		return phi;
	};
	auto bernoulli = [&](Real r, Real z) {
		auto const value = problems::RotatingStar::profile(r, z);
		return problems::RotatingStar::gamma*value.internalEnergy/value.density + potential(r, z)
			- .5*problems::RotatingStar::spin*problems::RotatingStar::spin*r*r;
	};
	Real const center = bernoulli(0, 0);
	EXPECT_NEAR(bernoulli(.325, .215), center, 2e-3);
	EXPECT_NEAR(bernoulli(.525, .115), center, 2e-3);
}

TEST(RotatingStar, IniSelectsProblemAndScalesDomainWithoutDependingOnWorkingDirectory) {
	test::TemporaryDirectory directory;
	auto const path = directory.path / "rotating.ini";
	{ std::ofstream out(path); out << "problem.name=rotatingStar\nstar.radius=2e9\nstar.centralDensity=2e6\nstar.center.x=3e8\n"; }
	auto const c = octotigerII::parseConfig({"--config=" + path.string()});
	EXPECT_EQ(c.problem, "rotatingStar");
	EXPECT_EQ(c.mesh.lower, units::Length::from_value(-4e9));
	EXPECT_EQ(c.mesh.upper, units::Length::from_value(4e9));
	EXPECT_EQ(c.star.center[0], units::Length::from_value(3e8));
	EXPECT_EQ(c.star.center[1], units::Length{});
	EXPECT_EQ(c.hydro.gamma, problems::RotatingStar::gamma);
	EXPECT_EQ(c.star.atmosphereFraction, 1e-10);
	EXPECT_NO_THROW(initialSnapshot(c, {0, {}}));
}

TEST(RotatingStar, StartupProbeFindsOffCenterCoreAndFinalStateUsesPhysicalScale) {
	auto const c = test::parseConfig({"--mesh.cells=4", "--mesh.level=0", "--amr.enabled=on", "--amr.maxLevel=4",
		"--amr.bufferCells=0", "--star.radius=2.2e8", "--mesh.lower=-1e9", "--mesh.upper=1e9",
		"--star.center.x=1.3e8", "--star.center.y=-1.7e8", "--runtime.stopTime=0", "--output.enabled=off"});
	Runtime runtime(c);
	problems::RotatingStar const star(c.star.radius, c.star.centralDensity, c.star.atmosphereFraction);
	hydro::HydroSystem const gas(c.hydro);
	units::Density peak{};
	bool finest = false;
	for (auto const& block : runtime.snapshots()) {
		finest = finest || block.location.level == c.amr.maxLevel;
		block.layout.forEachInterior([&](auto const& cell, std::size_t i) {
			auto relative = block.layout.cellCenter(block.lower, block.cellWidth, cell);
			for (int d = 0; d < 3; ++d) relative[d] -= c.star.center[d];
			test::expectStateNear(block.hydro.values()[i], gas.conservedState(star(relative)));
			peak = std::max(peak, block.hydro.values()[i].density());
		});
	}
	EXPECT_TRUE(finest);
	EXPECT_GT(peak, .5*c.star.centralDensity);
}

TEST(RotatingStar, StellarSpinIsIndependentOfGridRotation) {
	auto const inertial = test::parseConfig({"--mesh.level=0", "--frame.omega=0"});
	auto const rotating = test::parseConfig({"--mesh.level=0", "--frame.omega=0.13319147161856015"});
	auto const a = initialSnapshot(inertial, {0, {}}), b = initialSnapshot(rotating, {0, {}});
	for (std::size_t i = 0; i < a.hydro.values().size(); ++i)
		test::expectStateNear(a.hydro.values()[i], b.hydro.values()[i], 0);
}

TEST(RotatingStar, RejectsModelChangesThatWouldInvalidateTheFixedEquilibrium) {
	EXPECT_THROW(test::parseConfig({"--star.polytropicIndex=3"}), std::invalid_argument);
	EXPECT_THROW(test::parseConfig({"--hydro.gamma=1.4"}), std::invalid_argument);
	EXPECT_THROW(test::parseConfig({"--mesh.boundary.xLower=periodic", "--mesh.boundary.xUpper=periodic"}), std::invalid_argument);
	EXPECT_THROW(test::parseConfig({"--star.radius=0"}), std::invalid_argument);
}
