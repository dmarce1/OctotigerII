#include "testSupport.hpp"
#include "octotigerII/output.hpp"
#include "octotigerII/radiation/couplingDiagnostics.hpp"
#include <fstream>
#include <map>
#include <memory>
#include <silo.h>
#include <sstream>

using namespace octotigerII;

namespace {
void near(Real actual, Real expected) {
	EXPECT_NEAR(actual, expected, 3e-13 * std::max({Real(1), std::abs(actual), std::abs(expected)}));
}

Config configuration() {
	return parseConfig({"--problem.name=sod", "--radiation.enabled=on", "--mesh.cells=4", "--mesh.level=0",
		"--output.enabled=off", "--verification.analytic=off"});
}

hydro::ConservedState gasState(Config const& c) {
	hydro::PrimitiveState p;
	p.density() = units::Density::from_value(2);
	p.pressure() = units::Pressure::from_value(10);
	p.velocity(0) = units::Velocity::from_value(-2);
	return hydro::HydroSystem(c.hydro).conservedState(p);
}

radiation::RadiationSystem::State radiationState() {
	radiation::RadiationSystem::State result;
	result.energy() = units::EnergyDensity::from_value(30 * units::value(constants::c));
	result.radiativeFlux(0) = Real(0.2) * constants::c * result.energy();
	return result;
}

Snapshot snapshot(Config const& c) {
	auto block = initialSnapshot(c, {});
	for (auto& gas : block.hydro.values()) { gas = gasState(c); }
	for (auto& rad : block.radiation.values()) { rad = radiationState(); }
	return block;
}

std::vector<std::string> split(std::string const& line) {
	std::vector<std::string> values;
	std::istringstream input(line);
	std::string value;
	while (std::getline(input, value, ',')) { values.push_back(value); }
	return values;
}
}

TEST(RadiationConservation, PhysicalAndWeightedIntegralsUsePhysicalFluxAndInertialPosition) {
	for (Real ratio : {Real(1), Real(0.125)}) {
		auto c = configuration();
		c.radiation.lightSpeedRatio = ratio;
		c.radiation.opacity = 3;
		c.radiation.diagnosticLength = 5;
		if constexpr (ndim >= 2) c.frame.omega = units::InverseTime::from_value(0.4);
		auto block = snapshot(c);
		block.time = units::Time::from_value(0.7);
		// Shift the geometry so angular momentum does not cancel by symmetry.
		if constexpr (ndim >= 2) block.lower[1] += units::Length::from_value(3);
		auto const d = diagnose({block}, c);
		auto const volume = block.layout.cellMeasure(block.cellWidth);
		Real const measure = units::value(volume) * block.layout.interiorCellCount();
		Real const gasEnergy = units::value(gasState(c).totalEnergy()) * measure;
		Real const radiationEnergy = units::value(radiationState().energy()) * measure;
		near(units::value(d.physicalTotalEnergy), gasEnergy + radiationEnergy);
		near(units::value(d.rslaTotalEnergy), gasEnergy + radiationEnergy / ratio);
		near(units::value(d.physicalTotalEnergyNorm), gasEnergy + radiationEnergy);
		near(units::value(d.rslaTotalEnergyNorm), gasEnergy + radiationEnergy / ratio);
		near(units::value(d.physicalTotalMomentum[0]), 2 * measure); // -4 gas + 6 radiation
		near(units::value(d.rslaTotalMomentum[0]), (-4 + 6 / ratio) * measure);
		near(units::value(d.physicalTotalMomentumNorm[0]), 10 * measure);
		near(units::value(d.rslaTotalMomentumNorm[0]), (4 + 6 / ratio) * measure);
		if constexpr (ndim >= 2) {
			Real gasAngular = 0, radiationAngular = 0;
			block.layout.forEachInterior([&](auto const& cell, std::size_t) {
				auto const x = physics::RotatingFrame(c.frame.omega).toInertial(
					block.layout.cellCenter(block.lower, block.cellWidth, cell), block.time);
				gasAngular += units::value(volume) * units::value(x[1]) * 4;
				radiationAngular -= units::value(volume) * units::value(x[1]) * 6;
			});
			near(units::value(d.radiationAngularMomentumZ), radiationAngular);
			near(units::value(d.physicalTotalAngularMomentumZ), gasAngular + radiationAngular);
			near(units::value(d.rslaTotalAngularMomentumZ), gasAngular + radiationAngular / ratio);
		}
		auto const estimate = radiation::couplingDiagnostics(gasState(c), radiationState(), block.cellWidth, c);
		near(d.maximumCellOpticalDepth, estimate.cellOpticalDepth);
		near(d.maximumTrappingParameter, estimate.trappingParameter);
		near(d.maximumRslaCriterion, estimate.rslaCriterion);
	}
}

TEST(RadiationConservation, CombinedCsvUsesBothBoundaryLedgersWithoutChangingLegacyColumns) {
	test::TemporaryDirectory directory;
	auto c = configuration();
	c.radiation.lightSpeedRatio = 0.25;
	c.output.directory = directory.path.string();
	auto block = snapshot(c);
	auto d = diagnose({block}, c);
	d.boundary.inward.gasEnergy = units::Energy::from_value(2);
	d.boundary.outward.gasEnergy = units::Energy::from_value(3);
	d.boundary.inward.potentialEnergy = units::Energy::from_value(5);
	d.boundary.outward.potentialEnergy = units::Energy::from_value(7);
	d.boundary.inward.radiationEnergy = units::Energy::from_value(11);
	d.boundary.outward.radiationEnergy = units::Energy::from_value(13);
	d.boundary.inward.momentum[0] = units::Momentum::from_value(17);
	d.boundary.outward.momentum[0] = units::Momentum::from_value(19);
	d.boundary.inward.radiationFlux[0] = units::Momentum::from_value(23) * constants::c * constants::c;
	d.boundary.outward.radiationFlux[0] = units::Momentum::from_value(29) * constants::c * constants::c;
	{
		Output output(c);
		output({block}, 0, d);
	}
	std::ifstream file(directory.path / "conservation.csv");
	std::string header, row;
	ASSERT_TRUE(std::getline(file, header));
	ASSERT_TRUE(std::getline(file, row));
	auto const names = split(header), values = split(row);
	ASSERT_EQ(names.size(), values.size());
	std::map<std::string, Real> fields;
	for (std::size_t i = 0; i < names.size(); ++i) { ASSERT_TRUE(fields.emplace(names[i], std::stod(values[i])).second); }
	EXPECT_LT(header.find("radiation_energy_erg_grid"), header.find("physical_total_energy_erg_grid"));
	near(fields.at("physical_total_energy_erg_in"), 18);
	near(fields.at("physical_total_energy_erg_out"), 23);
	near(fields.at("physical_total_energy_erg_corrected"), units::value(d.physicalTotalEnergy) + 5);
	near(fields.at("rsla_total_energy_erg_in"), 51);
	near(fields.at("rsla_total_energy_erg_out"), 62);
	near(fields.at("rsla_total_energy_erg_corrected"), units::value(d.rslaTotalEnergy) + 11);
	near(fields.at("physical_total_momentum_x_g_cm_s_corrected"), units::value(d.physicalTotalMomentum[0]) + 8);
	near(fields.at("rsla_total_momentum_x_g_cm_s_corrected"), units::value(d.rslaTotalMomentum[0]) + 26);
}

TEST(RadiationCouplingDiagnostics, PhysicalLengthDoesNotShrinkUnderRefinementOrMeshRotation) {
	auto c = configuration();
	c.radiation.opacity = 3;
	c.radiation.diagnosticLength = 5;
	c.radiation.lightSpeedRatio = 0.2;
	auto const gas = gasState(c);
	auto const rad = radiationState();
	auto const a = radiation::couplingDiagnostics(gas, rad, units::Length::from_value(0.5), c);
	if constexpr (ndim >= 2) c.frame.omega = units::InverseTime::from_value(100);
	auto const b = radiation::couplingDiagnostics(gas, rad, units::Length::from_value(0.25), c);
	near(a.cellOpticalDepth, 3);
	near(b.cellOpticalDepth, 1.5);
	near(a.trappingParameter, 60 / units::value(constants::c));
	near(a.trappingParameter, b.trappingParameter);
	near(a.rslaCriterion, b.rslaCriterion);
	Real const acoustic = std::sqrt((c.hydro.gamma * 10 + 4 * units::value(rad.energy()) / 9) / 2);
	near(a.rslaCriterion, (2 + acoustic) * 30 / (0.2 * units::value(constants::c)));
}

TEST(RadiationCouplingDiagnostics, DomainDefaultAndLargeValuesRemainNonfatal) {
	auto c = configuration();
	c.radiation.opacity = 3;
	c.mesh.lower = units::Length::from_value(-4);
	c.mesh.upper = units::Length::from_value(6);
	c.radiation.diagnosticLength = 0;
	auto const gas = gasState(c);
	auto const rad = radiationState();
	auto const domain = radiation::couplingDiagnostics(gas, rad, units::Length::from_value(0.5), c);
	c.radiation.diagnosticLength = 10;
	auto const explicitLength = radiation::couplingDiagnostics(gas, rad, units::Length::from_value(0.5), c);
	near(domain.rslaCriterion, explicitLength.rslaCriterion);
	c.radiation.opacity = 1e16;
	c.radiation.lightSpeedRatio = 1;
	auto const thick = radiation::couplingDiagnostics(gas, rad, units::Length::from_value(0.5), c);
	EXPECT_GT(thick.trappingParameter, 1);
	EXPECT_GT(thick.rslaCriterion, 1);
	EXPECT_TRUE(std::isfinite(thick.rslaCriterion));
	c.radiation.opacity = 0;
	auto const transparent = radiation::couplingDiagnostics(gas, rad, units::Length::from_value(0.5), c);
	EXPECT_EQ(transparent.cellOpticalDepth, 0);
	EXPECT_EQ(transparent.trappingParameter, 0);
	EXPECT_GT(transparent.rslaCriterion, 0);
}

TEST(RadiationCouplingDiagnostics, SiloFieldsRoundTripAtNonfatalLargeValues) {
	test::TemporaryDirectory directory;
	auto c = configuration();
	c.output.enabled = true;
	c.output.directory = directory.path.string();
	c.radiation.opacity = 1e16;
	c.radiation.diagnosticLength = 5;
	c.radiation.lightSpeedRatio = 0.1;
	auto block = snapshot(c);
	auto const expected = radiation::couplingDiagnostics(gasState(c), radiationState(), block.cellWidth, c);
	ASSERT_GT(expected.rslaCriterion, 1);
	{ Output output(c); output({block}, 0, diagnose({block}, c)); }
	std::unique_ptr<DBfile, decltype(&DBClose)> file(
		DBOpen((directory.path / "frame_000000.silo").c_str(), DB_HDF5, DB_READ), &DBClose);
	ASSERT_TRUE(file);
	for (auto const& [name, value] : std::map<std::string, Real>{
		{"cellOpticalDepth", expected.cellOpticalDepth},
		{"radiationTrappingParameter", expected.trappingParameter},
		{"rslaCriterion", expected.rslaCriterion}}) {
		std::unique_ptr<DBmultivar, decltype(&DBFreeMultivar)> multi(DBGetMultivar(file.get(), name.c_str()), &DBFreeMultivar);
		ASSERT_TRUE(multi);
		ASSERT_EQ(multi->nvars, 1);
		std::unique_ptr<DBquadvar, decltype(&DBFreeQuadvar)> field(DBGetQuadvar(file.get(), multi->varnames[0]), &DBFreeQuadvar);
		ASSERT_TRUE(field);
		ASSERT_EQ(field->datatype, DB_DOUBLE);
		ASSERT_EQ(field->centering, DB_ZONECENT);
		ASSERT_EQ(field->nels, int(block.layout.interiorCellCount()));
		ASSERT_EQ(field->nvals, 1);
		for (int i = 0; i < field->nels; ++i) { near(static_cast<double const*>(field->vals[0])[i], value); }
	}
}
