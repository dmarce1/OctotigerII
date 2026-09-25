#include "testSupport.hpp"
#include "octotigerII/composition/transport.hpp"
#include "octotigerII/storage/registry.hpp"
#include "octotigerII/runtime.hpp"
#include "octotigerII/simulation.hpp"
#include "octotigerII/amr/hierarchy.hpp"
#include <cctype>
#include <fstream>
#include <set>
using namespace octotigerII;

namespace {
Config config() {
	return parseConfig({"--problem.name=sod", "--mesh.cells=4", "--mesh.level=1", "--output.enabled=off",
		"--massFractions.enabled=on", "--massFractions.species=fuel:0.7:He=70%,O=30%;ash:0.3:A=56,Z=26;dye:0.2:A=0,Z=0"});
}
std::vector<units::Mass> masses(std::vector<Snapshot> const& blocks) {
	std::vector<units::Mass> sums(blocks.front().species.size());
	for (auto const& b : blocks) for (std::size_t s = 0; s < sums.size(); ++s)
		for (auto q : b.species[s].values()) sums[s] += q * b.layout.cellMeasure(b.cellWidth);
	return sums;
}
void check(std::vector<Snapshot> const& blocks, Config const& c) {
	for (auto const& b : blocks) for (std::size_t i = 0; i < b.hydro.values().size(); ++i) {
		units::Density rho{};
		for (std::size_t s = 0; s < b.species.size(); ++s) {
			auto q = b.species[s].values()[i];
			EXPECT_GE(q, units::Density{});
			if (!c.massFractions.species[s].tracer()) rho += q;
		}
		EXPECT_EQ(rho, b.hydro.values()[i].density());
		for (std::size_t s = 0; s < b.species.size(); ++s)
			EXPECT_NEAR(Real(b.species[s].values()[i] / rho), c.massFractions.species[s].initialFraction, 5e-13);
	}
}
}
TEST(Species, AllElementsAliasesAndMixtureMoments) {
	auto const& table = composition::periodicTable();
	ASSERT_EQ(table.size(), 118u);
	for (int z = 1; z <= 118; ++z) {
		auto const& e = table[z - 1];
		EXPECT_EQ(e.atomicNumber, z);
		std::string symbol(e.symbol), name(e.name);
		for (char& c : symbol) c = std::tolower(c);
		for (char& c : name) c = std::toupper(c);
		EXPECT_EQ(composition::element(symbol).atomicNumber, z);
		EXPECT_EQ(composition::element(name).atomicNumber, z);
		EXPECT_GT(e.atomicMass, Real(z));
	}
	EXPECT_EQ(composition::element("AlUmInUm").atomicNumber, 13);
	EXPECT_EQ(composition::element("CESIUM").atomicNumber, 55);
	EXPECT_EQ(composition::element("sulphur").atomicNumber, 16);
	auto c = config();
	auto const& s = c.massFractions.species[0];
	Real inverse = 0.7 / composition::element("He").atomicMass + 0.3 / composition::element("O").atomicMass;
	EXPECT_NEAR(1 / s.atomicMass, inverse, 1e-15);
	EXPECT_NEAR(s.atomicNumber / s.atomicMass, 0.7 * 2 / composition::element("He").atomicMass + 0.3 * 8 / composition::element("O").atomicMass, 1e-15);
	EXPECT_TRUE(c.massFractions.species.back().tracer());
}
TEST(Species, InvalidDefinitionsFail) {
	for (auto const& text : {"a:1:A=0,Z=2", "a:1:A=4", "a:1:He=20%,O=30%", "a:1:He=50%,helium=50%",
		"a:1:Unobtainium", "a:0.5:He", "a:1:A=nan,Z=2", "a:1:A=0,Z=0", "a:0.5:He;A:0.5:O", "a:-1:He;b:2:O", "a:1:A=4,Z=2,O=0"})
		EXPECT_THROW(composition::parseSpecies(text), std::exception) << text;
}
TEST(Species, HydroDensityHasNoAllocatedColumn) {
	auto c = config();
	storage::Layout layout({4}, 1);
	std::vector<storage::Locality> localities;
#ifdef OCTOTIGERII_WITH_HPX
	localities = {hpx::find_here()};
#else
	localities = {0};
#endif
	FieldRepository repository(c, layout, localities);
	auto const& fields = repository.directory();
	auto const& density = std::get<0>(fields.hydro.fields);
	EXPECT_EQ(density.id, 0u);
	ASSERT_EQ(density.sumSources.size(), 2u);
	for (std::size_t s = 0; s < fields.species.size(); ++s) {
		auto out = fields.species[s].output(layout.ranges()[0], 0);
		std::fill_n(out.data(), 4, units::Density::from_value(s == 2 ? 100 : s + 1));
		fields.species[s].commit(layout.ranges()[0], 0, out);
	}
	auto values = density.read(layout.ranges()[0], 0).get();
	for (std::size_t i = 0; i < 4; ++i) EXPECT_EQ(units::value(values.data()[i]), 3);
}
TEST(Species, NonuniformContactFluxesConserveEverySpecies) {
	auto c = config();
	mesh::MeshLayout layout(16);
	std::vector<std::vector<units::Density>> old(3, std::vector<units::Density>(layout.interiorCellCount()));
	layout.forEachInterior([&](auto const& cell, std::size_t i) {
		Real x = cell[0] < 8 ? 0.9 : 0.1;
		old[0][i] = units::Density::from_value(x);
		old[1][i] = units::Density::from_value(1 - x);
		old[2][i] = units::Density::from_value(cell[0] < 8 ? 1 : 0);
	});
	auto read = [&](std::size_t s, auto cell) {
		for (auto& x : cell) x = (x + 16) % 16;
		return old[s][layout.index(cell)];
	};
	auto flux = composition::fluxes(c.massFractions, layout, read,
		[](int d, auto const&) { return units::MassFlux::from_value(d == 0 ? 0.3 : 0); });
	std::array<Real, 3> before{}, after{};
	layout.forEachInterior([&](auto const& cell, std::size_t i) {
		Real rho = 0;
		for (std::size_t s = 0; s < 3; ++s) {
			auto next = composition::update(old[s][i], flux[s], layout, cell, units::TimePerLength::from_value(0.4));
			before[s] += units::value(old[s][i]); after[s] += units::value(next);
			EXPECT_GE(next, units::Density{});
			if (s < 2) rho += units::value(next);
		}
		EXPECT_NEAR(rho, 1, 1e-15);
	});
	for (int s = 0; s < 3; ++s) EXPECT_NEAR(before[s], after[s], 1e-11);
}
TEST(Species, RuntimeAndAMRPreserveCompositionAndMaterialMass) {
	for (bool adaptive : {false, true}) {
		auto c = config(); c.mesh.boundary = physics::BoundaryConditions::periodic();
		c.amr.enabled = adaptive; c.amr.maxLevel = 2; c.amr.shadowTolerance = 0; c.amr.bufferCells = 0;
		refinement::Criteria criteria{[](refinement::CellView const& cell) { return cell.center[0] < units::Length::from_value(0.2) && cell.level < 2 ? Real(2) : Real(0); }};
		Runtime runtime(c, criteria);
		auto before = masses(runtime.snapshots());
		for (int step = 0; step < 4; ++step) runtime.advance(0.2 * runtime.stableTimestep());
		if (adaptive) runtime.regrid(units::Time{}, true);
		auto blocks = runtime.snapshots(); check(blocks, c);
		auto after = masses(blocks);
		for (std::size_t s = 0; s < before.size(); ++s) EXPECT_NEAR(units::value(before[s]), units::value(after[s]), 2e-12);
	}
}

TEST(Species, IniMixturesAndCommandLineOverride) {
	test::TemporaryDirectory directory;
	auto path = directory.path / "species.ini";
	{ std::ofstream file(path); file << "[massFractions]\nenabled=on\nspecies=gas:1:helium=70%,O=30%;dye:2:A=0,Z=0\n"; }
	auto c = parseConfig({"--problem.name=sod", "--config=" + path.string()});
	ASSERT_EQ(c.massFractions.species.size(), 2u);
	EXPECT_EQ(c.massFractions.species[0].mixture.size(), 2u);
	EXPECT_TRUE(c.massFractions.species[1].tracer());
	c = parseConfig({"--problem.name=sod", "--config=" + path.string(), "--massFractions.species=iron:1:A=56,Z=26"});
	ASSERT_EQ(c.massFractions.species.size(), 1u);
	EXPECT_EQ(c.massFractions.species[0].atomicMass, 56);
}
TEST(Species, NonuniformFractionsSurviveConservativeSplitAndMerge) {
	auto c = config();
	auto parent = initialSnapshot(c, {0, {}});
	parent.layout.forEachInterior([&](auto const& cell, std::size_t i) {
		auto rho = parent.hydro.values()[i].density();
		Real fraction = cell[0] < 2 ? 0.9 : 0.1;
		parent.species[0].values()[i] = fraction * rho;
		parent.species[1].values()[i] = (1 - fraction) * rho;
		parent.species[2].values()[i] = cell[0] < 2 ? rho : units::Density{};
	});
	amr::Hierarchy source(c, {parent});
	std::vector<Snapshot> children;
	for (int slot = 0; slot < (1 << ndim); ++slot) children.push_back(source.transfer(parent.location.child(slot)));
	auto before = masses({parent}), after = masses(children);
	for (std::size_t s = 0; s < before.size(); ++s) EXPECT_NEAR(units::value(before[s]), units::value(after[s]), 2e-13);
	auto merged = amr::Hierarchy(c, children).transfer(parent.location);
	for (std::size_t s = 0; s < parent.species.size(); ++s)
		for (std::size_t i = 0; i < parent.hydro.values().size(); ++i)
			EXPECT_NEAR(units::value(merged.species[s].values()[i]), units::value(parent.species[s].values()[i]), 2e-13);
}
