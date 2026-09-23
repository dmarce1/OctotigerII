#include <fstream>
#include <gtest/gtest.h>
#include <limits>
#include "octotigerII/gravity/boundary.hpp"
#include "octotigerII/problems.hpp"
#include "octotigerII/runtime.hpp"
#include "octotigerII/storage/registry.hpp"
#include "octotigerII/subgrid/view.hpp"
#include "octotigerII/verification/analytic.hpp"
#include "testSupport.hpp"
#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/runtime_distributed/find_all_localities.hpp>
#endif

using namespace octotigerII;

namespace {

using Rule = physics::BoundaryCondition;
using Boundaries = physics::BoundaryConditions;

TEST(BoundaryOptions, EveryActiveFaceAndPeriodicPairAreValidated) {
	for (int axis = 0; axis < ndim; ++axis) {
		auto const prefix = std::string("--mesh.boundary.") + "xyz"[axis];
		for (auto side : {"Lower", "Upper"}) {
			EXPECT_THROW(parseConfig({"--mesh.periodic=off", prefix + side + "=periodic"}), std::invalid_argument);
			EXPECT_THROW(parseConfig({prefix + side + "=invalid"}), std::invalid_argument);
			for (auto rule : {"reflecting", "outflow", "inflow", "analytic"}) {
				auto args = std::vector<std::string>{"--mesh.periodic=off", prefix + side + "=" + rule};
				bool const available = std::string(rule) == "outflow" || std::string(rule) == "inflow" || std::string(rule) == "reflecting" ||
					(!build::gravity &&
						(std::string(rule) == "reflecting" || std::string(build::problem) == "sod" || std::string(build::problem) == "streaming"));
				if (available)
					EXPECT_NO_THROW(parseConfig(args));
				else
					EXPECT_THROW(parseConfig(args), std::invalid_argument);
			}
		}
		auto args = std::vector<std::string>{"--mesh.periodic=off", prefix + "Lower=periodic", prefix + "Upper=periodic"};
		EXPECT_TRUE(parseConfig(args).mesh.boundary.periodic(axis));
	}
	for (int axis = ndim; axis < 3; ++axis)
		EXPECT_THROW(parseConfig({std::string("--mesh.boundary.") + "xyz"[axis] + "Lower=outflow"}), std::exception);
}

TEST(BoundaryOptions, IniCliAndLegacyPrecedence) {
	test::TemporaryDirectory temporary;
	auto const path = (temporary.path / "boundaries.ini").string();
	{
		std::ofstream out(path);
		out << "[mesh]\nperiodic=off\n[mesh.boundary]\nxLower=reflecting\nxUpper=outflow\n";
	}
	auto c = parseConfig({"--config=" + path, "--mesh.boundary.xUpper=reflecting"});
	EXPECT_EQ(c.mesh.boundary.lower[0], Rule::Reflecting);
	EXPECT_EQ(c.mesh.boundary.upper[0], Rule::Reflecting);
	EXPECT_TRUE(parseConfig({"--config=" + path, "--mesh.periodic=on"}).mesh.boundary.all(Rule::Periodic));
	c = parseConfig({"--mesh.boundary.xLower=reflecting", "--mesh.periodic=off"});
	EXPECT_EQ(c.mesh.boundary.lower[0], Rule::Reflecting);
	EXPECT_EQ(c.mesh.boundary.upper[0], Rule::Outflow);
	EXPECT_THROW(parseConfig({"--mesh.periodic=on", "--mesh.boundary.xLower=outflow"}), std::invalid_argument);
	{
		std::ofstream out(path);
		out << "mesh.periodic=off\nmesh.boundary.xLower=periodic\n";
	}
	EXPECT_TRUE(parseConfig({"--config=" + path, "--mesh.boundary.xUpper=periodic"}).mesh.boundary.periodic(0));
}

TEST(BoundaryMapping, ReflectionMirrorsBothLayersAndAllCornerComponents) {
	for (int normal = 0; normal < ndim; ++normal) {
		Boundaries bc;
		bc.lower[normal] = bc.upper[normal] = Rule::Reflecting;
		bc.validate();
		for (int x : {-2, -1, 8, 9}) {
			auto cell = mesh::filledCoordinates(3);
			cell[normal] = x;
			auto const mapped = bc.map(cell, 8);
			EXPECT_EQ(mapped.source[normal], x < 0 ? -x - 1 : 15 - x);
			EXPECT_EQ(mapped.reflectionMask, 1u << normal);
			EXPECT_FALSE(mapped.analytic);
		}
	}
	auto bc = Boundaries::uniform(Rule::Reflecting);
	auto corner = bc.map(mesh::filledCoordinates(-2), 8);
	EXPECT_EQ(corner.source, mesh::filledCoordinates(1));
	EXPECT_EQ(corner.reflectionMask, (1u << ndim) - 1);
	// More than one domain width is handled without sampling another ghost.
	EXPECT_EQ(bc.map(mesh::filledCoordinates(-17), 8).source, mesh::filledCoordinates(0));
}

TEST(BoundaryMapping, AsymmetricFacesAndInvalidEnums) {
	Boundaries bc;
	bc.lower[0] = Rule::Reflecting;
	auto cell = mesh::filledCoordinates(2);
	cell[0] = -2;
	EXPECT_EQ(bc.map(cell, 8).source[0], 1);
	cell[0] = 9;
	EXPECT_EQ(bc.map(cell, 8).source[0], 7);
	EXPECT_EQ(bc.map(cell, 8).reflectionMask, 0u);
	bc.lower[0] = static_cast<Rule>(100);
	EXPECT_THROW(bc.validate(), std::invalid_argument);
}

TEST(GravityBoundaries, OutflowPeriodicAndReflectingAreSupported) {
	EXPECT_NO_THROW(gravity::validateBoundaries(Boundaries{}));
	for (int axis = 0; axis < ndim; ++axis) {
		for (bool lower : {true, false}) {
			for (auto rule : {Rule::Reflecting, Rule::Inflow, Rule::Analytic}) {
				Boundaries bc;
				(lower ? bc.lower : bc.upper)[axis] = rule;
				if (bc.contains(Rule::Analytic))
					EXPECT_THROW(gravity::validateBoundaries(bc), std::invalid_argument);
				else
					EXPECT_NO_THROW(gravity::validateBoundaries(bc));
			}
		}
		Boundaries bc;
		bc.lower[axis] = bc.upper[axis] = Rule::Periodic;
		EXPECT_NO_THROW(gravity::validateBoundaries(bc));
	}
}

TEST(AnalyticBoundaries, StreamingWrapsOnlyPeriodicAxes) {
	Config c;
	c.radiation.lightSpeedRatio = 0.25;
	c.mesh.boundary = Boundaries::uniform(Rule::Analytic);
	mesh::PhysicalCoordinates position{};
	auto const time = units::Time::from_value(0.5 / units::value(constants::c));
	auto const displacement = c.radiation.lightSpeedRatio * constants::c * time / std::sqrt(Real(ndim));
	position.fill(units::Length::from_value(0.25) + displacement);
	EXPECT_NEAR(units::value(verification::streamingState(c, position, time).radiation.energy()), 1.000001, 1e-14);
	position[0] += units::Length::from_value(1);
	EXPECT_NEAR(units::value(verification::streamingState(c, position, time).radiation.energy()), 1e-6, 1e-14);
	c.mesh.boundary.lower[0] = c.mesh.boundary.upper[0] = Rule::Periodic;
	EXPECT_NEAR(units::value(verification::streamingState(c, position, time).radiation.energy()), 1.000001, 1e-14);
	EXPECT_TRUE(verification::streamingReference(c).evaluate);
	c.mesh.boundary.upper[0] = c.mesh.boundary.lower[0] = Rule::Reflecting;
	EXPECT_FALSE(verification::streamingReference(c).evaluate);
}

TEST(AnalyticBoundaries, ProblemReferencePolicyTracksBoundaryChoice) {
	if (std::string(build::problem) == "sod") {
		auto c = parseConfig({"--mesh.boundary.xLower=analytic", "--mesh.boundary.xUpper=analytic"});
		EXPECT_TRUE(problemBoundary(c));
		EXPECT_TRUE(std::isinf(units::value(problemReference(c).validUntil)));
		c.mesh.boundary = Boundaries::periodic();
		EXPECT_FALSE(problemReference(c).evaluate);
	}
}

#if (OCTOTIGERII_HYDRO || OCTOTIGERII_RADIATION) && !OCTOTIGERII_GRAVITY
#if OCTOTIGERII_HYDRO
using Systems = ::testing::Types<hydro::HydroSystem>;
#else
using Systems = ::testing::Types<radiation::RadiationSystem>;
#endif

template <typename System>
class BoundaryTransport : public ::testing::Test {
protected:
	using State = typename System::State;
	System system = [] {
		if constexpr (std::is_same_v<System, hydro::HydroSystem>)
			return System(1.4);
		else
			return System(0.25 * constants::c);
	}();

	State state(Real scale) const {
		if constexpr (std::is_same_v<System, hydro::HydroSystem>) {
			hydro::PrimitiveState p;
			p.density() = units::Density::from_value(scale);
			p.pressure() = units::Pressure::from_value(2);
			for (int axis = 0; axis < ndim; ++axis)
				p.velocity(axis) = units::Velocity::from_value(0.1 * (axis + 1) * (int(scale) % 2 ? -1 : 1));
			return system.conservedState(p);
		} else {
			State value;
			value.energy() = units::EnergyDensity::from_value(scale);
			for (int axis = 0; axis < ndim; ++axis)
				value.radiativeFlux(axis) = (0.1 * (axis + 1) / ndim) * (int(scale) % 2 ? -1.0 : 1.0) * constants::c * value.energy();
			return value;
		}
	}
};

TYPED_TEST_SUITE(BoundaryTransport, Systems);

TYPED_TEST(BoundaryTransport, OutflowClampsOnlyInwardNormalAndInflowCopiesBothSigns) {
	using State = typename TypeParam::State;
	for (auto rule : {Rule::Outflow, Rule::Inflow})
		for (int axis = 0; axis < ndim; ++axis)
			for (bool lower : {true, false})
				for (Real scale : {1.0, 2.0}) {
					mesh::PatchData<State> patch(mesh::MeshLayout(4, 2), units::Length::from_value(0.25));
					auto const donor = this->state(scale);
					patch.layout().forEachInterior([&](auto cell, auto) { patch.atInterior(cell) = donor; });
					auto boundaries = Boundaries::uniform(Rule::Inflow);
					(lower ? boundaries.lower : boundaries.upper)[axis] = rule;
					physics::fillGhostCells(patch, boundaries, this->system);
					auto expected = donor;
					bool const inward = lower ? scale == 2 : scale == 1;
					if (rule == Rule::Outflow && inward) {
						if constexpr (std::is_same_v<TypeParam, hydro::HydroSystem>)
							expected.momentum(axis) = {};
						else
							expected.radiativeFlux(axis) = {};
					}
					for (int layer = 0; layer < 2; ++layer) {
						auto ghost = mesh::filledCoordinates(3);
						ghost[axis] = lower ? layer : 6 + layer;
						test::expectStateNear(patch.atStorage(ghost), expected, 0);
					}
				}
}

TYPED_TEST(BoundaryTransport, DistributedHaloMatchesIndependentDonorValuesOnEveryFace) {
	using State = typename TypeParam::State;
	std::vector<storage::Locality> owners;
#ifdef OCTOTIGERII_WITH_HPX
	owners = hpx::find_all_localities();
#else
	owners = {0};
#endif
	// Partitions intentionally cut across block ownership and range boundaries.
	owners.push_back(owners.front());
	storage::PartitionSet store(owners);
	auto c = parseConfig({"--mesh.periodic=off", "--mesh.cells=4", "--mesh.level=1"});
	c.mesh.lower = units::Length{};
	c.mesh.upper = units::Length::from_value(1);
	CartesianTopology topology(c, owners.size());
	storage::ColumnFields<State> fields(topology.storageLayout(), store, "boundary-test");
	for (auto const& block : topology.blocks()) {
		auto output = fields.handle().output(block.interior, 0);
		block.layout.forEachInterior([&](auto cell, auto i) {
			for (int d = 0; d < ndim; ++d)
				cell[d] += 4 * block.location.coordinates[d];
			output.put(i, this->state(1 + mesh::linearIndex(cell, mesh::filledCoordinates(8))));
		});
		fields.handle().commit(block.interior, 0, output);
	}
	// Rotate which axis is periodic, asymmetric, or analytic to cover all 6 faces.
	for (int normal = 0; normal < ndim; ++normal) {
		for (auto upperRule : {Rule::Outflow, Rule::Inflow, Rule::Analytic}) {
			c.mesh.boundary = Boundaries::uniform(Rule::Reflecting);
			c.mesh.boundary.upper[normal] = upperRule;
			if (upperRule == Rule::Inflow) c.mesh.boundary.lower[normal] = Rule::Outflow;
			if (ndim > 1) c.mesh.boundary.lower[(normal + 1) % ndim] = c.mesh.boundary.upper[(normal + 1) % ndim] = Rule::Periodic;
			physics::AnalyticBoundary<State> analytic = [&](auto const& x, auto time) {
				Real scale = 50 + units::value(time);
				for (int d = 0; d < ndim; ++d)
					scale += (d + 1) * units::value(x[d]);
				return this->state(scale);
			};
			for (auto const& block : topology.blocks()) {
				auto const plan = makeHaloPlan(c, topology.blocks(), block.id);
				mesh::MeshLayout padded(4, 2);
				for (Real seconds : {0.0, 0.7}) {
					auto const time = units::Time::from_value(seconds);
					std::vector<State> ghosts;
					readHalo(fields.handle(), plan, 0, ghosts);
					applyHaloBoundaries(plan, ghosts, this->system, time, analytic);
					mesh::forEachCoordinate(padded.extents(), [&](auto storageCell) {
						if (padded.isInterior(storageCell)) return;
						mesh::Coordinates donor{};
						mesh::PhysicalCoordinates position{};
						unsigned reflections = 0;
						unsigned outflowLower = 0, outflowUpper = 0;
						bool prescribed = false;
						for (int d = 0; d < ndim; ++d) {
							int const x = 4 * block.location.coordinates[d] + storageCell[d] - 2;
							position[d] = c.mesh.lower + (Real(x) + 0.5) * block.cellWidth;
							donor[d] = x;
							if (x >= 0 && x < 8) continue;
							auto const rule = x < 0 ? c.mesh.boundary.lower[d] : c.mesh.boundary.upper[d];
							if (rule == Rule::Periodic) donor[d] = (x + 8) % 8;
							if (rule == Rule::Outflow || rule == Rule::Inflow) donor[d] = std::clamp(x, 0, 7);
							if (rule == Rule::Outflow) (x < 0 ? outflowLower : outflowUpper) |= 1u << d;
							if (rule == Rule::Reflecting) {
								donor[d] = x < 0 ? -1 - x : 15 - x;
								reflections |= 1u << d;
							}
							if (rule == Rule::Analytic) prescribed = true;
						}
						State expected;
						if (prescribed)
							expected = analytic(position, time);
						else {
							expected = this->state(1 + mesh::linearIndex(donor, mesh::filledCoordinates(8)));
							for (int d = 0; d < ndim; ++d) {
								auto& normal = [&]() -> auto& {
									if constexpr (std::is_same_v<TypeParam, hydro::HydroSystem>)
										return expected.momentum(d);
									else
										return expected.radiativeFlux(d);
								}();
								using Quantity = std::remove_reference_t<decltype(normal)>;
								if (((outflowLower & (1u << d)) && normal > Quantity{}) || ((outflowUpper & (1u << d)) && normal < Quantity{})) normal = {};
								if (reflections & (1u << d)) expected = this->system.reflected(expected, d);
							}
						}
						test::expectStateNear(ghosts.at(plan.ghostIndices.at(padded.index(storageCell))), expected, 0);
					});
				}
			}
		}
	}
}

TYPED_TEST(BoundaryTransport, AnalyticPatchUsesCurrentTimeAndRejectsInvalidData) {
	using State = typename TypeParam::State;
	mesh::PatchData<State> patch(mesh::MeshLayout(4, 2), units::Length::from_value(0.25));
	patch.layout().forEachInterior([&](auto cell, auto) { patch.atInterior(cell) = this->state(1); });
	auto const boundaries = Boundaries::uniform(Rule::Analytic);
	EXPECT_THROW(physics::fillGhostCells(patch, boundaries, this->system), std::invalid_argument);
	physics::AnalyticBoundary<State> invalid = [&](auto const&, auto) { return Real(-1) * this->state(1); };
	EXPECT_THROW(physics::fillGhostCells(patch, boundaries, this->system, invalid), std::runtime_error);
	std::vector<units::Time> times;
	physics::AnalyticBoundary<State> evaluator = [&](auto const&, auto time) {
		times.push_back(time);
		return this->state(1);
	};
	physics::MusclHancock<TypeParam> solver(this->system);
	auto const dt = solver.stableTimestep(patch, 0.25);
	solver.advance(patch, dt, boundaries, evaluator);
	EXPECT_EQ(times.front(), units::Time{});
	EXPECT_EQ(times.back(), dt);
	EXPECT_EQ(std::count(times.begin(), times.end(), dt), times.size() / 2);
}

TYPED_TEST(BoundaryTransport, ReflectingBoxConservesMassAndEnergy) {
	using State = typename TypeParam::State;
	mesh::PatchData<State> patch(mesh::MeshLayout(8, 2), units::Length::from_value(0.25));
	patch.layout().forEachInterior([&](auto cell, auto) { patch.atInterior(cell) = this->state(1 + 0.02 * cell[0]); });
	auto total = [&] {
		State sum;
		patch.layout().forEachInterior([&](auto cell, auto) { sum += patch.atInterior(cell); });
		return sum;
	};
	auto const before = total();
	physics::MusclHancock<TypeParam> solver(this->system);
	for (int i = 0; i < 5; ++i)
		solver.advance(patch, solver.stableTimestep(patch, 0.2), Boundaries::uniform(Rule::Reflecting));
	auto const after = total();
	if constexpr (std::is_same_v<TypeParam, hydro::HydroSystem>) {
		EXPECT_NEAR(Real(after.density() / before.density()), 1, 2e-13);
		EXPECT_NEAR(Real(after.totalEnergy() / before.totalEnergy()), 1, 2e-13);
	} else
		EXPECT_NEAR(Real(after.energy() / before.energy()), 1, 2e-13);
}

TEST(BoundaryRuntime, MixedAndAnalyticBoundariesAreIndependentOfDecomposition) {
#if OCTOTIGERII_HYDRO
	using State = hydro::ConservedState;
#else
	using State = radiation::RadiationSystem::State;
#endif
	// Run through the real executor and work-stealing path, including changing stage times.
	for (auto rule : {Rule::Reflecting, Rule::Inflow, Rule::Outflow, Rule::Analytic}) {
		auto fine = parseConfig({"--mesh.periodic=off", "--mesh.cells=4", "--mesh.level=1", "--output.enabled=off"});
		if (rule == Rule::Analytic && !problemBoundary(fine)) continue;
		fine.mesh.boundary.lower[0] = rule;
		if (ndim > 1) fine.mesh.boundary.lower[1] = fine.mesh.boundary.upper[1] = Rule::Periodic;
		fine.runtime.workerTasks = 2;
		fine.verification.analytic = "off";
		auto single = fine;
		single.mesh.cells = 8;
		single.mesh.level = 0;
		Runtime a(fine), b(single);
		for (int step = 0; step < 3; ++step) {
			auto const dt = std::min(a.stableTimestep(), b.stableTimestep());
			a.advance(dt);
			b.advance(dt);
		}
		auto flatten = [](Runtime const& runtime) {
			std::vector<State> states(mesh::MeshLayout(8).interiorCellCount());
			for (auto const& snapshot : runtime.snapshots()) {
				snapshot.layout.forEachInterior([&](auto cell, auto i) {
					for (int d = 0; d < ndim; ++d)
						cell[d] += snapshot.location.coordinates[d] * snapshot.layout.cellsPerActiveDimension();
#if OCTOTIGERII_HYDRO
					states[mesh::linearIndex(cell, mesh::filledCoordinates(8))] = snapshot.hydro.values()[i];
#else
					states[mesh::linearIndex(cell, mesh::filledCoordinates(8))] = snapshot.radiation.values()[i];
#endif
				});
			}
			return states;
		};
		auto actual = flatten(a), expected = flatten(b);
		for (std::size_t i = 0; i < actual.size(); ++i)
			test::expectStateNear(actual[i], expected[i], 3e-12);
	}
}
#endif

}	 // namespace
