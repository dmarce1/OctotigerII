#include "testSupport.hpp"
#include "octotigerII/subgrid/view.hpp"
#include "octotigerII/runtime.hpp"
#include "octotigerII/simulation.hpp"
#include "octotigerII/radiation/diffusionFlux.hpp"
#include "octotigerII/radiation/opacity.hpp"
#include "../src/runtime/haloCache.hpp"
#include <algorithm>
#include <iostream>
#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/runtime_distributed/find_all_localities.hpp>
#endif

using namespace octotigerII;

#ifdef OCTOTIGERII_WITH_HPX
HaloReadStatistics haloTrafficSnapshot() { return haloReadStatistics(); }
HPX_PLAIN_ACTION(haloTrafficSnapshot, HaloTrafficSnapshotAction)
#endif

namespace {
std::vector<storage::Locality> owners() {
#ifdef OCTOTIGERII_WITH_HPX
	auto result = hpx::find_all_localities();
	result.push_back(hpx::find_here()); // A different partition can still be local.
	return result;
#else
	return {0, 0};
#endif
}

std::uint64_t remotePartitions(std::vector<storage::Locality> const& placement) {
#ifdef OCTOTIGERII_WITH_HPX
	return std::count_if(placement.begin(), placement.end(), [](auto const& owner) { return owner != hpx::find_here(); });
#else
	(void) placement;
	return 0;
#endif
}

HaloPlan plan(storage::Layout const& layout) {
	HaloPlan result;
	for (auto const& range : layout.ranges()) {
		// Request four scalar values but use only two. Traffic follows the
		// transmitted range, not the number of ghost copies or unique donors.
		result.reads.push_back({range.slice(1, 4), {{0, result.ghostCount, 1}, {3, result.ghostCount + 1, 1}}, 0});
		result.ghostCount += 2;
	}
	return result;
}

void expectCounts(HaloReadStatistics const& counts, std::uint64_t fields, std::uint64_t remote) {
	EXPECT_EQ(counts.haloCalls, 1u);
	EXPECT_EQ(counts.fieldReads, fields);
	EXPECT_EQ(counts.values, 4 * fields);
	EXPECT_EQ(counts.remoteFieldReads, remote);
	EXPECT_EQ(counts.remoteValues, 4 * remote);
	EXPECT_EQ(counts.remotePayloadBytes, 4 * remote * sizeof(Real));
}

HaloReadStatistics allLocalityHaloStatistics() {
#ifdef OCTOTIGERII_WITH_HPX
	std::vector<hpx::future<HaloReadStatistics>> pending;
	for (auto const& locality : hpx::find_all_localities()) pending.push_back(hpx::async<HaloTrafficSnapshotAction>(locality));
	HaloReadStatistics result;
	for (auto& snapshot : pending) result += snapshot.get();
	return result;
#else
	return haloReadStatistics();
#endif
}
}

TEST(HaloTraffic, CountsRequestedRangesAndTimeEndpoints) {
	auto const placement = owners();
	storage::Layout layout(std::vector<std::size_t>(placement.size(), 6), placement.size());
	storage::PartitionSet store(placement);
	storage::Field<units::Density> field(layout, store, 3);
	for (unsigned bank = 0; bank < 3; ++bank) for (auto const& range : layout.ranges()) {
		auto output = field.handle().output(range, bank);
		for (std::size_t i = 0; i < range.count; ++i) output.data()[i] = units::Density::from_value(10 * bank + i);
		field.handle().commit(range, bank, output);
	}
	auto const halo = plan(layout);
	std::vector<units::Density> ghosts;
	auto before = haloReadStatistics();
	readHalo(field.handle(), halo, 0, ghosts);
	expectCounts(haloReadStatistics() - before, placement.size(), remotePartitions(placement));
	for (std::size_t i = 0; i < placement.size(); ++i) {
		EXPECT_EQ(ghosts[2 * i], units::Density::from_value(1));
		EXPECT_EQ(ghosts[2 * i + 1], units::Density::from_value(4));
	}

	before = haloReadStatistics();
	readHalo(field.handle(), halo, 0, ghosts, {{0, 0.25, 2}});
	auto const interpolated = haloReadStatistics() - before;
	expectCounts(interpolated, 2 * placement.size(), 2 * remotePartitions(placement));
	for (std::size_t i = 0; i < placement.size(); ++i) {
		EXPECT_EQ(ghosts[2 * i], units::Density::from_value(6));
		EXPECT_EQ(ghosts[2 * i + 1], units::Density::from_value(9));
	}

	// A zero fraction must neither request nor count an unavailable endpoint.
	before = haloReadStatistics();
	readHalo(field.handle(), halo, 0, ghosts, {{0, 0, 99}});
	expectCounts(haloReadStatistics() - before, placement.size(), remotePartitions(placement));
	std::cout << "halo endpoint requests: values=" << interpolated.values
		<< " remote_values=" << interpolated.remoteValues << " remote_payload_bytes=" << interpolated.remotePayloadBytes << '\n';
}

TEST(HaloTraffic, ExpandsDerivedColumnsIntoTheirActualSourceReads) {
	using State = units::State<units::Density, units::EnergyDensity>;
	auto const placement = owners();
	storage::Layout layout(std::vector<std::size_t>(placement.size(), 6), placement.size());
	storage::PartitionSet store(placement);
	storage::ColumnFields<State> fields(layout, store, "halo.columns");
	auto handle = fields.handle();
	for (auto const& range : layout.ranges()) {
		auto output = handle.output(range, 0);
		for (std::size_t i = 0; i < range.count; ++i)
			output.put(i, State(units::Density::from_value(2), units::EnergyDensity::from_value(3)));
		handle.commit(range, 0, output);
	}
	// Shift ownership to exercise derived sources on different localities.
	auto otherPlacement = placement;
	std::rotate(otherPlacement.begin(), otherPlacement.begin() + 1, otherPlacement.end());
	storage::PartitionSet otherStore(otherPlacement);
	storage::Field<units::Density> other(layout, otherStore, 1);
	for (auto const& range : layout.ranges()) {
		auto output = other.handle().output(range, 0);
		std::fill_n(output.data(), range.count, units::Density::from_value(7));
		other.handle().commit(range, 0, output);
	}
	auto const first = std::get<0>(handle.fields);
	std::get<0>(handle.fields).sumSources = {first, other.handle()};
	std::get<0>(handle.fields).sumWeights = {2, 0.5};
	std::vector<State> ghosts;
	auto const before = haloReadStatistics();
	readHalo(handle, plan(layout), 0, ghosts);
	auto const counts = haloReadStatistics() - before;
	expectCounts(counts, 3 * placement.size(), 2 * remotePartitions(placement) + remotePartitions(otherPlacement));
	for (auto const& ghost : ghosts) {
		EXPECT_EQ(ghost.template get<0>(), units::Density::from_value(7.5));
		EXPECT_EQ(ghost.template get<1>(), units::EnergyDensity::from_value(3));
	}
	// Totals over every partition could hide assigning all derived sources to
	// the first source's locality. Isolate a range whose source owners differ.
	auto selected = plan(layout);
	selected.reads.resize(1);
	selected.ghostCount = 2;
	auto const partition = selected.reads.front().range.partition;
	auto const selectedBefore = haloReadStatistics();
	readHalo(handle, selected, 0, ghosts);
	expectCounts(haloReadStatistics() - selectedBefore, 3,
		remotePartitions({placement[partition], placement[partition], otherPlacement[partition]}));
}

TEST(HaloTraffic, AnalyticBoundariesRequestNoDonorValues) {
	// An empty plan must not touch this deliberately unallocated handle.
	storage::FieldHandle<units::Density> field;
	HaloPlan halo;
	halo.ghostCount = 2;
	std::vector<units::Density> ghosts;
	auto const before = haloReadStatistics();
	readHalo(field, halo, 0, ghosts);
	expectCounts(haloReadStatistics() - before, 0, 0);
	EXPECT_EQ(ghosts.size(), halo.ghostCount);
}

#if OCTOTIGERII_HYDRO && OCTOTIGERII_RADIATION
namespace {
template <typename System, typename Corrector>
void compareCompactMidpoint(Subgrid const& block, HaloPlan const& plan, System const& system,
	storage::ColumnHandle<typename System::State> const& fields,
	std::vector<typename System::State> const& initialGhosts, std::vector<typename System::State> const& midpointGhosts,
	units::Time time, units::Time dt, Corrector const& corrector) {
	using State = typename System::State;
	using Solver = finiteVolume::MusclHancock<System>;
	runtime_detail::InitialHalo<State> cache;
	cache.save(block, plan, initialGhosts, true);
	mesh::MeshLayout const padded(block.layout.cellsPerActiveDimension(), 2);
	auto const shell = mesh::MeshLayout(block.layout.cellsPerActiveDimension(), 1).cellCount() - block.layout.cellCount();
	ASSERT_EQ(cache.values.size(), shell);
	std::vector<State> compactGhosts;
	cache.restore(block, plan, compactGhosts);
	ASSERT_EQ(compactGhosts.size(), initialGhosts.size());
	// The compact representation contains already transformed first-layer cells;
	// outer cells and raw prolongation donors must not be transformed again.
	auto original = fields.read(block.interior, 0).get();
	PatchView<State> full(block, original, plan, initialGhosts);
	PatchView<State> compact(block, original, plan, compactGhosts);
	PatchView<State> midpoint(block, fields.read(block.interior, 1).get(), plan, midpointGhosts);
	auto predict = [&](State, mesh::Coordinates const& cell) { return midpoint.atStorage(cell); };
	typename Solver::Workspace fullWork, compactWork;
	std::vector<State> fullResult(block.interior.count), compactResult(block.interior.count);
	auto advance = [&](auto const& input, auto& work, auto& result) {
		Solver(system, finiteVolume::RotatingFrame{}, time).advanceInto(input, dt, work,
			[&](auto const& cell, State const& value) { result[block.layout.index(cell)] = value; }, predict, false, corrector);
	};
	advance(full, fullWork, fullResult);
	advance(compact, compactWork, compactResult);
	for (std::size_t i = 0; i < fullResult.size(); ++i) test::expectStateNear(compactResult[i], fullResult[i], 0);
	for (int axis = 0; axis < ndim; ++axis) {
		ASSERT_EQ(compactWork.fluxes[axis].size(), fullWork.fluxes[axis].size());
		for (std::size_t i = 0; i < fullWork.fluxes[axis].size(); ++i)
			test::expectStateNear(compactWork.fluxes[axis][i], fullWork.fluxes[axis][i], 0);
	}
	mesh::forEachCoordinate(padded.extents(), [&](auto const& cell) {
		if (padded.isInterior(cell)) return;
		bool firstLayer = true;
		for (auto coordinate : cell) firstLayer = firstLayer && coordinate > 0 && coordinate < padded.cellsPerActiveDimension() + 3;
		test::expectStateNear(compact.atStorage(cell), firstLayer ? full.atStorage(cell) : State{}, 0);
	});
	for (std::size_t i = plan.ghostCount; i < compactGhosts.size(); ++i)
		test::expectStateNear(compactGhosts[i], State{}, 0);
}

void coupledCacheCheck(bool openingKick) {
	auto config = parseConfig({"--problem.name=photon-source", "--mesh.cells=4", "--mesh.level=1", "--output.enabled=off"});
	config.massFractions.enabled = false;
	config.amr.enabled = false;
	if (openingKick) {
		config.hydro.acceleration[0] = units::Acceleration::from_value(1e10);
		config.radiation.lightSpeedRatio = 0.2;
	}
	config.validate();
	ASSERT_FALSE(config.gravityEnabled());
#ifdef OCTOTIGERII_WITH_HPX
	std::size_t const localities = hpx::find_all_localities().size();
#else
	std::size_t const localities = 1;
#endif
	CartesianTopology const topology(config, localities);
	std::uint64_t donorReads = 0, donorValues = 0, remoteReads = 0, remoteValues = 0;
	for (auto const& block : topology.blocks()) {
		auto const halo = topology.halo(block.id);
		auto const padded = mesh::MeshLayout(config.mesh.cells, 2);
		ASSERT_EQ(halo.ghostCount, padded.cellCount() - block.layout.interiorCellCount());
		for (auto const& read : halo.reads) {
			++donorReads;
			donorValues += read.range.count;
			if (read.range.partition != block.interior.partition) {
				++remoteReads;
				remoteValues += read.range.count;
			}
		}
	}
	constexpr auto gasFields = hydro::ConservedState::size();
	constexpr auto radiationFields = radiation::RadiationSystem::State::size();
	auto const fields = (openingKick ? 3 : 2) * gasFields + 2 * radiationFields;
	std::vector<Snapshot> reference;
	struct Schedule { int workers; bool stealing; };
	for (auto const schedule : {Schedule{1, false}, Schedule{2, false}, Schedule{2, true}}) {
		// HPX tasks may resume on another OS thread; gtest's thread-local
		// ScopedTrace cannot span the yielding runtime operations below.
		config.runtime.workerTasks = schedule.workers;
		config.runtime.workStealing = schedule.stealing;
		Runtime runtime(config);
		auto const dt = std::min(0.1 * runtime.stableTimestep(), units::Time::from_value(1e-5));
		auto const initial = openingKick ? diagnose(runtime.snapshots(), config) : Diagnostics{};
		for (int step = 0; step < 2; ++step) {
			auto const before = allLocalityHaloStatistics();
			runtime.advanceCoupled(dt);
			auto const counts = allLocalityHaloStatistics() - before;
			std::cout << "cached midpoint halo step=" << step << " dimensions=" << ndim << " localities=" << localities
				<< " workers=" << schedule.workers << " stealing=" << schedule.stealing << " opening_kick=" << openingKick
				<< " calls=" << counts.haloCalls << " field_reads=" << counts.fieldReads << " values=" << counts.values
				<< " remote_field_reads=" << counts.remoteFieldReads << " remote_values=" << counts.remoteValues
				<< " remote_payload_bytes=" << counts.remotePayloadBytes
				<< " initial_and_midpoint_expected_values=" << fields * donorValues
				<< " legacy_5gas_3rad_calculated_values=" << (5 * gasFields + 3 * radiationFields) * donorValues << '\n';
			// Cache reuse must survive a change of worker or execution locality.
			// A second accepted step must fetch its own fresh initial/midpoint data.
			EXPECT_EQ(counts.haloCalls, (openingKick ? 5 : 4) * runtime.size());
			EXPECT_EQ(counts.fieldReads, fields * donorReads);
			EXPECT_EQ(counts.values, fields * donorValues);
			if (!schedule.stealing) {
				EXPECT_EQ(counts.remoteFieldReads, fields * remoteReads);
				EXPECT_EQ(counts.remoteValues, fields * remoteValues);
			} else {
				// Scheduling can change the remote fraction while total requests
				// stay fixed. Enabling stealing does not guarantee a stolen task.
				EXPECT_LE(counts.remoteFieldReads, counts.fieldReads);
				EXPECT_LE(counts.remoteValues, counts.values);
			}
			EXPECT_EQ(counts.remotePayloadBytes, counts.remoteValues * sizeof(Real));
			if (openingKick) {
				// Periodic uniform forcing has no net boundary momentum flux.
				// Gas-radiation exchange cancels in Pgas + Frad/(c*chat), leaving
				// precisely the imposed mass times acceleration times interval.
				auto const snapshots = runtime.snapshots();
				auto const current = diagnose(snapshots, config);
				auto const interval = Real(step + 1) * dt;
				auto const impulseScale = initial.mass * config.hydro.acceleration[0] * interval;
				ASSERT_GT(impulseScale, units::Momentum{});
				for (int d = 0; d < ndim; ++d) {
					auto const expected = initial.mass * config.hydro.acceleration[d] * interval;
					auto const error = current.rslaTotalMomentum[d] - initial.rslaTotalMomentum[d] - expected;
					EXPECT_NEAR(Real(error / impulseScale), 0, 3e-12) << "momentum axis " << d;
				}
				// A stale pre-kick donor must not create a cell-boundary jump
				// in this spatially homogeneous evolving solution.
				ASSERT_FALSE(snapshots.empty());
				auto const& gas = snapshots.front().hydro.values().front();
				auto const& radiation = snapshots.front().radiation.values().front();
				for (auto const& block : snapshots)
					for (std::size_t i = 0; i < block.hydro.values().size(); ++i) {
						test::expectStateNear(block.hydro.values()[i], gas, 0);
						test::expectStateNear(block.radiation.values()[i], radiation, 0);
					}
			}
		}
		auto const result = runtime.snapshots();
		if (reference.empty()) reference = result;
		else {
			ASSERT_EQ(result.size(), reference.size());
			for (std::size_t b = 0; b < result.size(); ++b) {
				ASSERT_EQ(result[b].hydro.values().size(), reference[b].hydro.values().size());
				ASSERT_EQ(result[b].radiation.values().size(), reference[b].radiation.values().size());
				for (std::size_t i = 0; i < result[b].hydro.values().size(); ++i) {
					test::expectStateNear(result[b].hydro.values()[i], reference[b].hydro.values()[i], 0);
					test::expectStateNear(result[b].radiation.values()[i], reference[b].radiation.values()[i], 0);
				}
			}
		}
	}
}
}

TEST(HaloTraffic, CoupledStepReusesInitialAndMidpointHalosAcrossWorkerSchedules) {
	coupledCacheCheck(false);
}

TEST(HaloTraffic, OpeningKickRefreshesGasHaloAndPreservesCombinedMomentum) {
	coupledCacheCheck(true);
}

TEST(HaloTraffic, CompactInitialHaloMatchesFullMidpointCorrectorOnMixedLevelsAndPhysicalBoundaries) {
	// Use a problem that permits analytic faces; this fixture supplies its own
	// gas/radiation values and boundary evaluators below.
	auto config = parseConfig({"--problem.name=sod", "--mesh.cells=4", "--mesh.level=1", "--output.enabled=off"});
	config.mesh.lower = units::Length{};
	config.mesh.upper = units::Length::from_value(1);
	config.mesh.boundary = finiteVolume::BoundaryConditions::uniform(finiteVolume::BoundaryCondition::Reflecting);
	config.mesh.boundary.upper[0] = finiteVolume::BoundaryCondition::Analytic;
	config.amr.enabled = true;
	config.amr.maxLevel = 2;
	std::vector<mesh::BlockLocation> leaves;
	for (int slot = 0; slot < (1 << ndim); ++slot) {
		auto const child = mesh::BlockLocation{}.child(slot);
		if (slot == 0) for (int sub = 0; sub < (1 << ndim); ++sub) leaves.push_back(child.child(sub));
		else leaves.push_back(child);
	}
	auto const placement = owners();
	CartesianTopology const topology(config, placement.size(), leaves);
	storage::PartitionSet store(placement);
	storage::ColumnFields<hydro::ConservedState> gas(topology.storageLayout(), store, "compact-halo.gas");
	storage::ColumnFields<radiation::RadiationSystem::State> radiation(topology.storageLayout(), store, "compact-halo.radiation");
	hydro::HydroSystem const gasSystem(1.4);
	radiation::RadiationSystem const radiationSystem(Real(.2) * constants::c);
	auto const crossing = config.mesh.upper / constants::c;
	auto const time = Real(.2) * crossing;
	auto const dt = Real(.005) * crossing / Real(16 * ndim);
	auto profile = [&](mesh::PhysicalCoordinates const& position, units::Time at) {
		Real value = Real(.2) * (at / crossing);
		for (int d = 0; d < ndim; ++d) {
			Real const x = position[d] / config.mesh.upper;
			value += (d + 1) * (x + Real(.3) * x * x) / ndim;
		}
		return value;
	};
	finiteVolume::AnalyticBoundary<hydro::ConservedState> gasState = [&](auto const& position, auto at) {
		Real const value = profile(position, at);
		hydro::PrimitiveState state;
		state.density() = units::Density::from_value(2 + Real(.2) * value);
		state.pressure() = units::Pressure::from_value(2e16 * (1 + Real(.1) * value));
		for (int d = 0; d < ndim; ++d)
			state.velocity(d) = units::Velocity::from_value(1e7 * (d + 1) * (1 + Real(.1) * value) / ndim);
		return gasSystem.conservedState(state);
	};
	finiteVolume::AnalyticBoundary<radiation::RadiationSystem::State> radiationState = [&](auto const& position, auto at) {
		Real const value = profile(position, at);
		radiation::RadiationSystem::State state;
		state.energy() = units::EnergyDensity::from_value(2 + Real(.3) * value);
		for (int d = 0; d < ndim; ++d)
			state.radiativeFlux(d) = Real(.03) * (d + 1) / ndim * constants::c * state.energy();
		return state;
	};
	for (auto const& block : topology.blocks()) for (unsigned bank = 0; bank < 2; ++bank) {
		auto gasOutput = gas.handle().output(block.interior, bank);
		auto radiationOutput = radiation.handle().output(block.interior, bank);
		block.layout.forEachInterior([&](auto const& cell, std::size_t i) {
			auto const position = block.layout.cellCenter(block.lower, block.cellWidth, cell);
			gasOutput.put(i, gasState(position, time + Real(bank) * dt / 2.0));
			radiationOutput.put(i, radiationState(position, time + Real(bank) * dt / 2.0));
		});
		gas.handle().commit(block.interior, bank, gasOutput);
		radiation.handle().commit(block.interior, bank, radiationOutput);
	}
	std::size_t prolongations = 0, reflections = 0, analytic = 0;
	for (auto const& block : topology.blocks()) {
		auto const halo = topology.halo(block.id);
		prolongations += halo.prolongations.size();
		reflections += std::count_if(halo.reflectionMasks.begin(), halo.reflectionMasks.end(), [](auto mask) { return mask != 0; });
		analytic += halo.analyticGhosts.size();
		std::array<std::vector<hydro::ConservedState>, 2> gasGhosts;
		std::array<std::vector<radiation::RadiationSystem::State>, 2> radiationGhosts;
		for (unsigned bank = 0; bank < 2; ++bank) {
			readHalo(gas.handle(), halo, bank, gasGhosts[bank]);
			applyHaloBoundaries(halo, gasGhosts[bank], gasSystem, time + Real(bank) * dt / 2.0, gasState);
			readHalo(radiation.handle(), halo, bank, radiationGhosts[bank]);
			applyHaloBoundaries(halo, radiationGhosts[bank], radiationSystem, time + Real(bank) * dt / 2.0, radiationState);
		}
		compareCompactMidpoint(block, halo, gasSystem, gas.handle(), gasGhosts[0], gasGhosts[1], time, dt,
			[](auto const& flux, auto const&...) { return flux; });
		PatchView<hydro::ConservedState> material(block, gas.handle().read(block.interior, 1).get(), halo, gasGhosts[1]);
		for (Real amplification : {Real(1), Real(1e6)}) {
			// The second pass rejects both corrected candidates and exercises the
			// realizability limiter's fallback to the original face-adjacent cells.
			auto correction = [&](auto const& flux, auto const& centerLeft, auto const& centerRight,
				auto const& faceLeft, auto const& faceRight, auto const& left, auto const& right,
				int axis, units::Length width, units::Velocity speed) {
				auto const l = material.atStorage(left), r = material.atStorage(right);
				radiation::MaterialVelocity vl{}, vr{};
				for (int d = 0; d < ndim; ++d) { vl[d] = l.momentum(d) / l.density(); vr[d] = r.momentum(d) / r.density(); }
				auto const opacity = radiation::Opacity::from_value(30);
				return amplification * radiation::diffusionCorrectedFlux(radiationSystem, flux, centerLeft, centerRight,
					faceLeft, faceRight, l.density() * opacity, r.density() * opacity, vl, vr, axis, width, speed);
			};
			compareCompactMidpoint(block, halo, radiationSystem, radiation.handle(), radiationGhosts[0], radiationGhosts[1], time, dt, correction);
		}
	}
	EXPECT_GT(prolongations, 0u);
	EXPECT_GT(reflections, 0u);
	EXPECT_GT(analytic, 0u);
}
#endif
