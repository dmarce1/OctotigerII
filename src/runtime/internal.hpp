/** @file
 * @brief Private runtime state and stage declarations shared by the implementation.
 * Keep public interfaces in octotigerII/runtime.hpp. Types sent in HPX actions
 * use a named namespace so every translation unit sees the same wire types.
 */
#pragma once
#include "octotigerII/runtime.hpp"
#include "octotigerII/composition/transport.hpp"
#include <algorithm>
#include <exception>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <sstream>
#include "octotigerII/amr/hierarchy.hpp"
#include "octotigerII/problems.hpp"
#include "octotigerII/profiling.hpp"
#include "octotigerII/radiation/matterCoupling.hpp"
#include "octotigerII/radiation/diffusionFlux.hpp"
#include "octotigerII/storage/registry.hpp"
#include "octotigerII/subgrid/view.hpp"
#include "octotigerII/verification/analytic.hpp"
#include "haloCache.hpp"
#if OCTOTIGERII_GRAVITY
#include "octotigerII/gravity/fieldSolver.hpp"
#include "octotigerII/gravity/fluxWork.hpp"
#include "octotigerII/gravity/rotationWork.hpp"
#endif

#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/include/components.hpp>
#include <hpx/serialization/map.hpp>
#include <hpx/runtime_distributed/find_all_localities.hpp>
#include <hpx/runtime_local/get_os_thread_count.hpp>
#include <hpx/synchronization/mutex.hpp>
#endif

namespace octotigerII {
namespace runtime_detail {


// Drain every task even on failure: input storage cannot be recycled while
// another task or transfer still references it.
#ifdef OCTOTIGERII_WITH_HPX
	template <typename T>
	std::vector<T> collect(std::vector<hpx::future<T>>& pending) {
		std::exception_ptr error;
		std::vector<T> results;
		for (auto& future : pending) {
			try {
				results.push_back(future.get());
			} catch (...) {
				if (!error) error = std::current_exception();
			}
		}
		if (error) std::rethrow_exception(error);
		return results;
	}

	inline void finish(std::vector<hpx::future<void>>& pending) {
		std::exception_ptr error;
		for (auto& future : pending) {
			try {
				future.get();
			} catch (...) {
				if (!error) error = std::current_exception();
			}
		}
		if (error) std::rethrow_exception(error);
	}
#endif

	template <typename State>
	void exportFields(storage::ColumnHandle<State> const& field, storage::Range range, unsigned bank, mesh::PatchData<State>& patch) {
		auto input = field.read(range, bank).get();
		for (std::size_t i = 0; i < range.count; ++i)
			patch.values()[i] = input.at(i);
	}

	template <typename State>
	void initializeFields(storage::ColumnHandle<State> const& field, storage::Range range, mesh::PatchData<State> const& patch) {
		auto output = field.output(range, 0);
		for (std::size_t i = 0; i < range.count; ++i)
			output.put(i, patch.values()[i]);
		field.commit(range, 0, output);
	}

	template <typename State>
	void copyFields(storage::ColumnHandle<State> const& field, storage::Range range, unsigned bank, unsigned destination = ~0u) {
		if (destination == ~0u) destination = bank ^ 1;
		auto input = field.read(range, bank).get();
		auto output = field.output(range, destination);
		for (std::size_t i = 0; i < range.count; ++i)
			output.put(i, input.at(i));
		field.commit(range, destination, output);
	}

	template <typename T>
	void copyScalar(storage::FieldHandle<T> const& field, storage::Range range, unsigned bank, unsigned destination = ~0u) {
		if (destination == ~0u) destination = bank ^ 1;
		auto input = field.read(range, bank).get();
		auto output = field.output(range, destination);
		std::copy_n(input.data(), range.count, output.data());
		field.commit(range, destination, output);
	}
	inline void initializeSpecies(FieldDirectory const& fields, Subgrid const& block, Snapshot const& snapshot) {
		for (std::size_t s = 0; s < fields.species.size(); ++s) {
			auto output = fields.species[s].output(block.interior, 0);
			std::copy_n(snapshot.species.at(s).values().data(), block.interior.count, output.data());
			fields.species[s].commit(block.interior, 0, output);
		}
	}

	// The convex remainder is O(dt²). It keeps a cold, accelerating midpoint
	// physical without adding energy to the accepted conservative update.
	inline void predictorKineticRemainder(hydro::ConservedState const& initial, hydro::ConservedState& predicted) {
		if (!(initial.density() > units::Density{}) || !(predicted.density() > units::Density{})) return;
		for (int d = 0; d < ndim; ++d) {
			auto const residual = predicted.momentum(d) - initial.momentum(d) * (predicted.density() / initial.density());
			predicted.totalEnergy() += 0.5 * residual * residual / predicted.density();
		}
	}

	// Isolate the radiation source impulse from its accepted half-step
	// transport drive. Stored radiation fluxes and gas momentum are inertial.
	inline std::array<units::MomentumDensity, ndim> radiationMidpointImpulse(
		radiation::RadiationSystem::State const& initial, radiation::RadiationSystem::State const& midpoint,
		radiation::RadiationSystem::State const& transport, Real halfFraction, Real lightSpeedRatio) {
		std::array<units::MomentumDensity, ndim> result;
		for (int d = 0; d < ndim; ++d)
			result[d] = -(midpoint.radiativeFlux(d) - initial.radiativeFlux(d)
				- halfFraction * transport.radiativeFlux(d)) / (lightSpeedRatio * constants::c * constants::c);
		return result;
	}

	// This cross work becomes kinetic energy when the paired source impulse
	// is applied later. Until then it belongs to the raw drive's thermal
	// energy, including when thermodynamics selects the entropy auxiliary.
	inline void addRadiationForceWork(hydro::ConservedState& state, units::EnergyDensity work,
		hydro::HydroSystem const& system, bool dualEnergy) {
		if (work == units::EnergyDensity{}) return;
		if (dualEnergy) state.auxiliary() = system.auxiliaryFromInternalEnergy(state, system.internalEnergy(state) + work);
		state.totalEnergy() += work;
	}

	template <typename System>
	class SolverWorkspace {
	public:
		using State = typename System::State;
		std::vector<State> ghosts;
		std::vector<State> midpointGhosts, rawInitial, rateGhosts;
		std::vector<hydro::ConservedState> materialGhosts;
		using Solver = finiteVolume::MusclHancock<System>;
		typename Solver::Workspace work;

		void advance(Subgrid const& block, storage::ColumnHandle<State> const& fields, HaloPlan const& plan, System const& system, unsigned bank,
			units::Time dt, units::Time time, finiteVolume::AnalyticBoundary<State> const& analytic, storage::ColumnHandle<typename System::Flux> const& fluxFields,
			bool amr, std::vector<HaloTime> const& times = {},
			storage::ColumnHandle<hydro::ConservedState> const& increment = {}, units::Time referenceStep = {}, finiteVolume::RotatingFrame const& frame = finiteVolume::RotatingFrame{},
			storage::ColumnHandle<State> const& midpoint = {}, Config const* couplingConfig = nullptr,
			storage::ColumnHandle<hydro::ConservedState> const& material = {}, unsigned materialBank = 0,
			finiteVolume::AnalyticBoundary<hydro::ConservedState> const& materialBoundary = {}, units::Time limiterInterval = {},
			ProblemRadiationMaterial const& radiationMaterial = {}, InitialHalo<State> const* cachedInitial = nullptr,
			bool captureRawInitial = false, std::vector<hydro::ConservedState> const* sharedMaterialGhosts = nullptr) {
			auto interior = fields.read(block.interior, bank).get();
			if (cachedInitial) cachedInitial->restore(block, plan, ghosts);
			else readHalo(fields, plan, bank, ghosts, times);
			if (captureRawInitial) rawInitial = ghosts;
			midpointGhosts.clear();
			if constexpr (std::is_same_v<System, hydro::HydroSystem>) {
				if (!std::get<0>(midpoint.fields).id && referenceStep > units::Time{} && dt > units::Time{}) {
					if (cachedInitial && cachedInitial->compact)
						throw std::logic_error("Gravity rate prediction requires a raw initial halo");
					readHalo(increment, plan, 0, rateGhosts);
					midpointGhosts = ghosts;
					for (std::size_t i = 0; i < midpointGhosts.size(); ++i) {
						midpointGhosts[i] += (dt / (2.0 * referenceStep)) * rateGhosts[i];
						predictorKineticRemainder(ghosts[i], midpointGhosts[i]);
					}
					// Prolong the midpoint donors themselves: prolongation of a rate
					// alone is not the derivative of the nonlinear spatial limiter.
					applyHaloBoundaries(plan, midpointGhosts, system, time + dt / 2.0, analytic, frame);
				}
			}
			if (!cachedInitial || !cachedInitial->compact)
				applyHaloBoundaries(plan, ghosts, system, time, analytic, frame);
			PatchView<State> input(block, std::move(interior), plan, ghosts);
			auto output = fields.output(block.interior, bank ^ 1);
			{
				profiling::Region profile(std::is_same_v<System, hydro::HydroSystem> ? "hydro.advance" : "radiation.advance");
				auto writer = [&](mesh::Coordinates const& cell, State const& state) { output.put(block.layout.index(cell), state); };
				auto advanceWith = [&](auto&& predictor, bool hancock) {
					if constexpr (build::hydro && std::is_same_v<System, radiation::RadiationSystem>) {
						if (couplingConfig && radiationMaterial) {
							auto gasInterior = material.read(block.interior, materialBank).get();
							hydro::HydroSystem const gasSystem(*couplingConfig);
							auto const* gasGhosts = sharedMaterialGhosts;
							if (!gasGhosts) {
								readHalo(material, plan, materialBank, materialGhosts,
									std::get<0>(midpoint.fields).id && dt > units::Time{} ? std::vector<HaloTime>{} : times);
								applyHaloBoundaries(plan, materialGhosts, gasSystem, time + dt / 2.0, materialBoundary, frame);
								gasGhosts = &materialGhosts;
							}
							PatchView<hydro::ConservedState> gasInput(block, std::move(gasInterior), plan, *gasGhosts);
							Solver(system, frame, time).advanceInto(input, dt, work, writer, predictor, hancock,
								[&](auto const& flux, auto const& centerLeft, auto const& centerRight, auto const& faceLeft, auto const& faceRight,
									auto const& left, auto const& right, int axis, units::Length width, units::Velocity speed) {
									auto const l = frame.toGridState(gasInput.atStorage(left), time + dt / 2.0);
									auto const r = frame.toGridState(gasInput.atStorage(right), time + dt / 2.0);
									radiation::MaterialVelocity vl{}, vr{};
									for (int d = 0; d < ndim; ++d) { vl[d] = l.momentum(d) / l.density(); vr[d] = r.momentum(d) / r.density(); }
									auto extinctionAt = [&](mesh::Coordinates const& cell, hydro::ConservedState const& state) {
										auto position = block.lower;
										for (int d = 0; d < ndim; ++d)
											position[d] += (cell[d] - input.layout().ghostWidth() + Real(0.5)) * block.cellWidth;
										auto const opacity = checkedRadiationMaterial(radiationMaterial, position, time + dt / 2.0).opacity;
										return state.density() * radiation::opacityLaw(*couplingConfig, opacity).evaluate(state, gasSystem).fluxExtinction;
									};
									return radiation::diffusionCorrectedFlux(system, flux, centerLeft, centerRight, faceLeft, faceRight,
										extinctionAt(left, l), extinctionAt(right, r), vl, vr, axis, width, speed);
								}, limiterInterval);
							return;
						}
					}
					Solver(system, frame, time).advanceInto(input, dt, work, writer, predictor, hancock,
						[](auto const& flux, auto const&...) { return flux; }, limiterInterval);
				};
				if (std::get<0>(midpoint.fields).id && dt > units::Time{}) {
					auto centers = midpoint.read(block.interior, 1).get();
					readHalo(midpoint, plan, 1, midpointGhosts);
					applyHaloBoundaries(plan, midpointGhosts, system, time + dt / 2.0, analytic, frame);
					advanceWith([&](State, mesh::Coordinates const& cell) {
						return input.layout().isInterior(cell) ? centers.at(block.layout.index(input.layout().interiorCoordinates(cell))) :
							midpointGhosts.at(plan.ghostIndices.at(input.layout().index(cell)));
					}, false);
				} else if constexpr (std::is_same_v<System, hydro::HydroSystem>) {
					if (referenceStep > units::Time{} && dt > units::Time{}) {
						auto rates = increment.read(block.interior, 0).get();
						Solver(system, frame, time).advanceInto(input, dt, work, writer, [&](State value, mesh::Coordinates const& cell) {
							auto const initial = value;
							if (input.layout().isInterior(cell)) {
								value += (dt / (2.0 * referenceStep)) * rates.at(block.layout.index(input.layout().interiorCoordinates(cell)));
								predictorKineticRemainder(initial, value);
							} else value = midpointGhosts.at(plan.ghostIndices.at(input.layout().index(cell)));
							if (!system.admissible(value)) throw std::runtime_error("Gravity midpoint predictor is inadmissible");
							return value;
						}, false);
					} else advanceWith([](State value, mesh::Coordinates const&) { return value; }, true);
					if (referenceStep > units::Time{} && dt == units::Time{}) {
						auto rate = increment.output(block.interior, 0);
						block.layout.forEachInterior([&](auto const& cell, std::size_t i) {
							typename System::Flux rhs{};
							for (int d = 0; d < ndim; ++d) {
								auto upper = cell; ++upper[d];
								rhs += work.fluxes[d][block.layout.faceIndex(d, cell)] - work.fluxes[d][block.layout.faceIndex(d, upper)];
							}
							rate.put(i, (referenceStep / block.cellWidth) * rhs);
						});
						increment.commit(block.interior, 0, rate);
					}
				} else advanceWith([](State value, mesh::Coordinates const&) { return value; }, true);
			}
			fields.commit(block.interior, bank ^ 1, output);
			if (amr) {
				auto flux = fluxFields.output(block.boundaryFlux, 0);
				int const n = block.layout.cellsPerActiveDimension();
				for (int axis = 0; axis < ndim; ++axis)
					for (bool upper : {false, true}) {
						auto extents = mesh::filledCoordinates(n);
						extents[axis] = 1;
						mesh::forEachCoordinate(extents, [&](auto face) {
							face[axis] = upper ? n : 0;
							flux.put(boundaryFluxIndex(n, axis, upper, face), work.fluxes[axis][block.layout.faceIndex(axis, face)]);
						});
					}
				fluxFields.commit(block.boundaryFlux, 0, flux);
			}
		}

		BoundaryTransport boundaryTransport(Subgrid const& block, finiteVolume::BoundaryConditions const& boundaries, units::Time dt) const {
			BoundaryTransport result;
			int const n = block.layout.cellsPerActiveDimension();
			// cellMeasure also supplies the unit transverse area in 1D/2D.
			auto const measure = dt * block.layout.cellMeasure(block.cellWidth) / block.cellWidth;
			auto add = [](auto q, auto& inward, auto& outward) {
				if (q < decltype(q){}) inward -= q;
				else outward += q;
			};
			for (int axis = 0; axis < ndim; ++axis) {
				if (boundaries.periodic(axis)) continue;
				for (bool upper : {false, true}) {
					if (block.location.coordinates[axis] != (upper ? (1 << block.location.level) - 1 : 0)) continue;
					auto extents = mesh::filledCoordinates(n);
					extents[axis] = 1;
					mesh::forEachCoordinate(extents, [&](auto face) {
						face[axis] = upper ? n : 0;
						auto const q = (upper ? measure : -measure) * work.fluxes[axis][block.layout.faceIndex(axis, face)];
						if constexpr (std::is_same_v<System, hydro::HydroSystem>) {
							add(q.template get<0>(), result.inward.mass, result.outward.mass);
							add(q.totalEnergy(), result.inward.gasEnergy, result.outward.gasEnergy);
							for (int d = 0; d < ndim; ++d) add(q.momentum(d), result.inward.momentum[d], result.outward.momentum[d]);
						} else {
							add(q.template get<0>(), result.inward.radiationEnergy, result.outward.radiationEnergy);
							for (int d = 0; d < ndim; ++d) add(q.radiativeFlux(d), result.inward.radiationFlux[d], result.outward.radiationFlux[d]);
						}
					});
				}
			}
			return result;
		}
	};

	struct RadiationStepFields {
		// Bank zero is the original state; bank one is the source prediction.
		storage::ColumnHandle<hydro::ConservedState> gas;
		storage::ColumnHandle<radiation::RadiationSystem::State> radiation;
		storage::ColumnHandle<hydro::ConservedState> gasRate;
		storage::ColumnHandle<radiation::RadiationSystem::State> radiationRate;
		storage::ColumnHandle<gravity::State> gravity;
		std::vector<storage::FieldHandle<units::Density>> species;
		units::Time referenceStep{};
		units::Time limiterInterval{};
		bool gravityInRate = false;
		bool cacheInitialGasInvalid = false;
		template <typename Archive> void serialize(Archive& ar, unsigned) { ar & gas & radiation & gasRate & radiationRate & gravity & species & referenceStep & limiterInterval & gravityInRate & cacheInitialGasInvalid; }
	};

	// All source predictor fields use bank zero. increment stores the RHS integrated
	// over referenceStep, avoiding a dimensionless or untyped rate register.
	struct SourcePredictor {
		storage::ColumnHandle<hydro::ConservedState> increment;
		units::Time referenceStep{};
		storage::FieldHandle<units::EnergyDensity> rotationWork{};
		RadiationStepFields radiation;
		SourcePredictor(storage::ColumnHandle<hydro::ConservedState> value = {}, units::Time step = {})
		  : increment(std::move(value)), referenceStep(step) {}
		template <typename Archive> void serialize(Archive& ar, unsigned) { ar & increment & referenceStep & rotationWork & radiation; }
	};

	class Workspace {
	public:
		SolverWorkspace<hydro::HydroSystem> hydro;
		SolverWorkspace<radiation::RadiationSystem> radiation;
	};

	class PhaseResult {
	public:
		SchedulingStatistics tasks;
		BoundaryTransport boundary;
		units::Energy radiationSourceEnergy{};
		units::Time timestep = units::Time::from_value(std::numeric_limits<Real>::infinity());
		std::array<units::Velocity, ndim> signalSpeed{};
		std::map<int, units::Time> levelTimestep;
		std::vector<std::uint64_t> cachedBlocks;

		template <typename Archive>
		void serialize(Archive& archive, unsigned) {
			archive & tasks & timestep & signalSpeed & boundary & levelTimestep & radiationSourceEnergy & cachedBlocks;
		}
	};

	enum class Operation
	{
		Backup, Restore, Normalize, ResetFlux,
		Timestep,
		Advance,
		Probe,
		Reflux,
		RefluxTracked,
		PrepareGravityEnergy,
		FinishGravityEnergy,
		KickTracked,
		Kick,
		SaveRadiationStep,
		PredictRadiationStep,
		FinishRadiationStep,
		ForecastRadiationStep,
		RestoreRadiationStep,
		CoupleRadiation
	};

} // namespace runtime_detail

using namespace runtime_detail; // Private implementation header only.

class LocalExecutor
#ifdef OCTOTIGERII_WITH_HPX
  : public hpx::components::component_base<LocalExecutor>
#endif
{

public:
	LocalExecutor() = default;

	LocalExecutor(Config config, std::vector<Subgrid> blocks, FieldDirectory fields, std::size_t owner);

	void initialize();

	void begin(std::uint64_t generation, unsigned bank, units::Time time, int level, std::vector<HaloTime> times, Real fluxWeight,
		SourcePredictor source, Operation operation, std::vector<std::uint64_t> cachedWork);

	std::vector<std::uint64_t> claim(std::uint64_t generation);

	PhaseResult run(Operation operation, units::Time dt, std::uint64_t generation, std::vector<storage::Locality> const& executors);

	std::vector<Snapshot> snapshots(unsigned bank, mesh::TimeState time) const;

#ifdef OCTOTIGERII_WITH_HPX
	HPX_DEFINE_COMPONENT_ACTION(LocalExecutor, initialize, InitializeAction)
	HPX_DEFINE_COMPONENT_ACTION(LocalExecutor, begin, BeginAction)
	HPX_DEFINE_COMPONENT_ACTION(LocalExecutor, claim, ClaimAction)
	HPX_DEFINE_COMPONENT_ACTION(LocalExecutor, run, RunAction)
	HPX_DEFINE_COMPONENT_ACTION(LocalExecutor, snapshots, SnapshotsAction)
#endif

private:
	Config config_;
	finiteVolume::AnalyticBoundary<hydro::ConservedState> hydroBoundary_;
	finiteVolume::AnalyticBoundary<radiation::RadiationSystem::State> radiationBoundary_;
	ProblemRadiationMaterial radiationMaterial_;
	std::vector<Subgrid> blocks_;
	FieldDirectory fields_;
	std::size_t owner_ = 0;
	std::vector<std::uint64_t> owned_, active_;
	std::vector<HaloTime> times_;
	Real fluxWeight_ = 0;
	SourcePredictor source_;
	std::map<std::uint64_t, HaloPlan> plans_;
	std::map<std::uint64_t, std::vector<FluxCorrection>> refluxPlans_;
	std::map<std::uint64_t, std::vector<GravityWorkFace>> gravityWorkPlans_;
	std::vector<Workspace> workspaces_;
	// One slot per block; only its assigned task touches the slot during a phase.
	// Empty slots borrow reusable buffers from this locality's pool.
	std::vector<std::unique_ptr<BlockHaloCache>> blockCaches_, cachePool_;
	std::size_t gasCacheCapacity_ = 0, radiationCacheCapacity_ = 0;
	bool allowStealing_ = true;
	std::unique_ptr<BlockHaloCache> makeCache() const;
	std::unique_ptr<BlockHaloCache> acquireCache();
	void releaseCache(std::unique_ptr<BlockHaloCache>& cache);
#ifdef OCTOTIGERII_WITH_HPX
	hpx::mutex queueMutex_;
#else
	std::mutex queueMutex_;
#endif
	std::size_t next_ = 0;
	std::uint64_t generation_ = 0;
	unsigned bank_ = 0;
	units::Time time_{};

	PhaseResult worker(Operation operation, units::Time dt, std::uint64_t generation, std::vector<storage::Locality> const& executors, Workspace& workspace);

	template <typename State>
	void fluxRegister(storage::ColumnHandle<State> const& field, storage::Range range, bool reset);

	void fluxRegisters(Subgrid const& block, bool reset);
	void resetFlux(Subgrid const& block);
	void accumulateFlux(Subgrid const& block);

	void advanceSpecies(Subgrid const& block, HaloPlan const& plan, units::Time dt, std::vector<hydro::ConservedState> const& gasGhosts);

	void refluxSpecies(Subgrid const& block, std::vector<FluxCorrection> const& plan, units::Time dt);

	units::Time timestep(Subgrid const& block, std::array<units::Velocity, ndim>& maximumSpeed) const;

	template <typename System>
	void reflux(Subgrid const& block, std::vector<FluxCorrection> const& plan, storage::ColumnHandle<typename System::State> const& fields,
		storage::ColumnHandle<typename System::Flux> const& fluxFields, System const& system, units::Time dt, bool synchronize = true);

	BoundaryTransport finishGravityEnergy(Subgrid const& block, std::vector<GravityWorkFace> const& plan, units::Time dt);

	void kick(Subgrid const& block, units::Time dt, bool trackWork = false);
	units::Energy radiationSource(Subgrid const& block, Operation operation, units::Time dt);
};

class Runtime::Impl {
public:
	Config config;
	std::vector<storage::Locality> localities;
	std::unique_ptr<CartesianTopology> topology;
	std::unique_ptr<FieldRepository> fields;
#if OCTOTIGERII_GRAVITY
	std::unique_ptr<gravity::FieldSolver> gravitySolver;
	struct GlobalRotationWork {
		storage::PartitionSet store;
		gravity::RotationWorkWorkspace scratch;
		storage::Field<units::VelocitySquared> potential;
		storage::Field<units::Density> density;
		storage::Field<units::EnergyDensity> work;
		GlobalRotationWork(storage::Layout const& cells, std::vector<storage::Locality> const& owners)
		  : store(owners), scratch(cells, store), potential(cells, store, 2, "rotationWork.globalPotential"),
			density(cells, store, 1, "rotationWork.globalDensity"), work(cells, store, 1, "rotationWork.globalWork") {}
	};
	std::unique_ptr<GlobalRotationWork> globalRotation;
	gravity::Statistics rotationStatistics;
#endif
#ifdef OCTOTIGERII_WITH_HPX
	std::vector<hpx::id_type> executors;
#else
	std::unique_ptr<LocalExecutor> executor;
#endif
#ifdef OCTOTIGERII_WITH_HPX
	mutable hpx::mutex apiMutex;
#else
	mutable std::mutex apiMutex;
#endif
	unsigned bank = 0;
	std::uint64_t generation = 0;
	std::uint64_t dispatch = 0;
	std::vector<std::size_t> cacheOwners;
	mesh::TimeState time;
	units::Time gravityTime{};
	bool gravityReady = false, gravityEnergyActive = false;
	bool regridEnergyPending = false;
	units::Time gravityEnergyStart{};
	SchedulingStatistics statistics;
	BoundaryTransport boundary;
	units::Energy radiationSourceEnergy{};
	units::Energy eosFloorEnergy{};
	std::uint64_t eosFloorCells = 0;
	void applyEosFloor();
	refinement::Criteria criteria;
	std::unique_ptr<amr::Hierarchy> shadow;
	std::uint64_t lastRegridStep = 0;
	std::array<units::Velocity, ndim> signalSpeed{};
	std::array<units::Length, ndim> travel{}, travelBudget{};
	bool regridInitialized = false;
	struct LevelState {
		unsigned bank = 0;
		units::Time begin{}, end{};
		bool pending = false;
		unsigned predictorBank = ~0u;
	};
	std::vector<LevelState> levels;
	std::optional<RadiationStepFields> coupledStep;

	bool timeRefinement() const;
	int coarsestLevel() const;


	std::vector<Snapshot> exportSnapshots(std::optional<unsigned> requestedBank = {}, std::optional<mesh::TimeState> requestedTime = {}) const;

	void install(std::vector<mesh::BlockLocation> const& leaves, amr::Hierarchy const& source);

	// Operates only on the unpublished bank produced by solveGravity/setGravity.
	void recoverRegridEnergy(unsigned targetBank);

	/// Drain all locality results and require exactly one completed task per block.
	/// The caller may publish a bank only after this function succeeds.
	PhaseResult phase(Operation operation, units::Time dt, int level = -1, std::optional<units::Time> at = {}, Real fluxWeight = 0, SourcePredictor source = {});

#if OCTOTIGERII_GRAVITY
	// Physical gas states are used throughout. Provisional source impulses make
	// second-order hydro predictors; each HOLD shell replaces their sum with its
	// paired endpoint quadrature after the accepted mass fluxes have refluxed.
	struct GravityFrame {
		int level;
		units::Time begin, duration;
		struct Rotation {
			storage::Field<units::VelocitySquared> nested, shell, force;
			Rotation(storage::Layout const& cells, storage::PartitionSet const& store)
			  : nested(cells, store, 2, "rotationWork.nested"), shell(cells, store, 2, "rotationWork.shell"),
				force(cells, store, 2, "rotationWork.force") {}
		};
		std::unique_ptr<Rotation> rotation;
		storage::ColumnFields<gravity::State> nested, shell, force;
		storage::Field<units::Density> density;
		storage::ColumnFields<hydro::ConservedState> impulse;
		storage::Field<units::EnergyDensity> appliedWork;
		storage::Field<units::MassFlux> flux;
		GravityFrame(int l, units::Time at, units::Time h, storage::Layout const& cells, storage::Layout const& faces, storage::PartitionSet const& store, bool rotate)
		  : level(l), begin(at), duration(h), nested(cells, store, "timeGravity.nested"), shell(cells, store, "timeGravity.shell"),
			force(cells, store, "timeGravity.force"), density(cells, store, 1, "timeGravity.initialDensity"),
			impulse(cells, store, "timeGravity.impulse", false, 1), appliedWork(cells, store, 1, "timeGravity.appliedWork"),
			flux(faces, store, 1, "timeGravity.massFlux") {
			if (rotate) rotation = std::make_unique<Rotation>(cells, store);
		}
	};
	struct GravityInterval {
		storage::PartitionSet store;
		storage::Layout cells, faces;
		storage::ColumnFields<gravity::State> predictor;
		storage::ColumnFields<hydro::ConservedState> rate;
		storage::Field<units::EnergyDensity> deferred;
		std::vector<std::vector<GravityWorkFace>> plans;
		std::vector<std::vector<FluxCorrection>> reflux;
		std::vector<GravityFrame*> stack;
		std::map<int, std::unique_ptr<storage::ColumnFields<gravity::State>>> cachedNested;
		std::map<int, units::Time> cachedTime;
		std::unique_ptr<gravity::RotationWorkWorkspace> rotation;
		std::unique_ptr<storage::Field<units::VelocitySquared>> predictorRotation;
		std::map<int, std::unique_ptr<storage::Field<units::VelocitySquared>>> cachedRotation;
		units::Time sourceTime, duration;
		gravity::Statistics work;
		BoundaryTransport boundary;
		GravityInterval(Impl const& runtime, units::Time h)
		  : store(runtime.localities), cells(runtime.topology->storageLayout()),
			faces(std::vector<std::size_t>(runtime.topology->blocks().size(), allFaceCount(runtime.config.mesh.cells)), runtime.localities.size()),
			predictor(cells, store, "timeGravity.predictor", false, 1), rate(cells, store, "timeGravity.sourceRate", false, 1),
			deferred(cells, store, 1, "timeGravity.deferredWork"), sourceTime(runtime.time.time), duration(h) {
			if (runtime.config.frame.omega != units::InverseTime{} && runtime.config.gravity.energyTreatment == "mullen") {
				rotation = std::make_unique<gravity::RotationWorkWorkspace>(cells, store);
				predictorRotation = std::make_unique<storage::Field<units::VelocitySquared>>(cells, store, 1, "rotationWork.predictor");
			}
			for (auto const& b : runtime.topology->blocks()) {
				plans.push_back(makeGravityWorkPlan(runtime.config, runtime.topology->blocks(), b.id));
				reflux.push_back(makeRefluxPlan(runtime.config, runtime.topology->blocks(), b.id));
			}
		}
	};

	void addGravityWork(gravity::Statistics& total, gravity::Statistics const& add);

	void partialGravity(GravityInterval& interval, storage::ColumnHandle<gravity::State> const& output, unsigned outputBank,
		int minimumLevel, units::Time at, bool fullTargets = false,
		storage::FieldHandle<units::VelocitySquared> const& rotationOutput = {});

	storage::ColumnHandle<gravity::State> nestedGravity(GravityInterval& interval, int level, units::Time at);

	void shellGravity(GravityInterval& interval, GravityFrame& frame, unsigned endpoint, int nextLevel, units::Time at);

	void assemblePredictor(GravityInterval& interval, GravityFrame const& frame);

	void gravitySourceRate(GravityInterval& interval, int level);

	void provisionalGravity(GravityInterval& interval, GravityFrame& frame, units::Time step);

	void closeGravityFrame(GravityInterval& interval, GravityFrame& frame);

	void advanceGravityLevel(GravityInterval& interval, std::vector<int> const& occupied, std::size_t index,
		units::Time begin, units::Time duration, bool root = false);
#endif

	BoundaryTransport advanceLevel(std::vector<int> const& occupied, std::size_t index,
		units::Time begin, units::Time duration, bool root = false);
};

} // namespace octotigerII

#ifdef OCTOTIGERII_WITH_HPX
HPX_REGISTER_ACTION_DECLARATION(octotigerII::LocalExecutor::InitializeAction, octotigerII_initialize)
HPX_REGISTER_ACTION_DECLARATION(octotigerII::LocalExecutor::BeginAction, octotigerII_begin)
HPX_REGISTER_ACTION_DECLARATION(octotigerII::LocalExecutor::ClaimAction, octotigerII_claim)
HPX_REGISTER_ACTION_DECLARATION(octotigerII::LocalExecutor::RunAction, octotigerII_run)
HPX_REGISTER_ACTION_DECLARATION(octotigerII::LocalExecutor::SnapshotsAction, octotigerII_snapshots)
#endif
