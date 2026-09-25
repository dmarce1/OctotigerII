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
#include "octotigerII/storage/registry.hpp"
#include "octotigerII/subgrid/view.hpp"
#include "octotigerII/verification/analytic.hpp"
#if OCTOTIGERII_GRAVITY
#include "octotigerII/gravity/fieldSolver.hpp"
#include "octotigerII/gravity/fluxWork.hpp"
#endif

#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/include/components.hpp>
#include <hpx/serialization/map.hpp>
#include <hpx/runtime_distributed/find_all_localities.hpp>
#include <hpx/runtime_local/get_os_thread_count.hpp>
#include <hpx/synchronization/mutex.hpp>
#endif

namespace octotigerII {

namespace {

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

	void finish(std::vector<hpx::future<void>>& pending) {
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
	void initializeSpecies(FieldDirectory const& fields, Subgrid const& block, Snapshot const& snapshot) {
		for (std::size_t s = 0; s < fields.species.size(); ++s) {
			auto output = fields.species[s].output(block.interior, 0);
			std::copy_n(snapshot.species.at(s).values().data(), block.interior.count, output.data());
			fields.species[s].commit(block.interior, 0, output);
		}
	}

	// The convex remainder is O(dt²). It keeps a cold, accelerating midpoint
	// physical without adding energy to the accepted conservative update.
	void predictorKineticRemainder(hydro::ConservedState const& initial, hydro::ConservedState& predicted) {
		if (!(initial.density() > units::Density{}) || !(predicted.density() > units::Density{})) return;
		for (int d = 0; d < ndim; ++d) {
			auto const residual = predicted.momentum(d) - initial.momentum(d) * (predicted.density() / initial.density());
			predicted.totalEnergy() += 0.5 * residual * residual / predicted.density();
		}
	}

	template <typename System>
	class SolverWorkspace {
	public:
		using State = typename System::State;
		std::vector<State> ghosts;
		using Solver = physics::MusclHancock<System>;
		typename Solver::Workspace work;

		void advance(Subgrid const& block, storage::ColumnHandle<State> const& fields, HaloPlan const& plan, System const& system, unsigned bank,
			units::Time dt, units::Time time, physics::AnalyticBoundary<State> const& analytic, storage::ColumnHandle<typename System::Flux> const& fluxFields,
			bool amr, std::vector<HaloTime> const& times = {},
			storage::ColumnHandle<hydro::ConservedState> const& increment = {}, units::Time referenceStep = {}) {
			auto interior = fields.read(block.interior, bank).get();
			readHalo(fields, plan, bank, ghosts, times);
			std::vector<State> midpointGhosts;
			if constexpr (std::is_same_v<System, hydro::HydroSystem>) {
				if (referenceStep > units::Time{} && dt > units::Time{}) {
					std::vector<State> rates;
					readHalo(increment, plan, 0, rates);
					midpointGhosts = ghosts;
					for (std::size_t i = 0; i < midpointGhosts.size(); ++i) {
						midpointGhosts[i] += (dt / (2.0 * referenceStep)) * rates[i];
						predictorKineticRemainder(ghosts[i], midpointGhosts[i]);
					}
					// Prolong the midpoint donors themselves: prolongation of a rate
					// alone is not the derivative of the nonlinear spatial limiter.
					applyHaloBoundaries(plan, midpointGhosts, system, time + dt / 2.0, analytic);
				}
			}
			applyHaloBoundaries(plan, ghosts, system, time, analytic);
			PatchView<State> input(block, std::move(interior), plan, ghosts);
			auto output = fields.output(block.interior, bank ^ 1);
			{
				profiling::Region profile(std::is_same_v<System, hydro::HydroSystem> ? "hydro.advance" : "radiation.advance");
				auto writer = [&](mesh::Coordinates const& cell, State const& state) { output.put(block.layout.index(cell), state); };
				if constexpr (std::is_same_v<System, hydro::HydroSystem>) {
					if (referenceStep > units::Time{} && dt > units::Time{}) {
						auto rates = increment.read(block.interior, 0).get();
						Solver(system).advanceInto(input, dt, work, writer, [&](State value, mesh::Coordinates const& cell) {
							auto const initial = value;
							if (input.layout().isInterior(cell)) {
								value += (dt / (2.0 * referenceStep)) * rates.at(block.layout.index(input.layout().interiorCoordinates(cell)));
								predictorKineticRemainder(initial, value);
							} else value = midpointGhosts.at(plan.ghostIndices.at(input.layout().index(cell)));
							if (!system.admissible(value)) throw std::runtime_error("Gravity midpoint predictor is inadmissible");
							return value;
						}, false);
					} else Solver(system).advanceInto(input, dt, work, writer);
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
				} else Solver(system).advanceInto(input, dt, work, writer);
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

		BoundaryTransport boundaryTransport(Subgrid const& block, physics::BoundaryConditions const& boundaries, units::Time dt) const {
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

	// All source predictor fields use bank zero. increment stores the RHS integrated
	// over referenceStep, avoiding a dimensionless or untyped rate register.
	struct SourcePredictor {
		storage::ColumnHandle<hydro::ConservedState> increment;
		units::Time referenceStep{};
		template <typename Archive> void serialize(Archive& ar, unsigned) { ar & increment & referenceStep; }
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
		units::Time timestep = units::Time::from_value(std::numeric_limits<Real>::infinity());
		std::array<units::Velocity, ndim> signalSpeed{};
		std::map<int, units::Time> levelTimestep;

		template <typename Archive>
		void serialize(Archive& archive, unsigned) {
			archive & tasks & timestep & signalSpeed & boundary & levelTimestep;
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
		Kick
	};

}	 // namespace

// Metadata and ready-work queue per locality. Fields are only handles into
// the independent repository. Worker workspaces are pooled across stages.
/** @brief Per-locality work queue with bounded reusable worker storage.
 * Fields are handles into independent storage components. Queue claims hold a
 * short lock; no queue lock spans communication or numerical work. A dispatch
 * token rejects stale claims. Remote stealing transfers results to the owner.
 * @ingroup runtime
 */
class LocalExecutor
#ifdef OCTOTIGERII_WITH_HPX
  : public hpx::components::component_base<LocalExecutor>
#endif
{

public:
	LocalExecutor() = default;

	LocalExecutor(Config config, std::vector<Subgrid> blocks, FieldDirectory fields, std::size_t owner)
	  : config_(std::move(config))
	  , blocks_(std::move(blocks))
	  , fields_(std::move(fields))
	  , owner_(owner) {
		if (config_.mesh.boundary.contains(physics::BoundaryCondition::Analytic)) {
			auto evaluator = problemBoundary(config_);
			if (build::hydro && config_.hydroEnabled()) {
				hydroBoundary_ = [evaluator, gas = hydro::HydroSystem(config_.hydro)](
									 auto const& position, auto time) { return gas.conservedState(evaluator(position, time).hydro); };
			}
			if (build::radiation && config_.radiationEnabled()) {
				radiationBoundary_ = [evaluator](auto const& position, auto time) { return evaluator(position, time).radiation; };
			}
		}
		for (auto const& block : blocks_)
			if (block.interior.partition == owner_) owned_.push_back(block.id);
		for (auto id : owned_)
			plans_.emplace(id, makeHaloPlan(config_, blocks_, id));
		if (config_.hydroEnabled() && config_.gravityEnabled())
			for (auto id : owned_) gravityWorkPlans_.emplace(id, makeGravityWorkPlan(config_, blocks_, id));
		if (config_.amr.enabled)
			for (auto id : owned_)
				refluxPlans_.emplace(id, makeRefluxPlan(config_, blocks_, id));
		std::size_t workers = 1;
#ifdef OCTOTIGERII_WITH_HPX
		workers = hpx::get_os_thread_count();
#endif
		if (config_.runtime.workerTasks > 0) workers = std::size_t(config_.runtime.workerTasks);
		workers = std::max(std::size_t(1), std::min(workers, blocks_.size()));
		workspaces_.resize(workers);
	}

	void initialize() {
		profiling::Elapsed profile("runtime.initialize.wall_ns");
		for (auto id : owned_) {
			auto const& block = blocks_.at(id);
			auto const initial = initialSnapshot(config_, block.location);
			initializeSpecies(fields_, block, initial);
			if (build::hydro && config_.hydroEnabled()) initializeFields(fields_.hydro, block.interior, initial.hydro);
			if (build::radiation && config_.radiationEnabled()) initializeFields(fields_.radiation, block.interior, initial.radiation);
			if (build::gravity && config_.gravityEnabled()) {
				initializeFields(fields_.gravity, block.interior, initial.gravity);
				if (!config_.hydroEnabled()) {
					auto output = fields_.density.output(block.interior, 0);
					std::copy(initial.density.values().begin(), initial.density.values().end(), output.data());
					fields_.density.commit(block.interior, 0, output);
				}
			}
		}
	}

	void begin(std::uint64_t generation, unsigned bank, units::Time time, int level, std::vector<HaloTime> times, Real fluxWeight, SourcePredictor source = {}) {
		std::lock_guard guard(queueMutex_);
		generation_ = generation;
		bank_ = bank;
		time_ = time;
		times_ = std::move(times);
		fluxWeight_ = fluxWeight;
		source_ = std::move(source);
		active_.clear();
		for (auto id : owned_) if (level < 0 || blocks_[id].location.level == level) active_.push_back(id);
		next_ = 0;
	}

	std::vector<std::uint64_t> claim(std::uint64_t generation) {
		std::lock_guard guard(queueMutex_);
		if (generation != generation_) throw std::logic_error("Stale work request");
		if (next_ == active_.size()) return {};
		return {active_[next_++]};
	}

	PhaseResult run(Operation operation, units::Time dt, std::uint64_t generation, std::vector<storage::Locality> const& executors) {
		if (generation != generation_) throw std::logic_error("Stale stage dispatch");
#ifdef OCTOTIGERII_WITH_HPX
		std::vector<hpx::future<PhaseResult>> pending;
		for (std::size_t i = 0; i < workspaces_.size(); ++i)
			pending.push_back(hpx::async(profiling::annotated([&, i] { return worker(operation, dt, generation, executors, workspaces_[i]); },
				operation == Operation::Timestep	? "runtime.timestep.worker" :
					operation == Operation::Advance ? "runtime.advance.worker" :
					operation == Operation::Reflux	? "runtime.reflux.worker" :
													  "runtime.gravity_kick.worker")));
		auto results = collect(pending);
#else
		std::vector<PhaseResult> results{worker(operation, dt, generation, executors, workspaces_.front())};
#endif
		PhaseResult result;
		for (auto const& part : results) {
			result.boundary += part.boundary;
			result.tasks.localTasks += part.tasks.localTasks;
			result.tasks.stolenTasks += part.tasks.stolenTasks;
			result.timestep = std::min(result.timestep, part.timestep);
			for (auto const& [level, dt] : part.levelTimestep) {
				auto [where, inserted] = result.levelTimestep.emplace(level, dt);
				if (!inserted) where->second = std::min(where->second, dt);
			}
			for (int d = 0; d < ndim; ++d)
				result.signalSpeed[d] = std::max(result.signalSpeed[d], part.signalSpeed[d]);
		}
		return result;
	}

	std::vector<Snapshot> snapshots(unsigned bank, mesh::TimeState time) const {
		profiling::Elapsed profile("runtime.snapshots.local.wall_ns");
		std::vector<Snapshot> result;
		for (auto id : owned_) {
			auto const& block = blocks_.at(id);
			Snapshot snapshot;
			snapshot.location = block.location;
			snapshot.layout = block.layout;
			snapshot.lower = block.lower;
			snapshot.cellWidth = block.cellWidth;
			snapshot.time = time.time;
			snapshot.hydroEnabled = config_.hydroEnabled();
			snapshot.radiationEnabled = config_.radiationEnabled();
			snapshot.gravityEnabled = config_.gravityEnabled();
			auto exportOne = [&](auto const& field, auto& patch) {
				using Patch = std::remove_reference_t<decltype(patch)>;
				patch = Patch(block.layout, block.cellWidth, block.lower);
				patch.timeState() = time;
				exportFields(field, block.interior, bank, patch);
			};
			for (auto const& field : fields_.species) {
				snapshot.species.emplace_back(block.layout, block.cellWidth, block.lower);
				snapshot.species.back().timeState() = time;
				auto input = field.read(block.interior, bank).get();
				std::copy_n(input.data(), input.size(), snapshot.species.back().values().data());
			}
			if (build::hydro && config_.hydroEnabled()) exportOne(fields_.hydro, snapshot.hydro);
			if (build::radiation && config_.radiationEnabled()) exportOne(fields_.radiation, snapshot.radiation);
			if (build::gravity && config_.gravityEnabled()) {
				exportOne(fields_.gravity, snapshot.gravity);
				if (!snapshot.hydroEnabled) {
					snapshot.density = mesh::PatchData<units::Density>(block.layout, block.cellWidth, block.lower);
					snapshot.density.timeState() = time;
					auto input = fields_.density.read(block.interior, bank).get();
					std::copy_n(input.data(), input.size(), snapshot.density.values().begin());
				}
			}
			result.push_back(std::move(snapshot));
		}
		return result;
	}

#ifdef OCTOTIGERII_WITH_HPX
	HPX_DEFINE_COMPONENT_ACTION(LocalExecutor, initialize, InitializeAction)
	HPX_DEFINE_COMPONENT_ACTION(LocalExecutor, begin, BeginAction)
	HPX_DEFINE_COMPONENT_ACTION(LocalExecutor, claim, ClaimAction)
	HPX_DEFINE_COMPONENT_ACTION(LocalExecutor, run, RunAction)
	HPX_DEFINE_COMPONENT_ACTION(LocalExecutor, snapshots, SnapshotsAction)
#endif

private:
	Config config_;
	physics::AnalyticBoundary<hydro::ConservedState> hydroBoundary_;
	physics::AnalyticBoundary<radiation::RadiationSystem::State> radiationBoundary_;
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
#ifdef OCTOTIGERII_WITH_HPX
	hpx::mutex queueMutex_;
#else
	std::mutex queueMutex_;
#endif
	std::size_t next_ = 0;
	std::uint64_t generation_ = 0;
	unsigned bank_ = 0;
	units::Time time_{};

	PhaseResult worker(Operation operation, units::Time dt, std::uint64_t generation, std::vector<storage::Locality> const& executors, Workspace& workspace) {
		PhaseResult result;
		auto process = [&](std::uint64_t id, bool stolen) {
			auto const& block = blocks_.at(id);
			if (operation == Operation::Backup || operation == Operation::Restore || operation == Operation::Normalize) {
				unsigned const source = operation == Operation::Restore ? 2 :
					(operation == Operation::Normalize ? times_.at(block.location.level).bank : bank_);
				unsigned const destination = operation == Operation::Backup ? 2 :
					(operation == Operation::Restore ? bank_ : bank_ ^ 1);
				if (source != destination) {
					for (auto const& field : fields_.species) copyScalar(field, block.interior, source, destination);
					if (build::hydro && config_.hydroEnabled()) copyFields(fields_.hydro, block.interior, source, destination);
					if (build::radiation && config_.radiationEnabled()) copyFields(fields_.radiation, block.interior, source, destination);
					if (build::gravity && config_.gravityEnabled()) copyFields(fields_.gravity, block.interior, source, destination);
				}
			} else if (operation == Operation::ResetFlux) {
				resetFlux(block);
			} else if (operation == Operation::Timestep) {
				auto const dt = timestep(block, result.signalSpeed);
				result.timestep = std::min(result.timestep, dt);
				auto [where, inserted] = result.levelTimestep.emplace(block.location.level, dt);
				if (!inserted) where->second = std::min(where->second, dt);
			} else if (operation == Operation::Advance || operation == Operation::Probe) {
				// Owned plans are immutable. Stolen plans have only metadata and
				// are temporary, keeping persistent topology caches bounded.
				std::optional<HaloPlan> temporary;
				if (stolen) temporary = makeHaloPlan(config_, blocks_, id);
				auto const& plan = stolen ? *temporary : plans_.at(id);
				if constexpr (build::hydro) if (config_.hydroEnabled()) {
					workspace.hydro.advance(block, fields_.hydro, plan, hydro::HydroSystem(config_.hydro), bank_, dt, time_, hydroBoundary_,
						fields_.hydroFlux, config_.amr.enabled, times_, source_.increment, source_.referenceStep);
					if (config_.massFractions.enabled || config_.gravityEnabled()) {
						auto stored = fields_.massFlux.output(block.massFlux, 0);
						for (int d = 0; d < ndim; ++d)
							mesh::forEachCoordinate(block.layout.faceExtents(d), [&](auto const& face) {
								stored.data()[allFaceIndex(block.layout, d, face)] = workspace.hydro.work.fluxes[d][block.layout.faceIndex(d, face)].template get<0>();
							});
						fields_.massFlux.commit(block.massFlux, 0, stored);
					}
					if (config_.massFractions.enabled && operation != Operation::Probe) advanceSpecies(block, plan, dt, workspace.hydro.ghosts);
					result.boundary += workspace.hydro.boundaryTransport(block, config_.mesh.boundary, dt);
				}
				if constexpr (build::radiation) if (config_.radiationEnabled() && operation != Operation::Probe) {
					workspace.radiation.advance(block, fields_.radiation, plan, radiation::RadiationSystem(config_.radiation.lightSpeedRatio * constants::c),
						bank_, dt, time_, radiationBoundary_, fields_.radiationFlux, config_.amr.enabled, times_);
					result.boundary += workspace.radiation.boundaryTransport(block, config_.mesh.boundary, dt);
				}
				if (build::gravity && config_.gravityEnabled()) copyFields(fields_.gravity, block.interior, bank_);
				if (fluxWeight_ > 0 && operation != Operation::Probe) accumulateFlux(block);
			} else if (operation == Operation::Reflux || operation == Operation::RefluxTracked) {
				auto const plan = !config_.amr.enabled ? std::vector<FluxCorrection>{} :
					(stolen ? makeRefluxPlan(config_, blocks_, id) : refluxPlans_.at(id));
				if (config_.massFractions.enabled) refluxSpecies(block, plan, dt);
				if (build::hydro && config_.hydroEnabled()) reflux(block, plan, fields_.hydro, fields_.hydroFlux, hydro::HydroSystem(config_.hydro), dt, operation != Operation::RefluxTracked);
				if (build::radiation && config_.radiationEnabled())
					reflux(block, plan, fields_.radiation, fields_.radiationFlux, radiation::RadiationSystem(config_.radiation.lightSpeedRatio * constants::c),
						dt);
			} else if (operation == Operation::PrepareGravityEnergy) {
				auto current = fields_.gravity.read(block.interior, bank_).get();
				auto old = fields_.oldGravity.output(block.interior, 0);
				auto work = fields_.gravityKickWork.output(block.interior, 0);
				for (std::size_t i = 0; i < block.interior.count; ++i) { old.put(i, current.at(i)); work.data()[i] = {}; }
				fields_.oldGravity.commit(block.interior, 0, old);
				fields_.gravityKickWork.commit(block.interior, 0, work);
			} else if (operation == Operation::FinishGravityEnergy) {
				auto const plan = stolen ? makeGravityWorkPlan(config_, blocks_, id) : gravityWorkPlans_.at(id);
				result.boundary += finishGravityEnergy(block, plan, dt);
				copyFields(fields_.gravity, block.interior, bank_);
				if (build::radiation && config_.radiationEnabled()) copyFields(fields_.radiation, block.interior, bank_);
				for (auto const& field : fields_.species) copyScalar(field, block.interior, bank_);
			} else if (build::hydro && config_.hydroEnabled()) {
				kick(block, dt, operation == Operation::KickTracked);
				for (auto const& field : fields_.species) copyScalar(field, block.interior, bank_);
				if (build::gravity && config_.gravityEnabled()) copyFields(fields_.gravity, block.interior, bank_);
				if (build::radiation && config_.radiationEnabled()) copyFields(fields_.radiation, block.interior, bank_);
			} else {
				throw std::logic_error("Gravity kicks are not part of this executable");
			}
			if (stolen)
				++result.tasks.stolenTasks;
			else
				++result.tasks.localTasks;
		};
		for (;;) {
			auto work = claim(generation);
			if (work.empty()) break;
			process(work.front(), false);
		}
#ifdef OCTOTIGERII_WITH_HPX
		if (config_.runtime.workStealing)
			for (std::size_t distance = 1; distance < executors.size(); ++distance) {
				auto const donor = (owner_ + distance) % executors.size();
				for (;;) {
					auto work = hpx::async<ClaimAction>(executors[donor], generation).get();
					if (work.empty()) break;
					process(work.front(), true);
				}
			}
#else
		(void) executors;
#endif
		return result;
	}

	template <typename State>
	void fluxRegister(storage::ColumnHandle<State> const& field, storage::Range range, bool reset) {
		auto output = field.output(range, 1);
		if (reset) {
			for (std::size_t i = 0; i < range.count; ++i) output.put(i, State{});
		} else {
			auto old = field.read(range, 1).get(), raw = field.read(range, 0).get();
			for (std::size_t i = 0; i < range.count; ++i) output.put(i, old.at(i) + fluxWeight_ * raw.at(i));
		}
		field.commit(range, 1, output);
	}

	void fluxRegisters(Subgrid const& block, bool reset) {
		if (build::hydro && config_.hydroEnabled()) fluxRegister(fields_.hydroFlux, block.boundaryFlux, reset);
		if (build::radiation && config_.radiationEnabled()) fluxRegister(fields_.radiationFlux, block.boundaryFlux, reset);
		for (auto const& field : fields_.speciesFlux) {
			auto output = field.output(block.boundaryFlux, 1);
			if (reset) std::fill_n(output.data(), block.boundaryFlux.count, units::MassFlux{});
			else {
				auto old = field.read(block.boundaryFlux, 1).get(), raw = field.read(block.boundaryFlux, 0).get();
				for (std::size_t i = 0; i < block.boundaryFlux.count; ++i) output.data()[i] = old.data()[i] + fluxWeight_ * raw.data()[i];
			}
			field.commit(block.boundaryFlux, 1, output);
		}
	}
	void resetFlux(Subgrid const& block) { fluxRegisters(block, true); }
	void accumulateFlux(Subgrid const& block) { fluxRegisters(block, false); }

	void advanceSpecies(Subgrid const& block, HaloPlan const& plan, units::Time dt, std::vector<hydro::ConservedState> const& gasGhosts) {
		std::vector<storage::Buffer<units::Density>> interiors;
		std::vector<std::vector<units::Density>> ghosts(fields_.species.size());
		for (std::size_t s = 0; s < fields_.species.size(); ++s) {
			interiors.push_back(fields_.species[s].read(block.interior, bank_).get());
			readHalo(fields_.species[s], plan, bank_, ghosts[s], times_);
			// Passive concentrations use positive injection on coarse donor cells.
			// Normalize the material partial densities together at each face below.
			for (auto const& p : plan.prolongations) ghosts[s][p.destination] = ghosts[s][p.center];
			for (auto const& g : plan.analyticGhosts)
				ghosts[s][g.destination] = config_.massFractions.species[s].initialFraction * gasGhosts.at(g.destination).density();
		}
		mesh::MeshLayout const padded(config_.mesh.cells, 2);
		auto read = [&](std::size_t s, mesh::Coordinates const& cell) {
			auto storage = cell;
			for (auto& coordinate : storage) coordinate += padded.ghostWidth();
			if (padded.isInterior(storage)) return interiors[s].data()[block.layout.index(cell)];
			return ghosts[s].at(plan.ghostIndices.at(padded.index(storage)));
		};
		auto mass = fields_.massFlux.read(block.massFlux, 0).get();
		auto flux = composition::fluxes(config_.massFractions, block.layout, read,
			[&](int d, auto const& face) { return mass.data()[allFaceIndex(block.layout, d, face)]; });
		for (std::size_t s = 0; s < fields_.species.size(); ++s) {
			auto output = fields_.species[s].output(block.interior, bank_ ^ 1);
			block.layout.forEachInterior([&](auto const& cell, std::size_t i) {
				output.data()[i] = composition::update(interiors[s].data()[i], flux[s], block.layout, cell, dt / block.cellWidth);
			});
			fields_.species[s].commit(block.interior, bank_ ^ 1, output);
			if (config_.amr.enabled) {
				auto boundary = fields_.speciesFlux[s].output(block.boundaryFlux, 0);
				int const n = config_.mesh.cells;
				for (int d = 0; d < ndim; ++d) for (bool upper : {false, true}) {
					auto extents = mesh::filledCoordinates(n); extents[d] = 1;
					mesh::forEachCoordinate(extents, [&](auto face) {
						face[d] = upper ? n : 0;
						boundary.data()[boundaryFluxIndex(n, d, upper, face)] = flux[s][d][block.layout.faceIndex(d, face)];
					});
				}
				fields_.speciesFlux[s].commit(block.boundaryFlux, 0, boundary);
			}
		}
	}

	void refluxSpecies(Subgrid const& block, std::vector<FluxCorrection> const& plan, units::Time dt) {
		if (plan.empty()) return;
		for (std::size_t s = 0; s < fields_.species.size(); ++s) {
			auto coarse = fields_.speciesFlux[s].read(block.boundaryFlux, 0).get();
			auto current = fields_.species[s].read(block.interior, bank_ ^ 1).get();
			std::vector<units::Density> correction(block.interior.count);
			for (auto const& face : plan) {
				units::MassFlux fine{};
				for (auto const& range : face.fineFluxes) fine += fields_.speciesFlux[s].read(range, times_.empty() ? 0 : 1).get().data()[0];
				fine /= Real(face.fineFluxes.size());
				correction[face.cell] += (Real(face.sign) * dt / block.cellWidth) * (fine - coarse.data()[face.coarseFlux]);
			}
			auto output = fields_.species[s].output(block.interior, bank_ ^ 1);
			for (std::size_t i = 0; i < block.interior.count; ++i) {
				auto const value = current.data()[i] + correction[i];
				if (!units::finite(value) || value < units::Density{}) throw std::runtime_error("Species reflux produced a negative/nonfinite partial density");
				output.data()[i] = value;
			}
			fields_.species[s].commit(block.interior, bank_ ^ 1, output);
		}
	}

	units::Time timestep(Subgrid const& block, std::array<units::Velocity, ndim>& maximumSpeed) const {
		auto result = units::Time::from_value(std::numeric_limits<Real>::infinity());
		std::array<units::Velocity, ndim> blockSpeed{};
		std::array<units::Acceleration, ndim> maximumAcceleration{};
		[[maybe_unused]] auto fieldStep = [&](auto const& fields, auto const& system) {
			auto input = fields.read(block.interior, bank_).get();
			profiling::Region profile("transport.signal_speed");
			units::InverseTime rate{};
			std::array<units::Velocity, ndim> speed{};
			for (std::size_t i = 0; i < block.interior.count; ++i)
				for (int d = 0; d < ndim; ++d)
					speed[d] = std::max(speed[d], system.maximumSignalSpeed(input.at(i), d));
			for (auto value : speed)
				rate += value / block.cellWidth;
			for (int d = 0; d < ndim; ++d)
				blockSpeed[d] = std::max(blockSpeed[d], speed[d]);
			if (!(rate > units::InverseTime{}) || !units::finite(rate)) throw std::runtime_error("No finite positive signal speed");
			return config_.timestep.cfl / rate;
		};
		if (build::hydro && config_.hydroEnabled()) {
			for (int d = 0; d < ndim; ++d)
				maximumAcceleration[d] = units::abs(config_.hydro.acceleration[d]);
			result = fieldStep(fields_.hydro, hydro::HydroSystem(config_.hydro));
			units::Acceleration acceleration{};
			for (auto component : config_.hydro.acceleration)
				acceleration += units::abs(component);
			if (build::gravity && config_.gravityEnabled()) {
				auto input = fields_.gravity.read(block.interior, bank_).get();
				for (std::size_t i = 0; i < block.interior.count; ++i) {
					units::Acceleration norm{};
					for (int d = 0; d < ndim; ++d) {
						auto const value = units::abs(input.at(i).acceleration(d) + config_.hydro.acceleration[d]);
						norm += value;
						maximumAcceleration[d] = std::max(maximumAcceleration[d], value);
					}
					acceleration = std::max(acceleration, norm);
				}
			}
			if (acceleration > units::Acceleration{}) {
				auto const b = config_.timestep.cfl / result;
				auto const a = 0.5 * acceleration / block.cellWidth;
				result = 2 * config_.timestep.cfl / (b + units::sqrt(b * b + 4.0 * a * config_.timestep.cfl));
				result = std::min(result, 0.2 * units::sqrt(block.cellWidth / acceleration));
			}
		}
		if (build::radiation && config_.radiationEnabled())
			result = std::min(result, fieldStep(fields_.radiation, radiation::RadiationSystem(config_.radiation.lightSpeedRatio * constants::c)));
		if (build::radiation && config_.radiationEnabled())
			for (auto& speed : blockSpeed)
				speed = std::max(speed, config_.radiation.lightSpeedRatio * constants::c);
		for (int d = 0; d < ndim; ++d) {
			if (units::finite(result)) blockSpeed[d] += maximumAcceleration[d] * result;
			maximumSpeed[d] = std::max(maximumSpeed[d], blockSpeed[d]);
		}
		return result;
	}

	template <typename System>
	void reflux(Subgrid const& block, std::vector<FluxCorrection> const& plan, storage::ColumnHandle<typename System::State> const& fields,
		storage::ColumnHandle<typename System::Flux> const& fluxFields, System const& system, units::Time dt, bool synchronize = true) {
		if (plan.empty()) {
			if constexpr (!std::is_same_v<System, hydro::HydroSystem>) return;
		}
		using State = typename System::State;
		storage::Columns<typename System::Flux> coarseFlux;
		if (!plan.empty()) coarseFlux = fluxFields.read(block.boundaryFlux, 0).get();
		auto current = fields.read(block.interior, bank_ ^ 1).get();
		std::vector<State> corrections(block.interior.count);
		for (auto const& face : plan) {
			typename System::Flux fine{};
			std::vector<storage::PendingColumns<typename System::Flux>> pending;
			for (auto const& range : face.fineFluxes)
				pending.push_back(fluxFields.read(range, times_.empty() ? 0 : 1));
			std::exception_ptr error;
			for (auto& read : pending) {
				try {
					fine += read.get().at(0);
				} catch (...) {
					if (!error) error = std::current_exception();
				}
			}
			if (error) std::rethrow_exception(error);
			fine /= Real(face.fineFluxes.size());
			corrections[face.cell] += (Real(face.sign) * dt / block.cellWidth) * (fine - coarseFlux.at(face.coarseFlux));
		}
		auto output = fields.output(block.interior, bank_ ^ 1);
		for (std::size_t i = 0; i < block.interior.count; ++i) {
			if constexpr (std::is_same_v<System, hydro::HydroSystem>)
				if (config_.massFractions.enabled) corrections[i].density() = {};
			auto candidate = system.correctRoundoff(current.at(i) + corrections[i], componentAbs(corrections[i]));
			if constexpr (requires { system.synchronize(candidate); }) if (synchronize) system.synchronize(candidate);
			if (!system.admissible(candidate)) throw std::runtime_error("AMR reflux/dual-energy synchronization produced an inadmissible state");
			output.put(i, candidate);
		}
		fields.commit(block.interior, bank_ ^ 1, output);
	}

	BoundaryTransport finishGravityEnergy(Subgrid const& block, std::vector<GravityWorkFace> const& plan, units::Time dt) {
		auto const old = fields_.oldGravity.read(block.interior, 0).get();
		auto const now = fields_.gravity.read(block.interior, bank_).get();
		auto const mass = fields_.massFlux.read(block.massFlux, 0).get();
		auto const kicks = fields_.gravityKickWork.read(block.interior, 0).get();
		std::vector<units::VelocitySquared> potential(block.interior.count);
		std::vector<units::EnergyDensity> work(block.interior.count);
		for (std::size_t i = 0; i < potential.size(); ++i) potential[i] = 0.5 * (old.at(i).potential() + now.at(i).potential());
		for (int d = 0; d < ndim; ++d) {
			auto extents = block.layout.faceExtents(d); extents[d] -= 2;
			mesh::forEachCoordinate(extents, [&](auto face) {
				++face[d]; auto left = face; --left[d];
				auto const l = block.layout.index(left), r = block.layout.index(face);
				auto const amount = (dt / block.cellWidth) * mass.data()[allFaceIndex(block.layout, d, face)];
				auto const energy = -0.5 * amount * (potential[r] - potential[l]);
				work[l] += energy; work[r] += energy;
			});
		}
		BoundaryTransport boundary;
		auto const volume = block.layout.cellMeasure(block.cellWidth);
		for (auto const& face : plan) {
			auto const flux = fields_.massFlux.read(face.massFlux, 0).get().data()[0];
			auto const outward = Real(face.sign) * face.areaFraction * (dt / block.cellWidth) * flux;
			units::VelocitySquared difference{};
			if (face.physical) {
				auto const acceleration = 0.5 * (old.at(face.cell).acceleration(face.axis) + now.at(face.cell).acceleration(face.axis));
				difference = -Real(face.sign) * 0.5 * block.cellWidth * acceleration;
				auto const transported = volume * outward * (potential[face.cell] + difference);
				if (transported < units::Energy{}) boundary.inward.potentialEnergy -= transported;
				else boundary.outward.potentialEnergy += transported;
			} else {
				auto const a = fields_.oldPotential.read(face.neighbor, 0).get().data()[0];
				auto const b = std::get<0>(fields_.gravity.fields).read(face.neighbor, bank_).get().data()[0];
				difference = face.potentialFraction * (0.5 * (a + b) - potential[face.cell]);
			}
			work[face.cell] -= outward * difference;
		}
		auto input = fields_.hydro.read(block.interior, bank_).get();
		auto output = fields_.hydro.output(block.interior, bank_ ^ 1);
		hydro::HydroSystem const gas(config_.hydro);
		for (std::size_t i = 0; i < block.interior.count; ++i) {
			auto state = input.at(i);
			if (config_.gravity.energyTreatment == "mullen") state.totalEnergy() += work[i] - kicks.data()[i];
			gas.synchronize(state);
			if (!gas.admissible(state)) {
				std::ostringstream message;
				message << "Conservative gravity work produced an inadmissible state: cell=" << i
					<< " rho=" << units::value(state.density()) << " E=" << units::value(state.totalEnergy())
					<< " work=" << units::value(work[i]) << " kickWork=" << units::value(kicks.data()[i])
					<< " dt=" << units::value(dt);
				throw std::runtime_error(message.str());
			}
			output.put(i, state);
		}
		fields_.hydro.commit(block.interior, bank_ ^ 1, output);
		return boundary;
	}

	void kick(Subgrid const& block, units::Time dt, bool trackWork = false) {
		auto input = fields_.hydro.read(block.interior, bank_).get();
		std::optional<storage::Columns<gravity::State>> gravity;
		if (build::gravity && config_.gravityEnabled()) gravity = fields_.gravity.read(block.interior, bank_).get();
		auto output = fields_.hydro.output(block.interior, bank_ ^ 1);
		{
			profiling::Region profile("gravity.kick");
			storage::Buffer<units::EnergyDensity> kickWork;
			if (trackWork) {
				auto prior = fields_.gravityKickWork.read(block.interior, 0).get();
				kickWork = fields_.gravityKickWork.output(block.interior, 0);
				std::copy_n(prior.data(), block.interior.count, kickWork.data());
			}
			for (std::size_t i = 0; i < block.interior.count; ++i) {
				auto state = input.at(i);
				units::EnergyDensity work{};
				for (int d = 0; d < ndim; ++d) {
					auto const old = state.momentum(d);
					auto acceleration = config_.hydro.acceleration[d];
					if (build::gravity && config_.gravityEnabled()) acceleration += gravity->at(i).acceleration(d);
					auto const impulse = dt * state.density() * acceleration;
					state.momentum(d) += impulse;
					work += impulse * (old + 0.5 * impulse) / state.density();
					if (trackWork) kickWork.data()[i] += dt * gravity->at(i).acceleration(d) * (old + 0.5 * impulse);
				}
				state.totalEnergy() += work;
				if (!trackWork) hydro::HydroSystem(config_.hydro).synchronize(state);
				if (!hydro::HydroSystem(config_.hydro).admissible(state)) throw std::runtime_error("Invalid gravity kick state");
				output.put(i, state);
			}
			if (trackWork) fields_.gravityKickWork.commit(block.interior, 0, kickWork);
		}
		fields_.hydro.commit(block.interior, bank_ ^ 1, output);
	}
};

}	 // namespace octotigerII

#ifdef OCTOTIGERII_WITH_HPX
using LocalExecutorComponent = hpx::components::component<octotigerII::LocalExecutor>;
HPX_REGISTER_COMPONENT(LocalExecutorComponent, octotigerII_executor)
HPX_REGISTER_ACTION(octotigerII::LocalExecutor::InitializeAction, octotigerII_initialize)
HPX_REGISTER_ACTION(octotigerII::LocalExecutor::BeginAction, octotigerII_begin)
HPX_REGISTER_ACTION(octotigerII::LocalExecutor::ClaimAction, octotigerII_claim)
HPX_REGISTER_ACTION(octotigerII::LocalExecutor::RunAction, octotigerII_run)
HPX_REGISTER_ACTION(octotigerII::LocalExecutor::SnapshotsAction, octotigerII_snapshots)
#endif

namespace octotigerII {

class Runtime::Impl {
public:
	Config config;
	std::vector<storage::Locality> localities;
	std::unique_ptr<CartesianTopology> topology;
	std::unique_ptr<FieldRepository> fields;
#if OCTOTIGERII_GRAVITY
	std::unique_ptr<gravity::FieldSolver> gravitySolver;
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
	mesh::TimeState time;
	units::Time gravityTime{};
	bool gravityReady = false, gravityEnergyActive = false;
	bool regridEnergyPending = false;
	units::Time gravityEnergyStart{};
	SchedulingStatistics statistics;
	BoundaryTransport boundary;
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

	bool timeRefinement() const {
		if (!config.amr.enabled || !config.timestep.refinement || config.hasExternalAcceleration()) return false;
		auto const& blocks = topology->blocks();
		return std::any_of(blocks.begin(), blocks.end(), [&](auto const& b) { return b.location.level != blocks.front().location.level; });
	}
	int coarsestLevel() const {
		int result = config.amr.maxLevel;
		for (auto const& b : topology->blocks()) result = std::min(result, b.location.level);
		return result;
	}


	std::vector<Snapshot> exportSnapshots(std::optional<unsigned> requestedBank = {}, std::optional<mesh::TimeState> requestedTime = {}) const {
		std::vector<Snapshot> result;
		auto const inputBank = requestedBank.value_or(bank);
		auto const inputTime = requestedTime.value_or(time);
#ifdef OCTOTIGERII_WITH_HPX
		std::vector<hpx::future<std::vector<Snapshot>>> pending;
		for (auto const& id : executors)
			pending.push_back(hpx::async<LocalExecutor::SnapshotsAction>(id, inputBank, inputTime));
		for (auto& partition : collect(pending))
			for (auto& block : partition)
				result.push_back(std::move(block));
#else
		result = executor->snapshots(inputBank, inputTime);
#endif
		return result;
	}

	void install(std::vector<mesh::BlockLocation> const& leaves, amr::Hierarchy const& source) {
		bool const conserveEnergy = config.hydroEnabled() && config.gravityEnabled() && config.gravity.conserveRegridEnergy;
		auto nextTopology = std::make_unique<CartesianTopology>(config, localities.size(), leaves);
		auto nextFields = std::make_unique<FieldRepository>(config, nextTopology->storageLayout(), localities);
		auto const& directory = nextFields->directory();
		for (auto const& block : nextTopology->blocks()) {
			auto const snapshot = source.transfer(block.location);
			initializeSpecies(directory, block, snapshot);
			if (build::hydro && config.hydroEnabled()) initializeFields(directory.hydro, block.interior, snapshot.hydro);
			if (conserveEnergy) {
				// Reuse the inactive kick register for conservatively transferred E+rho*phi/2.
				// Hydro remains gas-only for mesh selection, EOS calls and shadow estimates.
				auto const energy = source.transferGasGravityEnergy(block.location);
				auto output = directory.gravityKickWork.output(block.interior, 0);
				std::copy(energy.begin(), energy.end(), output.data());
				directory.gravityKickWork.commit(block.interior, 0, output);
			}
			if (build::radiation && config.radiationEnabled()) initializeFields(directory.radiation, block.interior, snapshot.radiation);
			if (build::gravity && config.gravityEnabled()) {
				initializeFields(directory.gravity, block.interior, snapshot.gravity);
				if (!config.hydroEnabled()) {
					auto output = directory.density.output(block.interior, 0);
					std::copy(snapshot.density.values().begin(), snapshot.density.values().end(), output.data());
					directory.density.commit(block.interior, 0, output);
				}
			}
		}
#ifdef OCTOTIGERII_WITH_HPX
		std::vector<hpx::future<hpx::id_type>> pending;
		for (std::size_t i = 0; i < localities.size(); ++i)
			pending.push_back(hpx::new_<LocalExecutor>(localities[i], config, nextTopology->blocks(), directory, i));
		auto nextExecutors = collect(pending);
		executors.swap(nextExecutors);
#else
		auto nextExecutor = std::make_unique<LocalExecutor>(config, nextTopology->blocks(), directory, 0);
		executor.swap(nextExecutor);
#endif
#if OCTOTIGERII_GRAVITY
		gravitySolver.reset();
#endif
		topology.swap(nextTopology);
		fields.swap(nextFields);
		bank = 0;
		++generation;
		gravityReady = false;
		regridEnergyPending = conserveEnergy;
	}

	// Operates only on the unpublished bank produced by solveGravity/setGravity.
	void recoverRegridEnergy(unsigned targetBank) {
		if (!regridEnergyPending) return;
		auto const& directory = fields->directory();
		hydro::HydroSystem const gas(config.hydro);
		for (auto const& block : topology->blocks()) {
			auto const combined = directory.gravityKickWork.read(block.interior, 0).get();
			auto const gravity = directory.gravity.read(block.interior, targetBank).get();
			auto const input = directory.hydro.read(block.interior, targetBank).get();
			auto output = directory.hydro.output(block.interior, targetBank);
			for (std::size_t i = 0; i < block.interior.count; ++i) {
				auto state = input.at(i);
				state.totalEnergy() = combined.data()[i] - 0.5 * state.density() * gravity.at(i).potential();
				gas.synchronize(state);
				if (!gas.admissible(state)) throw std::runtime_error("Regrid gravity-energy recovery produced an inadmissible gas state");
				output.put(i, state);
			}
			directory.hydro.commit(block.interior, targetBank, output);
		}
	}

	/// Drain all locality results and require exactly one completed task per block.
	/// The caller may publish a bank only after this function succeeds.
	PhaseResult phase(Operation operation, units::Time dt, int level = -1, std::optional<units::Time> at = {}, Real fluxWeight = 0, SourcePredictor source = {}) {
		++dispatch;
		auto const stageTime = at.value_or(time.time);
		auto const stageBank = level >= 0 && !levels.empty() ? levels.at(level).bank : bank;
		std::vector<HaloTime> haloTimes;
		for (auto const& state : levels) {
			Real alpha = 0;
			if (state.pending) {
				alpha = (stageTime - state.begin) / (state.end - state.begin);
				auto const tolerance = 128 * epsilonR * std::max({units::abs(stageTime), units::abs(state.begin), units::abs(state.end)});
				if (stageTime < state.begin - tolerance || stageTime > state.end + tolerance) throw std::logic_error("Halo time outside coarse predictor interval");
				alpha = std::clamp(alpha, Real(0), Real(1));
			}
			haloTimes.push_back({state.bank, alpha, state.predictorBank});
		}
#ifdef OCTOTIGERII_WITH_HPX
		std::vector<hpx::future<void>> starts;
		for (auto const& id : executors)
			starts.push_back(hpx::async<LocalExecutor::BeginAction>(id, dispatch, stageBank, stageTime, level, haloTimes, fluxWeight, source));
		finish(starts);
		std::vector<hpx::future<PhaseResult>> pending;
		for (auto const& id : executors)
			pending.push_back(hpx::async<LocalExecutor::RunAction>(id, operation, dt, dispatch, executors));
		auto results = collect(pending);
#else
		executor->begin(dispatch, stageBank, stageTime, level, haloTimes, fluxWeight, source);
		std::vector<PhaseResult> results{executor->run(operation, dt, dispatch, {})};
#endif
		PhaseResult result;
		for (auto const& part : results) {
			result.boundary += part.boundary;
			result.tasks.localTasks += part.tasks.localTasks;
			result.tasks.stolenTasks += part.tasks.stolenTasks;
			result.timestep = std::min(result.timestep, part.timestep);
			for (auto const& [level, dt] : part.levelTimestep) {
				auto [where, inserted] = result.levelTimestep.emplace(level, dt);
				if (!inserted) where->second = std::min(where->second, dt);
			}
			for (int d = 0; d < ndim; ++d)
				result.signalSpeed[d] = std::max(result.signalSpeed[d], part.signalSpeed[d]);
		}
		auto const expected = std::count_if(topology->blocks().begin(), topology->blocks().end(), [&](auto const& b) { return level < 0 || b.location.level == level; });
		if (result.tasks.localTasks + result.tasks.stolenTasks != std::size_t(expected))
			throw std::logic_error("Stage did not complete every output range");
		statistics.localTasks += result.tasks.localTasks;
		statistics.stolenTasks += result.tasks.stolenTasks;
		return result;
	}

#if OCTOTIGERII_GRAVITY
	// Physical gas states are used throughout. Provisional source impulses make
	// second-order hydro predictors; each HOLD shell replaces their sum with its
	// paired endpoint quadrature after the accepted mass fluxes have refluxed.
	struct GravityFrame {
		int level;
		units::Time duration;
		storage::ColumnFields<gravity::State> nested, shell, force;
		storage::Field<units::Density> density;
		storage::ColumnFields<hydro::ConservedState> impulse;
		storage::Field<units::EnergyDensity> appliedWork;
		storage::Field<units::MassFlux> flux;
		GravityFrame(int l, units::Time h, storage::Layout const& cells, storage::Layout const& faces, storage::PartitionSet const& store)
		  : level(l), duration(h), nested(cells, store, "timeGravity.nested"), shell(cells, store, "timeGravity.shell"),
			force(cells, store, "timeGravity.force"), density(cells, store, 1, "timeGravity.initialDensity"),
			impulse(cells, store, "timeGravity.impulse", false, 1), appliedWork(cells, store, 1, "timeGravity.appliedWork"),
			flux(faces, store, 1, "timeGravity.massFlux") {}
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
		units::Time duration;
		gravity::Statistics work;
		BoundaryTransport boundary;
		GravityInterval(Impl const& runtime, units::Time h)
		  : store(runtime.localities), cells(runtime.topology->storageLayout()),
			faces(std::vector<std::size_t>(runtime.topology->blocks().size(), allFaceCount(runtime.config.mesh.cells)), runtime.localities.size()),
			predictor(cells, store, "timeGravity.predictor", false, 1), rate(cells, store, "timeGravity.sourceRate", false, 1),
			deferred(cells, store, 1, "timeGravity.deferredWork"), duration(h) {
			for (auto const& b : runtime.topology->blocks()) {
				plans.push_back(makeGravityWorkPlan(runtime.config, runtime.topology->blocks(), b.id));
				reflux.push_back(makeRefluxPlan(runtime.config, runtime.topology->blocks(), b.id));
			}
		}
	};

	void addGravityWork(gravity::Statistics& total, gravity::Statistics const& add) {
		total.multipolePairs += add.multipolePairs; total.directPairs += add.directPairs;
		total.workerTasks += add.workerTasks; total.ewaldPairs += add.ewaldPairs; total.reflectedPairs += add.reflectedPairs;
		if (total.localityCells.size() < add.localityCells.size()) total.localityCells.resize(add.localityCells.size());
		for (std::size_t i = 0; i < add.localityCells.size(); ++i) total.localityCells[i] += add.localityCells[i];
	}

	void partialGravity(GravityInterval& interval, storage::ColumnHandle<gravity::State> const& output, unsigned outputBank,
		int minimumLevel, units::Time at, bool fullTargets = false) {
		if (!gravitySolver) gravitySolver = std::make_unique<gravity::FieldSolver>(config, topology->blocks(), fields->directory(), localities);
		gravity::FieldSolveRequest request;
		request.output = output; request.outputBank = outputBank;
		for (auto const& b : topology->blocks()) {
			auto const& state = levels.at(b.location.level);
			Real alpha = 0;
			if (state.pending) alpha = std::clamp(Real((at - state.begin) / (state.end - state.begin)), Real(0), Real(1));
			bool const selected = b.location.level >= minimumLevel;
			request.sources.push_back({state.bank, state.predictorBank == ~0u ? (state.bank ^ 1) : state.predictorBank,
				selected ? 1 - alpha : 0, selected ? alpha : 0});
			request.targets.push_back(fullTargets || selected);
			if (!fullTargets && !selected) {
				auto out = output.output(b.interior, outputBank);
				for (std::size_t i = 0; i < b.interior.count; ++i) out.put(i, {});
				output.commit(b.interior, outputBank, out);
			}
		}
		addGravityWork(interval.work, gravitySolver->solve(request));
	}

	storage::ColumnHandle<gravity::State> nestedGravity(GravityInterval& interval, int level, units::Time at) {
		auto& cache = interval.cachedNested[level];
		if (!cache) cache = std::make_unique<storage::ColumnFields<gravity::State>>(interval.cells, interval.store, "timeGravity.cachedNested", false, 1);
		if (!interval.cachedTime.contains(level) || interval.cachedTime.at(level) != at) {
			partialGravity(interval, cache->handle(), 0, level, at);
			interval.cachedTime[level] = at;
		}
		return cache->handle();
	}

	void shellGravity(GravityInterval& interval, GravityFrame& frame, unsigned endpoint, int nextLevel, units::Time at) {
		auto const nested = nestedGravity(interval, frame.level, at);
		std::optional<storage::ColumnHandle<gravity::State>> fast;
		if (nextLevel >= 0) fast = nestedGravity(interval, nextLevel, at);
		for (auto const& b : topology->blocks()) {
			auto const all = nested.read(b.interior, 0).get();
			std::optional<storage::Columns<gravity::State>> subset;
			if (fast) subset = fast->read(b.interior, 0).get();
			auto n = frame.nested.handle().output(b.interior, endpoint);
			auto output = frame.shell.handle().output(b.interior, endpoint);
			for (std::size_t i = 0; i < b.interior.count; ++i) {
				n.put(i, all.at(i));
				output.put(i, all.at(i) - (subset ? subset->at(i) : gravity::State{}));
			}
			frame.nested.handle().commit(b.interior, endpoint, n);
			frame.shell.handle().commit(b.interior, endpoint, output);
		}
	}

	void assemblePredictor(GravityInterval& interval, GravityFrame const& frame) {
		bool const conventional = config.gravity.timeIntegration == "conventional";
		for (auto const& b : topology->blocks()) {
			auto output = interval.predictor.handle().output(b.interior, 0);
			auto current = (conventional ? frame.force.handle() : frame.nested.handle()).read(b.interior, 0).get();
			std::vector<storage::Columns<gravity::State>> ancestors;
			if (!conventional) for (auto const* f : interval.stack) if (f != &frame) ancestors.push_back(f->shell.handle().read(b.interior, 0).get());
			for (std::size_t i = 0; i < b.interior.count; ++i) {
				auto value = current.at(i);
				for (auto const& a : ancestors) value += a.at(i);
				output.put(i, value);
			}
			interval.predictor.handle().commit(b.interior, 0, output);
		}
	}

	void gravitySourceRate(GravityInterval& interval, int level) {
		for (auto const& b : topology->blocks()) if (level < 0 || b.location.level >= level) {
			auto work = gravity::fluxWork(b, interval.plans[b.id], interval.predictor.handle(), 0,
				fields->directory().massFlux, 0, interval.duration, true);
			auto prior = interval.rate.handle().read(b.interior, 0).get();
			auto out = interval.rate.handle().output(b.interior, 0);
			std::vector<hydro::ConservedState> correction(b.interior.count);
			auto coarse = fields->directory().hydroFlux.read(b.boundaryFlux, 0).get();
			for (auto const& face : interval.reflux[b.id]) {
				hydro::ConservedFlux fine{};
				for (auto const& range : face.fineFluxes) fine += fields->directory().hydroFlux.read(range, 0).get().at(0);
				fine /= Real(face.fineFluxes.size());
				correction[face.cell] += (Real(face.sign) * interval.duration / b.cellWidth) * (fine - coarse.at(face.coarseFlux));
			}
			auto gas = fields->directory().hydro.read(b.interior, levels.at(b.location.level).bank).get();
			auto force = interval.predictor.handle().read(b.interior, 0).get();
			for (std::size_t i = 0; i < b.interior.count; ++i) {
				auto value = prior.at(i) + correction[i];
				for (int d = 0; d < ndim; ++d) {
					value.momentum(d) += interval.duration * gas.at(i).density() * force.at(i).acceleration(d);
					if (config.gravity.energyTreatment == "naive") value.totalEnergy() += interval.duration * gas.at(i).momentum(d) * force.at(i).acceleration(d);
				}
				if (config.gravity.energyTreatment == "mullen") value.totalEnergy() += work.work[i];
				out.put(i, value);
			}
			interval.rate.handle().commit(b.interior, 0, out);
		}
	}

	void provisionalGravity(GravityInterval& interval, GravityFrame& frame, units::Time step) {
		auto const& directory = fields->directory();
		bool const conventional = config.gravity.timeIntegration == "conventional";
		for (auto const& b : topology->blocks()) if (b.location.level == frame.level) {
			auto const bank = levels.at(frame.level).bank;
			auto const old = directory.hydro.read(b.interior, bank).get();
			auto const next = directory.hydro.read(b.interior, bank ^ 1).get();
			auto const force = interval.predictor.handle().read(b.interior, 0).get();
			auto const mass = directory.massFlux.read(b.massFlux, 0).get();
			std::vector<units::EnergyDensity> heat(b.interior.count);
			for (auto* f : interval.stack) {
				// On the level being forecast, the nested field also forecasts the
				// as-yet-unopened descendants' energy transfer across its boundary.
				auto const field = f == &frame ? f->nested.handle() : f->shell.handle();
				auto const work = gravity::fluxWork(b, interval.plans[b.id], field, 0, directory.massFlux, 0, step, false);
				auto priorWork = f->appliedWork.handle().read(b.interior, 0).get();
				auto applied = f->appliedWork.handle().output(b.interior, 0);
				for (std::size_t i = 0; i < b.interior.count; ++i) { applied.data()[i] = priorWork.data()[i] + work.work[i]; heat[i] += work.work[i]; }
				f->appliedWork.handle().commit(b.interior, 0, applied);
				auto priorFlux = f->flux.handle().read(b.massFlux, 0).get();
				auto integrated = f->flux.handle().output(b.massFlux, 0);
				for (std::size_t i = 0; i < b.massFlux.count; ++i) integrated.data()[i] = priorFlux.data()[i] + (step / f->duration) * mass.data()[i];
				f->flux.handle().commit(b.massFlux, 0, integrated);
				if (!conventional || f == &frame) {
					auto g = (conventional ? f->force.handle() : f->shell.handle()).read(b.interior, 0).get();
					auto prior = f->impulse.handle().read(b.interior, 0).get();
					auto impulses = f->impulse.handle().output(b.interior, 0);
					for (std::size_t i = 0; i < b.interior.count; ++i) {
						auto p = prior.at(i);
						for (int d = 0; d < ndim; ++d) p.momentum(d) += (step / 2.0) * (old.at(i).density() + next.at(i).density()) * g.at(i).acceleration(d);
						impulses.put(i, p);
					}
					f->impulse.handle().commit(b.interior, 0, impulses);
				}
			}
			auto output = directory.hydro.output(b.interior, bank ^ 1);
			for (std::size_t i = 0; i < b.interior.count; ++i) {
				auto value = next.at(i);
				for (int d = 0; d < ndim; ++d) {
					value.momentum(d) += (step / 2.0) * (old.at(i).density() + value.density()) * force.at(i).acceleration(d);
					if (config.gravity.energyTreatment == "naive") value.totalEnergy() += (step / 2.0) * force.at(i).acceleration(d) * (old.at(i).momentum(d) + value.momentum(d));
				}
				if (config.gravity.energyTreatment == "mullen") value.totalEnergy() += heat[i];
				hydro::HydroSystem const system(config.hydro);
				system.synchronize(value);
				if (!system.admissible(value)) throw std::runtime_error("Provisional gravity source produced an inadmissible state");
				output.put(i, value);
			}
			directory.hydro.commit(b.interior, bank ^ 1, output);
		}
	}

	void closeGravityFrame(GravityInterval& interval, GravityFrame& frame) {
		auto const& directory = fields->directory();
		bool const conventional = config.gravity.timeIntegration == "conventional";
		for (auto const& b : topology->blocks()) {
			auto work = gravity::fluxWork(b, interval.plans[b.id], frame.shell.handle(), 0, frame.shell.handle(), 1,
				frame.flux.handle(), 0, frame.duration);
			interval.boundary += work.boundary;
			auto const deferred = interval.deferred.handle().read(b.interior, 0).get();
			if (b.location.level < frame.level) {
				auto out = interval.deferred.handle().output(b.interior, 0);
				for (std::size_t i = 0; i < b.interior.count; ++i) out.data()[i] = deferred.data()[i] + work.work[i];
				interval.deferred.handle().commit(b.interior, 0, out);
				continue;
			}
			auto const bank = levels.at(b.location.level).bank;
			auto input = directory.hydro.read(b.interior, bank).get();
			auto output = directory.hydro.output(b.interior, bank);
			auto const rho = frame.density.handle().read(b.interior, 0).get();
			auto const force = conventional ? frame.force.handle() : frame.shell.handle();
			auto const oldForce = force.read(b.interior, 0).get(), newForce = force.read(b.interior, 1).get();
			auto const impulse = frame.impulse.handle().read(b.interior, 0).get();
			auto const applied = frame.appliedWork.handle().read(b.interior, 0).get();
			for (std::size_t i = 0; i < b.interior.count; ++i) {
				auto value = input.at(i);
				if (!conventional || b.location.level == frame.level) for (int d = 0; d < ndim; ++d) {
					auto const dp = (frame.duration / 2.0) * (rho.data()[i] * oldForce.at(i).acceleration(d)
						+ value.density() * newForce.at(i).acceleration(d)) - impulse.at(i).momentum(d);
					if (config.gravity.energyTreatment == "naive") value.totalEnergy() += dp * (value.momentum(d) + 0.5 * dp) / value.density();
					value.momentum(d) += dp;
				}
				if (config.gravity.energyTreatment == "mullen") value.totalEnergy() += work.work[i] - applied.data()[i]
					+ (b.location.level == frame.level ? deferred.data()[i] : units::EnergyDensity{});
				hydro::HydroSystem const gas(config.hydro);
				gas.synchronize(value);
				if (!gas.admissible(value)) throw std::runtime_error("Gravity shell correction produced an inadmissible state");
				output.put(i, value);
			}
			directory.hydro.commit(b.interior, bank, output);
		}
	}

	void advanceGravityLevel(GravityInterval& interval, std::vector<int> const& occupied, std::size_t index,
		units::Time begin, units::Time duration, bool root = false) {
		int const level = occupied.at(index);
		phase(Operation::ResetFlux, {}, level, begin);
		Real fraction = root ? Real(1) : Real(1) / Real(std::uint64_t(1) << (level - occupied.at(index - 1)));
		Real elapsed = 0;
		while (elapsed < 1) {
			auto const now = begin + elapsed * duration;
			auto const limit = phase(Operation::Timestep, {}, level, now);
			for (int d = 0; d < ndim; ++d) signalSpeed[d] = std::max(signalSpeed[d], limit.signalSpeed[d]);
			while (fraction * duration > limit.timestep || fraction > 1 - elapsed) fraction /= 2;
			auto const step = fraction * duration;
			if (!(step > units::Time{}) || now + step == now) throw std::runtime_error("Gravity level timestep cannot advance time");
			GravityFrame frame(level, step, interval.cells, interval.faces, interval.store);
			int const nextLevel = index + 1 < occupied.size() ? occupied[index + 1] : -1;
			shellGravity(interval, frame, 0, nextLevel, now);
			if (config.gravity.timeIntegration == "conventional") partialGravity(interval, frame.force.handle(), 0, 0, now, true);
			for (auto const& b : topology->blocks()) {
				auto const old = fields->directory().hydro.read(b.interior, levels.at(b.location.level).bank).get();
				auto rho = frame.density.handle().output(b.interior, 0);
				auto impulse = frame.impulse.handle().output(b.interior, 0);
				auto work = frame.appliedWork.handle().output(b.interior, 0);
				auto flux = frame.flux.handle().output(b.massFlux, 0);
				for (std::size_t i = 0; i < b.interior.count; ++i) { rho.data()[i] = old.at(i).density(); impulse.put(i, {}); work.data()[i] = {}; }
				std::fill_n(flux.data(), b.massFlux.count, units::MassFlux{});
				frame.density.handle().commit(b.interior, 0, rho); frame.impulse.handle().commit(b.interior, 0, impulse);
				frame.appliedWork.handle().commit(b.interior, 0, work); frame.flux.handle().commit(b.massFlux, 0, flux);
				if (b.location.level == level) {
					auto pending = interval.deferred.handle().output(b.interior, 0);
					std::fill_n(pending.data(), b.interior.count, units::EnergyDensity{});
					interval.deferred.handle().commit(b.interior, 0, pending);
				}
			}
			interval.stack.push_back(&frame);
			assemblePredictor(interval, frame);
			for (std::size_t j = index; j < occupied.size(); ++j)
				phase(Operation::Probe, {}, occupied[j], now, 0, {interval.rate.handle(), interval.duration});
			gravitySourceRate(interval, level);
			auto& state = levels.at(level);
			state.begin = now; state.end = now + step;
			interval.boundary += phase(Operation::Advance, step, level, now, step / duration,
				{interval.rate.handle(), interval.duration}).boundary;
			provisionalGravity(interval, frame, step);
			// The halo forecast uses the accepted initial numerical RHS. The raw
			// coarse transport endpoint has not been refluxed and is only O(H)
			// accurate at a refinement boundary, even with a midpoint flux.
			for (auto const& b : topology->blocks()) if (b.location.level == level) {
				auto old = fields->directory().hydro.read(b.interior, state.bank).get();
				auto rate = interval.rate.handle().read(b.interior, 0).get();
				auto out = fields->directory().hydro.output(b.interior, 3);
				std::vector<units::Density> density(b.interior.count);
				for (std::size_t i = 0; i < b.interior.count; ++i) {
					auto value = old.at(i) + (step / interval.duration) * rate.at(i);
					predictorKineticRemainder(old.at(i), value);
					if (!hydro::HydroSystem(config.hydro).admissible(value)) throw std::runtime_error("Coarse gravity halo predictor is inadmissible");
					density[i] = value.density(); out.put(i, value);
				}
				for (auto const& species : fields->directory().species) {
					auto before = species.read(b.interior, state.bank).get();
					auto predicted = species.output(b.interior, 3);
					for (std::size_t i = 0; i < b.interior.count; ++i) predicted.data()[i] = before.data()[i] * (density[i] / old.at(i).density());
					species.commit(b.interior, 3, predicted);
				}
				fields->directory().hydro.commit(b.interior, 3, out);
				if (build::radiation && config.radiationEnabled()) copyFields(fields->directory().radiation, b.interior, state.bank ^ 1, 3);
			}
			state.predictorBank = 3;
			state.pending = true;
			if (nextLevel >= 0) advanceGravityLevel(interval, occupied, index + 1, now, step);
			phase(Operation::RefluxTracked, step, level, now + step);
			state.pending = false; state.bank ^= 1;
			shellGravity(interval, frame, 1, nextLevel, now + step);
			if (config.gravity.timeIntegration == "conventional") partialGravity(interval, frame.force.handle(), 1, 0, now + step, true);
			closeGravityFrame(interval, frame);
			interval.stack.pop_back();
			// A physical full field is retained for the next local CFL estimate.
			for (auto const& b : topology->blocks()) if (b.location.level >= level) {
				auto g = frame.nested.handle().read(b.interior, 1).get();
				std::vector<storage::Columns<gravity::State>> ancestors;
				for (auto const* f : interval.stack) ancestors.push_back(f->shell.handle().read(b.interior, 0).get());
				auto out = fields->directory().gravity.output(b.interior, levels.at(b.location.level).bank);
				for (std::size_t i = 0; i < b.interior.count; ++i) {
					auto value = g.at(i); for (auto const& a : ancestors) value += a.at(i); out.put(i, value);
				}
				fields->directory().gravity.commit(b.interior, levels.at(b.location.level).bank, out);
			}
			++statistics.levelSteps.at(level);
			elapsed += fraction;
		}
	}
#endif

	BoundaryTransport advanceLevel(std::vector<int> const& occupied, std::size_t index,
		units::Time begin, units::Time duration, bool root = false) {
		int const level = occupied.at(index);
		phase(Operation::ResetFlux, {}, level, begin);
		Real fraction = root ? Real(1) : Real(1) / Real(std::uint64_t(1) << (level - occupied.at(index - 1)));
		Real elapsed = 0;
		BoundaryTransport transported;
		while (elapsed < 1) {
			auto const now = begin + elapsed * duration;
			auto const limit = phase(Operation::Timestep, {}, level, now);
			for (int d = 0; d < ndim; ++d) signalSpeed[d] = std::max(signalSpeed[d], limit.signalSpeed[d]);
			while (fraction * duration > limit.timestep || fraction > 1 - elapsed) fraction /= 2;
			auto const step = fraction * duration;
			if (!(step > units::Time{}) || now + step == now) throw std::runtime_error("Level timestep cannot advance time");
			auto& state = levels.at(level);
			state.begin = now;
			state.end = now + step;
			transported += phase(Operation::Advance, step, level, now, step / duration).boundary;
			state.pending = true;
			if (index + 1 < occupied.size()) transported += advanceLevel(occupied, index + 1, now, step);
			// Fine registers now cover precisely this coarse step. Reflux before
			// publishing the coarse endpoint and resetting the child's registers.
			phase(Operation::Reflux, step, level, now + step);
			state.pending = false;
			state.bank ^= 1;
			++statistics.levelSteps.at(level);
			elapsed += fraction;
		}
		return transported;
	}
};

Runtime::Runtime(Config const& config, refinement::Criteria additionalCriteria)
  : impl_(std::make_unique<Impl>()) {
	config.validate();
	impl_->config = config;
	impl_->criteria = refinement::makeCriteria(config);
	for (auto& criterion : additionalCriteria)
		impl_->criteria.push_back(std::move(criterion));
#ifdef OCTOTIGERII_WITH_HPX
	impl_->localities = hpx::find_all_localities();
#else
	impl_->localities = {0};
#endif
	std::optional<amr::InitialMesh> startup;
	if (config.amr.enabled) {
		startup = amr::initializeMesh(config, impl_->criteria);
		impl_->topology = std::make_unique<CartesianTopology>(config, impl_->localities.size(), startup->leaves);
	} else
		impl_->topology = std::make_unique<CartesianTopology>(config, impl_->localities.size());
	impl_->fields = std::make_unique<FieldRepository>(config, impl_->topology->storageLayout(), impl_->localities);
#ifdef OCTOTIGERII_WITH_HPX
	std::vector<hpx::future<hpx::id_type>> pending;
	for (std::size_t i = 0; i < impl_->localities.size(); ++i)
		pending.push_back(hpx::new_<LocalExecutor>(impl_->localities[i], config, impl_->topology->blocks(), impl_->fields->directory(), i));
	impl_->executors = collect(pending);
	std::vector<hpx::future<void>> initialization;
	for (auto const& id : impl_->executors)
		initialization.push_back(hpx::async<LocalExecutor::InitializeAction>(id));
	finish(initialization);
#else
	impl_->executor = std::make_unique<LocalExecutor>(config, impl_->topology->blocks(), impl_->fields->directory(), 0);
	impl_->executor->initialize();
#endif
	if (startup) {
		impl_->shadow = std::make_unique<amr::Hierarchy>(config, impl_->exportSnapshots());
		impl_->regridInitialized = true;
		impl_->signalSpeed = startup->signalSpeed;
		for (int d = 0; d < ndim; ++d)
			impl_->travelBudget[d] = config.amr.signalBuffer * startup->signalSpeed[d] * startup->timestep * Real(config.amr.regridEvery);
	}
}

Runtime::~Runtime() = default;

std::size_t Runtime::size() const {
	return impl_->topology->blocks().size();
}

std::uint64_t Runtime::generation() const {
	std::lock_guard guard(impl_->apiMutex);
	return impl_->generation;
}

SchedulingStatistics Runtime::statistics() const {
	std::lock_guard guard(impl_->apiMutex);
	return impl_->statistics;
}

char const* Runtime::backend() {
#ifdef OCTOTIGERII_WITH_HPX
	return "HPX distributed fields";
#else
	return "serial fields";
#endif
}

std::vector<Snapshot> Runtime::snapshots() const {
	profiling::Elapsed profile("runtime.snapshots.wall_ns");
	std::lock_guard guard(impl_->apiMutex);
	auto result = impl_->exportSnapshots();
	// The current adapter and contiguous placement preserve directory order.
	if (result.size() != size()) throw std::logic_error("Incomplete snapshot directory");
	return result;
}

units::Time Runtime::stableTimestep() const {
	profiling::Elapsed profile("runtime.timestep.wall_ns");
	std::lock_guard guard(impl_->apiMutex);
	if (impl_->regridEnergyPending) throw std::logic_error("Solve gravity after regridding before computing a timestep");
	auto const result = impl_->phase(Operation::Timestep, {});
	impl_->signalSpeed = result.signalSpeed;
	return impl_->timeRefinement() ? result.levelTimestep.at(impl_->coarsestLevel()) : result.timestep;
}

void Runtime::advance(units::Time dt) {
	profiling::Elapsed profile("runtime.advance.wall_ns");
	std::lock_guard guard(impl_->apiMutex);
	if (impl_->regridEnergyPending) throw std::logic_error("Solve gravity after regridding before advancing");
	if (!(dt > units::Time{}) || !units::finite(dt) || impl_->time.time + dt == impl_->time.time) throw std::invalid_argument("Invalid step size");
	if (!impl_->config.hydroEnabled() && !impl_->config.radiationEnabled()) throw std::logic_error("No transport fields to advance");
	std::unique_ptr<amr::Hierarchy> nextShadow;
	if (impl_->shadow) {
		nextShadow = std::make_unique<amr::Hierarchy>(*impl_->shadow);
		nextShadow->advance(dt);
	}
	if (impl_->timeRefinement() && !impl_->config.gravityEnabled()) {
		auto const statistics = impl_->statistics;
		auto const speed = impl_->signalSpeed;
		std::vector<int> occupied;
		for (auto const& b : impl_->topology->blocks()) occupied.push_back(b.location.level);
		std::sort(occupied.begin(), occupied.end());
		occupied.erase(std::unique(occupied.begin(), occupied.end()), occupied.end());
		impl_->statistics.levelSteps.resize(occupied.back() + 1);
		impl_->phase(Operation::Backup, {});
		impl_->levels.assign(occupied.back() + 1, {impl_->bank, {}, {}, false});
		try {
			auto const transported = impl_->advanceLevel(occupied, 0, impl_->time.time, dt, true);
			impl_->phase(Operation::Normalize, {});
			auto nextTime = impl_->time;
			nextTime.completeStep(dt);
			if (nextShadow) nextShadow->refreshLeaves(impl_->exportSnapshots(impl_->bank ^ 1, nextTime));
			impl_->bank ^= 1;
			impl_->time = nextTime;
			++impl_->generation;
			impl_->boundary += transported;
			if (nextShadow) impl_->shadow = std::move(nextShadow);
			for (int d = 0; d < ndim; ++d) impl_->travel[d] += impl_->signalSpeed[d] * dt;
			impl_->levels.clear();
			return;
		} catch (...) {
			impl_->levels.clear();
			impl_->phase(Operation::Restore, {});
			impl_->statistics = statistics;
			impl_->signalSpeed = speed;
			throw;
		}
	}
	auto const transport = impl_->phase(Operation::Advance, dt);
	if (impl_->config.amr.enabled || impl_->config.hydroEnabled())
		impl_->phase(impl_->gravityEnergyActive && impl_->config.gravity.energyTreatment == "mullen" ? Operation::RefluxTracked : Operation::Reflux, dt);
	auto nextTime = impl_->time;
	nextTime.completeStep(dt);
	if (nextShadow) nextShadow->refreshLeaves(impl_->exportSnapshots(impl_->bank ^ 1, nextTime));
	impl_->bank ^= 1;
	++impl_->generation;
	impl_->time = nextTime;
	impl_->boundary += transport.boundary;
	if (nextShadow) impl_->shadow = std::move(nextShadow);
	for (int d = 0; d < ndim; ++d)
		impl_->travel[d] += impl_->signalSpeed[d] * dt;
}

gravity::Statistics Runtime::advanceGravity(units::Time dt) {
	if (!impl_->config.hydroEnabled() || !impl_->config.gravityEnabled()) throw std::logic_error("Coupled gravity requires self-gravitating gas");
	if (!(dt > units::Time{}) || !units::finite(dt) || impl_->time.time + dt == impl_->time.time) throw std::invalid_argument("Invalid step size");
	if (!impl_->config.amr.enabled || !impl_->config.timestep.refinement || impl_->config.hasExternalAcceleration()) {
		beginGravityEnergy(); kickGravity(dt / 2.0); advance(dt);
		auto work = solveGravity(); kickGravity(dt / 2.0); finishGravityEnergy(dt);
		return work;
	}
	std::lock_guard guard(impl_->apiMutex);
	if (impl_->regridEnergyPending || impl_->gravityEnergyActive || !impl_->gravityReady || impl_->gravityTime != impl_->time.time)
		throw std::logic_error("Coupled gravity requires a closed, synchronized gravity field");
#if OCTOTIGERII_GRAVITY
	profiling::Elapsed profile("runtime.gravity_time_integration.wall_ns");
	auto const savedStatistics = impl_->statistics;
	auto const savedSpeed = impl_->signalSpeed;
	std::vector<int> occupied;
	for (auto const& b : impl_->topology->blocks()) occupied.push_back(b.location.level);
	std::sort(occupied.begin(), occupied.end()); occupied.erase(std::unique(occupied.begin(), occupied.end()), occupied.end());
	impl_->statistics.levelSteps.resize(occupied.back() + 1);
	impl_->phase(Operation::Backup, {});
	impl_->levels.assign(occupied.back() + 1, {impl_->bank, {}, {}, false});
	try {
		Impl::GravityInterval interval(*impl_, dt);
		// Initialize all halo source rates before the first coarse forecast. Later
		// inactive rates are old by O(H), sufficient for a midpoint state to O(H²).
		for (auto const& b : impl_->topology->blocks()) {
			auto g = impl_->fields->directory().gravity.read(b.interior, impl_->bank).get();
			auto out = interval.predictor.handle().output(b.interior, 0);
			for (std::size_t i = 0; i < b.interior.count; ++i) out.put(i, g.at(i));
			interval.predictor.handle().commit(b.interior, 0, out);
		}
		impl_->phase(Operation::Probe, {}, -1, {}, 0, {interval.rate.handle(), interval.duration});
		impl_->gravitySourceRate(interval, -1);
		std::unique_ptr<amr::Hierarchy> nextShadow;
		if (impl_->shadow) {
			nextShadow = std::make_unique<amr::Hierarchy>(*impl_->shadow);
			nextShadow->advanceGravity(dt);
		}
		impl_->advanceGravityLevel(interval, occupied, 0, impl_->time.time, dt, true);
		impl_->phase(Operation::Normalize, {});
		auto nextTime = impl_->time; nextTime.completeStep(dt);
		if (nextShadow) {
			auto snapshots = impl_->exportSnapshots(impl_->bank ^ 1, nextTime);
			nextShadow->refreshLeaves(snapshots);
			nextShadow->refreshGravity(snapshots);
		}
		impl_->bank ^= 1; impl_->time = nextTime; impl_->gravityTime = nextTime.time;
		++impl_->generation;
		impl_->boundary += interval.boundary;
		if (nextShadow) impl_->shadow = std::move(nextShadow);
		for (int d = 0; d < ndim; ++d) impl_->travel[d] += impl_->signalSpeed[d] * dt;
		impl_->levels.clear();
		return interval.work;
	} catch (...) {
		impl_->levels.clear();
		impl_->phase(Operation::Restore, {});
		impl_->statistics = savedStatistics; impl_->signalSpeed = savedSpeed;
		throw;
	}
#else
	throw std::logic_error("Gravity is not part of this executable");
#endif
}

BoundaryTransport Runtime::boundaryTransport() const {
	std::lock_guard guard(impl_->apiMutex);
	return impl_->boundary;
}

std::size_t Runtime::shadowCellCount() const {
	std::lock_guard guard(impl_->apiMutex);
	return impl_->shadow ? impl_->shadow->size() : 0;
}

bool Runtime::regrid(units::Time nextStep, bool force) {
	std::lock_guard guard(impl_->apiMutex);
	auto const& config = impl_->config;
	if (!config.amr.enabled) return false;
	if (impl_->gravityEnergyActive) throw std::logic_error("Cannot regrid inside a gravity energy step");
	if (impl_->regridEnergyPending) throw std::logic_error("Solve gravity before regridding again");
	if (!(nextStep >= units::Time{}) || !units::finite(nextStep)) throw std::invalid_argument("Invalid AMR lookahead timestep");
	bool due = force || !impl_->regridInitialized || impl_->time.step - impl_->lastRegridStep >= std::uint64_t(config.amr.regridEvery);
	for (int d = 0; d < ndim; ++d)
		if (impl_->travel[d] + impl_->signalSpeed[d] * nextStep > impl_->travelBudget[d] * (1 + 64 * epsilonR)) due = true;
	if (!due) return false;
	profiling::Elapsed profile("amr.regrid.wall_ns");
	bool changed = false;
	auto const horizon = Real(config.amr.regridEvery) * nextStep;
	auto const original = impl_->exportSnapshots();
	amr::Hierarchy const source(config, original);
	auto candidate = original;
	std::vector<mesh::BlockLocation> leaves;
	std::unique_ptr<amr::Hierarchy> nextShadow;
	std::array<units::Velocity, ndim> bufferSpeed{};
	// Every transfer reads the same immutable pre-regrid state. Build the full
	// Morton-ordered destination before replacing any published field handles.
	for (int pass = 0; pass <= config.amr.maxLevel + 1; ++pass) {
		auto shadow = std::make_unique<amr::Hierarchy>(config, candidate);
		auto const& estimate = pass == 0 && impl_->shadow ? *impl_->shadow : *shadow;
		auto const selected = amr::selectMesh(config, candidate, estimate, horizon, impl_->criteria, pass == 0 && impl_->regridInitialized, shadow.get());
		if (pass == 0) bufferSpeed = selected.signalSpeed;
		if (!selected.changed) {
			nextShadow = std::move(shadow);
			break;
		}
		leaves = selected.leaves;
		candidate.clear();
		for (auto const& leaf : leaves)
			candidate.push_back(source.transfer(leaf));
		changed = true;
		if (pass == config.amr.maxLevel + 1) throw std::runtime_error("AMR refinement failed to settle");
	}
	if (changed) {
		if (config.hydroEnabled() && config.gravityEnabled() && config.gravity.conserveRegridEnergy &&
			(!impl_->gravityReady || impl_->gravityTime != impl_->time.time))
			throw std::logic_error("Energy-conserving regrid requires synchronized gravity on the old mesh");
		impl_->install(leaves, source);
	}
	impl_->shadow = std::move(nextShadow);
	impl_->lastRegridStep = impl_->time.step;
	impl_->regridInitialized = true;
	impl_->travel.fill({});
	for (int d = 0; d < ndim; ++d)
		impl_->travelBudget[d] = config.amr.signalBuffer * bufferSpeed[d] * horizon;
	return changed;
}

void Runtime::kickGravity(units::Time dt) {
	profiling::Elapsed profile("runtime.gravity_kick.wall_ns");
	std::lock_guard guard(impl_->apiMutex);
	if (!impl_->config.hydroEnabled() || (!impl_->config.gravityEnabled() && !impl_->config.hasExternalAcceleration()) ||
		(impl_->config.gravityEnabled() && (!impl_->gravityReady || impl_->gravityTime != impl_->time.time)) || !(dt > units::Time{}) || !units::finite(dt))
		throw std::logic_error("Gravity kick requires synchronized gas and gravity");
	std::unique_ptr<amr::Hierarchy> nextShadow;
	if (impl_->shadow) {
		nextShadow = std::make_unique<amr::Hierarchy>(*impl_->shadow);
		nextShadow->kick(dt);
	}
	impl_->phase(impl_->gravityEnergyActive && impl_->config.gravity.energyTreatment == "mullen" ? Operation::KickTracked : Operation::Kick, dt);
	impl_->bank ^= 1;
	++impl_->generation;
	if (nextShadow) impl_->shadow = std::move(nextShadow);
}

void Runtime::beginGravityEnergy() {
	std::lock_guard guard(impl_->apiMutex);
	if (impl_->gravityEnergyActive || !impl_->config.hydroEnabled() || !impl_->config.gravityEnabled() ||
		!impl_->gravityReady || impl_->gravityTime != impl_->time.time)
		throw std::logic_error("Conservative gravity step needs synchronized gas and gravity");
	impl_->phase(Operation::PrepareGravityEnergy, {});
	impl_->gravityEnergyStart = impl_->time.time;
	impl_->gravityEnergyActive = true;
}

void Runtime::finishGravityEnergy(units::Time dt) {
	std::lock_guard guard(impl_->apiMutex);
	if (!impl_->gravityEnergyActive || !(dt > units::Time{}) || !units::finite(dt) ||
		impl_->gravityEnergyStart + dt != impl_->time.time || impl_->gravityTime != impl_->time.time)
		throw std::logic_error("Conservative gravity work requires the completed transport interval and endpoint gravity");
	auto result = impl_->phase(Operation::FinishGravityEnergy, dt);
	std::unique_ptr<amr::Hierarchy> nextShadow;
	if (impl_->shadow) {
		nextShadow = std::make_unique<amr::Hierarchy>(*impl_->shadow);
		nextShadow->refreshLeaves(impl_->exportSnapshots(impl_->bank ^ 1));
	}
	impl_->bank ^= 1;
	++impl_->generation;
	impl_->boundary += result.boundary;
	impl_->gravityEnergyActive = false;
	if (nextShadow) impl_->shadow = std::move(nextShadow);
}

gravity::Statistics Runtime::solveGravity() {
	if (!impl_->config.gravityEnabled()) throw std::logic_error("Selected problem does not enable gravity");
	profiling::Elapsed profile("runtime.gravity.wall_ns");
	std::lock_guard guard(impl_->apiMutex);
#if OCTOTIGERII_GRAVITY
	if (!impl_->gravitySolver)
		impl_->gravitySolver = std::make_unique<gravity::FieldSolver>(impl_->config, impl_->topology->blocks(), impl_->fields->directory(), impl_->localities);
	auto result = impl_->gravitySolver->solve(impl_->bank);
	impl_->recoverRegridEnergy(impl_->bank ^ 1);
	std::unique_ptr<amr::Hierarchy> nextShadow;
	if (impl_->shadow) {
		if (impl_->regridEnergyPending)
			nextShadow = std::make_unique<amr::Hierarchy>(impl_->config, impl_->exportSnapshots(impl_->bank ^ 1));
		else {
			nextShadow = std::make_unique<amr::Hierarchy>(*impl_->shadow);
			nextShadow->refreshGravity(impl_->exportSnapshots(impl_->bank ^ 1));
		}
	}
	impl_->bank ^= 1;
	++impl_->generation;
	impl_->gravityTime = impl_->time.time;
	impl_->gravityReady = true;
	impl_->regridEnergyPending = false;
	if (nextShadow) impl_->shadow = std::move(nextShadow);
	return result;
#else
	throw std::logic_error("Gravity is not part of this executable");
#endif
}

void Runtime::setGravity(std::vector<std::vector<gravity::State>> const& fields) {
	std::lock_guard guard(impl_->apiMutex);
	if (!impl_->config.gravityEnabled() || fields.size() != size()) throw std::invalid_argument("Gravity directory size/field mismatch");
	for (std::size_t b = 0; b < size(); ++b) {
		if (fields[b].size() != impl_->topology->blocks()[b].interior.count) throw std::invalid_argument("Gravity block size mismatch");
		for (auto const& value : fields[b])
			if (!finite(value)) throw std::invalid_argument("Nonfinite gravity assignment");
	}
	// This is a synchronized source-field publication, with no live readers.
	// Fill the other bank first so a failed transfer cannot corrupt published data.
	auto const& directory = impl_->fields->directory();
	for (std::size_t b = 0; b < size(); ++b) {
		auto const range = impl_->topology->blocks()[b].interior;
		if (impl_->config.hydroEnabled())
			copyFields(directory.hydro, range, impl_->bank);
		else {
			auto input = directory.density.read(range, impl_->bank).get();
			auto output = directory.density.output(range, impl_->bank ^ 1);
			std::copy_n(input.data(), input.size(), output.data());
			directory.density.commit(range, impl_->bank ^ 1, output);
		}
		if (build::radiation && impl_->config.radiationEnabled()) copyFields(directory.radiation, range, impl_->bank);
		for (auto const& field : directory.species) copyScalar(field, range, impl_->bank);
		auto output = directory.gravity.output(range, impl_->bank ^ 1);
		for (std::size_t i = 0; i < range.count; ++i)
			output.put(i, fields[b][i]);
		directory.gravity.commit(range, impl_->bank ^ 1, output);
	}
	impl_->recoverRegridEnergy(impl_->bank ^ 1);
	std::unique_ptr<amr::Hierarchy> nextShadow;
	if (impl_->shadow) {
		if (impl_->regridEnergyPending)
			nextShadow = std::make_unique<amr::Hierarchy>(impl_->config, impl_->exportSnapshots(impl_->bank ^ 1));
		else {
			nextShadow = std::make_unique<amr::Hierarchy>(*impl_->shadow);
			nextShadow->refreshGravity(impl_->exportSnapshots(impl_->bank ^ 1));
		}
	}
	impl_->bank ^= 1;
	++impl_->generation;
	impl_->gravityTime = impl_->time.time;
	impl_->gravityReady = true;
	impl_->regridEnergyPending = false;
	if (nextShadow) impl_->shadow = std::move(nextShadow);
}

}	 // namespace octotigerII
