#include "octotigerII/runtime.hpp"
#include <algorithm>
#include <exception>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include "octotigerII/amr/hierarchy.hpp"
#include "octotigerII/problems.hpp"
#include "octotigerII/profiling.hpp"
#include "octotigerII/storage/registry.hpp"
#include "octotigerII/subgrid/view.hpp"
#include "octotigerII/verification/analytic.hpp"
#if OCTOTIGERII_GRAVITY
#include "octotigerII/gravity/fieldSolver.hpp"
#endif

#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/include/components.hpp>
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
	void copyFields(storage::ColumnHandle<State> const& field, storage::Range range, unsigned bank) {
		auto input = field.read(range, bank).get();
		auto output = field.output(range, bank ^ 1);
		for (std::size_t i = 0; i < range.count; ++i)
			output.put(i, input.at(i));
		field.commit(range, bank ^ 1, output);
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
			bool amr) {
			auto interior = fields.read(block.interior, bank).get();
			readHalo(fields, plan, bank, ghosts);
			applyHaloBoundaries(plan, ghosts, system, time, analytic);
			PatchView<State> input(block, std::move(interior), plan, ghosts);
			auto output = fields.output(block.interior, bank ^ 1);
			{
				profiling::Region profile(std::is_same_v<System, hydro::HydroSystem> ? "hydro.advance" : "radiation.advance");
				Solver(system).advanceInto(
					input, dt, work, [&](mesh::Coordinates const& cell, State const& state) { output.put(block.layout.index(cell), state); });
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
	};

	class Workspace {
	public:
		SolverWorkspace<hydro::HydroSystem> hydro;
		SolverWorkspace<radiation::RadiationSystem> radiation;
	};

	class PhaseResult {
	public:
		SchedulingStatistics tasks;
		units::Time timestep = units::Time::from_value(std::numeric_limits<Real>::infinity());
		std::array<units::Velocity, ndim> signalSpeed{};

		template <typename Archive>
		void serialize(Archive& archive, unsigned) {
			archive & tasks & timestep & signalSpeed;
		}
	};

	enum class Operation
	{
		Timestep,
		Advance,
		Reflux,
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
			if constexpr (build::hydro) {
				hydroBoundary_ = [evaluator, gas = hydro::HydroSystem(config_.hydro.gamma)](
									 auto const& position, auto time) { return gas.conservedState(evaluator(position, time).hydro); };
			}
			if constexpr (build::radiation) {
				radiationBoundary_ = [evaluator](auto const& position, auto time) { return evaluator(position, time).radiation; };
			}
		}
		for (auto const& block : blocks_)
			if (block.interior.partition == owner_) owned_.push_back(block.id);
		for (auto id : owned_)
			plans_.emplace(id, makeHaloPlan(config_, blocks_, id));
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
			if constexpr (build::hydro) initializeFields(fields_.hydro, block.interior, initial.hydro);
			if constexpr (build::radiation) initializeFields(fields_.radiation, block.interior, initial.radiation);
			if constexpr (build::gravity) {
				initializeFields(fields_.gravity, block.interior, initial.gravity);
				if (!config_.hydroEnabled()) {
					auto output = fields_.density.output(block.interior, 0);
					std::copy(initial.density.values().begin(), initial.density.values().end(), output.data());
					fields_.density.commit(block.interior, 0, output);
				}
			}
		}
	}

	void begin(std::uint64_t generation, unsigned bank, units::Time time) {
		std::lock_guard guard(queueMutex_);
		generation_ = generation;
		bank_ = bank;
		time_ = time;
		next_ = 0;
	}

	std::vector<std::uint64_t> claim(std::uint64_t generation) {
		std::lock_guard guard(queueMutex_);
		if (generation != generation_) throw std::logic_error("Stale work request");
		if (next_ == owned_.size()) return {};
		return {owned_[next_++]};
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
			result.tasks.localTasks += part.tasks.localTasks;
			result.tasks.stolenTasks += part.tasks.stolenTasks;
			result.timestep = std::min(result.timestep, part.timestep);
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
			if constexpr (build::hydro) exportOne(fields_.hydro, snapshot.hydro);
			if constexpr (build::radiation) exportOne(fields_.radiation, snapshot.radiation);
			if constexpr (build::gravity) {
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
	std::vector<std::uint64_t> owned_;
	std::map<std::uint64_t, HaloPlan> plans_;
	std::map<std::uint64_t, std::vector<FluxCorrection>> refluxPlans_;
	std::vector<Workspace> workspaces_;
	std::mutex queueMutex_;
	std::size_t next_ = 0;
	std::uint64_t generation_ = 0;
	unsigned bank_ = 0;
	units::Time time_{};

	PhaseResult worker(Operation operation, units::Time dt, std::uint64_t generation, std::vector<storage::Locality> const& executors, Workspace& workspace) {
		PhaseResult result;
		auto process = [&](std::uint64_t id, bool stolen) {
			auto const& block = blocks_.at(id);
			if (operation == Operation::Timestep)
				result.timestep = std::min(result.timestep, timestep(block, result.signalSpeed));
			else if (operation == Operation::Advance) {
				// Owned plans are immutable. Stolen plans have only metadata and
				// are temporary, keeping persistent topology caches bounded.
				std::optional<HaloPlan> temporary;
				if (stolen) temporary = makeHaloPlan(config_, blocks_, id);
				auto const& plan = stolen ? *temporary : plans_.at(id);
				if constexpr (build::hydro)
					workspace.hydro.advance(block, fields_.hydro, plan, hydro::HydroSystem(config_.hydro.gamma), bank_, dt, time_, hydroBoundary_,
						fields_.hydroFlux, config_.amr.enabled);
				if constexpr (build::radiation)
					workspace.radiation.advance(block, fields_.radiation, plan, radiation::RadiationSystem(config_.radiation.lightSpeedRatio * constants::c),
						bank_, dt, time_, radiationBoundary_, fields_.radiationFlux, config_.amr.enabled);
				if constexpr (build::gravity) copyFields(fields_.gravity, block.interior, bank_);
			} else if (operation == Operation::Reflux) {
				auto const plan = stolen ? makeRefluxPlan(config_, blocks_, id) : refluxPlans_.at(id);
				if constexpr (build::hydro) reflux(block, plan, fields_.hydro, fields_.hydroFlux, hydro::HydroSystem(config_.hydro.gamma), dt);
				if constexpr (build::radiation)
					reflux(block, plan, fields_.radiation, fields_.radiationFlux, radiation::RadiationSystem(config_.radiation.lightSpeedRatio * constants::c),
						dt);
			} else if constexpr (build::hydro) {
				kick(block, dt);
				if constexpr (build::gravity) copyFields(fields_.gravity, block.interior, bank_);
				if constexpr (build::radiation) copyFields(fields_.radiation, block.interior, bank_);
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
		if constexpr (build::hydro) {
			for (int d = 0; d < ndim; ++d)
				maximumAcceleration[d] = units::abs(config_.hydro.acceleration[d]);
			result = fieldStep(fields_.hydro, hydro::HydroSystem(config_.hydro.gamma));
			units::Acceleration acceleration{};
			for (auto component : config_.hydro.acceleration)
				acceleration += units::abs(component);
			if constexpr (build::gravity) {
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
		if constexpr (build::radiation)
			result = std::min(result, fieldStep(fields_.radiation, radiation::RadiationSystem(config_.radiation.lightSpeedRatio * constants::c)));
		if constexpr (build::radiation)
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
		storage::ColumnHandle<typename System::Flux> const& fluxFields, System const& system, units::Time dt) {
		if (plan.empty()) return;
		using State = typename System::State;
		auto coarseFlux = fluxFields.read(block.boundaryFlux, 0).get();
		auto current = fields.read(block.interior, bank_ ^ 1).get();
		std::vector<State> corrections(block.interior.count);
		for (auto const& face : plan) {
			typename System::Flux fine{};
			std::vector<storage::PendingColumns<typename System::Flux>> pending;
			for (auto const& range : face.fineFluxes)
				pending.push_back(fluxFields.read(range, 0));
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
			auto candidate = system.correctRoundoff(current.at(i) + corrections[i], componentAbs(corrections[i]));
			if (!system.admissible(candidate)) throw std::runtime_error("AMR reflux produced an inadmissible state");
			output.put(i, candidate);
		}
		fields.commit(block.interior, bank_ ^ 1, output);
	}

	void kick(Subgrid const& block, units::Time dt) {
		auto input = fields_.hydro.read(block.interior, bank_).get();
		std::optional<storage::Columns<gravity::State>> gravity;
		if constexpr (build::gravity) gravity = fields_.gravity.read(block.interior, bank_).get();
		auto output = fields_.hydro.output(block.interior, bank_ ^ 1);
		{
			profiling::Region profile("gravity.kick");
			for (std::size_t i = 0; i < block.interior.count; ++i) {
				auto state = input.at(i);
				units::EnergyDensity work{};
				for (int d = 0; d < ndim; ++d) {
					auto const old = state.momentum(d);
					auto acceleration = config_.hydro.acceleration[d];
					if constexpr (build::gravity) acceleration += gravity->at(i).acceleration(d);
					auto const impulse = dt * state.density() * acceleration;
					state.momentum(d) += impulse;
					work += impulse * (old + 0.5 * impulse) / state.density();
				}
				state.totalEnergy() += work;
				if (!hydro::HydroSystem(config_.hydro.gamma).admissible(state)) throw std::runtime_error("Invalid gravity kick state");
				output.put(i, state);
			}
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
	bool gravityReady = false;
	SchedulingStatistics statistics;
	refinement::Criteria criteria;
	std::unique_ptr<amr::Hierarchy> shadow;
	std::uint64_t lastRegridStep = 0;
	std::array<units::Velocity, ndim> signalSpeed{};
	std::array<units::Length, ndim> travel{}, travelBudget{};
	bool regridInitialized = false;

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
		auto nextTopology = std::make_unique<CartesianTopology>(config, localities.size(), leaves);
		auto nextFields = std::make_unique<FieldRepository>(config, nextTopology->storageLayout(), localities);
		auto const& directory = nextFields->directory();
		for (auto const& block : nextTopology->blocks()) {
			auto const snapshot = source.transfer(block.location);
			if constexpr (build::hydro) initializeFields(directory.hydro, block.interior, snapshot.hydro);
			if constexpr (build::radiation) initializeFields(directory.radiation, block.interior, snapshot.radiation);
			if constexpr (build::gravity) {
				initializeFields(directory.gravity, block.interior, snapshot.gravity);
				if constexpr (!build::hydro) {
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
	}

	/// Drain all locality results and require exactly one completed task per block.
	/// The caller may publish a bank only after this function succeeds.
	PhaseResult phase(Operation operation, units::Time dt) {
		++dispatch;
#ifdef OCTOTIGERII_WITH_HPX
		std::vector<hpx::future<void>> starts;
		for (auto const& id : executors)
			starts.push_back(hpx::async<LocalExecutor::BeginAction>(id, dispatch, bank, time.time));
		finish(starts);
		std::vector<hpx::future<PhaseResult>> pending;
		for (auto const& id : executors)
			pending.push_back(hpx::async<LocalExecutor::RunAction>(id, operation, dt, dispatch, executors));
		auto results = collect(pending);
#else
		executor->begin(dispatch, bank, time.time);
		std::vector<PhaseResult> results{executor->run(operation, dt, dispatch, {})};
#endif
		PhaseResult result;
		for (auto const& part : results) {
			result.tasks.localTasks += part.tasks.localTasks;
			result.tasks.stolenTasks += part.tasks.stolenTasks;
			result.timestep = std::min(result.timestep, part.timestep);
			for (int d = 0; d < ndim; ++d)
				result.signalSpeed[d] = std::max(result.signalSpeed[d], part.signalSpeed[d]);
		}
		if (result.tasks.localTasks + result.tasks.stolenTasks != topology->blocks().size())
			throw std::logic_error("Stage did not complete every output range");
		statistics.localTasks += result.tasks.localTasks;
		statistics.stolenTasks += result.tasks.stolenTasks;
		return result;
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
	if (config.amr.enabled) {
		auto dt = stableTimestep();
		if (!units::finite(dt)) dt = {};
		regrid(dt, true);
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
	auto const result = impl_->phase(Operation::Timestep, {});
	impl_->signalSpeed = result.signalSpeed;
	return result.timestep;
}

void Runtime::advance(units::Time dt) {
	profiling::Elapsed profile("runtime.advance.wall_ns");
	std::lock_guard guard(impl_->apiMutex);
	if (!(dt > units::Time{}) || !units::finite(dt) || impl_->time.time + dt == impl_->time.time) throw std::invalid_argument("Invalid step size");
	if (!impl_->config.hydroEnabled() && !impl_->config.radiationEnabled()) throw std::logic_error("No transport fields to advance");
	std::unique_ptr<amr::Hierarchy> nextShadow;
	if (impl_->shadow) {
		nextShadow = std::make_unique<amr::Hierarchy>(*impl_->shadow);
		nextShadow->advance(dt);
	}
	impl_->phase(Operation::Advance, dt);
	if (impl_->config.amr.enabled) impl_->phase(Operation::Reflux, dt);
	auto nextTime = impl_->time;
	nextTime.completeStep(dt);
	if (nextShadow) nextShadow->refreshLeaves(impl_->exportSnapshots(impl_->bank ^ 1, nextTime));
	impl_->bank ^= 1;
	++impl_->generation;
	impl_->time = nextTime;
	if (nextShadow) impl_->shadow = std::move(nextShadow);
	for (int d = 0; d < ndim; ++d)
		impl_->travel[d] += impl_->signalSpeed[d] * dt;
}

std::size_t Runtime::shadowCellCount() const {
	std::lock_guard guard(impl_->apiMutex);
	return impl_->shadow ? impl_->shadow->size() : 0;
}

bool Runtime::regrid(units::Time nextStep, bool force) {
	std::lock_guard guard(impl_->apiMutex);
	auto const& config = impl_->config;
	if (!config.amr.enabled) return false;
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
	if (changed) impl_->install(leaves, source);
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
	if (!impl_->config.hydroEnabled() || (!build::gravity && !impl_->config.hasExternalAcceleration()) ||
		(build::gravity && (!impl_->gravityReady || impl_->gravityTime != impl_->time.time)) || !(dt > units::Time{}) || !units::finite(dt))
		throw std::logic_error("Gravity kick requires synchronized gas and gravity");
	std::unique_ptr<amr::Hierarchy> nextShadow;
	if (impl_->shadow) {
		nextShadow = std::make_unique<amr::Hierarchy>(*impl_->shadow);
		nextShadow->kick(dt);
	}
	impl_->phase(Operation::Kick, dt);
	impl_->bank ^= 1;
	++impl_->generation;
	if (nextShadow) impl_->shadow = std::move(nextShadow);
}

gravity::Statistics Runtime::solveGravity() {
	profiling::Elapsed profile("runtime.gravity.wall_ns");
	std::lock_guard guard(impl_->apiMutex);
#if OCTOTIGERII_GRAVITY
	if (!impl_->gravitySolver)
		impl_->gravitySolver = std::make_unique<gravity::FieldSolver>(impl_->config, impl_->topology->blocks(), impl_->fields->directory(), impl_->localities);
	auto result = impl_->gravitySolver->solve(impl_->bank);
	std::unique_ptr<amr::Hierarchy> nextShadow;
	if (impl_->shadow) {
		nextShadow = std::make_unique<amr::Hierarchy>(*impl_->shadow);
		nextShadow->refreshGravity(impl_->exportSnapshots(impl_->bank ^ 1));
	}
	impl_->bank ^= 1;
	++impl_->generation;
	impl_->gravityTime = impl_->time.time;
	impl_->gravityReady = true;
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
		if constexpr (build::radiation) copyFields(directory.radiation, range, impl_->bank);
		auto output = directory.gravity.output(range, impl_->bank ^ 1);
		for (std::size_t i = 0; i < range.count; ++i)
			output.put(i, fields[b][i]);
		directory.gravity.commit(range, impl_->bank ^ 1, output);
	}
	std::unique_ptr<amr::Hierarchy> nextShadow;
	if (impl_->shadow) {
		nextShadow = std::make_unique<amr::Hierarchy>(*impl_->shadow);
		nextShadow->refreshGravity(impl_->exportSnapshots(impl_->bank ^ 1));
	}
	impl_->bank ^= 1;
	++impl_->generation;
	impl_->gravityTime = impl_->time.time;
	impl_->gravityReady = true;
	if (nextShadow) impl_->shadow = std::move(nextShadow);
}

}	 // namespace octotigerII
