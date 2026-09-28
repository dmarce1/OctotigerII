/** @file
 * @brief Per-locality worker queues, halo assembly, and block operations.
 */
#include "internal.hpp"

namespace octotigerII {

template <typename State>
void LocalExecutor::fluxRegister(storage::ColumnHandle<State> const& field, storage::Range range, bool reset) {
	auto output = field.output(range, 1);
	if (reset) {
		for (std::size_t i = 0; i < range.count; ++i) output.put(i, State{});
	} else {
		auto old = field.read(range, 1).get(), raw = field.read(range, 0).get();
		for (std::size_t i = 0; i < range.count; ++i) output.put(i, old.at(i) + fluxWeight_ * raw.at(i));
	}
	field.commit(range, 1, output);
}

template <typename System>
void LocalExecutor::reflux(Subgrid const& block, std::vector<FluxCorrection> const& plan, storage::ColumnHandle<typename System::State> const& fields,
	storage::ColumnHandle<typename System::Flux> const& fluxFields, System const& system, units::Time dt, bool synchronize) {
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

LocalExecutor::LocalExecutor(Config config, std::vector<Subgrid> blocks, FieldDirectory fields, std::size_t owner)
  : config_(std::move(config))
  , blocks_(std::move(blocks))
  , fields_(std::move(fields))
  , owner_(owner) {
	radiationMaterial_ = problemRadiationMaterial(config_);
	if (!radiationMaterial_) throw std::invalid_argument("Problem radiation material factory returned no evaluator");
	if (config_.mesh.boundary.contains(physics::BoundaryCondition::Analytic)) {
		auto evaluator = problemBoundary(config_);
		if (build::hydro && config_.hydroEnabled()) {
			hydroBoundary_ = [evaluator, gas = hydro::HydroSystem(config_)](
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
	blockCaches_.resize(blocks_.size());
	if (config_.hydroEnabled() && (config_.gravityEnabled() ||
		radiation::radiationCouplingEnabled(config_) || problemHasRadiationMaterial(config_))) {
		bool const coupled = config_.radiationEnabled() &&
			(radiation::radiationCouplingEnabled(config_) || problemHasRadiationMaterial(config_));
		auto const shell = mesh::MeshLayout(config_.mesh.cells, 1).cellCount() -
			mesh::MeshLayout(config_.mesh.cells, 0).cellCount();
		if (coupled) gasCacheCapacity_ = radiationCacheCapacity_ = shell;
		else for (auto const& [id, plan] : plans_)
			gasCacheCapacity_ = std::max(gasCacheCapacity_, plan.valueCount ? plan.valueCount : plan.ghostCount);
		cachePool_.reserve(owned_.size() + workers);
		for (std::size_t i = 0; i < owned_.size() + workers; ++i) cachePool_.push_back(makeCache());
	}
}

std::unique_ptr<BlockHaloCache> LocalExecutor::makeCache() const {
	auto cache = std::make_unique<BlockHaloCache>();
	cache->gas.values.reserve(gasCacheCapacity_);
	cache->radiation.values.reserve(radiationCacheCapacity_);
	cache->donors.reserve(config_.amr.maxLevel + 1);
	return cache;
}

std::unique_ptr<BlockHaloCache> LocalExecutor::acquireCache() {
	{
		std::lock_guard guard(queueMutex_);
		if (!cachePool_.empty()) {
			auto cache = std::move(cachePool_.back());
			cachePool_.pop_back();
			return cache;
		}
	}
	// A locality can steal more probes than its initial reserve. Retain the
	// additional buffer after use so subsequent steps reuse that capacity.
	return makeCache();
}

void LocalExecutor::releaseCache(std::unique_ptr<BlockHaloCache>& cache) {
	if (!cache) return;
	std::lock_guard guard(queueMutex_);
	cachePool_.push_back(std::move(cache));
}

void LocalExecutor::initialize() {
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

void LocalExecutor::begin(std::uint64_t generation, unsigned bank, units::Time time, int level, std::vector<HaloTime> times, Real fluxWeight,
	SourcePredictor source, Operation operation, std::vector<std::uint64_t> cachedWork) {
	std::lock_guard guard(queueMutex_);
	generation_ = generation;
	bank_ = bank;
	time_ = time;
	times_ = std::move(times);
	fluxWeight_ = fluxWeight;
	source_ = std::move(source);
	bool const cachedAdvance = operation == Operation::Advance && source_.referenceStep > units::Time{};
	allowStealing_ = !cachedAdvance;
	active_.clear();
	if (cachedAdvance) active_ = std::move(cachedWork);
	else for (auto id : owned_) if (level < 0 || blocks_[id].location.level == level) active_.push_back(id);
	// Phases have drained before begin. Invalidate on rollback and before a
	// replacement probe, including caches held here for remotely owned blocks.
	if (operation == Operation::Probe || operation == Operation::Backup || operation == Operation::Restore ||
		operation == Operation::RestoreRadiationStep || operation == Operation::Normalize)
		for (auto const& block : blocks_) if (level < 0 || block.location.level == level)
			if (blockCaches_[block.id]) cachePool_.push_back(std::move(blockCaches_[block.id]));
	next_ = 0;
}

std::vector<std::uint64_t> LocalExecutor::claim(std::uint64_t generation) {
	std::lock_guard guard(queueMutex_);
	if (generation != generation_) throw std::logic_error("Stale work request");
	if (next_ == active_.size()) return {};
	return {active_[next_++]};
}

PhaseResult LocalExecutor::run(Operation operation, units::Time dt, std::uint64_t generation, std::vector<storage::Locality> const& executors) {
	if (generation != generation_) throw std::logic_error("Stale stage dispatch");
#ifdef OCTOTIGERII_WITH_HPX
	std::vector<hpx::future<PhaseResult>> pending;
	for (std::size_t i = 0; i < workspaces_.size(); ++i)
		pending.push_back(hpx::async(profiling::annotated([&, i] { return worker(operation, dt, generation, executors, workspaces_[i]); },
			operation == Operation::Timestep	? "runtime.timestep.worker" :
				operation == Operation::Advance ? "runtime.advance.worker" :
				operation == Operation::PredictRadiationStep || operation == Operation::FinishRadiationStep ||
					operation == Operation::ForecastRadiationStep || operation == Operation::CoupleRadiation ? "runtime.radiation_source.worker" :
				operation == Operation::Reflux	? "runtime.reflux.worker" :
												  "runtime.gravity_kick.worker")));
	auto results = collect(pending);
#else
	std::vector<PhaseResult> results{worker(operation, dt, generation, executors, workspaces_.front())};
#endif
	PhaseResult result;
	for (auto const& part : results) {
		result.cachedBlocks.insert(result.cachedBlocks.end(), part.cachedBlocks.begin(), part.cachedBlocks.end());
		result.boundary += part.boundary;
		result.radiationSourceEnergy += part.radiationSourceEnergy;
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

std::vector<Snapshot> LocalExecutor::snapshots(unsigned bank, mesh::TimeState time) const {
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

PhaseResult LocalExecutor::worker(Operation operation, units::Time dt, std::uint64_t generation, std::vector<storage::Locality> const& executors, Workspace& workspace) {
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
		} else if (operation == Operation::SaveRadiationStep || operation == Operation::PredictRadiationStep ||
			operation == Operation::FinishRadiationStep || operation == Operation::ForecastRadiationStep ||
			operation == Operation::RestoreRadiationStep || operation == Operation::CoupleRadiation) {
			result.radiationSourceEnergy += radiationSource(block, operation, dt);
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
			auto const ownedPlan = plans_.find(id);
			if (ownedPlan == plans_.end()) temporary = makeHaloPlan(config_, blocks_, id);
			auto const& plan = temporary ? *temporary : ownedPlan->second;
			bool const probe = operation == Operation::Probe;
			bool const coupled = std::get<0>(source_.radiation.gas.fields).id != 0;
			auto& cache = blockCaches_[id];
			bool const reuse = !probe && cache && cache->matches(time_, source_.referenceStep,
				std::get<0>(source_.increment.fields).id, bank_, times_, source_.radiation.cacheInitialGasInvalid);
			if constexpr (build::hydro) if (config_.hydroEnabled()) {
				workspace.hydro.advance(block, fields_.hydro, plan, hydro::HydroSystem(config_), bank_, dt, time_, hydroBoundary_,
					fields_.hydroFlux, config_.amr.enabled, times_, source_.increment, source_.referenceStep, physics::RotatingFrame(config_.frame.omega),
					source_.radiation.gas, nullptr, {}, 0, {}, source_.radiation.limiterInterval, {},
					reuse && !source_.radiation.cacheInitialGasInvalid ? &cache->gas : nullptr, probe && !coupled);
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
			if constexpr (build::radiation) if (config_.radiationEnabled() && (operation != Operation::Probe || coupled)) {
				workspace.radiation.advance(block, fields_.radiation, plan, radiation::RadiationSystem(config_.radiation.lightSpeedRatio * constants::c,
					{config_.radiation.closedBoundary, config_.mesh.lower, config_.mesh.upper}),
					bank_, dt, time_, radiationBoundary_, fields_.radiationFlux, config_.amr.enabled, times_, {}, {}, physics::RotatingFrame(config_.frame.omega),
					source_.radiation.radiation, coupled ? &config_ : nullptr,
					coupled && dt > units::Time{} ? source_.radiation.gas : fields_.hydro,
					coupled && dt > units::Time{} ? 1 : bank_, hydroBoundary_, source_.radiation.limiterInterval, radiationMaterial_,
					reuse && coupled ? &cache->radiation : nullptr, false,
					coupled ? (probe ? &workspace.hydro.ghosts : &workspace.hydro.midpointGhosts) : nullptr);
				result.boundary += workspace.radiation.boundaryTransport(block, config_.mesh.boundary, dt);
				if (operation == Operation::Probe && coupled) {
					auto rate = source_.radiation.radiationRate.output(block.interior, 0);
					block.layout.forEachInterior([&](auto const& cell, std::size_t i) {
						radiation::RadiationSystem::Flux rhs{};
						for (int d = 0; d < ndim; ++d) {
							auto upper = cell; ++upper[d];
							rhs += workspace.radiation.work.fluxes[d][block.layout.faceIndex(d, cell)] -
								workspace.radiation.work.fluxes[d][block.layout.faceIndex(d, upper)];
						}
						rate.put(i, (source_.referenceStep / block.cellWidth) * rhs);
					});
					source_.radiation.radiationRate.commit(block.interior, 0, rate);
				}
			}
			if (build::gravity && config_.gravityEnabled()) copyFields(fields_.gravity, block.interior, bank_);
			if (fluxWeight_ > 0 && operation != Operation::Probe) accumulateFlux(block);
			if (probe) {
				cache = acquireCache();
				cache->time = time_; cache->referenceStep = source_.referenceStep;
				cache->rateId = std::get<0>(source_.increment.fields).id;
				cache->bank = bank_; cache->donors = times_;
				cache->gas.save(block, plan, coupled ? workspace.hydro.ghosts : workspace.hydro.rawInitial, coupled);
				if (coupled) cache->radiation.save(block, plan, workspace.radiation.ghosts, true);
				result.cachedBlocks.push_back(id);
			} else releaseCache(cache);
		} else if (operation == Operation::Reflux || operation == Operation::RefluxTracked) {
			auto const plan = !config_.amr.enabled ? std::vector<FluxCorrection>{} :
				(stolen ? makeRefluxPlan(config_, blocks_, id) : refluxPlans_.at(id));
			if (config_.massFractions.enabled) refluxSpecies(block, plan, dt);
			if (build::hydro && config_.hydroEnabled()) reflux(block, plan, fields_.hydro, fields_.hydroFlux, hydro::HydroSystem(config_), dt, operation != Operation::RefluxTracked);
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
		// A cached corrector can be assigned here after its probe was stolen.
		// Statistics describe execution locality, including that affinity work.
		if (block.interior.partition != owner_)
			++result.tasks.stolenTasks;
		else
			++result.tasks.localTasks;
	};
	auto processWithContext = [&](std::uint64_t id, bool stolen) {
		try {
			process(id, stolen);
		} catch (std::runtime_error const& error) {
			std::ostringstream context;
			context.precision(17);
			context << "Runtime block " << id << " level " << blocks_.at(id).location.level
				<< " operation " << static_cast<int>(operation) << " time " << time_.value()
				<< " dt " << dt.value() << " bank " << bank_ << ": " << error.what();
			throw std::runtime_error(context.str());
		}
	};
	for (;;) {
		auto work = claim(generation);
		if (work.empty()) break;
		processWithContext(work.front(), false);
	}
#ifdef OCTOTIGERII_WITH_HPX
	if (config_.runtime.workStealing && allowStealing_)
		for (std::size_t distance = 1; distance < executors.size(); ++distance) {
			auto const donor = (owner_ + distance) % executors.size();
			for (;;) {
				auto work = hpx::async<ClaimAction>(executors[donor], generation).get();
				if (work.empty()) break;
				processWithContext(work.front(), true);
			}
		}
#else
	(void) executors;
#endif
	return result;
}

void LocalExecutor::fluxRegisters(Subgrid const& block, bool reset) {
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

void LocalExecutor::resetFlux(Subgrid const& block) { fluxRegisters(block, true); }

void LocalExecutor::accumulateFlux(Subgrid const& block) { fluxRegisters(block, false); }

void LocalExecutor::advanceSpecies(Subgrid const& block, HaloPlan const& plan, units::Time dt, std::vector<hydro::ConservedState> const& gasGhosts) {
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

void LocalExecutor::refluxSpecies(Subgrid const& block, std::vector<FluxCorrection> const& plan, units::Time dt) {
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

units::Time LocalExecutor::timestep(Subgrid const& block, std::array<units::Velocity, ndim>& maximumSpeed) const {
	auto result = units::Time::from_value(std::numeric_limits<Real>::infinity());
	std::array<units::Velocity, ndim> blockSpeed{};
	std::array<units::Acceleration, ndim> maximumAcceleration{};
	physics::RotatingFrame const frame(config_.frame.omega);
	auto const external = frame.toGrid(config_.hydro.acceleration, time_);
	[[maybe_unused]] auto fieldStep = [&](auto const& fields, auto const& system) {
		auto input = fields.read(block.interior, bank_).get();
		profiling::Region profile("transport.signal_speed");
		units::InverseTime rate{};
		std::array<units::Velocity, ndim> speed{};
		block.layout.forEachInterior([&](auto const& cell, std::size_t i) {
			if (!frame.active()) {
				for (int d = 0; d < ndim; ++d) speed[d] = std::max(speed[d], system.maximumSignalSpeed(input.at(i), d));
				return;
			}
			auto const state = frame.toGridState(input.at(i), time_);
			mesh::PhysicalCoordinates position;
			for (int d = 0; d < ndim; ++d) position[d] = block.lower[d] + (Real(cell[d]) + 0.5) * block.cellWidth;
			for (int d = 0; d < ndim; ++d) {
				auto value = system.maximumSignalSpeed(state, d, frame.normalSpeed(position, d));
				if (d < 2) value += 0.5 * units::abs(config_.frame.omega) * block.cellWidth;
				speed[d] = std::max(speed[d], value);
			}
		});
		for (auto value : speed)
			rate += value / block.cellWidth;
		for (int d = 0; d < ndim; ++d)
			blockSpeed[d] = std::max(blockSpeed[d], speed[d]);
		if (!(rate > units::InverseTime{}) || !units::finite(rate)) throw std::runtime_error("No finite positive signal speed");
		return config_.timestep.cfl / rate;
	};
	if (build::hydro && config_.hydroEnabled()) {
		for (int d = 0; d < ndim; ++d)
			maximumAcceleration[d] = units::abs(external[d]);
		result = fieldStep(fields_.hydro, hydro::HydroSystem(config_));
		units::Acceleration acceleration{};
		for (auto component : external)
			acceleration += units::abs(component);
		if (build::gravity && config_.gravityEnabled()) {
			auto input = fields_.gravity.read(block.interior, bank_).get();
			for (std::size_t i = 0; i < block.interior.count; ++i) {
				units::Acceleration norm{};
				for (int d = 0; d < ndim; ++d) {
					auto const value = units::abs(input.at(i).acceleration(d) + external[d]);
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
	result = std::min(result, frame.maximumTimestep());
	for (int d = 0; d < ndim; ++d) {
		if (units::finite(result)) blockSpeed[d] += maximumAcceleration[d] * result;
		maximumSpeed[d] = std::max(maximumSpeed[d], blockSpeed[d]);
	}
	return result;
}

BoundaryTransport LocalExecutor::finishGravityEnergy(Subgrid const& block, std::vector<GravityWorkFace> const& plan, units::Time dt) {
#if OCTOTIGERII_HYDRO && OCTOTIGERII_GRAVITY
	auto const old = fields_.oldGravity.read(block.interior, 0).get();
	auto const now = fields_.gravity.read(block.interior, bank_).get();
	auto const mass = fields_.massFlux.read(block.massFlux, 0).get();
	auto const kicks = fields_.gravityKickWork.read(block.interior, 0).get();
	std::optional<storage::Buffer<units::EnergyDensity>> rotation;
	if (source_.rotationWork.id) rotation = source_.rotationWork.read(block.interior, 0).get();
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
	hydro::HydroSystem const gas(config_);
	for (std::size_t i = 0; i < block.interior.count; ++i) {
		auto state = input.at(i);
		if (config_.gravity.energyTreatment == "mullen") state.totalEnergy() += work[i] - kicks.data()[i] + (rotation ? rotation->data()[i] : units::EnergyDensity{});
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
#else
	(void) block; (void) plan; (void) dt;
	throw std::logic_error("Gravity energy work requires hydro and gravity in this build");
#endif
}

void LocalExecutor::kick(Subgrid const& block, units::Time dt, bool trackWork) {
#if OCTOTIGERII_HYDRO
	auto input = fields_.hydro.read(block.interior, bank_).get();
	std::optional<storage::Columns<radiation::RadiationSystem::State>> radiationInitial, radiationMidpoint, radiationRate;
	if (std::get<0>(source_.radiation.gas.fields).id &&
		(!trackWork || config_.hasExternalAcceleration())) {
		radiationInitial = source_.radiation.radiation.read(block.interior, 0).get();
		radiationMidpoint = source_.radiation.radiation.read(block.interior, 1).get();
		radiationRate = source_.radiation.radiationRate.read(block.interior, 0).get();
	}
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
			auto const sourceImpulse = radiationInitial ? radiationMidpointImpulse(radiationInitial->at(i), radiationMidpoint->at(i),
				radiationRate->at(i), dt / source_.radiation.referenceStep, config_.radiation.lightSpeedRatio) :
				std::array<units::MomentumDensity, ndim>{};
			std::array<units::Acceleration, ndim> self{};
			if (gravity) {
				for (int d = 0; d < ndim; ++d) self[d] = gravity->at(i).acceleration(d);
				self = physics::RotatingFrame(config_.frame.omega).toInertial(self, time_);
			}
			units::EnergyDensity work{}, sourceWork{};
			for (int d = 0; d < ndim; ++d) {
				auto const old = state.momentum(d);
				auto acceleration = config_.hydro.acceleration[d];
				if (build::gravity && config_.gravityEnabled()) acceleration += self[d];
				auto const impulse = dt * state.density() * acceleration;
				state.momentum(d) += impulse;
				work += impulse * (old + 0.5 * impulse) / state.density();
				// Mullen replaces self-gravity kick work with midpoint mass-flux
				// work, which already contains this source momentum. Only the
				// external-force contribution remains additional in that mode.
				sourceWork += dt * (trackWork ? config_.hydro.acceleration[d] : acceleration) * sourceImpulse[d];
				if (trackWork) kickWork.data()[i] += dt * self[d] * (old + 0.5 * impulse);
			}
			state.totalEnergy() += work;
			addRadiationForceWork(state, sourceWork, hydro::HydroSystem(config_), config_.hydro.dualEnergy.enabled);
			if (!trackWork) hydro::HydroSystem(config_).synchronize(state);
			if (!hydro::HydroSystem(config_).admissible(state)) throw std::runtime_error("Invalid gravity kick state");
			output.put(i, state);
		}
		if (trackWork) fields_.gravityKickWork.commit(block.interior, 0, kickWork);
	}
	fields_.hydro.commit(block.interior, bank_ ^ 1, output);
#else
	(void) block; (void) dt; (void) trackWork;
	throw std::logic_error("Gravity kicks require hydro in this build");
#endif
}

} // namespace octotigerII

#ifdef OCTOTIGERII_WITH_HPX
using LocalExecutorComponent = hpx::components::component<octotigerII::LocalExecutor>;
HPX_REGISTER_COMPONENT(LocalExecutorComponent, octotigerII_executor)
HPX_REGISTER_ACTION(octotigerII::LocalExecutor::InitializeAction, octotigerII_initialize)
HPX_REGISTER_ACTION(octotigerII::LocalExecutor::BeginAction, octotigerII_begin)
HPX_REGISTER_ACTION(octotigerII::LocalExecutor::ClaimAction, octotigerII_claim)
HPX_REGISTER_ACTION(octotigerII::LocalExecutor::RunAction, octotigerII_run)
HPX_REGISTER_ACTION(octotigerII::LocalExecutor::SnapshotsAction, octotigerII_snapshots)
#endif
