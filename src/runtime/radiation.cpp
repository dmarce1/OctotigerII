/** @file
 * @brief Paired radiation exchange stages and complete coupled interval rollback.
 */
#include "internal.hpp"

namespace octotigerII {
namespace {
template <typename State>
void copyBetween(storage::ColumnHandle<State> const& from, unsigned source,
	storage::ColumnHandle<State> const& to, unsigned destination, storage::Range range) {
	auto input = from.read(range, source).get();
	auto output = to.output(range, destination);
	for (std::size_t i = 0; i < range.count; ++i) output.put(i, input.at(i));
	to.commit(range, destination, output);
}

struct CouplingStorage {
	storage::PartitionSet store;
	storage::ColumnFields<hydro::ConservedState> gas, gasRate;
	storage::ColumnFields<radiation::RadiationSystem::State> radiation, radiationRate;
	storage::ColumnFields<gravity::State> gravity;
	std::vector<std::unique_ptr<storage::Field<units::Density>>> species;
	CouplingStorage(storage::Layout const& cells, std::vector<storage::Locality> const& localities, std::size_t speciesCount)
	  : store(localities), gas(cells, store, "coupling.gas", false, 2), gasRate(cells, store, "coupling.gasRate", false, 1),
		radiation(cells, store, "coupling.radiation", false, 2), radiationRate(cells, store, "coupling.radiationRate", false, 1),
		gravity(cells, store, "coupling.gravity", false, 1) {
		for (std::size_t i = 0; i < speciesCount; ++i)
			species.push_back(std::make_unique<storage::Field<units::Density>>(cells, store, 1, "coupling.species"));
	}
	RadiationStepFields handles(units::Time dt) const {
		RadiationStepFields result;
		result.gas = gas.handle(); result.gasRate = gasRate.handle();
		result.radiation = radiation.handle(); result.radiationRate = radiationRate.handle();
		result.gravity = gravity.handle(); result.referenceStep = dt; result.limiterInterval = dt / 2.0;
		for (auto const& field : species) result.species.push_back(field->handle());
		return result;
	}
};
}

units::Energy LocalExecutor::radiationSource(Subgrid const& block, Operation operation, units::Time dt) {
#if OCTOTIGERII_HYDRO && OCTOTIGERII_RADIATION
	auto const& scratch = source_.radiation;
	auto copySpecies = [&](bool save, unsigned destination) {
		for (std::size_t s = 0; s < fields_.species.size(); ++s) {
			auto const& from = save ? fields_.species[s] : scratch.species[s];
			auto const& to = save ? scratch.species[s] : fields_.species[s];
			auto input = from.read(block.interior, save ? bank_ : 0).get();
			auto output = to.output(block.interior, destination);
			std::copy_n(input.data(), block.interior.count, output.data());
			to.commit(block.interior, destination, output);
		}
	};
	if (operation == Operation::SaveRadiationStep) {
		copyBetween(fields_.hydro, bank_, scratch.gas, 0, block.interior);
		copyBetween(fields_.radiation, bank_, scratch.radiation, 0, block.interior);
		if (config_.gravityEnabled()) copyBetween(fields_.gravity, bank_, scratch.gravity, 0, block.interior);
		copySpecies(true, 0);
		return {};
	}
	if (operation == Operation::RestoreRadiationStep) {
		copyBetween(scratch.gas, 0, fields_.hydro, bank_, block.interior);
		copyBetween(scratch.radiation, 0, fields_.radiation, bank_, block.interior);
		if (config_.gravityEnabled()) copyBetween(scratch.gravity, 0, fields_.gravity, bank_, block.interior);
		copySpecies(false, bank_);
		return {};
	}
	hydro::HydroSystem const system(config_);
	// These positions are the problem's reference grid coordinates, the same
	// coordinates used for its initialized material and photon source profile.
	auto evaluationTime = time_;
	if (operation == Operation::PredictRadiationStep) evaluationTime += dt / 4.0;
	else if (operation == Operation::ForecastRadiationStep) evaluationTime += dt / 2.0;
	else if (operation == Operation::FinishRadiationStep) evaluationTime -= dt / 2.0;
	std::vector<RadiationMaterial> material(block.interior.count);
	block.layout.forEachInterior([&](auto const& cell, std::size_t i) {
		material[i] = checkedRadiationMaterial(radiationMaterial_, block.layout.cellCenter(block.lower, block.cellWidth, cell), evaluationTime);
	});
	if (operation == Operation::PredictRadiationStep) {
		unsigned const inputBank = times_.empty() ? bank_ : times_.at(block.location.level).bank;
		auto gas = fields_.hydro.read(block.interior, inputBank).get();
		auto rad = fields_.radiation.read(block.interior, inputBank).get();
		auto gasRate = source_.increment.read(block.interior, 0).get();
		auto radRate = scratch.radiationRate.read(block.interior, 0).get();
		std::optional<storage::Columns<hydro::ConservedState>> gasEndpoint;
		std::optional<storage::Columns<radiation::RadiationSystem::State>> radEndpoint;
		Real alpha = 0;
		if (!times_.empty()) {
			auto const& interval = times_.at(block.location.level);
			alpha = interval.fraction;
			if (alpha > 0) {
				unsigned const endpoint = interval.nextBank == ~0u ? inputBank ^ 1 : interval.nextBank;
				gasEndpoint = fields_.hydro.read(block.interior, endpoint).get();
				radEndpoint = fields_.radiation.read(block.interior, endpoint).get();
			}
		}
		std::optional<storage::Columns<gravity::State>> gravity;
		if (config_.gravityEnabled() && !scratch.gravityInRate) gravity = fields_.gravity.read(block.interior, inputBank).get();
		auto gasOutput = scratch.gas.output(block.interior, 1);
		auto radOutput = scratch.radiation.output(block.interior, 1);
		for (std::size_t i = 0; i < block.interior.count; ++i) {
			auto g = gas.at(i); auto r = rad.at(i);
			if (alpha > 0) { g = (1 - alpha) * g + alpha * gasEndpoint->at(i); r = (1 - alpha) * r + alpha * radEndpoint->at(i); }
			auto dg = (dt / (2.0 * source_.referenceStep)) * gasRate.at(i);
			auto dr = (dt / (2.0 * source_.referenceStep)) * radRate.at(i);
			dr.energy() += config_.radiation.lightSpeedRatio * material[i].photonPower * (dt / 2.0);
			if (!scratch.gravityInRate) {
				auto acceleration = config_.hydro.acceleration;
				if (gravity) {
					std::array<units::Acceleration, ndim> self{};
					for (int d = 0; d < ndim; ++d) self[d] = gravity->at(i).acceleration(d);
					self = physics::RotatingFrame(config_.frame.omega).toInertial(self, time_);
					for (int d = 0; d < ndim; ++d) acceleration[d] += self[d];
				}
				for (int d = 0; d < ndim; ++d) {
					dg.momentum(d) += (dt / 2.0) * g.density() * acceleration[d];
					dg.totalEnergy() += (dt / 2.0) * g.momentum(d) * acceleration[d];
				}
			}
			radiation::coupleForcedWithOpacityLaw(g, r, dg, dr, system, radiation::opacityLaw(config_, material[i].opacity),
				config_.radiation.lightSpeedRatio, dt / 2.0);
			gasOutput.put(i, g); radOutput.put(i, r);
		}
		scratch.gas.commit(block.interior, 1, gasOutput);
		scratch.radiation.commit(block.interior, 1, radOutput);
		return {};
	}
	unsigned const inputBank = operation == Operation::ForecastRadiationStep ? bank_ ^ 1 : bank_;
	unsigned const outputBank = operation == Operation::ForecastRadiationStep ? 3 : bank_ ^ 1;
	auto currentGas = fields_.hydro.read(block.interior, inputBank).get();
	auto currentRad = fields_.radiation.read(block.interior, inputBank).get();
	std::optional<storage::Columns<hydro::ConservedState>> initialGas;
	std::optional<storage::Columns<radiation::RadiationSystem::State>> initialRad;
	if (operation != Operation::CoupleRadiation) {
		initialGas = scratch.gas.read(block.interior, 0).get();
		initialRad = scratch.radiation.read(block.interior, 0).get();
	}
	auto gasOutput = fields_.hydro.output(block.interior, outputBank);
	auto radOutput = fields_.radiation.output(block.interior, outputBank);
	units::Energy injected{}, compensation{};
	auto const volume = block.layout.cellMeasure(block.cellWidth);
	for (std::size_t i = 0; i < block.interior.count; ++i) {
		auto gas = initialGas ? initialGas->at(i) : currentGas.at(i);
		auto rad = initialRad ? initialRad->at(i) : currentRad.at(i);
		if (initialGas) {
			auto radiationDrive = radiation::RadiationSystem::State(currentRad.at(i) - rad);
			auto const photons = config_.radiation.lightSpeedRatio * material[i].photonPower * dt;
			radiationDrive.energy() += photons;
			radiation::coupleForcedWithOpacityLaw(gas, rad, currentGas.at(i) - gas, radiationDrive,
				system, radiation::opacityLaw(config_, material[i].opacity), config_.radiation.lightSpeedRatio, dt);
			if (operation == Operation::FinishRadiationStep) {
				auto const corrected = volume * photons - compensation;
				auto const sum = injected + corrected;
				compensation = (sum - injected) - corrected;
				injected = sum;
			}
		} else radiation::coupleWithOpacityLaw(gas, rad, system, radiation::opacityLaw(config_, material[i].opacity), config_.radiation.lightSpeedRatio, dt);
		gasOutput.put(i, gas); radOutput.put(i, rad);
	}
	fields_.hydro.commit(block.interior, outputBank, gasOutput);
	fields_.radiation.commit(block.interior, outputBank, radOutput);
	if (config_.gravityEnabled()) copyFields(fields_.gravity, block.interior, inputBank, outputBank);
	for (auto const& species : fields_.species) copyScalar(species, block.interior, inputBank, outputBank);
	return injected;
#else
	(void) block; (void) operation; (void) dt;
	throw std::logic_error("Radiation-matter exchange requires hydro and radiation in this build");
#endif
}

void Runtime::coupleRadiation(units::Time dt) {
	std::lock_guard guard(impl_->apiMutex);
	if (!impl_->config.hydroEnabled() || !impl_->config.radiationEnabled() || dt < units::Time{} || !units::finite(dt))
		throw std::invalid_argument("Local radiation exchange requires hydro, radiation, and a finite nonnegative interval");
	if (impl_->gravityEnergyActive || impl_->regridEnergyPending || !impl_->levels.empty())
		throw std::logic_error("Local radiation exchange requires a synchronized closed interval");
	if (dt == units::Time{} || (!radiation::radiationCouplingEnabled(impl_->config) && !problemHasRadiationMaterial(impl_->config))) return;
	impl_->phase(Operation::CoupleRadiation, dt);
	std::unique_ptr<amr::Hierarchy> nextShadow;
	if (impl_->shadow) {
		nextShadow = std::make_unique<amr::Hierarchy>(*impl_->shadow);
		nextShadow->refreshLeaves(impl_->exportSnapshots(impl_->bank ^ 1));
	}
	impl_->bank ^= 1;
	++impl_->generation;
	if (nextShadow) impl_->shadow = std::move(nextShadow);
	impl_->applyEosFloor();
}

gravity::Statistics Runtime::advanceCoupled(units::Time dt) {
	std::lock_guard guard(impl_->apiMutex);
	auto const& config = impl_->config;
	if (!config.hydroEnabled() || !config.radiationEnabled() || !(dt > units::Time{}) || !units::finite(dt) || impl_->time.time + dt == impl_->time.time)
		throw std::invalid_argument("Coupled advance requires hydro, radiation, and a positive finite interval");
	if (impl_->gravityEnergyActive || impl_->regridEnergyPending || !impl_->levels.empty() || impl_->coupledStep ||
		(config.gravityEnabled() && (!impl_->gravityReady || impl_->gravityTime != impl_->time.time)))
		throw std::logic_error("Coupled advance requires a synchronized closed interval");
	if (dt > physics::RotatingFrame(config.frame.omega).maximumTimestep() * (1 + 64 * epsilonR))
		throw std::invalid_argument("Rotating-grid step exceeds the angular-phase limit");
	if (!radiation::radiationCouplingEnabled(config) && !problemHasRadiationMaterial(config)) {
		if (config.gravityEnabled()) { auto result=advanceGravityUnlocked(dt); impl_->applyEosFloor(); return result; }
		if (config.hasExternalAcceleration()) kickGravityUnlocked(dt / 2.0);
		advanceUnlocked(dt);
		if (config.hasExternalAcceleration()) kickGravityUnlocked(dt / 2.0);
		impl_->applyEosFloor();
		return {};
	}
	profiling::Elapsed profile("runtime.radiation_coupled.wall_ns");
	unsigned const savedBank = impl_->bank;
	auto const savedTime = impl_->time;
	auto const savedGeneration = impl_->generation;
	auto const savedStatistics = impl_->statistics;
	auto const savedBoundary = impl_->boundary;
	auto const savedSourceEnergy = impl_->radiationSourceEnergy;
	auto const savedFloorEnergy=impl_->eosFloorEnergy;
	auto const savedFloorCells=impl_->eosFloorCells;
	auto const savedSpeed = impl_->signalSpeed;
	auto const savedTravel = impl_->travel;
	auto const savedGravityTime = impl_->gravityTime;
	bool const savedGravityReady = impl_->gravityReady;
	auto const savedGravityStart = impl_->gravityEnergyStart;
#if OCTOTIGERII_GRAVITY
	auto const savedRotationStatistics = impl_->rotationStatistics;
#endif
	std::unique_ptr<amr::Hierarchy> savedShadow;
	if (impl_->shadow) savedShadow = std::make_unique<amr::Hierarchy>(*impl_->shadow);
	CouplingStorage backup(impl_->topology->storageLayout(), impl_->localities, impl_->fields->directory().species.size());
	SourcePredictor backupSource; backupSource.radiation = backup.handles(dt);
	impl_->phase(Operation::SaveRadiationStep, {}, -1, {}, 0, backupSource);
	CouplingStorage step(impl_->topology->storageLayout(), impl_->localities, impl_->fields->directory().species.size());
	impl_->coupledStep = step.handles(dt);
	try {
		gravity::Statistics gravityWork;
		bool const refinedGravity = config.gravityEnabled() && config.amr.enabled && config.timestep.refinement && !config.hasExternalAcceleration();
		bool const refinedTransport = !config.gravityEnabled() && impl_->timeRefinement();
		if (refinedGravity) gravityWork = advanceGravityUnlocked(dt);
		else if (refinedTransport) {
			impl_->phase(Operation::Probe, {});
			advanceUnlocked(dt);
		} else {
			impl_->phase(Operation::SaveRadiationStep, {});
			impl_->phase(Operation::Probe, {});
			impl_->phase(Operation::PredictRadiationStep, dt);
			if (config.gravityEnabled()) beginGravityEnergyUnlocked();
			bool const kick = config.gravityEnabled() || config.hasExternalAcceleration();
			if (kick) {
				kickGravityUnlocked(dt / 2.0);
				// The predictor probe predates the accepted opening gas kick.
				// Radiation is unchanged, but the conservative gas base needs fresh halos.
				impl_->coupledStep->cacheInitialGasInvalid = true;
			}
			advanceUnlocked(dt);
			if (config.gravityEnabled()) gravityWork = solveGravityUnlocked();
			if (kick) kickGravityUnlocked(dt / 2.0);
			if (config.gravityEnabled()) {
				finishGravityEnergyUnlocked(dt);
#if OCTOTIGERII_GRAVITY
				impl_->addGravityWork(gravityWork, impl_->rotationStatistics);
#endif
			}
			impl_->phase(Operation::FinishRadiationStep, dt);
			impl_->bank ^= 1;
			++impl_->generation;
		}
		impl_->applyEosFloor();
		if (impl_->shadow) impl_->shadow->refreshLeaves(impl_->exportSnapshots());
		impl_->coupledStep.reset();
		return gravityWork;
	} catch (...) {
		auto failure = std::current_exception();
		impl_->levels.clear();
		impl_->coupledStep.reset();
		impl_->bank = savedBank;
		impl_->time = savedTime;
		impl_->phase(Operation::RestoreRadiationStep, {}, -1, {}, 0, backupSource);
		impl_->generation = savedGeneration;
		impl_->statistics = savedStatistics;
		impl_->boundary = savedBoundary;
		impl_->radiationSourceEnergy = savedSourceEnergy;
		impl_->eosFloorEnergy=savedFloorEnergy; impl_->eosFloorCells=savedFloorCells;
		impl_->signalSpeed = savedSpeed;
		impl_->travel = savedTravel;
		impl_->gravityTime = savedGravityTime;
		impl_->gravityReady = savedGravityReady;
		impl_->gravityEnergyStart = savedGravityStart;
		impl_->gravityEnergyActive = false;
#if OCTOTIGERII_GRAVITY
		impl_->globalRotation.reset();
		impl_->rotationStatistics = savedRotationStatistics;
#endif
		impl_->shadow = std::move(savedShadow);
		std::rethrow_exception(failure);
	}
}

} // namespace octotigerII
