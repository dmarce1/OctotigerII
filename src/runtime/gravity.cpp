/** @file
 * @brief Gravity rung integration and endpoint energy accounting.
 */
#include "internal.hpp"

namespace octotigerII {

#if OCTOTIGERII_GRAVITY
void Runtime::Impl::addGravityWork(gravity::Statistics& total, gravity::Statistics const& add) {
	total.multipolePairs += add.multipolePairs; total.directPairs += add.directPairs;
	total.workerTasks += add.workerTasks; total.ewaldPairs += add.ewaldPairs; total.reflectedPairs += add.reflectedPairs;
	if (total.localityCells.size() < add.localityCells.size()) total.localityCells.resize(add.localityCells.size());
	for (std::size_t i = 0; i < add.localityCells.size(); ++i) total.localityCells[i] += add.localityCells[i];
}

void Runtime::Impl::partialGravity(GravityInterval& interval, storage::ColumnHandle<gravity::State> const& output, unsigned outputBank,
	int minimumLevel, units::Time at, bool fullTargets,
	storage::FieldHandle<units::VelocitySquared> const& rotationOutput) {
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
	if (rotationOutput.id) addGravityWork(interval.work, gravity::rotationWorkPotential(*gravitySolver, request, topology->blocks(),
		std::get<0>(fields->directory().hydro.fields), output, outputBank, *interval.rotation,
		rotationOutput, outputBank, config.mesh.upper - config.mesh.lower));
}

storage::ColumnHandle<gravity::State> Runtime::Impl::nestedGravity(GravityInterval& interval, int level, units::Time at) {
	auto& cache = interval.cachedNested[level];
	if (!cache) cache = std::make_unique<storage::ColumnFields<gravity::State>>(interval.cells, interval.store, "timeGravity.cachedNested", false, 1);
	if (!interval.cachedTime.contains(level) || interval.cachedTime.at(level) != at) {
		storage::FieldHandle<units::VelocitySquared> rotation;
		if (interval.rotation) {
			auto& value = interval.cachedRotation[level];
			if (!value) value = std::make_unique<storage::Field<units::VelocitySquared>>(interval.cells, interval.store, 1, "rotationWork.cachedNested");
			rotation = value->handle();
		}
		partialGravity(interval, cache->handle(), 0, level, at, false, rotation);
		interval.cachedTime[level] = at;
	}
	return cache->handle();
}

void Runtime::Impl::shellGravity(GravityInterval& interval, GravityFrame& frame, unsigned endpoint, int nextLevel, units::Time at) {
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
		if (frame.rotation) {
			auto allRotation = interval.cachedRotation.at(frame.level)->handle().read(b.interior, 0).get();
			std::optional<storage::Buffer<units::VelocitySquared>> fastRotation;
			if (nextLevel >= 0) fastRotation = interval.cachedRotation.at(nextLevel)->handle().read(b.interior, 0).get();
			auto n = frame.rotation->nested.handle().output(b.interior, endpoint);
			auto s = frame.rotation->shell.handle().output(b.interior, endpoint);
			for (std::size_t i = 0; i < b.interior.count; ++i) {
				n.data()[i] = allRotation.data()[i];
				s.data()[i] = allRotation.data()[i] - (fastRotation ? fastRotation->data()[i] : units::VelocitySquared{});
			}
			frame.rotation->nested.handle().commit(b.interior, endpoint, n);
			frame.rotation->shell.handle().commit(b.interior, endpoint, s);
		}
	}
}

void Runtime::Impl::assemblePredictor(GravityInterval& interval, GravityFrame const& frame) {
	interval.sourceTime = frame.begin;
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
		if (frame.rotation) {
			auto out = interval.predictorRotation->handle().output(b.interior, 0);
			auto base = (conventional ? frame.rotation->force : frame.rotation->nested).handle().read(b.interior, 0).get();
			std::vector<storage::Buffer<units::VelocitySquared>> parents;
			if (!conventional) for (auto const* f : interval.stack) if (f != &frame)
				parents.push_back(f->rotation->shell.handle().read(b.interior, 0).get());
			for (std::size_t i = 0; i < b.interior.count; ++i) {
				out.data()[i] = base.data()[i];
				for (auto const& p : parents) out.data()[i] += p.data()[i];
			}
			interval.predictorRotation->handle().commit(b.interior, 0, out);
		}
	}
}

void Runtime::Impl::gravitySourceRate(GravityInterval& interval, int level) {
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
		std::optional<storage::Buffer<units::VelocitySquared>> rotation;
		if (interval.rotation) rotation = interval.predictorRotation->handle().read(b.interior, 0).get();
		for (std::size_t i = 0; i < b.interior.count; ++i) {
			auto value = prior.at(i) + correction[i];
			std::array<units::Acceleration, ndim> acceleration;
			for (int d = 0; d < ndim; ++d) acceleration[d] = force.at(i).acceleration(d);
			acceleration = physics::RotatingFrame(config.frame.omega).toInertial(acceleration, interval.sourceTime);
			for (int d = 0; d < ndim; ++d) {
				value.momentum(d) += interval.duration * gas.at(i).density() * acceleration[d];
				if (config.gravity.energyTreatment == "naive") value.totalEnergy() += interval.duration * gas.at(i).momentum(d) * acceleration[d];
			}
			if (config.gravity.energyTreatment == "mullen") {
				value.totalEnergy() += work.work[i];
				if (rotation) value.totalEnergy() += (interval.duration * config.frame.omega) * gas.at(i).density() * rotation->data()[i];
			}
			out.put(i, value);
		}
		interval.rate.handle().commit(b.interior, 0, out);
	}
}

void Runtime::Impl::provisionalGravity(GravityInterval& interval, GravityFrame& frame, units::Time step) {
	auto const& directory = fields->directory();
	bool const conventional = config.gravity.timeIntegration == "conventional";
	for (auto const& b : topology->blocks()) if (b.location.level == frame.level) {
		auto const bank = levels.at(frame.level).bank;
		auto const old = directory.hydro.read(b.interior, bank).get();
		auto const next = directory.hydro.read(b.interior, bank ^ 1).get();
		auto const force = interval.predictor.handle().read(b.interior, 0).get();
		auto const mass = directory.massFlux.read(b.massFlux, 0).get();
		std::optional<storage::Columns<radiation::RadiationSystem::State>> radiationInitial, radiationMidpoint, radiationRate;
		if (coupledStep && config.gravity.energyTreatment == "naive") {
			radiationInitial = coupledStep->radiation.read(b.interior, 0).get();
			radiationMidpoint = coupledStep->radiation.read(b.interior, 1).get();
			radiationRate = coupledStep->radiationRate.read(b.interior, 0).get();
		}
		std::vector<units::EnergyDensity> heat(b.interior.count);
		for (auto* f : interval.stack) {
			// On the level being forecast, the nested field also forecasts the
			// as-yet-unopened descendants' energy transfer across its boundary.
			auto const field = f == &frame ? f->nested.handle() : f->shell.handle();
			auto work = gravity::fluxWork(b, interval.plans[b.id], field, 0, directory.massFlux, 0, step, false);
			if (f->rotation) {
				auto rotation = (f == &frame ? f->rotation->nested : f->rotation->shell).handle().read(b.interior, 0).get();
				for (std::size_t i = 0; i < b.interior.count; ++i)
					work.work[i] += (step * config.frame.omega / 2.0) * (old.at(i).density() + next.at(i).density()) * rotation.data()[i];
			}
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
					std::array<units::Acceleration, ndim> acceleration;
					for (int d = 0; d < ndim; ++d) acceleration[d] = g.at(i).acceleration(d);
					acceleration = physics::RotatingFrame(config.frame.omega).toInertial(acceleration, frame.begin + step / 2.0);
					for (int d = 0; d < ndim; ++d) p.momentum(d) += (step / 2.0) * (old.at(i).density() + next.at(i).density()) * acceleration[d];
					impulses.put(i, p);
				}
				f->impulse.handle().commit(b.interior, 0, impulses);
			}
		}
		auto output = directory.hydro.output(b.interior, bank ^ 1);
		for (std::size_t i = 0; i < b.interior.count; ++i) {
			auto value = next.at(i);
			auto const sourceImpulse = radiationInitial ? radiationMidpointImpulse(radiationInitial->at(i), radiationMidpoint->at(i),
				radiationRate->at(i), step / (2.0 * coupledStep->referenceStep), config.radiation.lightSpeedRatio) :
				std::array<units::MomentumDensity, ndim>{};
			units::EnergyDensity sourceWork{};
			std::array<units::Acceleration, ndim> acceleration;
			for (int d = 0; d < ndim; ++d) acceleration[d] = force.at(i).acceleration(d);
			acceleration = physics::RotatingFrame(config.frame.omega).toInertial(acceleration, frame.begin + step / 2.0);
			for (int d = 0; d < ndim; ++d) {
				value.momentum(d) += (step / 2.0) * (old.at(i).density() + value.density()) * acceleration[d];
				if (config.gravity.energyTreatment == "naive") value.totalEnergy() += (step / 2.0) * acceleration[d] * (old.at(i).momentum(d) + value.momentum(d));
				sourceWork += step * acceleration[d] * sourceImpulse[d];
			}
			if (config.gravity.energyTreatment == "mullen") value.totalEnergy() += heat[i];
			hydro::HydroSystem const system(config.hydro);
			addRadiationForceWork(value, sourceWork, system, config.hydro.dualEnergy.enabled);
			system.synchronize(value);
			if (!system.admissible(value)) throw std::runtime_error("Provisional gravity source produced an inadmissible state");
			output.put(i, value);
		}
		directory.hydro.commit(b.interior, bank ^ 1, output);
	}
}

void Runtime::Impl::closeGravityFrame(GravityInterval& interval, GravityFrame& frame) {
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
		std::optional<storage::Buffer<units::VelocitySquared>> oldRotation, newRotation;
		if (frame.rotation) {
			oldRotation = frame.rotation->shell.handle().read(b.interior, 0).get();
			newRotation = frame.rotation->shell.handle().read(b.interior, 1).get();
		}
		for (std::size_t i = 0; i < b.interior.count; ++i) {
			auto value = input.at(i);
			std::array<units::Acceleration, ndim> ga, gb;
			for (int d = 0; d < ndim; ++d) { ga[d] = oldForce.at(i).acceleration(d); gb[d] = newForce.at(i).acceleration(d); }
			ga = physics::RotatingFrame(config.frame.omega).toInertial(ga, frame.begin);
			gb = physics::RotatingFrame(config.frame.omega).toInertial(gb, frame.begin + frame.duration);
			if (frame.rotation) work.work[i] += (frame.duration * config.frame.omega / 2.0)
				* (rho.data()[i] * oldRotation->data()[i] + value.density() * newRotation->data()[i]);
			if (!conventional || b.location.level == frame.level) for (int d = 0; d < ndim; ++d) {
				auto const dp = (frame.duration / 2.0) * (rho.data()[i] * ga[d]
					+ value.density() * gb[d]) - impulse.at(i).momentum(d);
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

void Runtime::Impl::advanceGravityLevel(GravityInterval& interval, std::vector<int> const& occupied, std::size_t index,
	units::Time begin, units::Time duration, bool root) {
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
		GravityFrame frame(level, now, step, interval.cells, interval.faces, interval.store, bool(interval.rotation));
		int const nextLevel = index + 1 < occupied.size() ? occupied[index + 1] : -1;
		shellGravity(interval, frame, 0, nextLevel, now);
		if (config.gravity.timeIntegration == "conventional") partialGravity(interval, frame.force.handle(), 0, 0, now, true,
			frame.rotation ? frame.rotation->force.handle() : storage::FieldHandle<units::VelocitySquared>{});
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
		if (coupledStep) {
			coupledStep->limiterInterval = step / 2.0;
			phase(Operation::SaveRadiationStep, {}, level, now);
		}
		for (std::size_t j = index; j < occupied.size(); ++j)
			phase(Operation::Probe, {}, occupied[j], now, 0, {interval.rate.handle(), interval.duration});
		gravitySourceRate(interval, level);
		if (coupledStep) phase(Operation::PredictRadiationStep, step, -1, now, 0, {interval.rate.handle(), interval.duration});
		auto& state = levels.at(level);
		state.begin = now; state.end = now + step;
		interval.boundary += phase(Operation::Advance, step, level, now, step / duration,
			{interval.rate.handle(), interval.duration}).boundary;
		provisionalGravity(interval, frame, step);
		// The halo forecast uses the accepted initial numerical RHS. The raw
		// coarse transport endpoint has not been refluxed and is only O(H)
		// accurate at a refinement boundary, even with a midpoint flux.
		for (auto const& b : topology->blocks()) if (!coupledStep && b.location.level == level) {
			auto old = fields->directory().hydro.read(b.interior, state.bank).get();
			auto rate = interval.rate.handle().read(b.interior, 0).get();
			auto out = fields->directory().hydro.output(b.interior, 3);
			std::vector<units::Density> density(b.interior.count);
			for (std::size_t i = 0; i < b.interior.count; ++i) {
				auto value = old.at(i) + (step / interval.duration) * rate.at(i);
				if (!coupledStep) predictorKineticRemainder(old.at(i), value);
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
		if (coupledStep) phase(Operation::ForecastRadiationStep, step, level, now);
		state.predictorBank = 3;
		state.pending = true;
		if (nextLevel >= 0) advanceGravityLevel(interval, occupied, index + 1, now, step);
		phase(Operation::RefluxTracked, step, level, now + step);
		state.pending = false; state.bank ^= 1;
		shellGravity(interval, frame, 1, nextLevel, now + step);
		if (config.gravity.timeIntegration == "conventional") partialGravity(interval, frame.force.handle(), 1, 0, now + step, true,
			frame.rotation ? frame.rotation->force.handle() : storage::FieldHandle<units::VelocitySquared>{});
		closeGravityFrame(interval, frame);
		if (coupledStep) {
			phase(Operation::FinishRadiationStep, step, level, now + step);
			state.bank ^= 1;
		}
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

gravity::Statistics Runtime::advanceGravity(units::Time dt) {
	std::lock_guard guard(impl_->apiMutex);
	return advanceGravityUnlocked(dt);
}

gravity::Statistics Runtime::advanceGravityUnlocked(units::Time dt) {
	if (!impl_->config.hydroEnabled() || !impl_->config.gravityEnabled()) throw std::logic_error("Coupled gravity requires self-gravitating gas");
	if (!(dt > units::Time{}) || !units::finite(dt) || impl_->time.time + dt == impl_->time.time) throw std::invalid_argument("Invalid step size");
	if (dt > physics::RotatingFrame(impl_->config.frame.omega).maximumTimestep() * (1 + 64 * epsilonR))
		throw std::invalid_argument("Rotating-grid step exceeds the angular-phase limit; use stableTimestep()");
	if (!impl_->config.amr.enabled || !impl_->config.timestep.refinement || impl_->config.hasExternalAcceleration()) {
		beginGravityEnergyUnlocked(); kickGravityUnlocked(dt / 2.0); advanceUnlocked(dt);
		auto work = solveGravityUnlocked(); kickGravityUnlocked(dt / 2.0); finishGravityEnergyUnlocked(dt);
#if OCTOTIGERII_GRAVITY
		impl_->addGravityWork(work, impl_->rotationStatistics);
#endif
		return work;
	}
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
		if (interval.rotation) {
			if (!impl_->gravitySolver) impl_->gravitySolver = std::make_unique<gravity::FieldSolver>(impl_->config,
				impl_->topology->blocks(), impl_->fields->directory(), impl_->localities);
			gravity::FieldSolveRequest request;
			request.sourceBank = impl_->bank;
			impl_->addGravityWork(interval.work, gravity::rotationWorkPotential(*impl_->gravitySolver, request, impl_->topology->blocks(),
				std::get<0>(impl_->fields->directory().hydro.fields), interval.predictor.handle(), 0, *interval.rotation,
				interval.predictorRotation->handle(), 0, impl_->config.mesh.upper - impl_->config.mesh.lower));
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

void Runtime::kickGravity(units::Time dt) {
	std::lock_guard guard(impl_->apiMutex);
	kickGravityUnlocked(dt);
}

void Runtime::kickGravityUnlocked(units::Time dt) {
	profiling::Elapsed profile("runtime.gravity_kick.wall_ns");
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
	beginGravityEnergyUnlocked();
}

void Runtime::beginGravityEnergyUnlocked() {
	if (impl_->gravityEnergyActive || !impl_->config.hydroEnabled() || !impl_->config.gravityEnabled() ||
		!impl_->gravityReady || impl_->gravityTime != impl_->time.time)
		throw std::logic_error("Conservative gravity step needs synchronized gas and gravity");
	impl_->phase(Operation::PrepareGravityEnergy, {});
#if OCTOTIGERII_GRAVITY
	impl_->rotationStatistics = {};
	if (impl_->config.frame.omega != units::InverseTime{} && impl_->config.gravity.energyTreatment == "mullen") {
		if (!impl_->gravitySolver) impl_->gravitySolver = std::make_unique<gravity::FieldSolver>(impl_->config,
			impl_->topology->blocks(), impl_->fields->directory(), impl_->localities);
		auto rotation = std::make_unique<Impl::GlobalRotationWork>(impl_->topology->storageLayout(), impl_->localities);
		auto const& fields = impl_->fields->directory();
		gravity::FieldSolveRequest request; request.sourceBank = impl_->bank;
		impl_->addGravityWork(impl_->rotationStatistics, gravity::rotationWorkPotential(*impl_->gravitySolver, request,
			impl_->topology->blocks(), std::get<0>(fields.hydro.fields), fields.gravity, impl_->bank,
			rotation->scratch, rotation->potential.handle(), 0, impl_->config.mesh.upper - impl_->config.mesh.lower));
		for (auto const& b : impl_->topology->blocks()) {
			auto input = fields.hydro.read(b.interior, impl_->bank).get();
			auto density = rotation->density.handle().output(b.interior, 0);
			for (std::size_t i = 0; i < b.interior.count; ++i) density.data()[i] = input.at(i).density();
			rotation->density.handle().commit(b.interior, 0, density);
		}
		impl_->globalRotation = std::move(rotation);
	}
#endif
	impl_->gravityEnergyStart = impl_->time.time;
	impl_->gravityEnergyActive = true;
}

void Runtime::finishGravityEnergy(units::Time dt) {
	std::lock_guard guard(impl_->apiMutex);
	finishGravityEnergyUnlocked(dt);
}

void Runtime::finishGravityEnergyUnlocked(units::Time dt) {
	if (!impl_->gravityEnergyActive || !(dt > units::Time{}) || !units::finite(dt) ||
		impl_->gravityEnergyStart + dt != impl_->time.time || impl_->gravityTime != impl_->time.time)
		throw std::logic_error("Conservative gravity work requires the completed transport interval and endpoint gravity");
	SourcePredictor source;
#if OCTOTIGERII_GRAVITY
	if (impl_->globalRotation) {
		auto& rotation = *impl_->globalRotation;
		auto const& fields = impl_->fields->directory();
		gravity::FieldSolveRequest request; request.sourceBank = impl_->bank;
		impl_->addGravityWork(impl_->rotationStatistics, gravity::rotationWorkPotential(*impl_->gravitySolver, request,
			impl_->topology->blocks(), std::get<0>(fields.hydro.fields), fields.gravity, impl_->bank,
			rotation.scratch, rotation.potential.handle(), 1, impl_->config.mesh.upper - impl_->config.mesh.lower));
		for (auto const& b : impl_->topology->blocks()) {
			auto density = rotation.density.handle().read(b.interior, 0).get();
			auto gas = fields.hydro.read(b.interior, impl_->bank).get();
			auto old = rotation.potential.handle().read(b.interior, 0).get();
			auto now = rotation.potential.handle().read(b.interior, 1).get();
			auto work = rotation.work.handle().output(b.interior, 0);
			for (std::size_t i = 0; i < b.interior.count; ++i) work.data()[i] = (dt * impl_->config.frame.omega / 2.0)
				* (density.data()[i] * old.data()[i] + gas.at(i).density() * now.data()[i]);
			rotation.work.handle().commit(b.interior, 0, work);
		}
		source.rotationWork = rotation.work.handle();
	}
#endif
	auto result = impl_->phase(Operation::FinishGravityEnergy, dt, -1, {}, 0, source);
	std::unique_ptr<amr::Hierarchy> nextShadow;
	if (impl_->shadow) {
		nextShadow = std::make_unique<amr::Hierarchy>(*impl_->shadow);
		nextShadow->refreshLeaves(impl_->exportSnapshots(impl_->bank ^ 1));
	}
	impl_->bank ^= 1;
	++impl_->generation;
	impl_->boundary += result.boundary;
	impl_->gravityEnergyActive = false;
#if OCTOTIGERII_GRAVITY
	impl_->globalRotation.reset();
#endif
	if (nextShadow) impl_->shadow = std::move(nextShadow);
}

gravity::Statistics Runtime::solveGravity() {
	std::lock_guard guard(impl_->apiMutex);
	return solveGravityUnlocked();
}

gravity::Statistics Runtime::solveGravityUnlocked() {
	if (!impl_->config.gravityEnabled()) throw std::logic_error("Selected problem does not enable gravity");
	profiling::Elapsed profile("runtime.gravity.wall_ns");
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

} // namespace octotigerII
