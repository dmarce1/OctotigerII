/** @file
 * @brief Transport timestep selection and AMR subcycling.
 */
#include "internal.hpp"

namespace octotigerII {

bool Runtime::Impl::timeRefinement() const {
	if (!config.amr.enabled || !config.timestep.refinement || config.hasExternalAcceleration()) return false;
	auto const& blocks = topology->blocks();
	return std::any_of(blocks.begin(), blocks.end(), [&](auto const& b) { return b.location.level != blocks.front().location.level; });
}

int Runtime::Impl::coarsestLevel() const {
	int result = config.amr.maxLevel;
	for (auto const& b : topology->blocks()) result = std::min(result, b.location.level);
	return result;
}

BoundaryTransport Runtime::Impl::advanceLevel(std::vector<int> const& occupied, std::size_t index,
	units::Time begin, units::Time duration, bool root) {
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
		if (coupledStep) {
			coupledStep->limiterInterval = step / 2.0;
			phase(Operation::SaveRadiationStep, {}, level, now);
			phase(Operation::Probe, {}, level, now);
			phase(Operation::PredictRadiationStep, step, -1, now);
		}
		transported += phase(Operation::Advance, step, level, now, step / duration).boundary;
		if (coupledStep) {
			phase(Operation::ForecastRadiationStep, step, level, now);
			state.predictorBank = 3;
		}
		state.pending = true;
		if (index + 1 < occupied.size()) transported += advanceLevel(occupied, index + 1, now, step);
		// Fine registers now cover precisely this coarse step. Reflux before
		// publishing the coarse endpoint and resetting the child's registers.
		phase(Operation::Reflux, step, level, now + step);
		state.pending = false;
		state.bank ^= 1;
		if (coupledStep) {
			phase(Operation::FinishRadiationStep, step, level, now + step);
			state.bank ^= 1;
		}
		++statistics.levelSteps.at(level);
		elapsed += fraction;
	}
	return transported;
}

units::Time Runtime::stableTimestep() const {
	profiling::Elapsed profile("runtime.timestep.wall_ns");
	std::lock_guard guard(impl_->apiMutex);
	if (impl_->regridEnergyPending) throw std::logic_error("Solve gravity after regridding before computing a timestep");
	auto const result = impl_->phase(Operation::Timestep, {});
	impl_->signalSpeed = result.signalSpeed;
	if (!impl_->timeRefinement()) return result.timestep;
	auto interval = result.levelTimestep.at(impl_->coarsestLevel());
	if (impl_->config.hydroEnabled() && impl_->config.radiationEnabled() &&
		(radiation::radiationCouplingEnabled(impl_->config) || problemHasRadiationMaterial(impl_->config))) {
		// Every leaf supplies source-aware midpoint data for coarse ghost
		// averages, including fine leaves beyond immediate face neighbors.
		// Their half-interval transport drives must respect their own CFL.
		interval = std::min(interval, 2.0 * result.timestep);
	}
	return interval;
}

void Runtime::advance(units::Time dt) {
	std::lock_guard guard(impl_->apiMutex);
	advanceUnlocked(dt);
	if (!impl_->gravityEnergyActive && !impl_->coupledStep) impl_->applyEosFloor();
}

void Runtime::advanceUnlocked(units::Time dt) {
	profiling::Elapsed profile("runtime.advance.wall_ns");
	if (impl_->regridEnergyPending) throw std::logic_error("Solve gravity after regridding before advancing");
	if (!(dt > units::Time{}) || !units::finite(dt) || impl_->time.time + dt == impl_->time.time) throw std::invalid_argument("Invalid step size");
	if (dt > finiteVolume::RotatingFrame(impl_->config.frame.omega).maximumTimestep() * (1 + 64 * epsilonR))
		throw std::invalid_argument("Rotating-grid step exceeds the angular-phase limit; use stableTimestep()");
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

} // namespace octotigerII
