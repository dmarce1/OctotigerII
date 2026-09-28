/** @file
 * @brief Stage dispatch and completion before field publication.
 */
#include "internal.hpp"

namespace octotigerII {

PhaseResult Runtime::Impl::phase(Operation operation, units::Time dt, int level, std::optional<units::Time> at, Real fluxWeight, SourcePredictor source) {
	if (coupledStep && !std::get<0>(source.radiation.gas.fields).id) {
		source.radiation = *coupledStep;
		source.radiation.gravityInRate = config.gravityEnabled() && std::get<0>(source.increment.fields).id;
		if (!std::get<0>(source.increment.fields).id) {
			source.increment = coupledStep->gasRate;
			source.referenceStep = coupledStep->referenceStep;
		}
	}
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
	if (cacheOwners.size() != topology->blocks().size()) {
		cacheOwners.resize(topology->blocks().size());
		for (auto const& block : topology->blocks()) cacheOwners[block.id] = block.interior.partition;
	}
	std::vector<std::vector<std::uint64_t>> cachedWork(localities.size());
	if (operation == Operation::Advance && source.referenceStep > units::Time{})
		for (auto const& block : topology->blocks()) if (level < 0 || block.location.level == level)
			cachedWork.at(cacheOwners[block.id]).push_back(block.id);
#ifdef OCTOTIGERII_WITH_HPX
	std::vector<hpx::future<void>> starts;
	for (std::size_t i = 0; i < executors.size(); ++i)
		starts.push_back(hpx::async<LocalExecutor::BeginAction>(executors[i], dispatch, stageBank, stageTime, level, haloTimes, fluxWeight, source, operation, std::move(cachedWork[i])));
	finish(starts);
	std::vector<hpx::future<PhaseResult>> pending;
	for (auto const& id : executors)
		pending.push_back(hpx::async<LocalExecutor::RunAction>(id, operation, dt, dispatch, executors));
	auto results = collect(pending);
#else
	executor->begin(dispatch, stageBank, stageTime, level, haloTimes, fluxWeight, source, operation, std::move(cachedWork[0]));
	std::vector<PhaseResult> results{executor->run(operation, dt, dispatch, {})};
#endif
	PhaseResult result;
	for (std::size_t i = 0; i < results.size(); ++i) {
		auto const& part = results[i];
		for (auto id : part.cachedBlocks) cacheOwners.at(id) = i;
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
	auto const expected = std::count_if(topology->blocks().begin(), topology->blocks().end(), [&](auto const& b) { return level < 0 || b.location.level == level; });
	if (result.tasks.localTasks + result.tasks.stolenTasks != std::size_t(expected))
		throw std::logic_error("Stage did not complete every output range");
	if (operation == Operation::FinishRadiationStep) radiationSourceEnergy += result.radiationSourceEnergy;
	statistics.localTasks += result.tasks.localTasks;
	statistics.stolenTasks += result.tasks.stolenTasks;
	return result;
}

} // namespace octotigerII
