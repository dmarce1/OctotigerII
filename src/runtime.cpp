/** @file
 * @brief Runtime lifecycle, snapshots, and public bookkeeping.
 */
#include "runtime/internal.hpp"

namespace octotigerII {

std::vector<Snapshot> Runtime::Impl::exportSnapshots(std::optional<unsigned> requestedBank, std::optional<mesh::TimeState> requestedTime) const {
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

BoundaryTransport Runtime::boundaryTransport() const {
	std::lock_guard guard(impl_->apiMutex);
	return impl_->boundary;
}

units::Energy Runtime::radiationSourceEnergy() const {
	std::lock_guard guard(impl_->apiMutex);
	return impl_->radiationSourceEnergy;
}

std::size_t Runtime::shadowCellCount() const {
	std::lock_guard guard(impl_->apiMutex);
	return impl_->shadow ? impl_->shadow->size() : 0;
}

} // namespace octotigerII
