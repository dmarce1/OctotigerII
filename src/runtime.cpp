#include "octotigerII/runtime.hpp"
#include <algorithm>
#include <limits>
#include <stdexcept>

#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/include/actions.hpp>
#include <hpx/include/async.hpp>
#include <hpx/include/components.hpp>
#include <hpx/include/serialization.hpp>
#include <hpx/runtime_distributed/find_all_localities.hpp>

namespace octotigerII {
class SubgridComponent : public hpx::components::component_base<SubgridComponent> {
  public:
	SubgridComponent() = default;
	SubgridComponent(Config c, mesh::BlockLocation location) : subgrid_(c, location) {}
	Snapshot snapshot() const {
		return subgrid_.snapshot();
	}
	Real stableTimestep() const {
		return subgrid_.stableTimestep();
	}
	void advance(std::vector<HydroSnapshot> const& h, std::vector<RadiationSnapshot> const& r,
				 Real dt) {
		subgrid_.advance(h, r, dt);
	}
	void setGravity(std::vector<gravity::State> const& fields) {
		subgrid_.setGravity(fields);
	}
	void kickGravity(Real dt) {
		subgrid_.kickGravity(dt);
	}
	HPX_DEFINE_COMPONENT_ACTION(SubgridComponent, snapshot, SnapshotAction)
	HPX_DEFINE_COMPONENT_ACTION(SubgridComponent, stableTimestep, TimestepAction)
	HPX_DEFINE_COMPONENT_ACTION(SubgridComponent, advance, AdvanceAction)
	HPX_DEFINE_COMPONENT_ACTION(SubgridComponent, setGravity, GravityAction)
	HPX_DEFINE_COMPONENT_ACTION(SubgridComponent, kickGravity, KickAction)

  private:
	Subgrid subgrid_;
};
} // namespace octotigerII
using SubgridComponentType = hpx::components::component<octotigerII::SubgridComponent>;
HPX_REGISTER_COMPONENT(SubgridComponentType, octotigerII_subgrid);
HPX_REGISTER_ACTION(octotigerII::SubgridComponent::SnapshotAction, octotigerII_snapshot);
HPX_REGISTER_ACTION(octotigerII::SubgridComponent::TimestepAction, octotigerII_timestep);
HPX_REGISTER_ACTION(octotigerII::SubgridComponent::AdvanceAction, octotigerII_advance);
HPX_REGISTER_ACTION(octotigerII::SubgridComponent::GravityAction, octotigerII_gravity);
HPX_REGISTER_ACTION(octotigerII::SubgridComponent::KickAction, octotigerII_kick);
#endif

namespace octotigerII {
struct Runtime::Impl {
#ifdef OCTOTIGERII_WITH_HPX
	std::vector<hpx::id_type> blocks;
#else
	std::vector<Subgrid> blocks;
#endif
};
Runtime::Runtime(Config const& c) : impl_(std::make_unique<Impl>()) {
	c.validate();
	int const n = 1 << c.level;
#ifdef OCTOTIGERII_WITH_HPX
	auto const localities = hpx::find_all_localities();
	if (localities.empty())
		throw std::runtime_error("No HPX localities");
#endif
	for (int z = 0; z < (c.dimensions > 2 ? n : 1); ++z)
		for (int y = 0; y < (c.dimensions > 1 ? n : 1); ++y)
			for (int x = 0; x < n; ++x) {
				mesh::BlockLocation const location{c.level, {x, y, z}, c.dimensions};
#ifdef OCTOTIGERII_WITH_HPX
				impl_->blocks.push_back(
					hpx::new_<SubgridComponent>(
						localities[impl_->blocks.size() % localities.size()], c, location)
						.get());
#else
				impl_->blocks.emplace_back(c, location);
#endif
			}
}
Runtime::~Runtime() = default;
std::size_t Runtime::size() const {
	return impl_->blocks.size();
}
char const* Runtime::backend() {
#ifdef OCTOTIGERII_WITH_HPX
	return "HPX";
#else
	return "serial CPU";
#endif
}
std::vector<Snapshot> Runtime::snapshots() const {
	std::vector<Snapshot> result;
#ifdef OCTOTIGERII_WITH_HPX
	std::vector<hpx::future<Snapshot>> pending;
	for (auto const& id : impl_->blocks)
		pending.push_back(hpx::async<SubgridComponent::SnapshotAction>(id));
	for (auto& future : pending)
		result.push_back(future.get());
#else
	for (auto const& block : impl_->blocks)
		result.push_back(block.snapshot());
#endif
	return result;
}
Real Runtime::stableTimestep() const {
	Real dt = std::numeric_limits<Real>::infinity();
#ifdef OCTOTIGERII_WITH_HPX
	std::vector<hpx::future<Real>> pending;
	for (auto const& id : impl_->blocks)
		pending.push_back(hpx::async<SubgridComponent::TimestepAction>(id));
	for (auto& future : pending)
		dt = std::min(dt, future.get());
#else
	for (auto const& block : impl_->blocks)
		dt = std::min(dt, block.stableTimestep());
#endif
	return dt;
}
void Runtime::advance(std::vector<Snapshot> const& old, Real dt) {
	std::vector<HydroSnapshot> hydro;
	std::vector<RadiationSnapshot> radiation;
	for (auto const& patch : old) {
		if (patch.hydroEnabled)
			hydro.push_back({patch.location, patch.hydro});
		if (patch.radiationEnabled)
			radiation.push_back({patch.location, patch.radiation});
	}
#ifdef OCTOTIGERII_WITH_HPX
	std::vector<hpx::future<void>> pending;
	for (auto const& id : impl_->blocks)
		pending.push_back(hpx::async<SubgridComponent::AdvanceAction>(id, hydro, radiation, dt));
	for (auto& future : pending)
		future.get();
#else
	for (auto& block : impl_->blocks)
		block.advance(hydro, radiation, dt);
#endif
}
void Runtime::setGravity(std::vector<std::vector<gravity::State>> const& fields) {
	if (fields.size() != size())
		throw std::invalid_argument("Gravity directory size mismatch");
#ifdef OCTOTIGERII_WITH_HPX
	std::vector<hpx::future<void>> pending;
	for (std::size_t i = 0; i < size(); ++i)
		pending.push_back(hpx::async<SubgridComponent::GravityAction>(impl_->blocks[i], fields[i]));
	for (auto& future : pending)
		future.get();
#else
	for (std::size_t i = 0; i < size(); ++i)
		impl_->blocks[i].setGravity(fields[i]);
#endif
}
void Runtime::kickGravity(Real dt) {
#ifdef OCTOTIGERII_WITH_HPX
	std::vector<hpx::future<void>> pending;
	for (auto const& id : impl_->blocks)
		pending.push_back(hpx::async<SubgridComponent::KickAction>(id, dt));
	for (auto& future : pending)
		future.get();
#else
	for (auto& block : impl_->blocks)
		block.kickGravity(dt);
#endif
}
} // namespace octotigerII
