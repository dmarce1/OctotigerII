// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#include "octotigerII/gravity/adaptiveFieldSolver.hpp"
#include <algorithm>
#include <atomic>
#include <bit>
#include <map>
#include <numeric>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include "octotigerII/gravity/ewald.hpp"
#include "octotigerII/gravity/images.hpp"
#include "octotigerII/profiling.hpp"
#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/include/components.hpp>
#include <hpx/runtime_local/get_os_thread_count.hpp>
#endif

namespace octotigerII::gravity {
namespace {
	using diagonal::Coefficients;
	using Location = mesh::BlockLocation;
	constexpr std::size_t absent = std::numeric_limits<std::size_t>::max();
	class Node {
	public:
		Location location;
		std::size_t parent = absent;
		std::array<std::size_t, 8> children{};
		storage::Range field;
		bool leaf() const {
			return field.count != 0;
		}
	};
	class Tree {
	public:
		Tree(Config const& config, std::vector<Subgrid> const& blocks) {
			static_assert(ndim == 3);
			int const bits = std::countr_zero(unsigned(config.mesh.cells));
			std::unordered_map<Location, storage::Range, mesh::BlockLocationHash> locations;
			for (auto const& block : blocks)
				block.layout.forEachInterior([&](auto const& cell, std::size_t i) {
					Location location{block.location.level + bits, {}};
					for (int d = 0; d < 3; ++d)
						location.coordinates[d] = block.location.coordinates[d] * config.mesh.cells + cell[d];
					if (!locations.emplace(location, block.interior.slice(i, 1)).second) throw std::logic_error("Overlapping gravity leaves");
					while (!location.isRoot()) {
						location = location.parent();
						locations.try_emplace(location);
					}
				});
			for (auto const& [location, range] : locations)
				nodes.push_back({location, absent, {}, range});
			std::sort(nodes.begin(), nodes.end(), [](auto const& a, auto const& b) {
				auto const& x = a.location;
				auto const& y = b.location;
				return std::tie(x.level, x.coordinates[2], x.coordinates[1], x.coordinates[0]) <
					std::tie(y.level, y.coordinates[2], y.coordinates[1], y.coordinates[0]);
			});
			std::unordered_map<Location, std::size_t, mesh::BlockLocationHash> ids;
			for (std::size_t i = 0; i < nodes.size(); ++i)
				ids.emplace(nodes[i].location, i);
			levels.resize(nodes.back().location.level + 2, nodes.size());
			for (std::size_t i = 0; i < nodes.size(); ++i) {
				auto& node = nodes[i];
				levels[node.location.level] = std::min(levels[node.location.level], i);
				if (!node.location.isRoot()) node.parent = ids.at(node.location.parent());
				if (!node.leaf())
					for (int slot = 0; slot < 8; ++slot)
						node.children[slot] = ids.at(node.location.child(slot));
			}
		}
		std::vector<Node> nodes;
		std::vector<std::size_t> levels;
	};
	class Interaction {
	public:
		std::size_t source = 0;
		diagonal::Offset separation{};
		int count = 1;
		unsigned reflection = 0;
		bool correction = false, direct = false;
	};
	class Segment {
	public:
		storage::Range range;
		std::vector<std::size_t> nodes;
	};
	enum class AdaptiveStage
	{
		Initialize,
		Upward,
		Downward,
		Publish
	};
	void add(Coefficients& to, Coefficients const& from) {
		for (std::size_t i = 0; i < to.size(); ++i)
			to[i] += from[i];
	}
	void accumulate(Statistics& a, Statistics const& b) {
		a.multipolePairs += b.multipolePairs;
		a.directPairs += b.directPairs;
		a.workerTasks += b.workerTasks;
		a.ewaldPairs += b.ewaldPairs;
		a.reflectedPairs += b.reflectedPairs;
		a.localityCells.insert(a.localityCells.end(), b.localityCells.begin(), b.localityCells.end());
	}
	Coefficients rescale(Coefficients value, Real ratio, int order) {
		auto const& indices = diagonal::indices(order);
		std::vector<Real> powers(order + 1, 1);
		for (int i = 1; i <= order; ++i)
			powers[i] = powers[i - 1] * ratio;
		for (std::size_t i = 0; i < std::min(indices.size(), value.size()); ++i)
			value[i] *= powers[indices[i].degree()];
		if (value.size() > indices.size()) value.back() *= ratio * ratio;
		return value;
	}
#ifdef OCTOTIGERII_WITH_HPX
	template <typename T>
	std::vector<T> collect(std::vector<hpx::future<T>>& pending) {
		std::exception_ptr error;
		std::vector<T> result;
		for (auto& future : pending) {
			try {
				result.push_back(future.get());
			} catch (...) {
				if (!error) error = std::current_exception();
			}
		}
		if (error) std::rethrow_exception(error);
		return result;
	}
#endif
}	 // namespace

class AdaptiveFmmPartition
#ifdef OCTOTIGERII_WITH_HPX
  : public hpx::components::component_base<AdaptiveFmmPartition>
#endif
{
public:
	AdaptiveFmmPartition() = default;
	AdaptiveFmmPartition(Config const& config, std::vector<Subgrid> const& blocks, FieldDirectory fields, std::size_t owner, std::size_t partitions)
	  : config_(config)
	  , images_(config.mesh.boundary)
	  , fields_(std::move(fields))
	  , tree_(std::make_unique<Tree>(config, blocks))
	  , owner_(owner)
	  , partitions_(partitions) {
		auto const n = tree_->nodes.size();
		begin_ = (n * owner + partitions - 1) / partitions;
		end_ = (n * (owner + 1) + partitions - 1) / partitions;
		order_ = config.gravity.multipoleOrder;
		coefficients_ = diagonal::coefficientCount(order_) + (images_.active() ? 1 : 0);
		length_ = units::value(config.mesh.upper - config.mesh.lower);
		moments_.resize(end_ - begin_);
		locals_.assign(end_ - begin_, Coefficients(coefficients_));
		interactions_.resize(end_ - begin_);
		for (auto i = begin_; i < end_; ++i)
			moments_[i - begin_].resize(tree_->nodes[i].leaf() ? 1 : coefficients_);
#ifdef OCTOTIGERII_WITH_HPX
		workers_ = config.runtime.workerTasks > 0 ? config.runtime.workerTasks : hpx::get_os_thread_count();
#endif
		std::vector<std::pair<std::pair<std::size_t, std::size_t>, std::size_t>> addresses;
		for (auto i = begin_; i < end_; ++i)
			if (tree_->nodes[i].leaf()) {
				auto const range = tree_->nodes[i].field;
				addresses.push_back({{range.partition, range.offset}, i});
			}
		std::sort(addresses.begin(), addresses.end());
		for (auto const& [address, id] : addresses) {
			auto const [partition, offset] = address;
			if (segments_.empty() || segments_.back().range.partition != partition || segments_.back().range.offset + segments_.back().range.count != offset)
				segments_.push_back({{partition, offset, 0}, {}});
			++segments_.back().range.count;
			segments_.back().nodes.push_back(id);
		}
		for (auto const& image : images_.images()) {
			walk(0, 0, image, false);
			if (images_.periodic()) walk(0, 0, image, true);
		}
	}
	Statistics execute(AdaptiveStage stage, unsigned depth, unsigned bank, std::vector<storage::Locality> const& peers) {
		switch (stage) {
		case AdaptiveStage::Initialize:
			return initialize(bank);
		case AdaptiveStage::Upward:
			return upward(depth, peers);
		case AdaptiveStage::Downward:
			return downward(depth, peers);
		case AdaptiveStage::Publish:
			return publish(bank);
		}
		throw std::logic_error("Unknown adaptive gravity stage");
	}
	std::vector<Coefficients> read(bool local, std::vector<std::size_t> const& ids) const {
		std::vector<Coefficients> result;
		for (auto id : ids) {
			if (id < begin_ || id >= end_) throw std::out_of_range("Wrong adaptive FMM owner");
			result.push_back((local ? locals_ : moments_)[id - begin_]);
		}
		return result;
	}
#ifdef OCTOTIGERII_WITH_HPX
	HPX_DEFINE_COMPONENT_ACTION(AdaptiveFmmPartition, execute, ExecuteAction)
	HPX_DEFINE_COMPONENT_ACTION(AdaptiveFmmPartition, read, ReadAction)
#endif
private:
	Config config_;
	ImageGeometry images_;
	FieldDirectory fields_;
	std::unique_ptr<Tree> tree_;
	std::size_t owner_ = 0, partitions_ = 1, begin_ = 0, end_ = 0, workers_ = 1;
	int order_ = 1, coefficients_ = 4;
	Real length_ = 1;
	std::vector<Coefficients> moments_, locals_;
	std::vector<std::vector<Interaction>> interactions_;
	std::vector<Segment> segments_;

	Real width(std::size_t id) const {
		return length_ / Real(1 << tree_->nodes[id].location.level);
	}
	diagonal::Vector childOffset(std::size_t id) const {
		auto const& c = tree_->nodes[id].location.coordinates;
		return {(c[0] & 1) ? 0.25 : -0.25, (c[1] & 1) ? 0.25 : -0.25, (c[2] & 1) ? 0.25 : -0.25};
	}
	Interaction geometry(std::size_t a, std::size_t b, SourceImage const& image, bool correction) const {
		auto const& x = tree_->nodes[a].location;
		auto const& y = tree_->nodes[b].location;
		int const level = std::max(x.level, y.level);
		int count = 1 << (level + 1);
		diagonal::Offset r{};
		auto const periods = images_.periods(count);
		for (int d = 0; d < 3; ++d) {
			int const target = (2 * x.coordinates[d] + 1) * (1 << (level - x.level));
			int source = (2 * y.coordinates[d] + 1) * (1 << (level - y.level));
			if (image.mask & (1u << d)) source = -source + image.translation[d] * count;
			r[d] = target - source;
			if (periods[d]) {
				r[d] %= periods[d];
				if (2 * r[d] > periods[d]) r[d] -= periods[d];
				if (2 * r[d] < -periods[d]) r[d] += periods[d];
			}
		}
		if (x.level == y.level) {
			count /= 2;
			for (auto& value : r)
				value /= 2;
		}
		return {b, r, count, image.mask, correction, tree_->nodes[a].leaf() && tree_->nodes[b].leaf()};
	}
	bool acceptable(std::size_t target, Interaction const& pair) const {
		using std::abs;
		using std::hypot;
		auto const& source = tree_->nodes[pair.source];
		auto const& sink = tree_->nodes[target];
		Real const a = Real(pair.count) / Real(1 << sink.location.level);
		Real const b = Real(pair.count) / Real(1 << source.location.level);
		auto const periods = images_.periods(pair.count);
		for (int d = 0; d < 3; ++d)
			if (periods[d] && 2 * abs(pair.separation[d]) + a + b > periods[d]) return false;
		Real const distance = pair.correction ? images_.nextImageDistance(pair.separation, pair.count) :
												hypot(Real(pair.separation[0]), Real(pair.separation[1]), Real(pair.separation[2]));
		return distance * config_.gravity.openingAngle > std::max(a, b);
	}
	void walk(std::size_t target, std::size_t source, SourceImage const& image, bool correction) {
		auto const& a = tree_->nodes[target];
		auto const& b = tree_->nodes[source];
		auto const pair = geometry(target, source, image, correction);
		if (pair.direct || acceptable(target, pair)) {
			if (!correction && pair.separation == diagonal::Offset{}) return;
			if (target >= begin_ && target < end_) interactions_[target - begin_].push_back(pair);
			return;
		}
		if (!a.leaf() && (b.leaf() || a.location.level <= b.location.level)) {
			for (auto child : a.children)
				walk(child, source, image, correction);
		} else {
			for (auto child : b.children)
				walk(target, child, image, correction);
		}
	}
	template <typename F>
	Statistics parallel(std::size_t count, char const* name, F const& f) const {
		if (!count) return {};
		std::atomic<std::size_t> next{0};
		auto work = [&] {
			Statistics result;
			result.workerTasks = 1;
			for (;;) {
				auto const start = next.fetch_add(16, std::memory_order_relaxed);
				if (start >= count) break;
				for (auto i = start; i < std::min(count, start + 16); ++i)
					f(i, result);
			}
			return result;
		};
#ifdef OCTOTIGERII_WITH_HPX
		std::vector<hpx::future<Statistics>> pending;
		for (std::size_t i = 0; i < std::min(workers_, (count + 15) / 16); ++i)
			pending.push_back(hpx::async(profiling::annotated(work, name)));
		Statistics result;
		for (auto const& part : collect(pending))
			accumulate(result, part);
		return result;
#else
		(void) name;
		return work();
#endif
	}
	class Sources {
	public:
		AdaptiveFmmPartition const& owner;
		bool local;
		std::unordered_map<std::size_t, Coefficients> remote;
		Coefficients const& at(std::size_t id) const {
			if (id >= owner.begin_ && id < owner.end_) return (local ? owner.locals_ : owner.moments_)[id - owner.begin_];
			return remote.at(id);
		}
	};
	Sources fetch(std::unordered_set<std::size_t> const& requests, bool local, std::vector<storage::Locality> const& peers) const {
		Sources result{*this, local, {}};
#ifdef OCTOTIGERII_WITH_HPX
		std::vector<std::vector<std::size_t>> byOwner(partitions_);
		for (auto id : requests)
			if (id < begin_ || id >= end_) byOwner[id * partitions_ / tree_->nodes.size()].push_back(id);
		std::vector<std::vector<std::size_t>> batches;
		std::vector<hpx::future<std::vector<Coefficients>>> pending;
		auto const batchSize = std::max(std::size_t(1), std::size_t(65536) / coefficients_);
		for (std::size_t owner = 0; owner < partitions_; ++owner) {
			auto& ids = byOwner[owner];
			std::sort(ids.begin(), ids.end());
			for (std::size_t i = 0; i < ids.size(); i += batchSize) {
				batches.emplace_back(ids.begin() + i, ids.begin() + std::min(ids.size(), i + batchSize));
				pending.push_back(hpx::async<ReadAction>(peers[owner], local, batches.back()));
			}
		}
		auto values = collect(pending);
		for (std::size_t i = 0; i < values.size(); ++i)
			for (std::size_t j = 0; j < values[i].size(); ++j)
				result.remote.emplace(batches[i][j], std::move(values[i][j]));
#else
		(void) requests;
		(void) peers;
#endif
		return result;
	}
	std::pair<std::size_t, std::size_t> levelRange(unsigned depth) const {
		auto const start = std::max(begin_, tree_->levels.at(depth));
		return {start, std::max(start, std::min(end_, tree_->levels.at(depth + 1)))};
	}
	Statistics initialize(unsigned bank) {
		using std::isfinite;
		for (auto& moment : moments_)
			std::fill(moment.begin(), moment.end(), 0);
		for (auto& local : locals_)
			std::fill(local.begin(), local.end(), 0);
		return parallel(segments_.size(), "gravity.adaptive.p2m", [&](std::size_t s, Statistics&) {
			auto const& segment = segments_[s];
			auto const handle = [&] {
				if (build::hydro && config_.hydroEnabled())
					return std::get<0>(fields_.hydro.fields);
				else
					return fields_.density;
			}();
			auto input = handle.read(segment.range, bank).get();
			for (std::size_t i = 0; i < segment.nodes.size(); ++i) {
				auto const id = segment.nodes[i];
				Real const h = width(id);
				Real const mass = units::value(input.data()[i]) * h * h * h;
				if (!(mass >= 0) || !isfinite(mass)) throw std::invalid_argument("Invalid adaptive gravity mass");
				moments_[id - begin_][0] = mass;
			}
		});
	}
	Statistics upward(unsigned depth, std::vector<storage::Locality> const& peers) {
		auto const [start, stop] = levelRange(depth);
		std::unordered_set<std::size_t> needs;
		for (auto id = start; id < stop; ++id)
			if (!tree_->nodes[id].leaf())
				for (auto child : tree_->nodes[id].children)
					needs.insert(child);
		auto const source = fetch(needs, false, peers);
		return parallel(stop - start, "gravity.adaptive.m2m", [&](std::size_t i, Statistics&) {
			auto const id = start + i;
			if (tree_->nodes[id].leaf()) return;
			for (auto child : tree_->nodes[id].children)
				add(moments_[id - begin_],
					(images_.active() ? imageShiftMultipole : diagonal::shiftMultipole)(source.at(child), childOffset(child), 0.5, order_));
		});
	}
	Statistics downward(unsigned depth, std::vector<storage::Locality> const& peers) {
		auto const [start, stop] = levelRange(depth);
		std::unordered_set<std::size_t> needs, parentIds;
		for (auto id = start; id < stop; ++id) {
			for (auto const& pair : interactions_[id - begin_])
				needs.insert(pair.source);
			if (tree_->nodes[id].parent != absent) parentIds.insert(tree_->nodes[id].parent);
		}
		auto const sources = fetch(needs, false, peers), parents = fetch(parentIds, true, peers);
		return parallel(stop - start, "gravity.adaptive.m2l_l2l", [&](std::size_t i, Statistics& stats) {
			auto const id = start + i;
			auto& local = locals_[id - begin_];
			for (auto const& pair : interactions_[id - begin_]) {
				auto const& moment = sources.at(pair.source);
				if (moment[0] == 0) continue;
				Real const unit = length_ / pair.count;
				Coefficients contribution(coefficients_);
				if (pair.direct) {
					if (pair.correction)
						ewald::addDirect(contribution, moment[0], pair.separation, images_.periods(pair.count), unit);
					else
						diagonal::addDirect(contribution, moment[0], pair.separation, unit);
					++stats.directPairs;
				} else {
					auto const normalized = rescale(reflectMultipole(moment, pair.reflection, order_), width(pair.source) / unit, order_);
					if (pair.correction)
						ewald::getOperator(order_, pair.separation, images_.periods(pair.count))->add(contribution, normalized, unit);
					else
						diagonal::getOperator(order_, pair.separation)->add(contribution, normalized, unit);
					++stats.multipolePairs;
				}
				add(local, rescale(std::move(contribution), width(id) / unit, order_));
				if (pair.correction) ++stats.ewaldPairs;
				if (pair.reflection) ++stats.reflectedPairs;
			}
			if (tree_->nodes[id].parent != absent)
				add(local, (images_.active() ? imageShiftLocal : diagonal::shiftLocal)(parents.at(tree_->nodes[id].parent), childOffset(id), 0.5, order_));
		});
	}
	template <typename StateType>
	void copy(storage::ColumnHandle<StateType> const& field, storage::Range range, unsigned bank) {
		auto input = field.read(range, bank).get();
		auto output = field.output(range, bank ^ 1);
		for (std::size_t i = 0; i < range.count; ++i)
			output.put(i, input.at(i));
		field.commit(range, bank ^ 1, output);
	}
	Statistics publish(unsigned bank) {
		auto result = parallel(segments_.size(), "gravity.adaptive.publish", [&](std::size_t s, Statistics&) {
			auto const& segment = segments_[s];
			auto output = fields_.gravity.output(segment.range, bank ^ 1);
			for (std::size_t i = 0; i < segment.nodes.size(); ++i) {
				auto const id = segment.nodes[i];
				auto const& local = locals_[id - begin_];
				auto const g = constants::G * units::Mass::from_value(1) / units::Length::from_value(1);
				State field;
				field.potential() = g * local[0];
				for (int d = 0; d < 3; ++d) {
					diagonal::Offset e{};
					e[d] = 1;
					field.acceleration(d) = -g * diagonal::derivative(local, e[0], e[1], e[2]) / units::Length::from_value(width(id));
				}
				if (!finite(field)) throw std::runtime_error("Nonfinite adaptive FMM solution");
				output.put(i, field);
			}
			fields_.gravity.commit(segment.range, bank ^ 1, output);
			if (build::hydro && config_.hydroEnabled())
				copy(fields_.hydro, segment.range, bank);
			else {
				auto input = fields_.density.read(segment.range, bank).get();
				auto next = fields_.density.output(segment.range, bank ^ 1);
				std::copy_n(input.data(), segment.range.count, next.data());
				fields_.density.commit(segment.range, bank ^ 1, next);
			}
			if (build::radiation && config_.radiationEnabled()) copy(fields_.radiation, segment.range, bank);
		});
		std::uint64_t leaves = 0;
		for (auto const& segment : segments_)
			leaves += segment.nodes.size();
		result.localityCells = {leaves};
		return result;
	}
};
}	 // namespace octotigerII::gravity

#ifdef OCTOTIGERII_WITH_HPX
using AdaptiveFmmComponent = hpx::components::component<octotigerII::gravity::AdaptiveFmmPartition>;
HPX_REGISTER_COMPONENT(AdaptiveFmmComponent, octotigerII_adaptive_fmm)
HPX_REGISTER_ACTION(octotigerII::gravity::AdaptiveFmmPartition::ExecuteAction, octotigerII_adaptive_fmm_execute)
HPX_REGISTER_ACTION(octotigerII::gravity::AdaptiveFmmPartition::ReadAction, octotigerII_adaptive_fmm_read)
#endif

namespace octotigerII::gravity {
class AdaptiveFieldSolver::Impl {
public:
	unsigned depths = 0;
#ifdef OCTOTIGERII_WITH_HPX
	std::vector<hpx::id_type> partitions;
#else
	std::unique_ptr<AdaptiveFmmPartition> partition;
#endif
	Statistics phase(AdaptiveStage stage, unsigned depth, unsigned bank) {
#ifdef OCTOTIGERII_WITH_HPX
		std::vector<hpx::future<Statistics>> pending;
		for (auto const& id : partitions)
			pending.push_back(hpx::async<AdaptiveFmmPartition::ExecuteAction>(id, stage, depth, bank, partitions));
		Statistics result;
		for (auto const& part : collect(pending))
			accumulate(result, part);
		return result;
#else
		return partition->execute(stage, depth, bank, {});
#endif
	}
};
AdaptiveFieldSolver::AdaptiveFieldSolver(
	Config const& config, std::vector<Subgrid> const& blocks, FieldDirectory const& fields, std::vector<storage::Locality> const& localities)
  : impl_(std::make_unique<Impl>()) {
	unsigned const bits = std::countr_zero(unsigned(config.mesh.cells));
	for (auto const& block : blocks)
		impl_->depths = std::max(impl_->depths, unsigned(block.location.level) + bits + 1);
#ifdef OCTOTIGERII_WITH_HPX
	std::vector<hpx::future<hpx::id_type>> pending;
	for (std::size_t i = 0; i < localities.size(); ++i)
		pending.push_back(hpx::new_<AdaptiveFmmPartition>(localities[i], config, blocks, fields, i, localities.size()));
	impl_->partitions = collect(pending);
#else
	(void) localities;
	impl_->partition = std::make_unique<AdaptiveFmmPartition>(config, blocks, fields, 0, 1);
#endif
}
AdaptiveFieldSolver::~AdaptiveFieldSolver() = default;
Statistics AdaptiveFieldSolver::solve(unsigned bank) {
	auto result = impl_->phase(AdaptiveStage::Initialize, 0, bank);
	for (unsigned depth = impl_->depths - 1; depth > 0; --depth)
		accumulate(result, impl_->phase(AdaptiveStage::Upward, depth - 1, bank));
	for (unsigned depth = 0; depth < impl_->depths; ++depth)
		accumulate(result, impl_->phase(AdaptiveStage::Downward, depth, bank));
	accumulate(result, impl_->phase(AdaptiveStage::Publish, 0, bank));
	return result;
}
}	 // namespace octotigerII::gravity
