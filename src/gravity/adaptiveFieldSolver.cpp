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
		std::size_t block = absent;
		bool leaf() const {
			return field.count != 0;
		}
	};
	class Tree {
	public:
		Tree(Config const& config, std::vector<Subgrid> const& blocks) {
			static_assert(ndim == 3);
			int const bits = std::countr_zero(unsigned(config.mesh.cells));
			std::unordered_map<Location, std::pair<storage::Range, std::size_t>, mesh::BlockLocationHash> locations;
			for (std::size_t b = 0; b < blocks.size(); ++b) {
				auto const& block = blocks[b];
				block.layout.forEachInterior([&](auto const& cell, std::size_t i) {
					Location location{block.location.level + bits, {}};
					for (int d = 0; d < 3; ++d)
						location.coordinates[d] = block.location.coordinates[d] * config.mesh.cells + cell[d];
					if (!locations.emplace(location, std::pair{block.interior.slice(i, 1), b}).second) throw std::logic_error("Overlapping gravity leaves");
					while (!location.isRoot()) {
						location = location.parent();
						locations.try_emplace(location);
					}
				});
			}
			for (auto const& [location, source] : locations)
				nodes.push_back({location, absent, {}, source.first, source.second});
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
		std::size_t block = 0;
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
	  , partitions_(partitions)
	  , blockCount_(blocks.size()) {
		auto const n = tree_->nodes.size();
		begin_ = (n * owner + partitions - 1) / partitions;
		end_ = (n * (owner + 1) + partitions - 1) / partitions;
		order_ = config.gravity.multipoleOrder;
		coefficients_ = diagonal::coefficientCount(order_) + (images_.active() ? 1 : 0);
		forceCoefficients_ = diagonal::coefficientCount(order_ + 1) + (images_.active() ? 1 : 0);
		length_ = units::value(config.mesh.upper - config.mesh.lower);
		moments_.resize(end_ - begin_);
		locals_.assign(end_ - begin_, Coefficients(coefficients_));
		forceLocals_.assign(end_ - begin_, Coefficients(forceCoefficients_));
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
			auto const block = tree_->nodes[id].block;
			if (segments_.empty() || segments_.back().block != block || segments_.back().range.partition != partition ||
				segments_.back().range.offset + segments_.back().range.count != offset)
				segments_.push_back({{partition, offset, 0}, {}, block});
			++segments_.back().range.count;
			segments_.back().nodes.push_back(id);
		}
		for (auto const& image : images_.images()) {
			walk(0, 0, image, false);
			if (images_.periodic()) walk(0, 0, image, true);
		}
	}
	void configure(FieldSolveRequest const& request, bool publishState) {
		if ((!request.sources.empty() && request.sources.size() != blockCount_) ||
			(!request.targets.empty() && request.targets.size() != blockCount_))
			throw std::invalid_argument("Gravity selection must contain one entry per block");
		for (auto const& source : request.sources)
			if (!std::isfinite(source.weight) || !std::isfinite(source.secondWeight))
				throw std::invalid_argument("Gravity density weights must be finite");
		request_ = request;
		publishState_ = publishState;
		activeTargets_.clear();
		if (!request.targets.empty()) {
			activeTargets_.assign(tree_->nodes.size(), 0);
			for (std::size_t i = tree_->nodes.size(); i-- > 0;) {
				auto const& node = tree_->nodes[i];
				if (node.leaf()) activeTargets_[i] = request.targets[node.block] != 0;
				if (activeTargets_[i] && node.parent != absent) activeTargets_[node.parent] = 1;
			}
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
	std::vector<Coefficients> read(unsigned field, std::vector<std::size_t> const& ids) const {
		std::vector<Coefficients> result;
		for (auto id : ids) {
			if (id < begin_ || id >= end_) throw std::out_of_range("Wrong adaptive FMM owner");
			result.push_back((field == 2 ? forceLocals_ : (field == 1 ? locals_ : moments_))[id - begin_]);
		}
		return result;
	}
#ifdef OCTOTIGERII_WITH_HPX
	HPX_DEFINE_COMPONENT_ACTION(AdaptiveFmmPartition, execute, ExecuteAction)
	HPX_DEFINE_COMPONENT_ACTION(AdaptiveFmmPartition, read, ReadAction)
	HPX_DEFINE_COMPONENT_ACTION(AdaptiveFmmPartition, configure, ConfigureAction)
#endif
private:
	Config config_;
	ImageGeometry images_;
	FieldDirectory fields_;
	std::unique_ptr<Tree> tree_;
	std::size_t owner_ = 0, partitions_ = 1, begin_ = 0, end_ = 0, workers_ = 1;
	std::size_t blockCount_ = 0;
	FieldSolveRequest request_;
	bool publishState_ = true;
	std::vector<unsigned char> activeTargets_;
	int order_ = 1, coefficients_ = 4, forceCoefficients_ = 9;
	Real length_ = 1;
	// Scalar locals preserve the reciprocal potential operator; the auxiliary
	// force locals retain degree order_+1 for symmetric degree-order_ gradients.
	std::vector<Coefficients> moments_, locals_, forceLocals_;
	std::vector<std::vector<Interaction>> interactions_;
	std::vector<Segment> segments_;
	bool targetActive(std::size_t id) const {
		return activeTargets_.empty() || activeTargets_[id];
	}

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
		if (!a.leaf() && !b.leaf() && a.location.level == b.location.level) {
			// Refine both equal-sized nodes together. Splitting only the target
			// can accept (child(a), b) before b is split, while the reversed walk
			// accepts (child(b), a). Those are different approximate kernels and
			// violate volume-weighted potential reciprocity on an adaptive mesh.
			for (auto targetChild : a.children)
				for (auto sourceChild : b.children)
					walk(targetChild, sourceChild, image, correction);
		} else if (!a.leaf() && (b.leaf() || a.location.level < b.location.level)) {
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
		unsigned field;
		std::unordered_map<std::size_t, Coefficients> remote;
		Coefficients const& at(std::size_t id) const {
			if (id >= owner.begin_ && id < owner.end_)
				return (field == 2 ? owner.forceLocals_ : (field == 1 ? owner.locals_ : owner.moments_))[id - owner.begin_];
			return remote.at(id);
		}
	};
	Sources fetch(std::unordered_set<std::size_t> const& requests, unsigned field, std::vector<storage::Locality> const& peers) const {
		Sources result{*this, field, {}};
#ifdef OCTOTIGERII_WITH_HPX
		std::vector<std::vector<std::size_t>> byOwner(partitions_);
		for (auto id : requests)
			if (id < begin_ || id >= end_) byOwner[id * partitions_ / tree_->nodes.size()].push_back(id);
		std::vector<std::vector<std::size_t>> batches;
		std::vector<hpx::future<std::vector<Coefficients>>> pending;
		auto const batchSize = std::max(std::size_t(1), std::size_t(65536) / (field == 2 ? forceCoefficients_ : coefficients_));
		for (std::size_t owner = 0; owner < partitions_; ++owner) {
			auto& ids = byOwner[owner];
			std::sort(ids.begin(), ids.end());
			for (std::size_t i = 0; i < ids.size(); i += batchSize) {
				batches.emplace_back(ids.begin() + i, ids.begin() + std::min(ids.size(), i + batchSize));
				pending.push_back(hpx::async<ReadAction>(peers[owner], field, batches.back()));
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
		for (auto& local : forceLocals_)
			std::fill(local.begin(), local.end(), 0);
		return parallel(segments_.size(), "gravity.adaptive.p2m", [&](std::size_t s, Statistics&) {
			auto const& segment = segments_[s];
			auto const& handle = request_.density.id || !request_.density.sumSources.empty() ? request_.density :
				(build::hydro && config_.hydroEnabled() ? std::get<0>(fields_.hydro.fields) : fields_.density);
			auto const source = request_.sources.empty() ? DensitySelection{bank, 0, 1, 0} : request_.sources[segment.block];
			storage::Buffer<units::Density> first, second;
			if (source.weight != 0) first = handle.read(segment.range, source.bank).get();
			if (source.secondWeight != 0) second = handle.read(segment.range, source.secondBank).get();
			for (std::size_t i = 0; i < segment.nodes.size(); ++i) {
				auto const id = segment.nodes[i];
				Real const h = width(id);
				auto const rho = (source.weight != 0 ? source.weight * first.data()[i] : units::Density{}) +
					(source.secondWeight != 0 ? source.secondWeight * second.data()[i] : units::Density{});
				Real const mass = units::value(rho) * h * h * h;
				if ((!request_.allowSignedDensity && !(mass >= 0)) || !isfinite(mass)) throw std::invalid_argument("Invalid adaptive gravity mass");
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
			if (!targetActive(id)) continue;
			for (auto const& pair : interactions_[id - begin_])
				needs.insert(pair.source);
			if (tree_->nodes[id].parent != absent) parentIds.insert(tree_->nodes[id].parent);
		}
		auto const sources = fetch(needs, false, peers), parents = fetch(parentIds, true, peers);
		auto const forceParents = fetch(parentIds, 2, peers);
		return parallel(stop - start, "gravity.adaptive.m2l_l2l", [&](std::size_t i, Statistics& stats) {
			auto const id = start + i;
			if (!targetActive(id)) return;
			auto& local = locals_[id - begin_];
			auto& forceLocal = forceLocals_[id - begin_];
			for (auto const& pair : interactions_[id - begin_]) {
				auto const& moment = sources.at(pair.source);
				if (std::all_of(moment.begin(), moment.end(), [](Real value) { return value == 0; })) continue;
				Real const unit = length_ / pair.count;
				Coefficients contribution(coefficients_);
				Coefficients forceContribution(forceCoefficients_);
				if (pair.direct) {
					if (pair.correction) {
						ewald::addDirect(contribution, moment[0], pair.separation, images_.periods(pair.count), unit);
						ewald::addDirect(forceContribution, moment[0], pair.separation, images_.periods(pair.count), unit);
					} else {
						diagonal::addDirect(contribution, moment[0], pair.separation, unit);
						diagonal::addDirect(forceContribution, moment[0], pair.separation, unit);
					}
					++stats.directPairs;
				} else {
					auto const normalized = rescale(reflectMultipole(moment, pair.reflection, order_), width(pair.source) / unit, order_);
					if (pair.correction) {
						ewald::getOperator(order_, pair.separation, images_.periods(pair.count))->add(contribution, normalized, unit);
						ewald::getOperator(order_, pair.separation, images_.periods(pair.count), true)->add(forceContribution, normalized, unit);
					} else {
						diagonal::getOperator(order_, pair.separation)->add(contribution, normalized, unit);
						diagonal::getOperator(order_, pair.separation, true)->add(forceContribution, normalized, unit);
					}
					++stats.multipolePairs;
				}
				add(local, rescale(std::move(contribution), width(id) / unit, order_));
				add(forceLocal, rescale(std::move(forceContribution), width(id) / unit, order_ + 1));
				if (pair.correction) ++stats.ewaldPairs;
				if (pair.reflection) ++stats.reflectedPairs;
			}
			if (tree_->nodes[id].parent != absent) {
				add(local, (images_.active() ? imageShiftLocal : diagonal::shiftLocal)(parents.at(tree_->nodes[id].parent), childOffset(id), 0.5, order_));
				add(forceLocal, (images_.active() ? imageShiftLocal : diagonal::shiftLocal)(forceParents.at(tree_->nodes[id].parent), childOffset(id), 0.5, order_ + 1));
			}
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
		auto const& destination = std::get<0>(request_.output.fields).id ? request_.output : fields_.gravity;
		auto const outputBank = publishState_ ? bank ^ 1 : request_.outputBank;
		auto result = parallel(segments_.size(), "gravity.adaptive.publish", [&](std::size_t s, Statistics&) {
			auto const& segment = segments_[s];
			if (!request_.targets.empty() && !request_.targets[segment.block]) return;
			auto output = destination.output(segment.range, outputBank);
			for (std::size_t i = 0; i < segment.nodes.size(); ++i) {
				auto const id = segment.nodes[i];
				auto const& local = locals_[id - begin_];
				auto const& forceLocal = forceLocals_[id - begin_];
				auto const g = constants::G * units::Mass::from_value(1) / units::Length::from_value(1);
				State field;
				field.potential() = g * local[0];
				for (int d = 0; d < 3; ++d) {
					diagonal::Offset e{};
					e[d] = 1;
					field.acceleration(d) = -g * diagonal::derivative(forceLocal, e[0], e[1], e[2]) / units::Length::from_value(width(id));
				}
				if (!finite(field)) throw std::runtime_error("Nonfinite adaptive FMM solution");
				output.put(i, field);
			}
			destination.commit(segment.range, outputBank, output);
			if (!publishState_) return;
			if (build::hydro && config_.hydroEnabled())
				copy(fields_.hydro, segment.range, bank);
			else {
				auto input = fields_.density.read(segment.range, bank).get();
				auto next = fields_.density.output(segment.range, bank ^ 1);
				std::copy_n(input.data(), segment.range.count, next.data());
				fields_.density.commit(segment.range, bank ^ 1, next);
			}
			if (build::radiation && config_.radiationEnabled()) copy(fields_.radiation, segment.range, bank);
			for (auto const& field : fields_.species) {
				auto input = field.read(segment.range, bank).get();
				auto output = field.output(segment.range, bank ^ 1);
				std::copy_n(input.data(), segment.range.count, output.data());
				field.commit(segment.range, bank ^ 1, output);
			}
		});
		std::uint64_t leaves = 0;
		for (auto const& segment : segments_)
			if (request_.targets.empty() || request_.targets[segment.block]) leaves += segment.nodes.size();
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
HPX_REGISTER_ACTION(octotigerII::gravity::AdaptiveFmmPartition::ConfigureAction, octotigerII_adaptive_fmm_configure)
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
	void configure(FieldSolveRequest const& request, bool publishState) {
#ifdef OCTOTIGERII_WITH_HPX
		std::vector<hpx::future<void>> pending;
		std::exception_ptr error;
		try {
			for (auto const& id : partitions)
				pending.push_back(hpx::async<AdaptiveFmmPartition::ConfigureAction>(id, request, publishState));
		} catch (...) { error = std::current_exception(); }
		for (auto& task : pending) {
			try { task.get(); } catch (...) { if (!error) error = std::current_exception(); }
		}
		if (error) std::rethrow_exception(error);
#else
		partition->configure(request, publishState);
#endif
	}
	Statistics solve(unsigned bank);
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
	impl_->configure({}, true);
	return impl_->solve(bank);
}
Statistics AdaptiveFieldSolver::solve(FieldSolveRequest const& request) {
	impl_->configure(request, false);
	return impl_->solve(request.sourceBank);
}
Statistics AdaptiveFieldSolver::Impl::solve(unsigned bank) {
	auto result = phase(AdaptiveStage::Initialize, 0, bank);
	for (unsigned depth = depths - 1; depth > 0; --depth)
		accumulate(result, phase(AdaptiveStage::Upward, depth - 1, bank));
	for (unsigned depth = 0; depth < depths; ++depth)
		accumulate(result, phase(AdaptiveStage::Downward, depth, bank));
	accumulate(result, phase(AdaptiveStage::Publish, 0, bank));
	return result;
}
}	 // namespace octotigerII::gravity
