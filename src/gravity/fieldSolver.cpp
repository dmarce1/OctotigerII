#include "octotigerII/gravity/fieldSolver.hpp"
#include "octotigerII/gravity/adaptiveFieldSolver.hpp"
#include "octotigerII/gravity/ewald.hpp"
#include "octotigerII/gravity/images.hpp"
#include "octotigerII/profiling.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <unordered_map>
#include <unordered_set>
#include "octotigerII/gravity/diagonal/fmm.hpp"
#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/include/components.hpp>
#include <hpx/runtime_local/get_os_thread_count.hpp>
#endif

namespace octotigerII::gravity {

namespace {

	using diagonal::Coefficients;
	using Coordinate = std::array<int, 3>;

	std::size_t cellCount(int n) {
		return std::size_t(n) * n * n;
	}

	std::size_t index(Coordinate const& c, int n) {
		return (std::size_t(c[2]) * n + c[1]) * n + c[0];
	}

	Coordinate coordinate(std::size_t id, int n) {
		return {int(id % n), int(id / n % n), int(id / (std::size_t(n) * n))};
	}

	Coordinate parent(Coordinate c) {
		for (auto& value : c)
			value /= 2;
		return c;
	}

	Coordinate child(Coordinate const& c, int slot) {
		return {2 * c[0] + (slot & 1), 2 * c[1] + ((slot >> 1) & 1), 2 * c[2] + ((slot >> 2) & 1)};
	}

	diagonal::Vector childOffset(int slot) {
		return {(slot & 1) ? 0.25 : -0.25, (slot & 2) ? 0.25 : -0.25, (slot & 4) ? 0.25 : -0.25};
	}

	void add(Coefficients& to, Coefficients const& from) {
		for (std::size_t i = 0; i < to.size(); ++i)
			to[i] += from[i];
	}

	void accumulate(Statistics& to, Statistics const& from) {
		to.multipolePairs += from.multipolePairs;
		to.directPairs += from.directPairs;
		to.workerTasks += from.workerTasks;
		to.ewaldPairs += from.ewaldPairs;
		to.reflectedPairs += from.reflectedPairs;
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

	class Level {
	public:
		Level(int count, Real width, int coefficients, int forceCoefficients, std::size_t owner, std::size_t partitions, bool leaf)
		  : count(count)
		  , width(width)
		  , begin((cellCount(count) * owner + partitions - 1) / partitions)
		  , end((cellCount(count) * (owner + 1) + partitions - 1) / partitions)
		  , moments(end - begin, Coefficients(leaf ? 1 : coefficients))
		  , locals(end - begin, Coefficients(coefficients))
		  , forceLocals(end - begin, Coefficients(forceCoefficients)) {}

		int count;
		Real width;
		std::size_t begin, end;
		// The auxiliary local has one extra degree: its gradient and the source
		// moments then have the same order under interaction reversal.
		std::vector<Coefficients> moments, locals, forceLocals;
	};

	class FieldRange {
	public:
		storage::Range range;
		std::size_t leaf = 0;
		std::size_t block = 0;
	};

	enum class Stage
	{
		Initialize,
		Upward,
		Downward,
		Publish
	};

}	 // namespace

/// Owns disjoint cells at every level. Read actions only see completed source
/// levels; numerical tasks exclusively own their target expansion.
class FmmPartition
#ifdef OCTOTIGERII_WITH_HPX
  : public hpx::components::component_base<FmmPartition>
#endif
{
public:
	FmmPartition() = default;

	FmmPartition(Config const& config, std::vector<Subgrid> const& blocks, FieldDirectory fields, std::size_t owner, std::size_t partitions)
	  : config_(config)
	  , images_(config.mesh.boundary)
	  , fields_(std::move(fields))
	  , owner_(owner)
	  , partitions_(partitions)
	  , blockCount_(blocks.size()) {
		static_assert(ndim == 3);
		config_.validate();
		n_ = config_.mesh.cells * (1 << config_.mesh.level);
		cellWidth_ = (config_.mesh.upper - config_.mesh.lower) / Real(n_);
		Real const h = units::value(cellWidth_);
		int const coefficients = diagonal::coefficientCount(config_.gravity.multipoleOrder) + (images_.active() ? 1 : 0);
		int const forceCoefficients = diagonal::coefficientCount(config_.gravity.multipoleOrder + 1) + (images_.active() ? 1 : 0);
		for (int count = 1; count <= n_; count *= 2)
			levels_.emplace_back(count, h * (n_ / count), coefficients, forceCoefficients, owner_, partitions_, count == n_);
#ifdef OCTOTIGERII_WITH_HPX
		workers_ = hpx::get_os_thread_count();
		if (config_.runtime.workerTasks > 0) workers_ = std::size_t(config_.runtime.workerTasks);
#endif
		workers_ = std::max(std::size_t(1), std::min(workers_, levels_.back().moments.size()));
		// Each range is contiguous in both the FMM leaf partition and a field
		// block. This also handles a partition boundary cutting through a block.
		auto const& leaves = levels_.back();
		for (std::size_t id = leaves.begin; id < leaves.end; ++id) {
			auto c = coordinate(id, n_);
			Coordinate b{};
			for (int d = 0; d < 3; ++d) {
				b[d] = c[d] / config_.mesh.cells;
				c[d] %= config_.mesh.cells;
			}
			auto const blockId = index(b, 1 << config_.mesh.level);
			auto const& block = blocks.at(blockId);
			auto const offset = block.interior.offset + block.layout.index(c);
			if (ranges_.empty() || ranges_.back().block != blockId || ranges_.back().range.partition != block.interior.partition ||
				ranges_.back().range.offset + ranges_.back().range.count != offset)
				ranges_.push_back({{block.interior.partition, offset, 0}, id - leaves.begin, blockId});
			++ranges_.back().range.count;
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
			for (auto const& level : levels_) activeTargets_.emplace_back(cellCount(level.count), 0);
			for (std::size_t id = 0; id < cellCount(n_); ++id) {
				auto c = coordinate(id, n_);
				for (auto& value : c) value /= config_.mesh.cells;
				activeTargets_.back()[id] = request.targets[index(c, 1 << config_.mesh.level)] != 0;
			}
			for (std::size_t depth = levels_.size() - 1; depth > 0; --depth)
				for (std::size_t id = 0; id < activeTargets_[depth].size(); ++id)
					if (activeTargets_[depth][id]) activeTargets_[depth - 1][index(parent(coordinate(id, levels_[depth].count)), levels_[depth - 1].count)] = 1;
		}
	}

	Statistics execute(Stage stage, unsigned depth, unsigned bank, std::vector<storage::Locality> const& peers) {
		switch (stage) {
		case Stage::Initialize:
			return initialize(bank);
		case Stage::Upward:
			return upward(depth, peers);
		case Stage::Downward:
			return downward(depth, peers);
		case Stage::Publish:
			return publish(bank);
		}
		throw std::logic_error("Unknown FMM stage");
	}

	std::vector<Coefficients> read(unsigned depth, unsigned field, std::vector<std::size_t> const& ids) const {
		profiling::Region profile("gravity.exchange.pack");
		auto const& level = levels_.at(depth);
		auto const& data = field == 2 ? level.forceLocals : (field == 1 ? level.locals : level.moments);
		std::vector<Coefficients> result;
		result.reserve(ids.size());
		for (auto id : ids) {
			if (id < level.begin || id >= level.end) throw std::out_of_range("FMM read sent to incorrect owner");
			result.push_back(data.at(id - level.begin));
		}
		return result;
	}

#ifdef OCTOTIGERII_WITH_HPX
	HPX_DEFINE_COMPONENT_ACTION(FmmPartition, execute, ExecuteAction)
	HPX_DEFINE_COMPONENT_ACTION(FmmPartition, read, ReadAction)
	HPX_DEFINE_COMPONENT_ACTION(FmmPartition, configure, ConfigureAction)
#endif

private:
	class Sources {
	public:
		Level const& level;
		unsigned field;
		std::unordered_map<std::size_t, Coefficients> remote;

		Coefficients const& at(std::size_t id) const {
			if (id >= level.begin && id < level.end)
				return (field == 2 ? level.forceLocals : (field == 1 ? level.locals : level.moments))[id - level.begin];
			return remote.at(id);
		}
	};

	Config config_;
	ImageGeometry images_;
	FieldDirectory fields_;
	std::size_t owner_ = 0, partitions_ = 1, workers_ = 1;
	std::size_t blockCount_ = 0;
	FieldSolveRequest request_;
	bool publishState_ = true;
	std::vector<std::vector<unsigned char>> activeTargets_;
	int n_ = 0;
	units::Length cellWidth_{};
	std::vector<Level> levels_;
	std::vector<FieldRange> ranges_;

	bool targetActive(unsigned depth, std::size_t id) const {
		return activeTargets_.empty() || activeTargets_[depth][id];
	}

	template <typename Function>
	Statistics parallel(std::size_t count, char const* name, Function const& function) const {
		if (count == 0) return {};
		constexpr std::size_t grain = 16;
		auto const workers = std::min(workers_, (count + grain - 1) / grain);
		std::atomic<std::size_t> next{0};
		auto work = [&](std::size_t worker) {
			Statistics result;
			result.workerTasks = 1;
			for (;;) {
				auto const begin = next.fetch_add(grain, std::memory_order_relaxed);
				if (begin >= count) break;
				for (auto i = begin; i < std::min(count, begin + grain); ++i)
					function(i, worker, result);
			}
			return result;
		};
		Statistics result;
#ifdef OCTOTIGERII_WITH_HPX
		std::vector<hpx::future<Statistics>> pending;
		pending.reserve(workers);
		try {
			for (std::size_t i = 0; i < workers; ++i)
				pending.push_back(hpx::async(profiling::annotated([&, i] { return work(i); }, name)));
		} catch (...) {
			auto error = std::current_exception();
			try {
				collect(pending);
			} catch (...) {
			}
			std::rethrow_exception(error);
		}
		for (auto const& part : collect(pending))
			accumulate(result, part);
#else
		(void) workers;
		(void) name;
		result = work(0);
#endif
		return result;
	}

	// A source belongs to the interaction list if its parent was not accepted.
	// With theta < 1/sqrt(3), an unaccepted parent implies all its ancestors
	// were unaccepted as well. Leaves perform the remaining direct interactions.
	template <typename Function>
	void interactions(unsigned depth, std::size_t target, Function const& function) const {
		using std::ceil;
		using std::hypot;
		if (!targetActive(depth, target)) return;

		auto const& level = levels_[depth];
		auto const a = coordinate(target, level.count);
		auto const p = parent(a);
		int const parents = level.count / 2;
		Real const theta = config_.gravity.openingAngle;
		int const radius = theta <= 1 / Real(parents) ? parents : int(ceil(1 / theta));
		bool const leaf = depth + 1 == levels_.size();
		for (int z = std::max(0, p[2] - radius); z <= std::min(parents - 1, p[2] + radius); ++z)
			for (int y = std::max(0, p[1] - radius); y <= std::min(parents - 1, p[1] + radius); ++y)
				for (int x = std::max(0, p[0] - radius); x <= std::min(parents - 1, p[0] + radius); ++x) {
					if (hypot(Real(p[0] - x), Real(p[1] - y), Real(p[2] - z)) * theta > 1) continue;
					for (int slot = 0; slot < 8; ++slot) {
						auto const b = child({x, y, z}, slot);
						auto const source = index(b, level.count);
						if (target == source) continue;
						diagonal::Offset const r{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
						if (leaf || hypot(Real(r[0]), Real(r[1]), Real(r[2])) * theta > 1) function(source, r);
					}
				}
	}

	// Plan on worker threads, deduplicate by owner, and launch all remote reads
	// before waiting. Only requested source cells are retained in the halo.
	template <typename Needs>
	Sources fetch(
		unsigned depth, unsigned field, std::size_t targets, Needs const& needs, std::vector<storage::Locality> const& peers, Statistics& statistics) const {
		profiling::Elapsed profile("gravity.exchange.wall_ns");
		auto const& level = levels_.at(depth);
		Sources sources{level, field, {}};
		if (partitions_ == 1) return sources;
		std::vector<std::unordered_set<std::size_t>> requests(workers_);
		accumulate(statistics, parallel(targets, "gravity.exchange.plan", [&](std::size_t i, std::size_t worker, Statistics&) {
			needs(i, [&](std::size_t id) {
				if (id < level.begin || id >= level.end) requests[worker].insert(id);
			});
		}));
		std::unordered_set<std::size_t> unique;
		for (auto& request : requests)
			unique.merge(request);
		std::vector<std::vector<std::size_t>> ids(partitions_);
		for (auto id : unique)
			ids.at(id * partitions_ / cellCount(level.count)).push_back(id);
		sources.remote.reserve(unique.size());
#ifdef OCTOTIGERII_WITH_HPX
		// Bound individual parcels, including the coefficient payload.
		std::size_t const batchSize = std::max(std::size_t(1), std::size_t(65536) / diagonal::coefficientCount(config_.gravity.multipoleOrder + (field == 2)));
		std::vector<std::vector<std::size_t>> batches;
		std::vector<hpx::future<std::vector<Coefficients>>> pending;
		try {
			for (std::size_t owner = 0; owner < partitions_; ++owner) {
				std::sort(ids[owner].begin(), ids[owner].end());
				for (std::size_t begin = 0; begin < ids[owner].size(); begin += batchSize) {
					batches.emplace_back(ids[owner].begin() + begin, ids[owner].begin() + std::min(ids[owner].size(), begin + batchSize));
					pending.push_back(hpx::async<ReadAction>(peers.at(owner), depth, field, batches.back()));
				}
			}
		} catch (...) {
			auto error = std::current_exception();
			try {
				collect(pending);
			} catch (...) {
			}
			std::rethrow_exception(error);
		}
		auto values = [&] {
			profiling::Elapsed wait("gravity.exchange.wait.wall_ns");
			return collect(pending);
		}();
		profiling::Region unpack("gravity.exchange.unpack");
		for (std::size_t b = 0; b < batches.size(); ++b) {
			if (values[b].size() != batches[b].size()) throw std::logic_error("Incomplete FMM source transfer");
			for (std::size_t i = 0; i < batches[b].size(); ++i)
				sources.remote.emplace(batches[b][i], std::move(values[b][i]));
		}
#else
		(void) peers;
#endif
		return sources;
	}

	Statistics initialize(unsigned bank) {
		profiling::Elapsed profile("gravity.initialize.wall_ns");
		Statistics result;
		for (auto& level : levels_)
			accumulate(result, parallel(level.moments.size(), "gravity.clear", [&](std::size_t i, std::size_t, Statistics&) {
				std::fill(level.moments[i].begin(), level.moments[i].end(), 0);
				std::fill(level.locals[i].begin(), level.locals[i].end(), 0);
				std::fill(level.forceLocals[i].begin(), level.forceLocals[i].end(), 0);
			}));
		auto const volume = cellWidth_ * cellWidth_ * cellWidth_;
		auto const gram = units::Mass::from_value(1);
		accumulate(result, parallel(ranges_.size(), "gravity.p2m", [&](std::size_t r, std::size_t, Statistics&) {
			auto const& segment = ranges_[r];
			auto fill = [&](auto density) {
				for (std::size_t i = 0; i < segment.range.count; ++i) {
					auto const rho = density(i);
					if ((!request_.allowSignedDensity && !(rho >= units::Density{})) || !units::finite(rho) || !units::finite(rho * volume))
						throw std::invalid_argument("Gravity density/mass must be finite and nonnegative unless signed sources are enabled");
					levels_.back().moments[segment.leaf + i][0] = rho * volume / gram;
				}
			};
			auto const& handle = request_.density.id || !request_.density.sumSources.empty() ? request_.density :
				(build::hydro && config_.hydroEnabled() ? std::get<0>(fields_.hydro.fields) : fields_.density);
			auto const source = request_.sources.empty() ? DensitySelection{bank, 0, 1, 0} : request_.sources[segment.block];
			storage::Buffer<units::Density> first, second;
			if (source.weight != 0) first = handle.read(segment.range, source.bank).get();
			if (source.secondWeight != 0) second = handle.read(segment.range, source.secondBank).get();
			fill([&](std::size_t i) {
				return (source.weight != 0 ? source.weight * first.data()[i] : units::Density{}) +
					(source.secondWeight != 0 ? source.secondWeight * second.data()[i] : units::Density{});
			});
		}));
		return result;
	}

	Statistics upward(unsigned depth, std::vector<storage::Locality> const& peers) {
		profiling::Elapsed profile("gravity.upward.wall_ns");
		auto& level = levels_.at(depth);
		auto const& fine = levels_.at(depth + 1);
		Statistics result;
		auto children = fetch(
			depth + 1, false, level.moments.size(),
			[&](std::size_t i, auto need) {
				auto const c = coordinate(level.begin + i, level.count);
				for (int slot = 0; slot < 8; ++slot)
					need(index(child(c, slot), fine.count));
			},
			peers, result);
		accumulate(result, parallel(level.moments.size(), "gravity.m2m", [&](std::size_t i, std::size_t, Statistics&) {
			auto const c = coordinate(level.begin + i, level.count);
			for (int slot = 0; slot < 8; ++slot)
				add(level.moments[i],
					(images_.active() ? imageShiftMultipole : diagonal::shiftMultipole)(
						children.at(index(child(c, slot), fine.count)), childOffset(slot), 0.5, config_.gravity.multipoleOrder));
		}));
		return result;
	}

	Statistics imageDownward(unsigned depth, std::vector<storage::Locality> const& peers) {
		auto& level = levels_.at(depth);
		Statistics result;
		bool const leaf = depth + 1 == levels_.size();
		int const order = config_.gravity.multipoleOrder;
		auto interactionsFor = [&](std::size_t target, auto const& f) {
			if (!targetActive(depth, target)) return;
			images_.interactions(depth, coordinate(target, level.count), leaf, config_.gravity.openingAngle,
				[&](Coordinate b, auto r, unsigned mask, bool correction) { f(index(b, level.count), r, mask, correction); });
		};
		auto sources = fetch(
			depth, false, level.locals.size(),
			[&](std::size_t i, auto need) { interactionsFor(level.begin + i, [&](std::size_t source, auto, unsigned, bool) { need(source); }); }, peers,
			result);
		// Root corrections are allowed; they have no parent local to fetch.
		auto parents = depth ?
			fetch(
				depth - 1, true, level.locals.size(),
				[&](std::size_t i, auto need) { if (targetActive(depth, level.begin + i)) need(index(parent(coordinate(level.begin + i, level.count)), level.count / 2)); }, peers, result) :
			Sources{level, true, {}};
		auto forceParents = depth ?
			fetch(depth - 1, 2, level.locals.size(),
				[&](std::size_t i, auto need) { if (targetActive(depth, level.begin + i)) need(index(parent(coordinate(level.begin + i, level.count)), level.count / 2)); }, peers, result) :
			Sources{level, 2, {}};
		accumulate(result, parallel(level.locals.size(), "gravity.images.m2l_l2l", [&](std::size_t i, std::size_t, Statistics& statistics) {
			if (!targetActive(depth, level.begin + i)) return;
			auto& local = level.locals[i];
			auto& forceLocal = level.forceLocals[i];
			interactionsFor(level.begin + i, [&](std::size_t source, diagonal::Offset r, unsigned mask, bool correction) {
				auto const& moment = sources.at(source);
				if (std::all_of(moment.begin(), moment.end(), [](Real value) { return value == 0; })) return;
				if (leaf) {
					if (correction) {
						ewald::addDirect(local, moment[0], r, images_.periods(level.count), level.width);
						ewald::addDirect(forceLocal, moment[0], r, images_.periods(level.count), level.width);
					} else {
						diagonal::addDirect(local, moment[0], r, level.width);
						diagonal::addDirect(forceLocal, moment[0], r, level.width);
					}
					++statistics.directPairs;
				} else {
					auto const reflected = reflectMultipole(moment, mask, order);
					if (correction) {
						ewald::getOperator(order, r, images_.periods(level.count))->add(local, reflected, level.width);
						ewald::getOperator(order, r, images_.periods(level.count), true)->add(forceLocal, reflected, level.width);
					} else {
						diagonal::getOperator(order, r)->add(local, reflected, level.width);
						diagonal::getOperator(order, r, true)->add(forceLocal, reflected, level.width);
					}
					++statistics.multipolePairs;
				}
				if (correction) ++statistics.ewaldPairs;
				if (mask) ++statistics.reflectedPairs;
			});
			if (depth) {
				auto const c = coordinate(level.begin + i, level.count);
				int const slot = (c[0] & 1) + 2 * (c[1] & 1) + 4 * (c[2] & 1);
				add(local, imageShiftLocal(parents.at(index(parent(c), level.count / 2)), childOffset(slot), 0.5, order));
				add(forceLocal, imageShiftLocal(forceParents.at(index(parent(c), level.count / 2)), childOffset(slot), 0.5, order + 1));
			}
		}));
		return result;
	}

	Statistics downward(unsigned depth, std::vector<storage::Locality> const& peers) {
		if (images_.active()) return imageDownward(depth, peers);
		if (depth == 0) return {};
		profiling::Elapsed profile("gravity.downward.wall_ns");
		auto& level = levels_.at(depth);
		Statistics result;
		auto sources = fetch(
			depth, false, level.locals.size(),
			[&](std::size_t i, auto need) { interactions(depth, level.begin + i, [&](std::size_t id, auto const&) { need(id); }); }, peers, result);
		auto parents = fetch(
			depth - 1, true, level.locals.size(),
			[&](std::size_t i, auto need) { if (targetActive(depth, level.begin + i)) need(index(parent(coordinate(level.begin + i, level.count)), level.count / 2)); }, peers, result);
		auto forceParents = fetch(
			depth - 1, 2, level.locals.size(),
			[&](std::size_t i, auto need) { if (targetActive(depth, level.begin + i)) need(index(parent(coordinate(level.begin + i, level.count)), level.count / 2)); }, peers, result);
		bool const leaf = depth + 1 == levels_.size();
		accumulate(result, parallel(level.locals.size(), leaf ? "gravity.p2p_l2l" : "gravity.m2l_l2l", [&](std::size_t i, std::size_t, Statistics& statistics) {
			auto const target = level.begin + i;
			if (!targetActive(depth, target)) return;
			interactions(depth, target, [&](std::size_t source, diagonal::Offset const& r) {
				auto const& moment = sources.at(source);
				if (!publishState_ && std::all_of(moment.begin(), moment.end(), [](Real value) { return value == 0; })) return;
				if (leaf) {
					diagonal::addDirect(level.locals[i], moment[0], r, level.width);
					diagonal::addDirect(level.forceLocals[i], moment[0], r, level.width);
					if (!publishState_ || target < source) ++statistics.directPairs;
				} else {
					diagonal::getOperator(config_.gravity.multipoleOrder, r)->add(level.locals[i], moment, level.width);
					diagonal::getOperator(config_.gravity.multipoleOrder, r, true)->add(level.forceLocals[i], moment, level.width);
					if (!publishState_ || target < source) ++statistics.multipolePairs;
				}
			});
			auto const c = coordinate(target, level.count);
			int const slot = (c[0] & 1) + 2 * (c[1] & 1) + 4 * (c[2] & 1);
			add(level.locals[i], diagonal::shiftLocal(parents.at(index(parent(c), level.count / 2)), childOffset(slot), 0.5, config_.gravity.multipoleOrder));
			add(level.forceLocals[i], diagonal::shiftLocal(forceParents.at(index(parent(c), level.count / 2)), childOffset(slot), 0.5, config_.gravity.multipoleOrder + 1));
		}));
		return result;
	}

	template <typename StateType>
	void copy(storage::ColumnHandle<StateType> const& field, storage::Range range, unsigned bank) const {
		auto input = field.read(range, bank).get();
		auto output = field.output(range, bank ^ 1);
		for (std::size_t i = 0; i < range.count; ++i)
			output.put(i, input.at(i));
		field.commit(range, bank ^ 1, output);
	}

	Statistics publish(unsigned bank) {
		profiling::Elapsed profile("gravity.publish.wall_ns");
		auto const gram = units::Mass::from_value(1);
		auto const cm = units::Length::from_value(1);
		auto const& destination = std::get<0>(request_.output.fields).id ? request_.output : fields_.gravity;
		auto const outputBank = publishState_ ? bank ^ 1 : request_.outputBank;
		auto result = parallel(ranges_.size(), "gravity.publish", [&](std::size_t r, std::size_t, Statistics&) {
			auto const& segment = ranges_[r];
			if (!request_.targets.empty() && !request_.targets[segment.block]) return;
			auto const range = segment.range;
			auto output = destination.output(range, outputBank);
			for (std::size_t i = 0; i < range.count; ++i) {
				auto const& local = levels_.back().locals[segment.leaf + i];
				auto const& forceLocal = levels_.back().forceLocals[segment.leaf + i];
				State field;
				field.potential() = constants::G * (gram / cm) * local[0];
				field.acceleration(0) = -constants::G * (gram / cm) * diagonal::derivative(forceLocal, 1, 0, 0) / cellWidth_;
				field.acceleration(1) = -constants::G * (gram / cm) * diagonal::derivative(forceLocal, 0, 1, 0) / cellWidth_;
				field.acceleration(2) = -constants::G * (gram / cm) * diagonal::derivative(forceLocal, 0, 0, 1) / cellWidth_;
				if (!finite(field)) throw std::runtime_error("Nonfinite gravity solution");
				output.put(i, field);
			}
			destination.commit(range, outputBank, output);
			if (!publishState_) return;
			if (build::hydro && config_.hydroEnabled())
				copy(fields_.hydro, range, bank);
			else {
				auto input = fields_.density.read(range, bank).get();
				auto density = fields_.density.output(range, bank ^ 1);
				std::copy_n(input.data(), range.count, density.data());
				fields_.density.commit(range, bank ^ 1, density);
			}
			if (build::radiation && config_.radiationEnabled()) copy(fields_.radiation, range, bank);
			for (auto const& field : fields_.species) {
				auto input = field.read(range, bank).get();
				auto output = field.output(range, bank ^ 1);
				std::copy_n(input.data(), range.count, output.data());
				field.commit(range, bank ^ 1, output);
			}
		});
		std::uint64_t targets = 0;
		for (auto const& segment : ranges_)
			if (request_.targets.empty() || request_.targets[segment.block]) targets += segment.range.count;
		result.localityCells = {targets};
		return result;
	}
};

}	 // namespace octotigerII::gravity

#ifdef OCTOTIGERII_WITH_HPX
using FmmComponent = hpx::components::component<octotigerII::gravity::FmmPartition>;
HPX_REGISTER_COMPONENT(FmmComponent, octotigerII_fmm)
HPX_REGISTER_ACTION(octotigerII::gravity::FmmPartition::ExecuteAction, octotigerII_fmm_execute)
HPX_REGISTER_ACTION(octotigerII::gravity::FmmPartition::ReadAction, octotigerII_fmm_read)
HPX_REGISTER_ACTION(octotigerII::gravity::FmmPartition::ConfigureAction, octotigerII_fmm_configure)
#endif

namespace octotigerII::gravity {

class FieldSolver::Impl {
public:
	std::unique_ptr<AdaptiveFieldSolver> adaptive;
	unsigned depths = 0;
#ifdef OCTOTIGERII_WITH_HPX
	std::vector<hpx::id_type> partitions;
#else
	std::unique_ptr<FmmPartition> partition;
#endif
	void configure(FieldSolveRequest const& request, bool publishState) {
#ifdef OCTOTIGERII_WITH_HPX
		std::vector<hpx::future<void>> pending;
		std::exception_ptr error;
		try {
			for (auto const& id : partitions)
				pending.push_back(hpx::async<FmmPartition::ConfigureAction>(id, request, publishState));
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

	Statistics phase(Stage stage, unsigned depth, unsigned bank) {
#ifdef OCTOTIGERII_WITH_HPX
		std::vector<hpx::future<Statistics>> pending;
		pending.reserve(partitions.size());
		try {
			for (auto const& id : partitions)
				pending.push_back(hpx::async<FmmPartition::ExecuteAction>(id, stage, depth, bank, partitions));
		} catch (...) {
			auto error = std::current_exception();
			try {
				collect(pending);
			} catch (...) {
			}
			std::rethrow_exception(error);
		}
		Statistics result;
		for (auto const& part : collect(pending)) {
			accumulate(result, part);
			result.localityCells.insert(result.localityCells.end(), part.localityCells.begin(), part.localityCells.end());
		}
		return result;
#else
		return partition->execute(stage, depth, bank, {});
#endif
	}
};

FieldSolver::FieldSolver(
	Config const& config, std::vector<Subgrid> const& blocks, FieldDirectory const& fields, std::vector<storage::Locality> const& localities)
  : impl_(std::make_unique<Impl>()) {
	profiling::Elapsed profile("gravity.setup.wall_ns");
	config.validate();
	if (localities.empty()) throw std::invalid_argument("FMM requires at least one locality");
	if (config.amr.enabled) {
		impl_->adaptive = std::make_unique<AdaptiveFieldSolver>(config, blocks, fields, localities);
		return;
	}
	int const n = config.mesh.cells * (1 << config.mesh.level);
	for (int count = 1; count <= n; count *= 2)
		++impl_->depths;
#ifdef OCTOTIGERII_WITH_HPX
	std::vector<hpx::future<hpx::id_type>> pending;
	pending.reserve(localities.size());
	try {
		for (std::size_t i = 0; i < localities.size(); ++i)
			pending.push_back(hpx::new_<FmmPartition>(localities[i], config, blocks, fields, i, localities.size()));
	} catch (...) {
		auto error = std::current_exception();
		try {
			collect(pending);
		} catch (...) {
		}
		std::rethrow_exception(error);
	}
	impl_->partitions = collect(pending);
#else
	(void) localities;
	impl_->partition = std::make_unique<FmmPartition>(config, blocks, fields, 0, 1);
#endif
}

FieldSolver::~FieldSolver() = default;

Statistics FieldSolver::solve(unsigned bank) {
	if (impl_->adaptive) return impl_->adaptive->solve(bank);
	impl_->configure({}, true);
	return impl_->solve(bank);
}

Statistics FieldSolver::solve(FieldSolveRequest const& request) {
	if (impl_->adaptive) return impl_->adaptive->solve(request);
	impl_->configure(request, false);
	return impl_->solve(request.sourceBank);
}

Statistics FieldSolver::Impl::solve(unsigned bank) {
	profiling::Elapsed profile("gravity.solve.wall_ns");
	auto result = phase(Stage::Initialize, 0, bank);
	for (unsigned depth = depths - 1; depth > 0; --depth)
		accumulate(result, phase(Stage::Upward, depth - 1, bank));
	for (unsigned depth = 0; depth < depths; ++depth)
		accumulate(result, phase(Stage::Downward, depth, bank));
	auto publication = phase(Stage::Publish, 0, bank);
	accumulate(result, publication);
	result.localityCells = std::move(publication.localityCells);
	profiling::sample("gravity.multipole_pairs", double(result.multipolePairs));
	profiling::sample("gravity.direct_pairs", double(result.directPairs));
	profiling::sample("gravity.worker_tasks", double(result.workerTasks));
	profiling::sample("gravity.ewald_pairs", double(result.ewaldPairs));
	profiling::sample("gravity.reflected_pairs", double(result.reflectedPairs));
	return result;
}

}	 // namespace octotigerII::gravity
