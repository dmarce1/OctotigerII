#include "octotigerII/gravity/fieldSolver.hpp"
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
		Level(int count, Real width, int coefficients, std::size_t owner, std::size_t partitions, bool leaf)
		  : count(count)
		  , width(width)
		  , begin((cellCount(count) * owner + partitions - 1) / partitions)
		  , end((cellCount(count) * (owner + 1) + partitions - 1) / partitions)
		  , moments(end - begin, Coefficients(leaf ? 1 : coefficients))
		  , locals(end - begin, Coefficients(coefficients)) {}

		int count;
		Real width;
		std::size_t begin, end;
		std::vector<Coefficients> moments, locals;
	};

	class FieldRange {
	public:
		storage::Range range;
		std::size_t leaf = 0;
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
	  , partitions_(partitions) {
		static_assert(ndim == 3);
		config_.validate();
		n_ = config_.mesh.cells * (1 << config_.mesh.level);
		cellWidth_ = (config_.mesh.upper - config_.mesh.lower) / Real(n_);
		Real const h = units::value(cellWidth_);
		int const coefficients = diagonal::coefficientCount(config_.gravity.multipoleOrder) + (images_.active() ? 1 : 0);
		for (int count = 1; count <= n_; count *= 2)
			levels_.emplace_back(count, h * (n_ / count), coefficients, owner_, partitions_, count == n_);
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
			auto const& block = blocks.at(index(b, 1 << config_.mesh.level));
			auto const offset = block.interior.offset + block.layout.index(c);
			if (ranges_.empty() || ranges_.back().range.partition != block.interior.partition ||
				ranges_.back().range.offset + ranges_.back().range.count != offset)
				ranges_.push_back({{block.interior.partition, offset, 0}, id - leaves.begin});
			++ranges_.back().range.count;
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

	std::vector<Coefficients> read(unsigned depth, bool local, std::vector<std::size_t> const& ids) const {
		profiling::Region profile("gravity.exchange.pack");
		auto const& level = levels_.at(depth);
		auto const& data = local ? level.locals : level.moments;
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
#endif

private:
	class Sources {
	public:
		Level const& level;
		bool local;
		std::unordered_map<std::size_t, Coefficients> remote;

		Coefficients const& at(std::size_t id) const {
			if (id >= level.begin && id < level.end) return (local ? level.locals : level.moments)[id - level.begin];
			return remote.at(id);
		}
	};

	Config config_;
	ImageGeometry images_;
	FieldDirectory fields_;
	std::size_t owner_ = 0, partitions_ = 1, workers_ = 1;
	int n_ = 0;
	units::Length cellWidth_{};
	std::vector<Level> levels_;
	std::vector<FieldRange> ranges_;

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
		unsigned depth, bool local, std::size_t targets, Needs const& needs, std::vector<storage::Locality> const& peers, Statistics& statistics) const {
		profiling::Elapsed profile("gravity.exchange.wall_ns");
		auto const& level = levels_.at(depth);
		Sources sources{level, local, {}};
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
		std::size_t const batchSize = std::max(std::size_t(1), std::size_t(65536) / diagonal::coefficientCount(config_.gravity.multipoleOrder));
		std::vector<std::vector<std::size_t>> batches;
		std::vector<hpx::future<std::vector<Coefficients>>> pending;
		try {
			for (std::size_t owner = 0; owner < partitions_; ++owner) {
				std::sort(ids[owner].begin(), ids[owner].end());
				for (std::size_t begin = 0; begin < ids[owner].size(); begin += batchSize) {
					batches.emplace_back(ids[owner].begin() + begin, ids[owner].begin() + std::min(ids[owner].size(), begin + batchSize));
					pending.push_back(hpx::async<ReadAction>(peers.at(owner), depth, local, batches.back()));
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
			}));
		auto const volume = cellWidth_ * cellWidth_ * cellWidth_;
		auto const gram = units::Mass::from_value(1);
		accumulate(result, parallel(ranges_.size(), "gravity.p2m", [&](std::size_t r, std::size_t, Statistics&) {
			auto const& segment = ranges_[r];
			auto fill = [&](auto density) {
				for (std::size_t i = 0; i < segment.range.count; ++i) {
					auto const rho = density(i);
					if (!(rho >= units::Density{}) || !units::finite(rho) || !units::finite(rho * volume))
						throw std::invalid_argument("Gravity density/mass must be finite and nonnegative");
					levels_.back().moments[segment.leaf + i][0] = rho * volume / gram;
				}
			};
			if constexpr (build::hydro) {
				auto input = std::get<0>(fields_.hydro.fields).read(segment.range, bank).get();
				fill([&](std::size_t i) { return input.data()[i]; });
			} else {
				auto input = fields_.density.read(segment.range, bank).get();
				fill([&](std::size_t i) { return input.data()[i]; });
			}
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
				[&](std::size_t i, auto need) { need(index(parent(coordinate(level.begin + i, level.count)), level.count / 2)); }, peers, result) :
			Sources{level, true, {}};
		accumulate(result, parallel(level.locals.size(), "gravity.images.m2l_l2l", [&](std::size_t i, std::size_t, Statistics& statistics) {
			auto& local = level.locals[i];
			interactionsFor(level.begin + i, [&](std::size_t source, diagonal::Offset r, unsigned mask, bool correction) {
				auto const& moment = sources.at(source);
				if (moment[0] == 0) return;
				if (leaf) {
					if (correction)
						ewald::addDirect(local, moment[0], r, images_.periods(level.count), level.width);
					else
						diagonal::addDirect(local, moment[0], r, level.width);
					++statistics.directPairs;
				} else {
					auto const reflected = reflectMultipole(moment, mask, order);
					if (correction)
						ewald::getOperator(order, r, images_.periods(level.count))->add(local, reflected, level.width);
					else
						diagonal::getOperator(order, r)->add(local, reflected, level.width);
					++statistics.multipolePairs;
				}
				if (correction) ++statistics.ewaldPairs;
				if (mask) ++statistics.reflectedPairs;
			});
			if (depth) {
				auto const c = coordinate(level.begin + i, level.count);
				int const slot = (c[0] & 1) + 2 * (c[1] & 1) + 4 * (c[2] & 1);
				add(local, imageShiftLocal(parents.at(index(parent(c), level.count / 2)), childOffset(slot), 0.5, order));
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
			[&](std::size_t i, auto need) { need(index(parent(coordinate(level.begin + i, level.count)), level.count / 2)); }, peers, result);
		bool const leaf = depth + 1 == levels_.size();
		accumulate(result, parallel(level.locals.size(), leaf ? "gravity.p2p_l2l" : "gravity.m2l_l2l", [&](std::size_t i, std::size_t, Statistics& statistics) {
			auto const target = level.begin + i;
			interactions(depth, target, [&](std::size_t source, diagonal::Offset const& r) {
				if (leaf) {
					diagonal::addDirect(level.locals[i], sources.at(source)[0], r, level.width);
					if (target < source) ++statistics.directPairs;
				} else {
					diagonal::getOperator(config_.gravity.multipoleOrder, r)->add(level.locals[i], sources.at(source), level.width);
					if (target < source) ++statistics.multipolePairs;
				}
			});
			auto const c = coordinate(target, level.count);
			int const slot = (c[0] & 1) + 2 * (c[1] & 1) + 4 * (c[2] & 1);
			add(level.locals[i], diagonal::shiftLocal(parents.at(index(parent(c), level.count / 2)), childOffset(slot), 0.5, config_.gravity.multipoleOrder));
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
		auto result = parallel(ranges_.size(), "gravity.publish", [&](std::size_t r, std::size_t, Statistics&) {
			auto const& segment = ranges_[r];
			auto const range = segment.range;
			auto output = fields_.gravity.output(range, bank ^ 1);
			for (std::size_t i = 0; i < range.count; ++i) {
				auto const& local = levels_.back().locals[segment.leaf + i];
				State field;
				field.potential() = constants::G * (gram / cm) * local[0];
				field.acceleration(0) = -constants::G * (gram / cm) * diagonal::derivative(local, 1, 0, 0) / cellWidth_;
				field.acceleration(1) = -constants::G * (gram / cm) * diagonal::derivative(local, 0, 1, 0) / cellWidth_;
				field.acceleration(2) = -constants::G * (gram / cm) * diagonal::derivative(local, 0, 0, 1) / cellWidth_;
				if (!finite(field)) throw std::runtime_error("Nonfinite gravity solution");
				output.put(i, field);
			}
			fields_.gravity.commit(range, bank ^ 1, output);
			if constexpr (build::hydro)
				copy(fields_.hydro, range, bank);
			else {
				auto input = fields_.density.read(range, bank).get();
				auto density = fields_.density.output(range, bank ^ 1);
				std::copy_n(input.data(), range.count, density.data());
				fields_.density.commit(range, bank ^ 1, density);
			}
			if constexpr (build::radiation) copy(fields_.radiation, range, bank);
		});
		result.localityCells = {levels_.back().moments.size()};
		return result;
	}
};

}	 // namespace octotigerII::gravity

#ifdef OCTOTIGERII_WITH_HPX
using FmmComponent = hpx::components::component<octotigerII::gravity::FmmPartition>;
HPX_REGISTER_COMPONENT(FmmComponent, octotigerII_fmm)
HPX_REGISTER_ACTION(octotigerII::gravity::FmmPartition::ExecuteAction, octotigerII_fmm_execute)
HPX_REGISTER_ACTION(octotigerII::gravity::FmmPartition::ReadAction, octotigerII_fmm_read)
#endif

namespace octotigerII::gravity {

class FieldSolver::Impl {
public:
	unsigned depths = 0;
#ifdef OCTOTIGERII_WITH_HPX
	std::vector<hpx::id_type> partitions;
#else
	std::unique_ptr<FmmPartition> partition;
#endif

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
	profiling::Elapsed profile("gravity.solve.wall_ns");
	auto result = impl_->phase(Stage::Initialize, 0, bank);
	for (unsigned depth = impl_->depths - 1; depth > 0; --depth)
		accumulate(result, impl_->phase(Stage::Upward, depth - 1, bank));
	for (unsigned depth = 0; depth < impl_->depths; ++depth)
		accumulate(result, impl_->phase(Stage::Downward, depth, bank));
	auto publication = impl_->phase(Stage::Publish, 0, bank);
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
