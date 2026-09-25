#include "testSupport.hpp"
#include "octotigerII/gravity/fieldSolver.hpp"
#include "octotigerII/gravity/ewald.hpp"
#include "octotigerII/gravity/images.hpp"
#include "octotigerII/subgrid/topology.hpp"
#include <iostream>
#include <numeric>
#ifdef OCTOTIGERII_WITH_HPX
#include <hpx/runtime_distributed/find_all_localities.hpp>
#endif

using namespace octotigerII;

namespace {
std::vector<storage::Locality> localities() {
#ifdef OCTOTIGERII_WITH_HPX
	return hpx::find_all_localities();
#else
	return {0};
#endif
}

Config config(bool adaptive, int boundary, int order = 3, Real openingAngle = 0.5) {
	auto c = test::parseConfig({"--mesh.cells=4", "--mesh.level=1", "--output.enabled=off"});
	c.amr.enabled = adaptive;
	c.gravity.multipoleOrder = order;
	c.gravity.openingAngle = openingAngle;
	if (boundary == 1) c.mesh.boundary = physics::BoundaryConditions::periodic();
	if (boundary == 2) c.mesh.boundary.lower[0] = physics::BoundaryCondition::Reflecting;
	return c;
}

std::vector<mesh::BlockLocation> leaves(bool adaptive) {
	std::vector<mesh::BlockLocation> result;
	for (int slot = 0; slot < 8; ++slot) {
		auto const child = mesh::BlockLocation{}.child(slot);
		if (adaptive && slot == 0)
			for (int sub = 0; sub < 8; ++sub) result.push_back(child.child(sub));
		else result.push_back(child);
	}
	return result;
}

class Fixture {
public:
	Fixture(bool adaptive, int boundary, int order = 3, Real openingAngle = 0.5)
	  : c(config(adaptive, boundary, order, openingAngle)), owners(localities()), topology(c, owners.size(), leaves(adaptive)),
		repository(c, topology.storageLayout(), owners), store(owners),
		density(topology.storageLayout(), store, 2, "partial.density"),
		output(topology.storageLayout(), store, "partial.output", false, 4),
		solver(c, topology.blocks(), repository.directory(), owners) {
		std::size_t index = 0;
		for (auto const& block : topology.blocks()) {
			for (unsigned bank = 0; bank < 2; ++bank) {
				auto rho = density.handle().output(block.interior, bank);
				auto source = repository.directory().density.output(block.interior, bank);
				auto field = repository.directory().gravity.output(block.interior, bank);
				block.layout.forEachInterior([&](auto const& cell, std::size_t i) {
					Real const base = 2 + Real(((index + i) * 17 + 31) % 127) / 53;
					rho.data()[i] = units::Density::from_value(base + (bank ? ((cell[0] & 1) ? Real(0.125) : Real(-0.125)) : 0));
					source.data()[i] = units::Density::from_value(base);
					field.put(i, {});
				});
				density.handle().commit(block.interior, bank, rho);
				repository.directory().density.commit(block.interior, bank, source);
				repository.directory().gravity.commit(block.interior, bank, field);
			}
			index += block.interior.count;
		}
	}

	gravity::FieldSolveRequest request(unsigned bank) const {
		gravity::FieldSolveRequest result;
		result.density = density.handle();
		result.output = output.handle();
		result.outputBank = bank;
		return result;
	}

	std::vector<gravity::State> read(storage::ColumnHandle<gravity::State> const& handle, unsigned bank) const {
		std::vector<gravity::State> result;
		for (auto const& block : topology.blocks()) {
			auto input = handle.read(block.interior, bank).get();
			for (std::size_t i = 0; i < block.interior.count; ++i) result.push_back(input.at(i));
		}
		return result;
	}

	Config c;
	std::vector<storage::Locality> owners;
	CartesianTopology topology;
	FieldRepository repository;
	storage::PartitionSet store;
	storage::Field<units::Density> density;
	storage::ColumnFields<gravity::State> output;
	gravity::FieldSolver solver;
};

void compare(std::vector<gravity::State> const& actual, std::vector<gravity::State> const& expected, Real tolerance = 1e-12,
	std::vector<gravity::State> const& referenceScale = {}) {
	ASSERT_EQ(actual.size(), expected.size());
	gravity::State{}.forEach([&](auto f, auto) {
		Real scale = 0, error = 0;
		bool valid = true;
		for (auto const& state : expected) scale = std::max(scale, std::abs(units::value(state.template get<f>())));
		for (auto const& state : referenceScale) scale = std::max(scale, std::abs(units::value(state.template get<f>())));
		for (std::size_t i = 0; i < actual.size(); ++i) {
			valid = valid && units::finite(actual[i].template get<f>()) && units::finite(expected[i].template get<f>());
			error = std::max(error, std::abs(units::value(actual[i].template get<f>() - expected[i].template get<f>())));
		}
		EXPECT_TRUE(valid);
		EXPECT_LE(error, tolerance * scale) << "field=" << int(f) << " normalized error=" << (scale ? error / scale : error);
	});
}
} // namespace

class PartialGravity : public ::testing::TestWithParam<std::tuple<bool, int>> {};

TEST_P(PartialGravity, SourceSplitIsLinearReciprocalAndIndependentOfPublication) {
	auto const [adaptive, boundary] = GetParam();
	Fixture f(adaptive, boundary);
	f.solver.solve(0);
	auto const baseline = f.read(f.repository.directory().gravity, 1);
	f.solver.solve(f.request(0));
	compare(f.read(f.output.handle(), 0), baseline, 0);
	for (unsigned group = 0; group < 2; ++group) {
		auto request = f.request(group + 1);
		for (std::size_t b = 0; b < f.topology.blocks().size(); ++b)
			request.sources.push_back({0, 0, b % 2 == group ? Real(1) : Real(0), 0});
		f.solver.solve(request);
	}
	auto const a = f.read(f.output.handle(), 1), b = f.read(f.output.handle(), 2);
	auto sum = a;
	for (std::size_t i = 0; i < sum.size(); ++i) sum[i] += b[i];
	compare(sum, baseline);
	long double ab = 0, ba = 0, norm = 0;
	std::size_t index = 0;
	for (std::size_t block = 0; block < f.topology.blocks().size(); ++block) {
		auto const& grid = f.topology.blocks()[block];
		auto rho = f.density.handle().read(grid.interior, 0).get();
		for (std::size_t i = 0; i < grid.interior.count; ++i, ++index) {
			long double const mass = units::value(rho.data()[i] * grid.layout.cellMeasure(grid.cellWidth));
			long double const work = mass * units::value((block % 2 ? a[index] : b[index]).potential());
			(block % 2 ? ba : ab) += work;
			norm += std::abs(work);
		}
	}
	EXPECT_LT(std::abs(ab - ba) / norm, 3e-14L);
	// Independent output must not publish either the hydro/density bank or the full field.
	compare(f.read(f.repository.directory().gravity, 1), baseline, 0);
	compare(f.read(f.repository.directory().gravity, 0), std::vector<gravity::State>(baseline.size()), 0);
	for (auto const& grid : f.topology.blocks()) {
		auto rho = f.density.handle().read(grid.interior, 0).get();
		for (unsigned bank = 0; bank < 2; ++bank) {
			auto original = f.repository.directory().density.read(grid.interior, bank).get();
			for (std::size_t i = 0; i < grid.interior.count; ++i) EXPECT_EQ(rho.data()[i], original.data()[i]);
		}
	}
	// A target mask retains the already computed values on excluded blocks.
	auto target = f.request(3);
	f.solver.solve(target);
	target.targets.resize(f.topology.blocks().size());
	for (std::size_t block = 0; block < target.targets.size(); ++block) {
		target.targets[block] = block % 2 == 0;
		target.sources.push_back({0, 0, block % 2 ? Real(1) : Real(0), 0});
	}
	auto const stats = f.solver.solve(target);
	auto expected = baseline;
	index = 0;
	std::uint64_t targets = 0;
	for (std::size_t block = 0; block < f.topology.blocks().size(); ++block)
		for (std::size_t i = 0; i < f.topology.blocks()[block].interior.count; ++i, ++index)
			if (target.targets[block]) { expected[index] = b[index]; ++targets; }
	compare(f.read(f.output.handle(), 3), expected, 0);
	EXPECT_EQ(std::accumulate(stats.localityCells.begin(), stats.localityCells.end(), std::uint64_t(0)), targets);
	f.solver.solve(0);
	compare(f.read(f.repository.directory().gravity, 1), baseline, 0);
}

TEST_P(PartialGravity, WeightedBanksRetainSignedHigherMomentsAndRecoverAfterRejection) {
	auto const [adaptive, boundary] = GetParam();
	Fixture f(adaptive, boundary);
	f.solver.solve(f.request(0));
	auto second = f.request(1);
	second.sourceBank = 1;
	f.solver.solve(second);
	auto const a = f.read(f.output.handle(), 0), b = f.read(f.output.handle(), 1);
	auto delta = f.request(2);
	delta.sources.assign(f.topology.blocks().size(), {1, 0, 1, -1});
	EXPECT_ANY_THROW(f.solver.solve(delta));
	delta.allowSignedDensity = true;
	f.solver.solve(delta);
	auto expected = b;
	for (std::size_t i = 0; i < expected.size(); ++i) expected[i] -= a[i];
	// Normalize to the input fields: symmetry can make an exact component zero,
	// leaving only subtraction roundoff in its endpoint difference. Each source
	// pair has zero monopole and a nonzero dipole.
	compare(f.read(f.output.handle(), 2), expected, 1e-12, a);
	compare(f.read(f.output.handle(), 0), a, 0);
	compare(f.read(f.output.handle(), 1), b, 0);
	// Block-specific bank choices and interpolation coefficients.
	auto mixed = f.request(3);
	for (std::size_t block = 0; block < f.topology.blocks().size(); ++block)
		mixed.sources.push_back(block % 2 ? gravity::DensitySelection{1, 0, Real(0.75), Real(0.25)} :
			gravity::DensitySelection{0, 1, Real(0.25), Real(0.75)});
	f.solver.solve(mixed);
	for (std::size_t i = 0; i < expected.size(); ++i) expected[i] = Real(0.25) * a[i] + Real(0.75) * b[i];
	compare(f.read(f.output.handle(), 3), expected);
}

INSTANTIATE_TEST_SUITE_P(Geometries, PartialGravity,
	::testing::Values(std::tuple{false, 0}, std::tuple{true, 0}, std::tuple{false, 1}, std::tuple{true, 1}, std::tuple{false, 2}, std::tuple{true, 2}));

TEST(PartialGravityMomentum, ShellForceBalance) {
	for (bool adaptive : {false, true}) for (Real theta : {Real(0.01), Real(0.5)}) for (int boundary : {0, 1}) {
		if (boundary && theta < 0.5) continue;
		Fixture f(adaptive, boundary, 5, theta);
		f.solver.solve(f.request(0));
		auto slow = f.request(1);
		for (std::size_t block = 0; block < f.topology.blocks().size(); ++block)
			slow.sources.push_back({0, 0, block % 2 ? Real(0) : Real(1), 0});
		f.solver.solve(slow);
		auto fast = f.request(2);
		for (std::size_t block = 0; block < f.topology.blocks().size(); ++block)
			fast.sources.push_back({0, 0, block % 2 ? Real(1) : Real(0), 0});
		f.solver.solve(fast);
		auto const full = f.read(f.output.handle(), 0);
		auto const slowField = f.read(f.output.handle(), 1);
		auto const fastField = f.read(f.output.handle(), 2);
		std::array<std::array<long double, ndim>, 3> force{}, norm{};
		std::size_t index = 0;
		for (std::size_t block = 0; block < f.topology.blocks().size(); ++block) {
			auto const& grid = f.topology.blocks()[block];
			auto rho = f.density.handle().read(grid.interior, 0).get();
			for (std::size_t i = 0; i < grid.interior.count; ++i, ++index) {
				long double const mass = units::value(rho.data()[i] * grid.layout.cellMeasure(grid.cellWidth));
				for (int d = 0; d < ndim; ++d) {
					// Full, slow shell (slow-slow plus both cross directions),
					// and cross-only action/reaction use the same physical masses.
					std::array<long double, 3> const acceleration{
						units::value(full[index].acceleration(d)),
						units::value((block % 2 ? slowField[index] : full[index]).acceleration(d)),
						units::value((block % 2 ? slowField[index] : fastField[index]).acceleration(d))};
					for (int kind = 0; kind < 3; ++kind) {
						auto const value = mass * acceleration[kind];
						force[kind][d] += value;
						norm[kind][d] += std::abs(value);
					}
				}
			}
		}
		for (int kind = 0; kind < 3; ++kind) {
			long double residual = 0;
			for (int d = 0; d < ndim; ++d) residual = std::max(residual, std::abs(force[kind][d]) / norm[kind][d]);
			std::cout << "gravity momentum adaptive=" << adaptive << " theta=" << theta << " periodic=" << boundary
				<< " kind=" << kind << " residual=" << residual << '\n';
			EXPECT_LT(residual, 3e-14L) << "adaptive=" << adaptive << " theta=" << theta << " periodic=" << boundary << " kind=" << kind;
		}
	}
}

TEST(PartialGravityMomentum, AuxiliaryForceLocalsAreMutualThroughMaximumOrder) {
	using gravity::diagonal::Coefficients;
	using gravity::diagonal::Vector;
	std::array<std::array<Vector, 3>, 2> const points{{
		{{{-0.3, 0.1, -0.2}, {0.4, -0.25, 0.35}, {-0.2, -0.15, 0.4}}},
		{{{0.35, -0.25, -0.1}, {-0.2, 0.3, 0.25}, {0.1, -0.4, -0.35}}}}};
	std::array<std::array<Real, 3>, 2> const masses{{{1.2, 0.8, 1.7}, {0.7, 1.3, 0.9}}};
	for (int order : {1, 5, 10}) for (bool periodicCorrection : {false, true}) {
		std::array<Coefficients, 2> moments, locals;
		for (int side = 0; side < 2; ++side) {
			moments[side].assign(gravity::diagonal::coefficientCount(order) + 1, 0);
			locals[side].assign(gravity::diagonal::coefficientCount(order + 1) + 1, 0);
			for (int point = 0; point < 3; ++point) {
				auto const contribution = gravity::imageShiftMultipole(Coefficients{masses[side][point]}, points[side][point], 1, order);
				for (std::size_t i = 0; i < contribution.size(); ++i) moments[side][i] += contribution[i];
			}
		}
		for (int side = 0; side < 2; ++side) {
			gravity::diagonal::Offset r{4, 1, -2};
			if (side) for (auto& coordinate : r) coordinate = -coordinate;
			if (periodicCorrection)
				gravity::ewald::getOperator(order, r, {16, 16, 16}, true)->add(locals[side], moments[side ^ 1], 1);
			else
				gravity::diagonal::getOperator(order, r, true)->add(locals[side], moments[side ^ 1], 1);
		}
		std::array<long double, ndim> force{}, norm{};
		for (int side = 0; side < 2; ++side) for (int point = 0; point < 3; ++point) {
			auto const local = gravity::imageShiftLocal(locals[side], points[side][point], 1, order + 1);
			for (int d = 0; d < ndim; ++d) {
				gravity::diagonal::Offset e{};
				e[d] = 1;
				long double const value = masses[side][point] * gravity::diagonal::derivative(local, e[0], e[1], e[2]);
				force[d] += value;
				norm[d] += std::abs(value);
			}
		}
		for (int d = 0; d < ndim; ++d)
			EXPECT_LT(std::abs(force[d]) / norm[d], 3e-14L) << "order=" << order << " periodic=" << periodicCorrection << " axis=" << d;
	}
}
