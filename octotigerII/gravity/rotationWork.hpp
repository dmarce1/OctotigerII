/** @file
 * @brief Conservative local gravitational work from rigid grid rotation.
 */
#pragma once
#include "octotigerII/gravity/fieldSolver.hpp"
#include <algorithm>
#include <stdexcept>

namespace octotigerII::gravity {

/// Scratch fields for the two signed, coordinate-weighted source solves.
struct RotationWorkWorkspace {
	storage::Field<units::Density> density;
	storage::ColumnFields<State> gravity;
	RotationWorkWorkspace(storage::Layout const& cells, storage::PartitionSet const& store)
	  : density(cells, store, 2, "rotationWork.weightedDensity"),
		gravity(cells, store, "rotationWork.weightedGravity", false, 2) {}
};

/// Compute psi_i = [x_i g_y-y_i g_x+L(g_y[x rho/L]-g_x[y rho/L])]/2.
/// Multiplication by Omega*rho gives the inertial energy source from grid motion.
/// Equal/opposite pair forces make its volume integral zero. For exact central
/// forces it is precisely rho*(Omega cross x).g; for an approximate mutual FMM
/// it differs only by a local force-approximation error, never a global correction.
/// The supplied request and physical source density must describe the ordinary
/// acceleration in ordinary, including identical masks and time interpolation.
inline Statistics rotationWorkPotential(FieldSolver& solver, FieldSolveRequest const& original,
	std::vector<Subgrid> const& blocks, storage::FieldHandle<units::Density> const& physicalDensity,
	storage::ColumnHandle<State> const& ordinary, unsigned ordinaryBank,
	RotationWorkWorkspace& scratch, storage::FieldHandle<units::VelocitySquared> const& output,
	unsigned outputBank, units::Length length) {
	if (!(length > units::Length{}) || !units::finite(length))
		throw std::invalid_argument("Rotation work requires a positive finite normalization length");
	for (std::size_t b = 0; b < blocks.size(); ++b) {
		auto const& block = blocks[b];
		auto const selection = original.sources.empty() ? DensitySelection{original.sourceBank, 0, 1, 0} : original.sources.at(b);
		storage::Buffer<units::Density> first, second;
		if (selection.weight != 0) first = physicalDensity.read(block.interior, selection.bank).get();
		if (selection.secondWeight != 0) second = physicalDensity.read(block.interior, selection.secondBank).get();
		for (int d = 0; d < 2; ++d) {
			auto weighted = scratch.density.handle().output(block.interior, d);
			block.layout.forEachInterior([&](auto const& cell, std::size_t i) {
				auto const rho = (selection.weight != 0 ? selection.weight * first.data()[i] : units::Density{})
					+ (selection.secondWeight != 0 ? selection.secondWeight * second.data()[i] : units::Density{});
				auto const position = block.lower[d] + (Real(cell[d]) + 0.5) * block.cellWidth;
				weighted.data()[i] = (position / length) * rho;
			});
			scratch.density.handle().commit(block.interior, d, weighted);
		}
	}
	Statistics total;
	for (unsigned d = 0; d < 2; ++d) {
		auto request = original;
		request.density = scratch.density.handle();
		request.sources.clear(); request.sourceBank = d;
		request.output = scratch.gravity.handle(); request.outputBank = d;
		request.allowSignedDensity = true;
		auto const result = solver.solve(request);
		total.multipolePairs += result.multipolePairs; total.directPairs += result.directPairs;
		total.workerTasks += result.workerTasks; total.ewaldPairs += result.ewaldPairs; total.reflectedPairs += result.reflectedPairs;
		if (total.localityCells.size() < result.localityCells.size()) total.localityCells.resize(result.localityCells.size());
		for (std::size_t i = 0; i < result.localityCells.size(); ++i) total.localityCells[i] += result.localityCells[i];
	}
	for (std::size_t b = 0; b < blocks.size(); ++b) {
		auto const& block = blocks[b];
		auto out = output.output(block.interior, outputBank);
		if (!original.targets.empty() && !original.targets.at(b)) {
			std::fill_n(out.data(), block.interior.count, units::VelocitySquared{});
		} else {
			auto const g = ordinary.read(block.interior, ordinaryBank).get();
			auto const gx = scratch.gravity.handle().read(block.interior, 0).get();
			auto const gy = scratch.gravity.handle().read(block.interior, 1).get();
			block.layout.forEachInterior([&](auto const& cell, std::size_t i) {
				auto const x = block.lower[0] + (Real(cell[0]) + 0.5) * block.cellWidth;
				auto const y = block.lower[1] + (Real(cell[1]) + 0.5) * block.cellWidth;
				out.data()[i] = 0.5 * (x * g.at(i).acceleration(1) - y * g.at(i).acceleration(0)
					+ length * (gx.at(i).acceleration(1) - gy.at(i).acceleration(0)));
			});
		}
		output.commit(block.interior, outputBank, out);
	}
	return total;
}

} // namespace octotigerII::gravity
