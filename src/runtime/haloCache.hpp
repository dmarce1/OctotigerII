#pragma once

#include "octotigerII/subgrid/view.hpp"

namespace octotigerII::runtime_detail {

// A coupled corrector only uses original states in the first ghost layer.
// The gravity-only rate path additionally needs raw prolongation donors to
// assemble its midpoint; retain the complete raw halo in that case.
template <typename State>
struct InitialHalo {
	std::vector<State> values;
	bool compact = false;

	void save(Subgrid const& block, HaloPlan const& plan, std::vector<State> const& input, bool firstLayer) {
		compact = firstLayer;
		if (!compact) { values.assign(input.begin(), input.end()); return; }
		values.clear();
		mesh::MeshLayout layout(block.layout.cellsPerActiveDimension(), 2);
		mesh::forEachCoordinate(mesh::filledCoordinates(1), mesh::filledCoordinates(layout.cellsPerActiveDimension() + 3),
			[&](auto const& cell) {
				if (!layout.isInterior(cell)) values.push_back(input.at(plan.ghostIndices.at(layout.index(cell))));
			});
	}

	void restore(Subgrid const& block, HaloPlan const& plan, std::vector<State>& output) const {
		if (!compact) { output.assign(values.begin(), values.end()); return; }
		output.assign(plan.valueCount ? plan.valueCount : plan.ghostCount, State{});
		mesh::MeshLayout layout(block.layout.cellsPerActiveDimension(), 2);
		std::size_t i = 0;
		mesh::forEachCoordinate(mesh::filledCoordinates(1), mesh::filledCoordinates(layout.cellsPerActiveDimension() + 3),
			[&](auto const& cell) {
				if (!layout.isInterior(cell)) output.at(plan.ghostIndices.at(layout.index(cell))) = values.at(i++);
			});
		if (i != values.size()) throw std::logic_error("Initial halo cache has the wrong shape");
	}
};

struct BlockHaloCache {
	InitialHalo<hydro::ConservedState> gas;
	InitialHalo<radiation::RadiationSystem::State> radiation;
	units::Time time{}, referenceStep{};
	std::uint64_t rateId = 0;
	unsigned bank = 0;
	std::vector<HaloTime> donors;

	bool matches(units::Time at, units::Time reference, std::uint64_t rate, unsigned inputBank,
		std::vector<HaloTime> const& times, bool unchangedBankCopy) const {
		if (time != at || referenceStep != reference || rateId != rate || (!unchangedBankCopy && bank != inputBank)
			|| donors.size() != times.size()) return false;
		for (std::size_t i = 0; i < times.size(); ++i)
			if (donors[i].bank != times[i].bank || donors[i].fraction != times[i].fraction || donors[i].nextBank != times[i].nextBank)
				return false;
		return true;
	}
};

} // namespace octotigerII::runtime_detail
