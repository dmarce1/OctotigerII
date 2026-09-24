// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#include "octotigerII/refinement/criteria.hpp"
#include <cmath>

namespace octotigerII::refinement {

Real ShadowCriterion::operator()(CellView const& cell) const {
	if (options_.shadowTolerance == 0) return 0;
	Real error = 0;
	auto inspect = [&](auto const& field) {
		field.value.forEach([&](auto i, auto value) {
			auto const shadow = field.shadow.template get<i>();
			auto const denominator = std::max(units::abs(value), units::abs(shadow)) + options_.shadowFloor * field.scale.template get<i>();
			if (denominator > decltype(denominator){}) error = std::max(error, Real(units::abs(value - shadow) / denominator));
		});
	};
	if (cell.hasHydro && options_.hydro) inspect(cell.hydro);
	if (cell.hasRadiation && options_.radiation) inspect(cell.radiation);
	return error / options_.shadowTolerance;
}

Criteria makeCriteria(Config const& config) {
	return {MassCriterion(config.amr.maxCellMass), DensityCriterion(config.amr.refineDensity), ShadowCriterion(config.amr)};
}

Real score(CellView const& cell, Criteria const& criteria) {
	using std::isfinite;
	Real result = 0;
	for (auto const& criterion : criteria) {
		auto const value = criterion(cell);
		if (!isfinite(value) || value < 0) throw std::runtime_error("Refinement criterion returned an invalid score");
		result = std::max(result, value);
	}
	return result;
}

}	 // namespace octotigerII::refinement
