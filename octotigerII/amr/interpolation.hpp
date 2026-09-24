// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#pragma once
#include "octotigerII/physics/finiteVolume.hpp"

namespace octotigerII::amr {

template <typename State>
State slope(State const& left, State const& center, State const& right) {
	State result{};
	center.forEach([&](auto i, auto value) {
		result.template get<i>() = physics::limitedSlope(value - left.template get<i>(), right.template get<i>() - value, physics::Limiter::Minmod);
	});
	return result;
}

/// One slope fraction for every child preserves the parent average and keeps
/// all children inside the convex admissible set of the conserved variables.
template <typename System>
typename System::State interpolate(
	typename System::State const& center, std::array<typename System::State, ndim> const& slopes, std::array<Real, ndim> const& offset, System const& system) {
	auto valid = [&](Real fraction) {
		for (int slot = 0; slot < (1 << ndim); ++slot) {
			auto candidate = center;
			for (int d = 0; d < ndim; ++d)
				candidate += (slot & (1 << d) ? 0.5 : -0.5) * fraction * slopes[d];
			if (!system.admissible(candidate)) return false;
		}
		return true;
	};
	Real fraction = 1;
	if (!valid(fraction)) {
		Real low = 0, high = 1;
		for (int i = 0; i < 48; ++i) {
			Real const middle = (low + high) / 2;
			if (valid(middle))
				low = middle;
			else
				high = middle;
		}
		fraction = low;
	}
	auto result = center;
	for (int d = 0; d < ndim; ++d)
		result += fraction * offset[d] * slopes[d];
	return result;
}

}	 // namespace octotigerII::amr
