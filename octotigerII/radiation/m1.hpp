/** @file
 * @brief M1 closure and HLL solver in temporary (E, F/c) calculation variables.
 * @ingroup numerics
 */
// Distributed under the Boost Software License, Version 1.0.
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include "octotigerII/math/FpeGuard.hpp"
#include "octotigerII/units/state.hpp"


namespace octotigerII::radiation {


/// M1 calculation kernel using energy-density variables (E, F/c).
/// These temporary quantities are distinct from the physical stored radiation state.
/// Closure: @ref ref_levermore1984 "Levermore (1984)"; wave speeds and transport:
/// @ref ref_skinner2013 "Skinner and Ostriker (2013)".
/// @ingroup numerics
class M1 {
public:

	using State = units::ScalarVectorState<units::EnergyDensity, units::EnergyDensity>;
	using Flux = units::ScalarVectorState<units::EnergyFlux, units::EnergyFlux>;
	static constexpr Real roundoff = 4096 * epsilonR;
	/// Normal flux and minimum/maximum characteristic velocities of the M1 system.
	/// @ingroup numerics
	class Waves {
	public:

		Flux flux;
		units::Velocity minus{}, plus{};
	};


	/// Return |F|/c from the temporary homogeneous energy-density state.
	static units::EnergyDensity magnitude(State const& u) {
		units::EnergyDensity result{};
		for (int axis = 0; axis < ndim; ++axis)
			result = units::hypot(result, u[axis + 1]);
		return result;
	}

	/// Check the closed realizable cone E≥0, |F|/c≤E, allowing scaled roundoff.
	static bool admissible(State const& u) {
		if (!finite(u) || u[0] < units::EnergyDensity{}) return false;
		auto norm = magnitude(u);
		return norm <= u[0] || (u[0] > units::EnergyDensity{} && norm - u[0] <= roundoff * u[0]);
	}

	/// Throw if the calculation state violates finiteness or the realizable cone.
	static void checkState(State const& u) {
		if (!admissible(u)) throw std::runtime_error("M1 requires finite E>=0 and |F|<=cE");
	}

	/// Validate and project only tolerated cone overshoot back to its boundary.
	static State canonical(State u) {
		checkState(u);
		auto norm = magnitude(u);
		if (norm > u[0])
			for (int i = 1; i < State::size(); ++i)
				u[i] *= Real(u[0] / norm);
		return u;
	}

	/// Return the closure flux and normal wave speeds for (E, F/c).
	/// The configured transport speed scales all characteristic velocities.
	static Waves physicalFlux(State const& u, int normal, units::Velocity chat) {
		using std::sqrt;

		FpeGuard guard;
		checkState(u);
		if (normal < 0 || normal >= ndim || (!(chat > units::Velocity{}) || !units::finite(chat))) throw std::invalid_argument("Invalid M1 normal or speed");
		std::array<Real, ndim> f{};
		if (u[0] > units::EnergyDensity{})
			for (int i = 0; i < ndim; ++i)
				f[i] = u[i + 1] / u[0];
		Real f2 = 0;
		for (auto component : f)
			f2 += component * component;
		if (f2 > 1) {
			for (auto& x : f)
				x /= sqrt(f2);
			f2 = 1;
		}
		Real fn = f[normal], mu2 = f2 > 0 ? std::min(Real(1), fn * fn / f2) : 0;
		Real s = sqrt(4 - 3 * f2);
		auto isotropic = u[0] * ((1 - f2) / (s + 1));
		auto directed = u[0] * (3 / (s + 2));
		Real radical = sqrt(std::max(Real(0), 2 * (1 - f2) * (s + mu2 * (s - 2)) / (s + 1)));
		Waves out;
		out.minus = chat * ((fn - radical) / s);
		out.plus = chat * ((fn + radical) / s);
		out.flux[0] = chat * u[0] * fn;
		for (int i = 0; i < ndim; ++i)
			out.flux[i + 1] = chat * (directed * fn * f[i] + (normal == i ? isotropic : units::Pressure{}));
		return out;
	}

	/// Two-wave HLL flux in homogeneous calculation variables;
	/// see @ref ref_harten1983 "Harten et al. (1983)".
	static Flux hll(State const& left, State const& right, int normal, units::Velocity chat) {
		FpeGuard guard;
		auto ul = canonical(left), ur = canonical(right);
		auto l = physicalFlux(ul, normal, chat), r = physicalFlux(ur, normal, chat);
		auto sm = std::min(l.minus, r.minus), sp = std::max(l.plus, r.plus);
		if (sm >= units::Velocity{}) return l.flux;
		if (sp <= units::Velocity{}) return r.flux;
		Real wr = sp / (sp - sm), wl = -sm / (sp - sm);
		return wr * (l.flux - sm * ul) + wl * (r.flux - sp * ur);
	}

	/// Use the component update scale to repair tiny negative energy or cone overshoot.
	static State roundoffState(State u, State const& updateScale) {
		auto scale = updateScale[0];
		for (int i = 1; i < State::size(); ++i)
			scale = std::max(scale, updateScale[i]);
		auto tolerance = roundoff * scale, norm = magnitude(u);
		if (u[0] < units::EnergyDensity{} && -u[0] <= tolerance && norm <= tolerance) return {};
		if (u[0] >= units::EnergyDensity{} && norm > u[0] && norm - u[0] <= tolerance)
			for (int i = 1; i < State::size(); ++i)
				u[i] *= Real(u[0] / norm);
		checkState(u);
		return u;
	}
};
}	 // namespace octotigerII::radiation
