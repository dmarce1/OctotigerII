/** @file
 * Image geometry and independent Newtonian/Ewald interaction lists.
 * Distributed under the Boost Software License, Version 1.0.
 */
#pragma once
#include <cmath>
#include <functional>
#include "octotigerII/gravity/boundary.hpp"
#include "octotigerII/gravity/diagonal/fmm.hpp"

namespace octotigerII::gravity {

/// A reflected fundamental source cell; paired walls also generate a lattice.
class SourceImage {
public:
	unsigned mask = 0;
	diagonal::Offset translation{};	   // in original domain lengths
};

/// Periods are 0 (open), 1 (periodic), or 2 (paired reflecting walls).
/// One reflecting wall gives one additional source copy. Orthogonal copies compose.
class ImageGeometry {
public:
	ImageGeometry(physics::BoundaryConditions const& boundaries = {});

	bool active() const {
		return images_.size() > 1 || periodic();
	}
	bool periodic() const {
		return periods_ != diagonal::Offset{};
	}
	diagonal::Offset periods(int count) const;
	std::vector<SourceImage> const& images() const {
		return images_;
	}
	diagonal::Offset separation(diagonal::Offset target, diagonal::Offset source, int count, SourceImage const& image) const;
	/// Integer physical centers on a common lattice, without the +1/2 cell offset.
	diagonal::Offset centerSeparation(diagonal::Offset target, diagonal::Offset source, int count, SourceImage const& image) const;
	bool acceptable(diagonal::Offset r, int count, double theta, bool correction) const;
	double nextImageDistance(diagonal::Offset r, int count) const;

	/// Independent walks: nearest-image Newtonian and all remaining-image Ewald.
	/// Refuse expansions crossing a nearest-image branch until it is uniform.
	/// This prevents differently accepted levels from subtracting different images.
	template <typename Function>
	void interactions(unsigned depth, diagonal::Offset target, bool leaf, double theta, Function const& function) const {
		for (auto const& image : images_)
			for (int kind = 0; kind < (periodic() ? 2 : 1); ++kind) {
				bool const correction = kind != 0;
				auto walk = [&](auto&& self, unsigned d, diagonal::Offset source) -> void {
					diagonal::Offset a{};
					for (int axis = 0; axis < 3; ++axis)
						a[axis] = target[axis] >> (depth - d);
					int const count = 1 << d;
					auto const r = separation(a, source, count, image);
					bool const accepted = acceptable(r, count, theta, correction);
					if (d == depth) {
						if (!(accepted || leaf)) return;
						if (!correction && r == diagonal::Offset{}) return;	   // only the physical self
						function(source, r, image.mask, correction);
					} else if (!accepted) {
						for (int slot = 0; slot < 8; ++slot)
							self(self, d + 1, {2 * source[0] + (slot & 1), 2 * source[1] + ((slot >> 1) & 1), 2 * source[2] + ((slot >> 2) & 1)});
					}
				};
				walk(walk, 0, {});
			}
	}

private:
	diagonal::Offset periods_{};
	std::vector<SourceImage> images_{{}};
};

/// The extra source coefficient is the raw z²/2 moment lost by harmonic reduction.
/// The extra local coefficient is its normalized Laplacian. They retain the 3P
/// neutralizing background exactly, including for expansion order one.
diagonal::Coefficients imageShiftMultipole(diagonal::Coefficients const&, diagonal::Vector const&, double ratio, int order);
diagonal::Coefficients imageShiftLocal(diagonal::Coefficients const&, diagonal::Vector const&, double ratio, int order);
diagonal::Coefficients reflectMultipole(diagonal::Coefficients const&, unsigned mask, int order);

}	 // namespace octotigerII::gravity
