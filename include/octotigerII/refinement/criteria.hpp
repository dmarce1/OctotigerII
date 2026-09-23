// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#pragma once
#include <functional>
#include "octotigerII/config.hpp"
#include "octotigerII/hydro/hydroSystem.hpp"
#include "octotigerII/radiation/radiationTransport.hpp"

namespace octotigerII::refinement {

/// The same typed interface serves every field and every refinement criterion.
/// Gradients have units of field/length. Shadow is the independently evolved
/// coarse field reconstructed at this cell, at the same physical time.
template <typename State>
class FieldSample {
public:
	using Gradient = decltype(State{} / units::Length::from_value(1));
	State value{}, shadow{}, scale{};
	std::array<Gradient, ndim> gradient{};
};

class CellView {
public:
	mesh::PhysicalCoordinates center{};
	units::Length width{};
	units::Volume volume{};
	units::Time time{};
	int level = 0;
	units::Mass mass{};
	std::array<units::Velocity, ndim> signalSpeed{};
	std::array<units::Acceleration, ndim> acceleration{};
	FieldSample<hydro::ConservedState> hydro;
	FieldSample<radiation::RadiationSystem::State> radiation;
	bool hasHydro = false, hasRadiation = false;
};

/// Score > 1 requests refinement. A smaller coarsening threshold supplies
/// hysteresis. Criteria only inspect fields; topology and buffering live in AMR.
using Criterion = std::function<Real(CellView const&)>;
using Criteria = std::vector<Criterion>;

class MassCriterion {
public:
	explicit MassCriterion(units::Mass maximum)
	  : maximum_(maximum) {}
	Real operator()(CellView const& cell) const {
		return maximum_ > units::Mass{} ? Real(cell.mass / maximum_) : 0;
	}

private:
	units::Mass maximum_;
};

class ShadowCriterion {
public:
	explicit ShadowCriterion(Config::AmrOptions const& options)
	  : options_(options) {}
	Real operator()(CellView const& cell) const;

private:
	Config::AmrOptions options_;
};

Criteria makeCriteria(Config const&);
Real score(CellView const&, Criteria const&);

}	 // namespace octotigerII::refinement
