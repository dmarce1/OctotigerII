#include "octotigerII/problems.hpp"
#include "octotigerII/problems/binaryScf.hpp"
#include "octotigerII/verification/analytic.hpp"
#include "octotigerII/subgrid/subgrid.hpp"

namespace octotigerII::binary_scf {
ProblemBoundary problemBoundary(Config const&) { return {}; }
verification::Reference problemReference(Config const&) {
	verification::Reference r;
	r.reason = "The binary SCF is a numerical equilibrium; use direct gravity and mesh convergence checks";
	return r;
}
void problemDefaults(Config& c) {
	c.mesh.cells = 8;
	c.mesh.level = 2;
	c.mesh.lower = -2.0 * c.scf.separation;
	c.mesh.upper = 2.0 * c.scf.separation;
	c.hydro.gamma = Real(5) / 3;
	c.runtime.stopTime = units::Time{};
	c.amr.shadowTolerance = 0;
	if constexpr (build::massFractions) {
		c.massFractions.enabled = true;
		c.massFractions.species = composition::parseSpecies("primary_core:0:He;primary_envelope:1:He;donor_core:0:He;donor_envelope:0:He;atmosphere:0:He");
	}
}
void validateProblem(Config const& c) {
	problems::BinaryScf::validate(c);
	if (c.massFractions.enabled) {
		std::array<std::string, 5> const names{"primary_core", "primary_envelope", "donor_core", "donor_envelope", "atmosphere"};
		if (c.massFractions.species.size() != names.size()) throw std::invalid_argument("binary-scf requires its five core/envelope/atmosphere material species");
		for (std::size_t s = 0; s < names.size(); ++s)
			if (c.massFractions.species[s].name != names[s] || c.massFractions.species[s].tracer())
				throw std::invalid_argument("binary-scf material species must be primary_core, primary_envelope, donor_core, donor_envelope, atmosphere in that order");
	}
}
void initializeProblem(Snapshot& data, Config const& c, bool) {
	auto const model = problems::BinaryScf::get(c);
	hydro::HydroSystem const gas(c);
	if (c.massFractions.enabled) for (int s = 0; s < 5; ++s) data.species.emplace_back(data.layout, data.cellWidth, data.lower);
	data.layout.forEachInterior([&](mesh::Coordinates const& cell, std::size_t i) {
		auto const x = data.layout.cellCenter(data.lower, data.cellWidth, cell);
		auto const value = model->average(x, data.cellWidth);
		auto primitive = model->state(value, x);
		if (c.massFractions.enabled) {
			units::Density rho{};
			for (int s = 0; s < 5; ++s) {
				auto const r = (s < 4 ? value.partial[s] : c.scf.atmosphereFraction) * model->densityUnit();
				data.species[s].values()[i] = r; rho += r;
			}
			primitive.density() = rho;
		}
		data.hydro.values()[i] = gas.conservedState(primitive);
	});
}
} // namespace octotigerII::binary_scf
