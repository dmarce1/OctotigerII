// Material-advection fixture for the Timmes Helmholtz closure.
#include <cmath>
#include "octotigerII/problems.hpp"
#include "octotigerII/verification/analytic.hpp"
namespace octotigerII::helmholtz_advection {
ProblemBoundary problemBoundary(Config const&) {
    return {};
}
verification::Reference problemReference(Config const&) {
    return {};
}
void problemDefaults(Config& c) {
    c.hydro.eos = "helmholtz";
    c.massFractions.enabled = true;
    c.massFractions.species = composition::parseSpecies("helium:0.5:A=4,Z=2;iron:0.5:A=56,Z=26");
    c.mesh.cells = 8;
    c.mesh.level = 0;
    c.mesh.lower = units::Length::from_value(-5e7);
    c.mesh.upper = units::Length::from_value(5e7);
    c.mesh.boundary = physics::BoundaryConditions::periodic();
    c.runtime.stopTime = units::Time::from_value(0.01);
    c.amr.shadowTolerance = 0;
}
void validateProblem(Config const& c) {
    if (c.hydro.eos != "helmholtz" || !c.massFractions.enabled || c.massFractions.species.size() != 2 || c.massFractions.species[0].tracer() ||
        c.massFractions.species[1].tracer())
        throw std::invalid_argument("helmholtz-advection requires Helmholtz and exactly two material species");
}
void initializeProblem(Snapshot& data, Config const& c, bool) {
    hydro::HydroSystem gas(c);
    for (int s = 0; s < 2; ++s) {
        data.species.emplace_back(data.layout, data.cellWidth, data.lower);
    }
    data.layout.forEachInterior([&](auto const& cell, std::size_t i) {
        auto const x = data.layout.cellCenter(data.lower, data.cellWidth, cell);
        Real const phase = 2 * piR * Real((x[0] - c.mesh.lower) / (c.mesh.upper - c.mesh.lower));
        auto const rho = units::Density::from_value(1e3);
        Real const fraction = 0.5 + 0.2 * std::sin(phase);
        std::vector<units::Density> species{fraction * rho, (1 - fraction) * rho};
        hydro::ConservedState composition;
        gas.setComposition(composition, species, c.massFractions);
        hydro::PrimitiveState primitive;
        primitive.density() = rho;
        primitive.pressure() = units::Pressure::from_value(1e18);
        primitive.velocity(0) = units::Velocity::from_value(1e6);
        primitive.nuclei() = composition.nuclei() / rho;
        primitive.electrons() = composition.electrons() / rho;
        data.hydro.values()[i] = gas.conservedState(primitive);
        for (int s = 0; s < 2; ++s) {
            data.species[s].values()[i] = species[s];
        }
    });
}
}    // namespace octotigerII::helmholtz_advection
