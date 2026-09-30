// Stationary Gresho-Chan vortex: Springel (2010), section 8.5, eqs. 121-122.
#include <cmath>
#include "octotigerII/problems.hpp"
#include "octotigerII/verification/analytic.hpp"

namespace octotigerII::gresho {
namespace {
hydro::PrimitiveState state(mesh::PhysicalCoordinates const& point, Config const& c) {
    auto const length = c.mesh.upper - c.mesh.lower;
    auto const center = (c.mesh.upper + c.mesh.lower) / 2.0;
    Real const x = (point[0] - center) / length, y = (point[1] - center) / length;
    Real const r = std::hypot(x, y);
    Real const v = r < 0.2 ? 5 * r : r < 0.4 ? 2 - 5 * r : 0;
    Real const p = r < 0.2 ? 5 + 12.5 * r * r : r < 0.4
        ? 9 + 12.5 * r * r - 20 * r + 4 * std::log(r / 0.2) : 3 + 4 * std::log(2.0);
    hydro::PrimitiveState q{};
    q.density() = units::Density::from_value(1);
    q.pressure() = units::Pressure::from_value(p);
    if (r > 0) {
        q.velocity(0) = units::Velocity::from_value(-v * y / r);
        q.velocity(1) = units::Velocity::from_value(v * x / r);
    }
    return q;
}
}
ProblemBoundary problemBoundary(Config const&) { return {}; }
verification::Reference problemReference(Config const& c) {
    verification::Reference ref;
    ref.name = "Stationary Gresho-Chan vortex";
    ref.evaluate = [c](auto const& point, units::Time) {
        verification::ExactState result;
        result.hydro = state(point, c);
        return result;
    };
    return ref;
}
void problemDefaults(Config& c) {
    c.mesh.lower = units::Length::from_value(-0.5);
    c.mesh.upper = units::Length::from_value(0.5);
    c.mesh.boundary = finiteVolume::BoundaryConditions::periodic();
    c.hydro.gamma = 5.0 / 3.0;
    c.runtime.stopTime = units::Time::from_value(3);
    c.amr.shadowTolerance = 0;
}
void validateProblem(Config const& c) {
    if (c.hydro.eos != "ideal" || c.frame.omega != units::InverseTime{} || c.radiation.enabled)
        throw std::invalid_argument("gresho requires ideal gas, an inertial frame, and radiation off");
    for (auto a : c.hydro.acceleration)
        if (a != units::Acceleration{}) throw std::invalid_argument("gresho requires zero external acceleration");
    if (!c.mesh.boundary.all(finiteVolume::BoundaryCondition::Periodic))
        throw std::invalid_argument("gresho requires periodic boundaries");
}
void initializeProblem(Snapshot& data, Config const& c, bool) {
    hydro::HydroSystem gas(c);
    data.layout.forEachInterior([&](auto const& cell, std::size_t i) {
        data.hydro.values()[i] = gas.conservedState(state(data.layout.cellCenter(data.lower, data.cellWidth, cell), c));
    });
}
}
