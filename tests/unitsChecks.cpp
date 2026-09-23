#include <boost/units/base_units/si/kilogram.hpp>
#include <boost/units/systems/si.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include "octotigerII/subgrid/subgrid.hpp"

using namespace octotigerII;
using namespace octotigerII::units;


namespace {
template <typename A, typename B>
concept Addable = requires(A a, B b) { a + b; };
static_assert(!std::is_convertible_v<Real, Density>);
static_assert(!std::is_convertible_v<Pressure, Density>);
static_assert(!Addable<Density, Pressure>);
#if OCTOTIGERII_HYDRO
static_assert(!Addable<hydro::ConservedState, hydro::ConservedFlux>);
#endif
static_assert(std::is_same_v<decltype(Density{} * Velocity{}), MomentumDensity>);
static_assert(std::is_same_v<decltype(Pressure{} * Velocity{}), EnergyFlux>);
#if OCTOTIGERII_HYDRO
static_assert(std::is_same_v<decltype(hydro::ConservedFlux{} * TimePerLength{}), hydro::ConservedState>);
#endif
static_assert(std::is_same_v<decltype(constants::G * Mass{} / (Length{} * Length{})), Acceleration>);
static_assert(std::is_same_v<decltype(constants::k_B * Temperature{}), Energy>);
static_assert(std::is_same_v<decltype(constants::a_r * boost::units::pow<4>(Temperature{})), EnergyDensity>);
#if OCTOTIGERII_RADIATION
static_assert(std::is_same_v<decltype(radiation::RadiationSystem::Flux{} * TimePerLength{}), radiation::RadiationSystem::State>);
static_assert(std::is_same_v<std::remove_cvref_t<decltype(radiation::RadiationSystem::State{}.radiativeFlux(0))>, EnergyFlux>);
#endif
static_assert(sizeof(Density) == sizeof(Real));
#if OCTOTIGERII_HYDRO
static_assert(sizeof(hydro::ConservedState) == (ndim + 2) * sizeof(Real));
#endif
#if OCTOTIGERII_RADIATION
static_assert(sizeof(radiation::RadiationSystem::State) == (ndim + 1) * sizeof(Real));
#endif

void require(bool ok, char const* message) {
	if (!ok) throw std::runtime_error(message);
}

template <typename Q>
void close(Q actual, Q expected, char const* message) {
	require(finite(actual) && finite(expected), "Nonfinite unit check");
	require(abs(actual - expected) <= 2e-14 * abs(expected), message);
}

void check() {
	using std::pow;

	namespace si = boost::units::si;
	close(Length(2.0 * si::meters), Length::from_value(200), "SI length to CGS");
	close(Mass(3.0 * si::kilograms), Mass::from_value(3000), "SI mass to CGS");
	close(Pressure(4.0 * si::pascals), Pressure::from_value(40), "SI pressure to CGS");
	close(Energy(5.0 * si::joules), Energy::from_value(5e7), "SI energy to CGS");
	close(constants::G, GravitationalConstant::from_value(6.67430e-8), "CGS G");
	close(constants::c, Velocity::from_value(2.99792458e10), "CGS c");
	close(constants::m_u, Mass::from_value(1.66053906892e-24), "CGS atomic mass unit");
	close(constants::k_B * Temperature::from_value(100), Energy::from_value(1.380649e-14), "k_B T");
	auto derived =
		(8 * pow(piR, 5) / 15) * boost::units::pow<4>(constants::k_B) / (boost::units::pow<3>(constants::planck) * boost::units::pow<3>(constants::c));
	close(constants::a_r, derived, "Radiation constant from k_B, h, c");
#if OCTOTIGERII_RADIATION
	auto e = EnergyDensity::from_value(12);
	std::array<EnergyFlux, ndim> f{};
	f[0] = e * constants::c;
	auto state = radiation::RadiationSystem::fromPhysical(e, f);
	auto restored = radiation::RadiationSystem::toPhysicalFlux(state);
	for (int axis = 0; axis < ndim; ++axis)
		require(state.radiativeFlux(axis) == f[axis], "Stored flux is not physical CGS flux");
	for (int axis = 0; axis < ndim; ++axis)
		require(restored[axis] == f[axis], "Physical radiation flux roundtrip");
	radiation::RadiationSystem reduced(0.25 * constants::c);
	close(reduced.physicalFlux(state, 0).energy(), f[0] * 0.25, "Reduced transport speed");
	radiation::RadiationSystem::State isotropic{};
	isotropic.energy() = e;
	auto flux = reduced.physicalFlux(isotropic, 0);
	close(flux.get<1>(), e * constants::c * reduced.reducedLightSpeed() / 3.0, "Isotropic radiation pressure");
	require(flux.get<0>() == EnergyFlux{}, "Isotropic energy flux");
	for (int axis = 1; axis < ndim; ++axis)
		require(flux.radiativeFlux(axis) == EnergyFluxTransport{}, "Isotropic flux symmetry");
	require(reduced.admissible({}), "Radiation vacuum");
	auto vacuumFlux = reduced.riemann({}, {}, 0);
	vacuumFlux.forEach([&](auto, auto q) { require(q == decltype(q){}, "Vacuum flux"); });
#endif
#if OCTOTIGERII_HYDRO
	hydro::HydroSystem gas(1.4);
	hydro::PrimitiveState primitive{};
	primitive.density() = Density::from_value(2);
	primitive.velocity(0) = Velocity::from_value(3);
	primitive.pressure() = Pressure::from_value(5);
	auto conserved = gas.conservedState(primitive);
	close(conserved.momentum(0), MomentumDensity::from_value(6), "Momentum dimensions");
	close(conserved.totalEnergy(), EnergyDensity::from_value(21.5), "Thermal plus kinetic energy");
	auto gasFlux = gas.physicalFlux(conserved, 0);
	close(gasFlux.mass(), MassFlux::from_value(6), "Mass flux");
	close(gasFlux.momentum(0), Pressure::from_value(23), "Momentum flux");
	close(gasFlux.energy(), EnergyFlux::from_value(79.5), "Energy flux");
#endif
	mesh::MeshLayout layout(4);
	close(layout.cellMeasure(Length::from_value(2)), Volume::from_value(pow(2, ndim)), "CGS unit transverse measure");
	static_assert(std::tuple_size_v<mesh::Coordinates> == ndim);
	require(layout.childCount() == (1 << ndim), "Child count must match compile-time dimension");
	std::size_t expectedIndex = 0;
	layout.forEachInterior([&](mesh::Coordinates const& cell, std::size_t index) {
		require(index == expectedIndex++, "ndim Cartesian traversal order");
		require(layout.index(cell) == index, "ndim Cartesian flattening");
	});
	require(expectedIndex == static_cast<std::size_t>(pow(4, ndim)), "Unexpected inactive-axis cells");
	bool rejected = false;
	try {
		layout.extent(ndim);
	} catch (std::out_of_range const&) {
		rejected = true;
	}
	require(rejected, "Inactive axis must not exist");
	mesh::BlockLocation parent{1, mesh::filledCoordinates(1)};
	for (int slot = 0; slot < (1 << ndim); ++slot) {
		auto child = parent.child(slot);
		require(child.parent() == parent && child.childSlot() == slot, "ndim child topology");
	}
}

}	 // namespace


int main() {
	try {
		check();
		std::cout << "CGS types, conversions, constants and selected physics passed\n";
		return 0;
	} catch (std::exception const& e) {
		std::cerr << e.what() << '\n';
		return 1;
	}
}
