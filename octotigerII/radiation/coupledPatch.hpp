/** @file
 * @brief Source-aware midpoint transport and conservative gray exchange on a patch.
 */
#pragma once
#include "octotigerII/radiation/diffusionFlux.hpp"
#include "octotigerII/radiation/matterCoupling.hpp"
#include "octotigerII/problems.hpp"
#include <functional>

namespace octotigerII::radiation {

struct CoupledPatchWorkspace {
	hydro::Solver::Workspace hydro;
	radiation::Solver::Workspace radiation;
};

/// Advance a coupled patch whose caller has supplied four valid ghost layers.
/// A forced half source/transport solve provides midpoint values throughout the
/// target's reconstruction stencil. The final accepted flux divergences drive
/// the full implicit source solve; pressure transport and stiff relaxation are
/// integrated together, rather than as successive independent operations.
/// write(interiorCell,gasState,radiationState) receives the coupled result.
/// Accepted face fluxes remain in workspace for reflux, species, and ledgers.
/// midpointBoundary can reapply physical exterior conditions after prediction;
/// the caller's periodic/internal ghost data already evolve consistently.
template <typename Writer>
void advanceCoupledPatch(hydro::Fields const& gas, radiation::Fields const& rad,
	hydro::HydroSystem const& gasSystem, RadiationSystem const& radSystem, Opacity opacity,
	units::Time dt, CoupledPatchWorkspace& workspace, Writer&& write,
	finiteVolume::RotatingFrame const& frame = finiteVolume::RotatingFrame{}, units::Time time = {},
	std::function<void(hydro::Fields&, radiation::Fields&, units::Time)> const& midpointBoundary = {},
	ProblemRadiationMaterial const& prescribedMaterial = {}, OpacityLaw const& configuredLaw = {}) {
	auto const& layout = gas.layout();
	if (layout.ghostWidth() != 4 || rad.layout().ghostWidth() != 4 ||
		layout.extents() != rad.layout().extents() || gas.cellWidth() != rad.cellWidth() || gas.lower() != rad.lower())
		throw std::invalid_argument("Coupled patch transport requires matching four-layer ghost patches");
	if (!(dt > units::Time{}) || !units::finite(dt)) throw std::invalid_argument("Coupled patch timestep must be positive and finite");
	Real const ratio = radSystem.reducedLightSpeed() / constants::c;
	auto extendedLower = gas.lower();
	for (auto& x : extendedLower) x -= 2.0 * gas.cellWidth();
	mesh::MeshLayout extendedLayout(layout.cellsPerActiveDimension() + 4, 2);
	hydro::Fields extendedGas(extendedLayout, gas.cellWidth(), extendedLower);
	radiation::Fields extendedRad(extendedLayout, rad.cellWidth(), extendedLower);
	// Both layouts have the same storage extent and physical storage positions.
	extendedGas.values() = gas.values();
	extendedRad.values() = rad.values();
	hydro::Fields midpointGas = gas;
	radiation::Fields midpointRad = rad;
	auto materialAt = [&](auto const& patch, mesh::Coordinates const& storageCell, units::Time at) {
		auto position = patch.lower();
		for (int d = 0; d < ndim; ++d)
			position[d] += (storageCell[d] - patch.layout().ghostWidth() + Real(0.5)) * patch.cellWidth();
		return prescribedMaterial ? checkedRadiationMaterial(prescribedMaterial, position, at) : RadiationMaterial{opacity, {}};
	};
	auto lawAt = [&](RadiationMaterial const& material) {
		auto law = configuredLaw;
		law.constantAbsorption = material.opacity;
		return law;
	};
	auto correction = [&](hydro::Fields const& material, units::Time at) {
		return [&, at, materialPtr = &material](RadiationSystem::Flux const& flux,
			RadiationSystem::State const& centerLeft, RadiationSystem::State const& centerRight,
			RadiationSystem::State const& faceLeft, RadiationSystem::State const& faceRight,
			mesh::Coordinates const& left, mesh::Coordinates const& right, int axis, units::Length width, units::Velocity speed) {
			auto const l = frame.toGridState(materialPtr->atStorage(left), at);
			auto const r = frame.toGridState(materialPtr->atStorage(right), at);
			MaterialVelocity lv{}, rv{};
			for (int d = 0; d < ndim; ++d) { lv[d] = l.momentum(d) / l.density(); rv[d] = r.momentum(d) / r.density(); }
			return diffusionCorrectedFlux(radSystem, flux, centerLeft, centerRight, faceLeft, faceRight,
				lawAt(materialAt(*materialPtr, left, at)).evaluate(l, gasSystem).fluxExtinction * l.density(),
				lawAt(materialAt(*materialPtr, right, at)).evaluate(r, gasSystem).fluxExtinction * r.density(), lv, rv, axis, width, speed);
		};
	};
	hydro::Solver::Workspace initialGasWork;
	radiation::Solver::Workspace initialRadWork;
	auto unchanged = [](auto state, mesh::Coordinates const&) { return state; };
	auto discard = [](auto const&, auto const&) {};
	// Fluxes are evaluated at t0, but their predictor increment must obey the
	// realizability limiter at its actual half-step interval. An unlimited
	// zero-dt AP pressure flux can leave a neighboring transparent cell's cone.
	hydro::Solver(gasSystem, frame, time).advanceInto(extendedGas, {}, initialGasWork, discard,
		unchanged, false, [](auto const& flux, auto const&...) { return flux; }, dt / 2.0);
	radiation::Solver(radSystem, frame, time).advanceInto(extendedRad, {}, initialRadWork, discard,
		unchanged, false, correction(extendedGas, time), dt / 2.0);
	extendedLayout.forEachInterior([&](auto const& cell, std::size_t i) {
		hydro::ConservedFlux gasDivergence{};
		RadiationSystem::Flux radDivergence{};
		for (int axis = 0; axis < ndim; ++axis) {
			auto upper = cell; ++upper[axis];
			gasDivergence += initialGasWork.fluxes[axis][extendedLayout.faceIndex(axis, cell)]
				- initialGasWork.fluxes[axis][extendedLayout.faceIndex(axis, upper)];
			radDivergence += initialRadWork.fluxes[axis][extendedLayout.faceIndex(axis, cell)]
				- initialRadWork.fluxes[axis][extendedLayout.faceIndex(axis, upper)];
		}
		auto g = extendedGas.atInterior(cell);
		auto r = extendedRad.atInterior(cell);
		auto const material = materialAt(extendedGas, extendedLayout.storageCoordinates(cell), time + dt / 4.0);
		auto radiationDrive = RadiationSystem::integratedFlux(radDivergence, dt / (2.0 * rad.cellWidth()));
		radiationDrive.energy() += ratio * material.photonPower * (dt / 2.0);
		coupleForcedWithOpacityLaw(g, r, hydro::HydroSystem::integratedFlux(gasDivergence, dt / (2.0 * gas.cellWidth())),
			radiationDrive, gasSystem, lawAt(material), ratio, dt / 2.0);
		midpointGas.values()[i] = g;
		midpointRad.values()[i] = r;
	});
	if (midpointBoundary) midpointBoundary(midpointGas, midpointRad, time + dt / 2.0);
	auto transportedGas = gas;
	auto transportedRad = rad;
	hydro::Solver(gasSystem, frame, time).advanceInto(gas, dt, workspace.hydro,
		[&](auto const& cell, auto const& value) { transportedGas.atInterior(cell) = value; },
		[&](auto const&, auto const& cell) { return midpointGas.atStorage(cell); }, false);
	radiation::Solver(radSystem, frame, time).advanceInto(rad, dt, workspace.radiation,
		[&](auto const& cell, auto const& value) { transportedRad.atInterior(cell) = value; },
		[&](auto const&, auto const& cell) { return midpointRad.atStorage(cell); }, false, correction(midpointGas, time + dt / 2.0));
	layout.forEachInterior([&](auto const& cell, std::size_t) {
		auto g = gas.atInterior(cell);
		auto r = rad.atInterior(cell);
		auto const dg = hydro::ConservedState(transportedGas.atInterior(cell) - g);
		auto dr = RadiationSystem::State(transportedRad.atInterior(cell) - r);
		auto const material = materialAt(gas, layout.storageCoordinates(cell), time + dt / 2.0);
		dr.energy() += ratio * material.photonPower * dt;
		coupleForcedWithOpacityLaw(g, r, dg, dr, gasSystem, lawAt(material), ratio, dt);
		write(cell, g, r);
	});
}

} // namespace octotigerII::radiation
