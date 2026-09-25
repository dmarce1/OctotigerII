/** @file
 * @brief Unsplit MUSCL-Hancock transport and slope limiting.
 * @ingroup numerics
 */
// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.
#pragma once

#include "octotigerII/mesh.hpp"
#include "octotigerII/physics/boundary.hpp"
#include "octotigerII/profiling.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>


namespace octotigerII::physics {

/// Return the implemented transport scheme name.
std::string_view finiteVolumeSchemeName();


enum class Limiter
{
	Minmod,
	VanLeer,
	MinmodTheta
};


/// Limit component differences without mixing their dimensions.
/// Supports minmod, van Leer harmonic mean, and generalized minmod with theta.
template <typename Q>
Q limitedSlope(Q leftDifference, Q rightDifference, Limiter limiter, Real theta = Real(1.5)) {
	if (!((leftDifference > Q{} && rightDifference > Q{}) || (leftDifference < Q{} && rightDifference < Q{}))) {
		return Q{};
	}
	Real const sign = (leftDifference > Q{} ? Real(1) : Real(-1));
	auto const left = boost::units::abs(leftDifference);
	auto const right = boost::units::abs(rightDifference);
	switch (limiter) {
	case Limiter::Minmod:
		return sign * std::min(left, right);
	case Limiter::VanLeer:
		return sign * (Real(2) * left * right / (left + right));
	case Limiter::MinmodTheta: {
		auto const centered = Real(0.5) * (left + right);
		return sign * std::min({theta * left, centered, theta * right});
	}
	}
	throw std::logic_error("Unknown finite-volume limiter");
}

/// Fill physical boundaries for an owning test patch, including corner cells.
/// The distributed runtime instead obtains its temporary halos from the mesh adapter.
template <typename System>
void fillGhostCells(mesh::PatchData<typename System::State>& patch, BoundaryConditions const& boundaries, System const& system,
	AnalyticBoundary<typename System::State> const& analytic = {}) {
	boundaries.validate();
	if (boundaries.contains(BoundaryCondition::Analytic) && !analytic)
		throw std::invalid_argument("Analytic boundary requires a problem evaluator");
	auto const& layout = patch.layout();
	mesh::forEachCoordinate(layout.extents(), [&](mesh::Coordinates const& destination) {
		if (layout.isInterior(destination)) return;
		auto cell = destination;
		for (int axis = 0; axis < ndim; ++axis)
			cell[axis] -= layout.ghostWidth();
		auto const mapped = boundaries.map(cell, layout.cellsPerActiveDimension());
		if (mapped.analytic) {
			mesh::PhysicalCoordinates position{};
			for (int axis = 0; axis < ndim; ++axis)
				position[axis] = patch.lower()[axis] + (Real(cell[axis]) + 0.5) * patch.cellWidth();
			patch.atStorage(destination) = evaluateBoundary(analytic, position, patch.timeState().time, system);
		} else {
			patch.atStorage(destination) =
				transformBoundary(patch.atInterior(mapped.source), mapped.reflectionMask, mapped.outflowLowerMask, mapped.outflowUpperMask, system);
		}
	});
}


// Dimensionally unsplit MUSCL-Hancock. Reconstruction is directional, the
// half-step predictor contains the divergence from every active direction,
// and the only intercell quantities are face-centered fluxes.
/// Dimensionally unsplit, second-order predictor/corrector for conservative systems.
/// Piecewise-linear reconstruction follows the MUSCL lineage of
/// @ref ref_vanleer1979 "van Leer (1979)". Face states are predicted by half of
/// the full multidimensional flux divergence, then a shared face flux updates cells.
/// The System adapter supplies states, fluxes, admissibility, and unit conversions.
/// @ingroup numerics
template <typename System>
class MusclHancock {
public:

	using State = typename System::State;
	using Reconstruction = typename System::Reconstruction;
	using Flux = typename System::Flux;
	using FaceFluxes = std::array<std::vector<Flux>, ndim>;

	/// Physical integration interval and the face fluxes produced by the step.
	/// @ingroup numerics
	class StepResult {
	public:

		mesh::TimeInterval timeInterval;
		FaceFluxes faceFluxes;
	};


	MusclHancock(System system, Limiter limiter = Limiter::VanLeer, Real theta = Real(1.5))
	  : system_(std::move(system))
	  , limiter_(limiter)
	  , theta_(theta) {
		if (!(theta_ >= 1 && theta_ <= 2)) {
			throw std::invalid_argument("minmod-theta parameter must lie in [1,2]");
		}
	}

	/// Return CFL divided by the sum of maximum directional speeds divided by cell width.
	/// Requires a positive finite speed and 0<CFL≤0.5.
	template <typename Patch>
	units::Time stableTimestep(Patch const& patch, Real courantNumber) const {
		if (!(courantNumber > 0 && courantNumber <= Real(0.5))) {
			throw std::invalid_argument("Invalid MUSCL-Hancock layout or Courant number");
		}
		std::array<units::Velocity, ndim> maximumSpeed{};
		patch.layout().forEachInterior([&](mesh::Coordinates const& cell, std::size_t) {
			State const state = patch.atInterior(cell);
			for (int axis = 0; axis < ndim; ++axis) {
				maximumSpeed[axis] = std::max(maximumSpeed[axis], system_.maximumSignalSpeed(state, axis));
			}
		});
		units::InverseTime inverseStep{};
		for (auto speed : maximumSpeed) {
			inverseStep += speed / patch.cellWidth();
		}
		if (!(inverseStep > units::InverseTime{}) || !units::finite(inverseStep)) {
			throw std::runtime_error("No finite positive signal speed in finite-volume patch");
		}
		return courantNumber / inverseStep;
	}

	/// Advance an owning test patch using the supplied physical boundary conditions.
	StepResult advance(mesh::PatchData<State>& patch, units::Time stepSize, BoundaryConditions const& boundaries,
		AnalyticBoundary<State> const& analytic = {}) const {
		return advanceWithBoundaryUpdater(
			patch, stepSize, [&](mesh::PatchData<State>& boundaryPatch, units::Time) { fillGhostCells(boundaryPatch, boundaries, system_, analytic); });
	}

	// AMR drivers can supply same-level, coarse/fine, or physical boundary
	// data at the requested physical time. The time-tagged interface avoids
	// baking a global-cycle assumption into the spatial integrator and can
	// later interpolate boundaries for level-dependent timesteps.
	/// Advance an owning patch with a caller-provided time-aware boundary fill.
	template <typename BoundaryUpdater>
	StepResult advanceWithBoundaryUpdater(mesh::PatchData<State>& patch, units::Time stepSize, BoundaryUpdater&& updateBoundaries) const {
		mesh::TimeInterval const interval{patch.timeState().time, patch.timeState().time + stepSize};
		updateBoundaries(patch, interval.begin);
		Workspace workspace;
		auto next = patch.values();
		advanceInto(patch, stepSize, workspace,
			[&](mesh::Coordinates const& cell, State value) {
				if constexpr (requires { system_.synchronize(value); }) system_.synchronize(value);
				next[patch.layout().index(patch.layout().storageCoordinates(cell))] = value;
			});
		patch.values().swap(next);
		patch.timeState().completeStep(stepSize);
		updateBoundaries(patch, interval.end);
		return StepResult{interval, std::move(workspace.fluxes)};
	}

	/// Reusable predictor and face-flux arrays for one concurrently executing worker.
	/// @ingroup numerics
	class Workspace {
	public:

		std::array<std::vector<State>, ndim> minus;
		std::array<std::vector<State>, ndim> plus;
		FaceFluxes fluxes;
	};


	// Input is immutable. The caller owns a disjoint output range and publishes
	// it only after every task in the stage succeeds. Workspace is reusable.
	/// Read an immutable patch and write each interior result through the supplied callback.
	/// Requires two ghost layers. Workspace is exclusive to this task; the caller
	/// publishes output only after all dependent work and transfers finish.
	template <typename Patch, typename Writer>
	void advanceInto(Patch const& patch, units::Time stepSize, Workspace& workspace, Writer&& write) const {
		mesh::MeshLayout const& layout = patch.layout();
		if (layout.ghostWidth() < 2) throw std::invalid_argument("MUSCL-Hancock needs two ghost cells");
		if (!(stepSize > units::Time{}) || !units::finite(stepSize)) throw std::invalid_argument("MUSCL-Hancock timestep must be positive and finite");
		auto& minus = workspace.minus;
		auto& plus = workspace.plus;
		for (int axis = 0; axis < ndim; ++axis) {
			minus[axis].resize(layout.cellCount());
			plus[axis].resize(layout.cellCount());
		}
		{
			profiling::Region profile("transport.reconstruct_predict");
			predictFaceStates(patch, stepSize, minus, plus);
		}

		auto& fluxes = workspace.fluxes;
		{
			profiling::Region profile("transport.fluxes");
			for (int axis = 0; axis < ndim; ++axis) {
				fluxes[axis].resize(layout.faceCount(axis));
				mesh::Coordinates const faceExtents = layout.faceExtents(axis);
				mesh::forEachCoordinate(faceExtents, [&](mesh::Coordinates const& face) {
					mesh::Coordinates left = face;
					for (int d = 0; d < ndim; ++d)
						left[d] = std::min(left[d], layout.interiorExtent(d) - 1);
					left = layout.storageCoordinates(left);
					mesh::Coordinates right = left;
					left[axis] = layout.ghostWidth() + face[axis] - 1;
					right[axis] = layout.ghostWidth() + face[axis];
					State const& leftState = plus[axis][layout.index(left)];
					State const& rightState = minus[axis][layout.index(right)];
					Flux highOrderFlux = system_.riemann(leftState, rightState, axis);
					highOrderFlux = system_.limitFlux(patch.atStorage(left), patch.atStorage(right), highOrderFlux, axis, stepSize / patch.cellWidth());
					fluxes[axis][layout.faceIndex(axis, face)] = highOrderFlux;
				});
			}
		}

		profiling::Region updateProfile("transport.update");
		layout.forEachInterior([&](mesh::Coordinates const& cell, std::size_t) {
			Flux update{};
			for (int axis = 0; axis < ndim; ++axis) {
				mesh::Coordinates lowerFace = cell;
				mesh::Coordinates upperFace = cell;
				++upperFace[axis];
				update += fluxes[axis][layout.faceIndex(axis, lowerFace)] - fluxes[axis][layout.faceIndex(axis, upperFace)];
			}
			State candidate = patch.atInterior(cell) + (stepSize / patch.cellWidth()) * update;
			candidate = system_.correctRoundoff(candidate, componentAbs((stepSize / patch.cellWidth()) * update));
			if (!system_.admissible(candidate)) {
				throw std::runtime_error("MUSCL-Hancock update produced an inadmissible state");
			}
			write(cell, candidate);
		});
	}

private:

	System system_;
	Limiter limiter_;
	Real theta_;

	/// Limit directional slopes, predict with the half-step unsplit divergence, and
	/// fall back to the cell center if a predicted face would be inadmissible.
	template <typename Patch>
	void predictFaceStates(
		Patch const& patch, units::Time stepSize, std::array<std::vector<State>, ndim>& minus, std::array<std::vector<State>, ndim>& plus) const {
		mesh::MeshLayout const& layout = patch.layout();
		int const ghostWidth = layout.ghostWidth();
		auto const begin = mesh::filledCoordinates(ghostWidth - 1);
		auto const end = mesh::filledCoordinates(ghostWidth + layout.cellsPerActiveDimension() + 1);

		mesh::forEachCoordinate(begin, end, [&](mesh::Coordinates const& cell) {
			std::size_t const cellIndex = layout.index(cell);
			State const centerState = patch.atStorage(cell);
			Reconstruction const center = system_.reconstructionVariables(centerState);
			std::array<Reconstruction, ndim> slopes{};
			for (int axis = 0; axis < ndim; ++axis) {
				mesh::Coordinates left = cell;
				mesh::Coordinates right = cell;
				--left[axis];
				++right[axis];
				Reconstruction const leftState = system_.reconstructionVariables(patch.atStorage(left));
				Reconstruction const rightState = system_.reconstructionVariables(patch.atStorage(right));
				slopes[axis].forEach([&](auto field, auto& slope) {
					slope = limitedSlope(center.template get<field>() - leftState.template get<field>(),
						rightState.template get<field>() - center.template get<field>(), limiter_, theta_);
				});
			}

			Real slopeFraction = 1;
			auto validFaces = [&](Real fraction) {
				for (int axis = 0; axis < ndim; ++axis) {
					State const lower = system_.conservedState(center - Real(0.5) * fraction * slopes[axis]);
					State const upper = system_.conservedState(center + Real(0.5) * fraction * slopes[axis]);
					if (!system_.admissible(lower) || !system_.admissible(upper)) {
						return false;
					}
				}
				return true;
			};
			if (!validFaces(slopeFraction)) {
				Real low = 0;
				Real high = 1;
				for (int iteration = 0; iteration < 48; ++iteration) {
					Real const fraction = Real(0.5) * (low + high);
					if (validFaces(fraction)) {
						low = fraction;
					} else {
						high = fraction;
					}
				}
				slopeFraction = low;
			}

			Flux predictorFlux{};
			for (int axis = 0; axis < ndim; ++axis) {
				minus[axis][cellIndex] = system_.conservedState(center - Real(0.5) * slopeFraction * slopes[axis]);
				plus[axis][cellIndex] = system_.conservedState(center + Real(0.5) * slopeFraction * slopes[axis]);
				predictorFlux += system_.physicalFlux(minus[axis][cellIndex], axis) - system_.physicalFlux(plus[axis][cellIndex], axis);
			}
			State const predictor = (Real(0.5) * stepSize / patch.cellWidth()) * predictorFlux;
			bool validPrediction = true;
			for (int axis = 0; axis < ndim; ++axis) {
				minus[axis][cellIndex] += predictor;
				plus[axis][cellIndex] += predictor;
				validPrediction = validPrediction && system_.admissible(minus[axis][cellIndex]) && system_.admissible(plus[axis][cellIndex]);
			}
			if (!validPrediction) {
				for (int axis = 0; axis < ndim; ++axis) {
					minus[axis][cellIndex] = centerState;
					plus[axis][cellIndex] = centerState;
				}
			}
		});
	}
};


}	 // namespace octotigerII::physics
